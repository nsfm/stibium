#include <algorithm>
#include <cstring>
#include <thread>
#include <unordered_map>

#include "fab/tree/triangulate.h"
#include "fab/tree/triangulate/mesher.h"

#include "fab/tree/tape.h"
#include "fab/tree/tree.h"
#include "fab/util/switches.h"

// Finds an array of vertices (as x,y,z float triplets).
// Sets *count to the number of vertices returned.
void triangulate(MathTree* tree, const Region r,
                 bool detect_edges, volatile int* halt,
                 float** const verts, unsigned* const count)
{
    Deck* deck = deck_from_tree(tree);
    Mesher t(deck, detect_edges, halt);

    // Top-level call to the recursive triangulation function.
    t.triangulate_region(r);

    // Copy data from the mesher to the output pointers.
    *verts = t.get_verts(count);
    deck_free(deck);
}

void triangulate_indexed(MathTree* tree, const Region r,
                         bool detect_edges, volatile int* halt,
                         std::vector<float>& verts,
                         std::vector<uint32_t>& indices,
                         std::atomic<uint64_t>* progress)
{
    Deck* deck = deck_from_tree(tree);
    Mesher t(deck, detect_edges, halt, progress);
    t.triangulate_region(r);
    t.get_mesh(verts, indices);
    deck_free(deck);
}

////////////////////////////////////////////////////////////////////////////////

namespace {

// FNV-1a over three packed floats (same keying as Mesher::intern)
struct VertHash {
    size_t operator()(const std::array<uint32_t, 3>& k) const {
        uint64_t h = 0xcbf29ce484222325ull;
        for (uint32_t w : k)
            h = (h ^ w) * 0x100000001b3ull;
        return h;
    }
};

uint64_t edge_key(uint32_t a, uint32_t b)
{
    return (uint64_t(a) << 32) | b;
}

/*
 *  Removes duplicate triangles (same three vertices in any winding),
 *  keeping the first copy - the global equivalent of
 *  Mesher::remove_dupes.
 */
void dedupe_indexed(std::vector<uint32_t>& idx)
{
    struct Key {
        uint32_t v[3];
        size_t seq;
    };
    const size_t nt = idx.size() / 3;
    std::vector<Key> keys;
    keys.reserve(nt);
    for (size_t t = 0; t < nt; ++t)
    {
        Key k = {{idx[t*3], idx[t*3 + 1], idx[t*3 + 2]}, t};
        std::sort(std::begin(k.v), std::end(k.v));
        keys.push_back(k);
    }

    std::sort(keys.begin(), keys.end(), [](const Key& a, const Key& b) {
        if (a.v[0] != b.v[0])   return a.v[0] < b.v[0];
        if (a.v[1] != b.v[1])   return a.v[1] < b.v[1];
        if (a.v[2] != b.v[2])   return a.v[2] < b.v[2];
        return a.seq < b.seq;
    });

    std::vector<char> dead(nt, 0);
    for (size_t i = 1; i < keys.size(); ++i)
        if (keys[i].v[0] == keys[i-1].v[0] &&
            keys[i].v[1] == keys[i-1].v[1] &&
            keys[i].v[2] == keys[i-1].v[2])
        {
            dead[keys[i].seq] = 1;
            keys[i] = keys[i-1];
        }

    size_t out = 0;
    for (size_t t = 0; t < nt; ++t)
    {
        if (dead[t])
            continue;
        idx[out*3]     = idx[t*3];
        idx[out*3 + 1] = idx[t*3 + 1];
        idx[out*3 + 2] = idx[t*3 + 2];
        out++;
    }
    idx.resize(out * 3);
}

/*
 *  Removes triangles with edges that aren't connected to the rest of
 *  the mesh - the global equivalent of Mesher::prune_flags (single
 *  snapshot pass: edges of pruned triangles still count as present).
 */
void prune_indexed(std::vector<uint32_t>& idx)
{
    const size_t nt = idx.size() / 3;
    std::vector<uint64_t> edges;
    edges.reserve(nt * 3);
    for (size_t t = 0; t < nt; ++t)
    {
        edges.push_back(edge_key(idx[t*3],     idx[t*3 + 1]));
        edges.push_back(edge_key(idx[t*3 + 1], idx[t*3 + 2]));
        edges.push_back(edge_key(idx[t*3 + 2], idx[t*3]));
    }
    std::sort(edges.begin(), edges.end());

    const auto has = [&edges](uint64_t e) {
        return std::binary_search(edges.begin(), edges.end(), e);
    };

    size_t out = 0;
    for (size_t t = 0; t < nt; ++t)
    {
        if (has(edge_key(idx[t*3 + 1], idx[t*3])) &&
            has(edge_key(idx[t*3 + 2], idx[t*3 + 1])) &&
            has(edge_key(idx[t*3],     idx[t*3 + 2])))
        {
            idx[out*3]     = idx[t*3];
            idx[out*3 + 1] = idx[t*3 + 1];
            idx[out*3 + 2] = idx[t*3 + 2];
            out++;
        }
    }
    idx.resize(out * 3);
}

/*
 *  Drops unreferenced vertices, renumbering the rest in first-use
 *  order (the numbering get_mesh would have produced).
 */
void compact_indexed(std::vector<float>& verts, std::vector<uint32_t>& idx)
{
    std::vector<uint32_t> remap(verts.size() / 3, UINT32_MAX);
    std::vector<float> packed;
    packed.reserve(verts.size());
    for (auto& i : idx)
    {
        if (remap[i] == UINT32_MAX)
        {
            remap[i] = packed.size() / 3;
            packed.insert(packed.end(),
                          verts.begin() + i*3, verts.begin() + i*3 + 3);
        }
        i = remap[i];
    }
    verts = std::move(packed);
}

}  // namespace

void triangulate_indexed_mt(MathTree* tree, const Region r,
                            bool detect_edges, volatile int* halt,
                            std::vector<float>& verts,
                            std::vector<uint32_t>& indices,
                            int threads,
                            std::atomic<uint64_t>* progress)
{
    if (threads <= 0)
    {
        const unsigned hw = std::thread::hardware_concurrency();
        threads = hw ? int(hw) : 1;
    }

    // Small regions aren't worth the thread + clone overhead.
    if (threads < 2 || r.voxels < uint64_t(MIN_VOLUME) * 64)
    {
        triangulate_indexed(tree, r, detect_edges, halt,
                            verts, indices, progress);
        return;
    }

    // Chunk the region: more chunks than threads, so a chunk that
    // lands on empty space doesn't leave its worker idle.
    std::vector<Region> chunks(size_t(threads) * 4);
    const int n_chunks = split(r, chunks.data(), chunks.size());
    chunks.resize(n_chunks);

    // One immutable deck serves every worker; each Mesher brings its
    // own evaluation workspace (no per-thread tree clones).
    Deck* deck = deck_from_tree(tree);

    struct ChunkMesh {
        std::vector<float> verts;
        std::vector<uint32_t> indices;
    };
    std::vector<ChunkMesh> results(threads);
    std::atomic<int> next(0);

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t)
    {
        pool.emplace_back([&, t]() {
            Mesher m(deck, detect_edges, halt, progress);
            int c;
            while ((c = next.fetch_add(1)) < n_chunks && !*halt)
                m.triangulate_region(chunks[c]);

            // Dedup / pruning must wait for the merged mesh: a seam
            // triangle's reverse edge lives in another worker's chunk.
            m.get_mesh(results[t].verts, results[t].indices, false);
        });
    }
    for (auto& th : pool)
        th.join();
    deck_free(deck);

    // Merge the per-worker meshes, re-interning the seam vertices
    // (workers that share a chunk boundary emit identical floats for
    // the shared vertices, so bit-pattern keying welds them).
    verts.clear();
    indices.clear();
    std::unordered_map<std::array<uint32_t, 3>, uint32_t, VertHash> ids;
    for (const auto& cm : results)
    {
        const size_t nv = cm.verts.size() / 3;
        std::vector<uint32_t> remap(nv);
        for (size_t v = 0; v < nv; ++v)
        {
            std::array<uint32_t, 3> key;
            memcpy(key.data(), &cm.verts[v*3], 3 * sizeof(float));

            auto found = ids.find(key);
            if (found != ids.end())
            {
                remap[v] = found->second;
            }
            else
            {
                const uint32_t id = verts.size() / 3;
                verts.insert(verts.end(),
                             cm.verts.begin() + v*3,
                             cm.verts.begin() + v*3 + 3);
                ids.emplace(key, id);
                remap[v] = id;
            }
        }
        for (auto i : cm.indices)
            indices.push_back(remap[i]);
    }

    // The global versions of the passes Mesher::finalize runs when
    // meshing single-threaded.
    if (detect_edges)
    {
        dedupe_indexed(indices);
        prune_indexed(indices);
        compact_indexed(verts, indices);
    }
}

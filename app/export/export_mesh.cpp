#include <Python.h>

#include <atomic>
#include <chrono>
#include <cstring>
#include <thread>
#include <vector>

#include <QFileDialog>
#include <QFileInfo>
#include <QRegularExpression>
#include <QStringList>
#include <QDateTime>

#include "app/settings.h"
#include <QMessageBox>

#include "export/export_mesh.h"
#include "vendor/meshoptimizer/meshoptimizer.h"

#include "dialog/resolution.h"

#include "fab/util/region.h"
#include "fab/tree/triangulate.h"
#include "fab/tree/triangulate/delaunay.h"
#include "fab/formats/stl.h"
#include "fab/formats/threemf.h"

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::run()
{
    // Sanity-check bounds
    if (std::isinf(bounds.xmin) || std::isinf(bounds.xmax) ||
        std::isinf(bounds.ymin) || std::isinf(bounds.ymax) ||
        std::isinf(bounds.zmin) || std::isinf(bounds.zmax))
    {
        QMessageBox::critical(NULL, "Export error",
                "<b>Export error:</b><br>"
                "Target shape has invalid (infinite) bounds");
        return;
    }

    // Get resolution, either hardcoded or from the user
    if (resolution == -1)
    {
        auto resolution_dialog = new ResolutionDialog(
                bounds, RESOLUTION_DIALOG_3D, UNITLESS, 1 << 22, NULL, 1);
        if (!resolution_dialog->exec())
            return;
        _resolution = resolution_dialog->getResolution();
        _detect_features = resolution_dialog->getDetectFeatures();
        _simplify = resolution_dialog->getSimplifyDeviation();
        _mesher = resolution_dialog->getMesher();
        _adv_autodense = resolution_dialog->getAutodense();
        _adv_density_cap = resolution_dialog->getDensityCap();
        _adv_decimate = resolution_dialog->getDecimate();
        _adv_snap = resolution_dialog->getSnap();
        _adv_stall = resolution_dialog->getStallPatience();
        _from_dialog = true;
        delete resolution_dialog;
    }
    else
    {
        _resolution = resolution;
        _detect_features = detect_features;
        _simplify = simplify;
        _mesher = getenv("STIBIUM_EXPORT_CLASSIC") ? 1 : 0;
    }

    if (_resolution == 0)
    {
        QMessageBox::critical(NULL, "Export error",
                "<b>Export error:</b><br>"
                "Resolution cannot be set to 0");
        return;
    }

    //  Get a target filename, either hardcoded or from the user
    if (filename.isEmpty())
    {
        QString filter = "3MF (*.3mf)";
        _filename = QFileDialog::getSaveFileName(
                NULL, "Export mesh",
                Settings::get("files/last_export_dir", "").toString(),
                "3MF (*.3mf);;STL (*.stl)", &filter);

        // If no recognized extension was typed, take it from the filter
        if (!_filename.isEmpty() &&
            !_filename.endsWith(".3mf", Qt::CaseInsensitive) &&
            !_filename.endsWith(".stl", Qt::CaseInsensitive))
        {
            _filename += filter.contains("*.stl") ? ".stl" : ".3mf";
        }
    }
    else
    {
        _filename = filename;
    }
    if (_filename.isEmpty())
        return;
    Settings::set("files/last_export_dir",
                  QFileInfo(_filename).absolutePath());

    if (checkWritable())
    {
        runAsync();
        /*  Post-export report (Nate's stats-dialog ask,
         *  2026-07-18): async() fills _stats; we are back on the
         *  GUI thread here.  Silent when cancelled or classic.  */
        if (!_stats.isEmpty() && !halt)
            QMessageBox::information(NULL, "Export complete",
                                     _stats);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool ExportMeshWorker::runHeadless(const QString& fname, float res,
                                   int detect)
{
    if (std::isinf(bounds.xmin) || std::isinf(bounds.xmax) ||
        std::isinf(bounds.ymin) || std::isinf(bounds.ymax) ||
        std::isinf(bounds.zmin) || std::isinf(bounds.zmax))
    {
        fprintf(stderr, "export: shape has infinite bounds\n");
        return false;
    }

    _filename = fname.isEmpty() ? filename : fname;
    _resolution = res > 0 ? res : resolution;
    _detect_features = detect < 0 ? detect_features : bool(detect);
    _simplify = simplify;
    /*  Headless default is Stibnite (integration 2026-07-18);
     *  STIBIUM_EXPORT_CLASSIC=1 routes to the classic sampler.
     *  STIBIUM_EXPORT_DMESH=1 remains accepted (now a no-op) so
     *  every documented harness command keeps working.  Knobs
     *  stay env-driven on this path.  */
    _mesher = getenv("STIBIUM_EXPORT_CLASSIC") ? 1 : 0;
    _from_dialog = false;

    if (_filename.isEmpty())
    {
        fprintf(stderr, "export: no filename (pass --export FILE or "
                        "set filename= in the script)\n");
        return false;
    }
    if (_resolution <= 0)
    {
        fprintf(stderr, "export: no resolution (pass --resolution R or "
                        "set resolution= in the script)\n");
        return false;
    }

    async();
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::async()
{
    Region r = {};
    r.ni = uint32_t((bounds.xmax - bounds.xmin) * _resolution);
    r.nj = uint32_t((bounds.ymax - bounds.ymin) * _resolution);
    r.nk = uint32_t((bounds.zmax - bounds.zmin) * _resolution);
    r.voxels = uint64_t(r.ni) * r.nj * r.nk;

    build_arrays(
            &r, bounds.xmin, bounds.ymin, bounds.zmin,
                bounds.xmax, bounds.ymax, bounds.zmax);

    // The mesher produces an indexed mesh directly (unique vertices +
    // three indices per triangle), so no welding pass is needed.
    progress_phase = PHASE_MESHING;
    progress_total = r.voxels;
    std::vector<float> verts;
    std::vector<uint32_t> indices;
    _stats.clear();
    /*  Stibnite (the adaptive-Delaunay mesher, doc/MESH-WAR.md) is
     *  the default; Classic is the dialog's other choice or
     *  STIBIUM_EXPORT_CLASSIC=1 headless.  */
    if (_mesher == 0)
    {
        /*  Dialog knobs -> process environment.  Safe: one export
         *  at a time is law, and the headless path never enters
         *  here with _from_dialog set.  */
        if (_from_dialog)
        {
            qputenv("STIBIUM_DMESH_AUTODENSE",
                    _adv_autodense ? "1" : "0");
            qputenv("STIBIUM_DMESH_AUTODENSE_MAX",
                    QByteArray::number(_adv_density_cap));
            qputenv("STIBIUM_DMESH_DECIMATE",
                    _adv_decimate ? "1" : "0");
            qputenv("STIBIUM_DMESH_SNAP", _adv_snap ? "1" : "0");
            qputenv("STIBIUM_DMESH_STALL",
                    QByteArray::number(_adv_stall));
        }

        /*  Mirror the mesher's progress sink into the dialog's
         *  bar until meshing returns: overall fraction, stage
         *  index, and a remaining-time estimate.
         *
         *  The raw estimate (elapsed x (1-f)/f) whipsaws while the
         *  early stages burn through their phase weights, so the
         *  display runs on a LOCKED DEADLINE: past 15% overall the
         *  smoothed estimate sets a finish time and the label just
         *  counts down toward it (smooth by construction).  The
         *  deadline only moves for sustained disagreement - pushed
         *  out when the smoothed estimate exceeds the countdown by
         *  25%, pulled in when it drops below half.  Overshooting
         *  the deadline shows "wrapping up" (eta sentinel -2)
         *  instead of a lie.  */
        progress_total = 1000;
        std::atomic<bool> mesh_done{ false };
        std::thread mirror([&]() {
            const auto m0 = std::chrono::steady_clock::now();
            double ema = -1;
            double deadline = -1;   // seconds since m0
            while (!mesh_done.load())
            {
                const float f = dmesh_progress()->overall.load();
                progress_done = uint64_t(f * 1000);
                progress_stage = dmesh_progress()->stage.load();
                const double el = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - m0)
                        .count();
                if (f > 0.15f && f < 1.f)
                {
                    const double raw = el * (1.f - f) / f;
                    ema = ema < 0 ? raw : 0.9 * ema + 0.1 * raw;
                    const double left0 = deadline - el;
                    if (deadline < 0 || ema > left0 * 1.25 ||
                        ema < left0 * 0.5)
                        deadline = el + ema;
                    const double left = deadline - el;
                    progress_eta_s = left >= 1 ? int(left) : -2;
                }
                else
                    progress_eta_s = -1;
                std::this_thread::sleep_for(
                        std::chrono::milliseconds(100));
            }
            progress_stage = -1;
            progress_eta_s = -1;
        });

        const auto t0 = std::chrono::steady_clock::now();
        Deck* deck = deck_from_tree(shape.tree.get());
        DMesh m;
        const bool ok = delaunay_mesh(deck, r, &halt, &m);
        deck_free(deck);
        mesh_done = true;
        mirror.join();
        const double wall = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();

        fprintf(stderr, "stibnite mesh report: %s, %zu tris, "
                "%llu open, "
                "%llu non-manifold, %llu constrained, %llu steiner, "
                "%llu repaired, %llu split, %llu snapped\n",
                ok ? "ok" : "FAILED (built without CGAL?)",
                m.tris.size() / 3,
                (unsigned long long)m.open_edges,
                (unsigned long long)m.nonmanifold_edges,
                (unsigned long long)m.constrained,
                (unsigned long long)m.steiner,
                (unsigned long long)m.repaired,
                (unsigned long long)m.split_verts,
                (unsigned long long)m.snapped);
        if (ok && !halt)
            _stats = QString(
                    "<b>Stibnite mesh report</b><br><br>"
                    "%1 triangles, %2 vertices<br>"
                    "open edges: %3<br>"
                    "feature edges constrained: %4<br>"
                    "repairs: %5 &nbsp; crease snaps: %6<br>"
                    "meshing time: %7 s")
                    .arg(m.tris.size() / 3)
                    .arg(m.verts.size() / 3)
                    .arg(m.open_edges)
                    .arg(m.constrained)
                    .arg(m.repaired)
                    .arg(m.snapped)
                    .arg(wall, 0, 'f', 1);
        verts = std::move(m.verts);
        indices = std::move(m.tris);
    }
    else
        triangulate_indexed_mt(shape.tree.get(), r, _detect_features,
                               &halt, verts, indices, -1,
                               &progress_done);

    // Simplification and file writing can't report granular progress;
    // drop the bar back to its busy animation instead of a full bar.
    progress_total = 0;

    /*  QEM on the dmesh path (quadric assessment 2026-07-17):
     *  STIBIUM_DMESH_SIMPLIFY=<mm> forces the vendored meshopt
     *  pass on dmesh exports - deviation-bounded, in model
     *  units.  Nate's standalone experiment (50% off, no visual
     *  defects, twice) green-lit QEM for this geometry; this
     *  delivers it from code already in the tree.  Phase 2
     *  (crease vertex_lock) is designed in doc/reviews/.  */
    float simplify = _simplify;
    if (_mesher == 0)
        if (const char* se = getenv("STIBIUM_DMESH_SIMPLIFY"))
            simplify = float(atof(se));
    if (simplify > 0 && indices.size() >= 3 && !halt)
    {
        progress_phase = PHASE_SIMPLIFYING;
        simplifyMesh(simplify, verts, indices);
        if (!_stats.isEmpty())
            _stats += QString("<br>after simplify: %1 triangles, "
                              "%2 vertices")
                    .arg(indices.size() / 3)
                    .arg(verts.size() / 3);
    }

    /*  Pinch split LAST (strictly after QEM): meshopt tears
     *  coincident sheet copies into real boundary holes if the
     *  split runs first, and its collapses mint fresh pinches of
     *  their own (bino 274 -> 325 nm edges, measured) - the tail
     *  position cures both.  Zero vertex motion; index-level
     *  only, which 3MF preserves.  */
    if (_mesher == 0 && indices.size() >= 3 && !halt)
    {
        const uint64_t copies = dmesh_split_pinches(verts, indices);
        if (copies && !_stats.isEmpty())
            _stats += QString("<br>pinch seams split: %1 sheet "
                              "copies").arg(copies);
    }

    progress_phase = PHASE_WRITING;

    /*  Provenance metadata (Nate's peace-of-mind ask,
     *  2026-07-20): mesher, resolution, date, and every
     *  STIBIUM_* environment override, stamped into the export
     *  so past meshes explain themselves.  3MF carries the full
     *  string; STL's 80-byte header gets the short form.  */
    QString prov = QString("mesher=%1 resolution=%2 detect=%3 "
                           "simplify=%4 date=%5")
            .arg(_mesher == 0 ? "stibnite" : "classic")
            .arg(_resolution)
            .arg(_detect_features ? 1 : 0)
            .arg(_simplify)
            .arg(QDateTime::currentDateTime().toString(Qt::ISODate));
    extern char** environ;
    for (char** e = environ; *e; ++e)
        if (strncmp(*e, "STIBIUM_", 8) == 0)
            prov += QString(" %1").arg(*e);
    {
        const QString plain =
                QString(_stats).replace("<br>", "; ")
                               .replace("&nbsp;", " ")
                               .replace(QRegularExpression("<[^>]*>"),
                                        "");
        if (!plain.isEmpty())
            prov += QString("; stats: %1").arg(plain);
    }

    // Format follows the file extension; STL is the fallback for
    // scripted exports with unrecognized names (the legacy behavior).
    if (_filename.endsWith(".3mf", Qt::CaseInsensitive))
        save_3mf_indexed(verts.data(), verts.size() / 3,
                         indices.data(), indices.size() / 3,
                         _filename.toStdString().c_str(),
                         prov.toStdString().c_str());
    else
    {
        const QString stamp = QString("Stibium stibnite r=%1 %2")
                .arg(_resolution)
                .arg(QDateTime::currentDateTime()
                             .toString("yyyy-MM-dd"));
        save_stl_indexed_stamped(verts.data(), indices.data(),
                                 indices.size() / 3,
                                 _filename.toStdString().c_str(),
                                 stamp.toStdString().c_str());
    }
    free_arrays(&r);
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::simplifyMesh(float deviation,
                                    std::vector<float>& verts,
                                    std::vector<uint32_t>& indices) const
{
    const size_t vertex_count = verts.size() / 3;

    // target_error is relative to mesh extent; deviation is in model units
    const float scale = meshopt_simplifyScale(
            verts.data(), vertex_count, 3 * sizeof(float));
    const float target_error = scale > 0 ? deviation / scale : deviation;

    std::vector<uint32_t> simplified(indices.size());
    const size_t out_index_count = meshopt_simplify(
            simplified.data(), indices.data(), indices.size(),
            verts.data(), vertex_count, 3 * sizeof(float),
            0, target_error, meshopt_SimplifyPrune, NULL);

    simplified.resize(out_index_count);
    indices = std::move(simplified);

    // Drop the vertices orphaned by simplification (renumbering the
    // rest in first-use order), so indexed formats don't carry dead
    // entries.
    std::vector<uint32_t> remap(vertex_count, UINT32_MAX);
    std::vector<float> packed;
    packed.reserve(verts.size());
    for (auto& i : indices)
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

////////////////////////////////////////////////////////////////////////////////

/*  --facedev support: load an EXISTING mesh (binary STL or 3MF)
 *  and sweep it against this worker's tape (dmesh_face_sweep) -
 *  the offline referee that lets past exports be ranked with the
 *  metric that tracks the eye.  */

#include <zlib.h>

#include "fab/tree/tape.h"

static bool load_stl_mesh(const QString& path,
                          std::vector<float>& verts,
                          std::vector<uint32_t>& tris)
{
    FILE* f = fopen(path.toLocal8Bit().constData(), "rb");
    if (!f)
        return false;
    char hdr[80];
    if (fread(hdr, 1, 80, f) != 80)
        { fclose(f); return false; }
    uint32_t n = 0;
    if (fread(&n, 4, 1, f) != 1)
        { fclose(f); return false; }
    verts.reserve(size_t(n) * 9);
    tris.reserve(size_t(n) * 3);
    for (uint32_t t = 0; t < n; ++t)
    {
        float rec[12];
        uint16_t attr;
        if (fread(rec, 4, 12, f) != 12 ||
            fread(&attr, 2, 1, f) != 1)
            { fclose(f); return false; }
        for (int v = 0; v < 3; ++v)
        {
            tris.push_back(uint32_t(verts.size() / 3));
            verts.push_back(rec[3 + v*3]);
            verts.push_back(rec[3 + v*3 + 1]);
            verts.push_back(rec[3 + v*3 + 2]);
        }
    }
    fclose(f);
    return true;
}

/*  Minimal central-directory zip walk + raw-deflate inflate:
 *  reads the 3D/*.model entry (the OPC part holding the mesh and
 *  our metadata) of any spec-conforming 3MF, including our own
 *  writer's, into `model`.  Shared by the geometry loader and the
 *  --provenance reader.  */
static bool read_3mf_model_xml(const QString& path, std::string& model)
{
    FILE* f = fopen(path.toLocal8Bit().constData(), "rb");
    if (!f)
        return false;
    fseek(f, 0, SEEK_END);
    const long fsz = ftell(f);
    const long tail = std::min(fsz, long(65557));
    std::vector<unsigned char> tb(tail);
    fseek(f, fsz - tail, SEEK_SET);
    if (fread(tb.data(), 1, tail, f) != size_t(tail))
        { fclose(f); return false; }
    long eocd = -1;
    for (long i = tail - 22; i >= 0; --i)
        if (tb[i] == 0x50 && tb[i+1] == 0x4b &&
            tb[i+2] == 0x05 && tb[i+3] == 0x06)
            { eocd = i; break; }
    if (eocd < 0)
        { fclose(f); return false; }
    const auto rd16 = [&](long o) {
        return uint32_t(tb[o]) | (uint32_t(tb[o+1]) << 8); };
    const auto rd32 = [&](long o) {
        return uint32_t(tb[o]) | (uint32_t(tb[o+1]) << 8) |
               (uint32_t(tb[o+2]) << 16) |
               (uint32_t(tb[o+3]) << 24); };
    const uint32_t nent = rd16(eocd + 10);
    uint32_t cdofs = rd32(eocd + 16);
    for (uint32_t e = 0; e < nent; ++e)
    {
        unsigned char ch[46];
        fseek(f, cdofs, SEEK_SET);
        if (fread(ch, 1, 46, f) != 46)
            break;
        const auto crd16 = [&](int o) {
            return uint32_t(ch[o]) | (uint32_t(ch[o+1]) << 8); };
        const auto crd32 = [&](int o) {
            return uint32_t(ch[o]) | (uint32_t(ch[o+1]) << 8) |
                   (uint32_t(ch[o+2]) << 16) |
                   (uint32_t(ch[o+3]) << 24); };
        const uint32_t method = crd16(10), csz = crd32(20),
                usz = crd32(24), nlen = crd16(28),
                xlen = crd16(30), clen = crd16(32),
                lofs = crd32(42);
        std::string name(nlen, 0);
        if (fread(name.data(), 1, nlen, f) != nlen)
            break;
        cdofs += 46 + nlen + xlen + clen;
        if (name.size() < 6 ||
            name.compare(name.size() - 6, 6, ".model") != 0)
            continue;
        /*  Local header: skip its (possibly different) name/extra
         *  lengths.  */
        unsigned char lh[30];
        fseek(f, lofs, SEEK_SET);
        if (fread(lh, 1, 30, f) != 30)
            break;
        const uint32_t lnlen = uint32_t(lh[26]) |
                (uint32_t(lh[27]) << 8);
        const uint32_t lxlen = uint32_t(lh[28]) |
                (uint32_t(lh[29]) << 8);
        fseek(f, lofs + 30 + lnlen + lxlen, SEEK_SET);
        std::vector<unsigned char> cbuf(csz);
        if (fread(cbuf.data(), 1, csz, f) != csz)
            break;
        model.resize(usz);
        if (method == 0)
            model.assign(cbuf.begin(), cbuf.end());
        else if (method == 8)
        {
            z_stream st = {};
            if (inflateInit2(&st, -MAX_WBITS) != Z_OK)
                break;
            st.next_in = cbuf.data();
            st.avail_in = csz;
            st.next_out =
                    reinterpret_cast<unsigned char*>(&model[0]);
            st.avail_out = usz;
            const int zr = inflate(&st, Z_FINISH);
            inflateEnd(&st);
            if (zr != Z_STREAM_END)
                model.clear();
        }
        break;
    }
    fclose(f);
    return !model.empty();
}

static bool load_3mf_mesh(const QString& path,
                          std::vector<float>& verts,
                          std::vector<uint32_t>& tris)
{
    std::string model;
    if (!read_3mf_model_xml(path, model))
        return false;
    /*  Light scan of the model XML: vertices in document order,
     *  triangles by index (matching our writer and the spec).  */
    const char* p = model.c_str();
    while ((p = strstr(p, "<vertex ")))
    {
        const char* px = strstr(p, "x=\"");
        const char* py = strstr(p, "y=\"");
        const char* pz = strstr(p, "z=\"");
        if (!px || !py || !pz)
            break;
        verts.push_back(strtof(px + 3, nullptr));
        verts.push_back(strtof(py + 3, nullptr));
        verts.push_back(strtof(pz + 3, nullptr));
        p += 8;
    }
    p = model.c_str();
    while ((p = strstr(p, "<triangle ")))
    {
        const char* v1 = strstr(p, "v1=\"");
        const char* v2 = strstr(p, "v2=\"");
        const char* v3 = strstr(p, "v3=\"");
        if (!v1 || !v2 || !v3)
            break;
        tris.push_back(uint32_t(strtoul(v1 + 4, nullptr, 10)));
        tris.push_back(uint32_t(strtoul(v2 + 4, nullptr, 10)));
        tris.push_back(uint32_t(strtoul(v3 + 4, nullptr, 10)));
        p += 10;
    }
    return !verts.empty() && !tris.empty();
}

int ExportMeshWorker::facedevHeadless(const QString& mesh,
                                      float res)
{
    const float r2 = res > 0 ? res : resolution;
    if (r2 <= 0)
    {
        fprintf(stderr, "facedev: no resolution (pass "
                        "--resolution R)\n");
        return 1;
    }
    std::vector<float> mverts;
    std::vector<uint32_t> mtris;
    const bool ok = mesh.endsWith(".3mf", Qt::CaseInsensitive)
            ? load_3mf_mesh(mesh, mverts, mtris)
            : load_stl_mesh(mesh, mverts, mtris);
    if (!ok)
    {
        fprintf(stderr, "facedev: could not read %s\n",
                mesh.toLocal8Bit().constData());
        return 1;
    }
    fprintf(stderr, "facedev: %zu tris from %s\n",
            mtris.size() / 3, mesh.toLocal8Bit().constData());
    Deck* deck = deck_from_tree(shape.tree.get());
    dmesh_face_sweep(deck, mverts, mtris, 1.f / r2,
                     getenv("STIBIUM_DMESH_FACE_DUMP"));
    deck_free(deck);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

/*  --provenance support: read only the mesh file's embedded stamp
 *  and print what it records about how the mesh was made.  No .sb,
 *  no Shape, no sampling - staged exports with ambiguous names
 *  ("r2mint" - revision or resolution?) explain themselves from the
 *  header alone.  3MF carries the full config in its Description
 *  metadata; a binary STL carries the short form in its 80-byte
 *  header.  */

static std::string xml_unescape(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); )
    {
        if (s[i] == '&')
        {
            if (!s.compare(i, 5, "&amp;"))  { out += '&';  i += 5; continue; }
            if (!s.compare(i, 4, "&lt;"))   { out += '<';  i += 4; continue; }
            if (!s.compare(i, 4, "&gt;"))   { out += '>';  i += 4; continue; }
            if (!s.compare(i, 6, "&quot;")) { out += '"';  i += 6; continue; }
            if (!s.compare(i, 6, "&apos;")) { out += '\''; i += 6; continue; }
        }
        out += s[i++];
    }
    return out;
}

int meshProvenance(const QString& mesh)
{
    const QByteArray path = mesh.toLocal8Bit();

    if (mesh.endsWith(".3mf", Qt::CaseInsensitive))
    {
        std::string model;
        if (!read_3mf_model_xml(mesh, model))
        {
            fprintf(stderr, "provenance: could not read %s\n",
                    path.constData());
            return 1;
        }
        const char* key = "<metadata name=\"Description\">";
        const size_t a = model.find(key);
        if (a == std::string::npos)
        {
            printf("%s: 3MF carries no Stibium provenance "
                   "(pre-stamp or foreign file)\n", path.constData());
            return 0;
        }
        const size_t s = a + strlen(key);
        const size_t e = model.find("</metadata>", s);
        const QString desc = QString::fromStdString(xml_unescape(
                model.substr(s, e == std::string::npos ? e : e - s)));

        printf("%s: 3MF provenance\n", path.constData());
        /*  One line per recorded field: the leading run of
         *  space-separated key=value fields, then each "; "-delimited
         *  mesh-report fragment.  A description with no delimiters
         *  falls through as a single line.  */
        const QStringList segs = desc.split("; ");
        for (int i = 0; i < segs.size(); ++i)
        {
            if (i == 0)
                for (const QString& kv :
                     segs[0].split(' ', Qt::SkipEmptyParts))
                    printf("  %s\n", kv.toLocal8Bit().constData());
            else if (!segs[i].trimmed().isEmpty())
                printf("  %s\n", segs[i].trimmed()
                                       .toLocal8Bit().constData());
        }
        return 0;
    }

    // Binary STL: the provenance is the 80-byte header.
    FILE* f = fopen(path.constData(), "rb");
    if (!f)
    {
        fprintf(stderr, "provenance: could not open %s\n",
                path.constData());
        return 1;
    }
    char hdr[80];
    const size_t got = fread(hdr, 1, sizeof(hdr), f);
    fclose(f);
    if (got != sizeof(hdr))
    {
        fprintf(stderr, "provenance: %s is too small to be a binary "
                        "STL\n", path.constData());
        return 1;
    }

    // Trim trailing NUL / space padding, then sanitize non-printables.
    size_t n = sizeof(hdr);
    while (n && (hdr[n-1] == '\0' || hdr[n-1] == ' '))
        --n;
    std::string stamp;
    for (size_t i = 0; i < n; ++i)
        stamp += (hdr[i] >= 0x20 && hdr[i] < 0x7f) ? hdr[i] : '.';

    if (strncmp(hdr, "Stibium", 7))
    {
        printf("%s: STL header carries no Stibium stamp "
               "(foreign or pre-stamp export)\n"
               "  header: \"%s\"\n", path.constData(), stamp.c_str());
        return 0;
    }
    printf("%s: STL provenance stamp\n  %s\n",
           path.constData(), stamp.c_str());
    return 0;
}

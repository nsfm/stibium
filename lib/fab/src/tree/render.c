#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#include "fab/tree/tape.h"
#include "fab/tree/tree.h"
#include "fab/tree/render.h"
#include "fab/tree/math/math_g.h"

#include "fab/util/switches.h"

/*  region8
 *
 *  Renders a tape pixel-by-pixel into a given region,
 *  using tape_eval_r to find an array of results in
 *  a single pass through the tape.
 *
 */
static
void region8(const Tape* tape, TapeCtx* ctx, Region region, uint8_t** img);

/*  region16
 *
 *  Renders a tape pixel-by-pixel into a given region,
 *  using tape_eval_r to find an array of results in
 *  a single pass through the tape.
 *
 */
static
void region16(const Tape* tape, TapeCtx* ctx, Region region, uint16_t** img);

/*  MPR-style tile fan-out (STIBIUM_TILE_RENDER): subdivide ambiguous
 *  regions into many children per tape push instead of two, so the
 *  hierarchy is shallow and each push specializes against a region
 *  a whole fan-out factor smaller.  =1 selects the default fan-out;
 *  an explicit number (8/16/64/128) tunes it; unset keeps the
 *  original binary bisection.  */
#define MAX_TILE_FANOUT 128

static
int tile_fanout(void)
{
    static int fanout = -1;
    if (fanout == -1) {
        const char* env = getenv("STIBIUM_TILE_RENDER");
        fanout = env ? atoi(env) : 0;
        if (fanout == 1)                fanout = 64;
        if (fanout < 0)                 fanout = 0;
        if (fanout > MAX_TILE_FANOUT)   fanout = MAX_TILE_FANOUT;
    }
    return fanout;
}

/*  Near-to-far ordering (descending z), so filled tiles occlude the
 *  work behind them via the existing already-light cull check.  */
static
int cmp_kmin_desc(const void* a, const void* b)
{
    const uint32_t ka = ((const Region*)a)->kmin;
    const uint32_t kb = ((const Region*)b)->kmin;
    return ka < kb ? 1 : ka > kb ? -1 : 0;
}

////////////////////////////////////////////////////////////////////////////////
void render8_tape(Tape* tape, TapeCtx* ctx, Region region,
                  uint8_t** img, volatile int* halt,
                  void (*callback)())
{
    // Special interrupt system, set asynchronously by on high
    if (*halt)  return;

    // Render pixel-by-pixel if we're below a certain size.
    if (region.voxels > 0 && region.voxels < MIN_VOLUME) {
        if (callback)   (*callback)();
        region8(tape, ctx, region, img);
        return;
    }


    // Pre-emptively halt evaluation if all the points in this
    // region are already light.
    uint8_t L = region.L[region.nk] >> 8;
    bool cull = true;
    for (unsigned row = region.jmin; cull && row < region.jmin + region.nj; ++row) {
        for (unsigned col = region.imin; cull && col < region.imin + region.ni; ++col) {
            if (L > img[row][col]) {
                cull = false;
                break;
            }
        }
    }
    if (cull) return;

    Interval X = {region.X[0], region.X[region.ni]},
             Y = {region.Y[0], region.Y[region.nj]},
             Z = {region.Z[0], region.Z[region.nk]};

    Interval result = tape_eval_i(tape, ctx, X, Y, Z);

    // If we're inside the object, fill with color.
    if (result.upper < 0) {
        for (unsigned row = region.jmin; row < region.jmin + region.nj; ++row) {
            for (unsigned col = region.imin; col < region.imin + region.ni; ++col) {
                if (L > img[row][col])  img[row][col] = L;
            }
        }
    }

    // In unambiguous cases, return immediately
    if (result.upper < 0 || result.lower >= 0)  return;

#if PRUNE
    Tape* sub = tape_push(tape, ctx, X, Y, Z, TAPE_PUSH_BINARY);
#else
    Tape* sub = tape_retain(tape);
#endif

    // Subdivide and recurse if we're not at voxel size.
    if (region.ni*region.nj*region.nk > 1) {
        const int fanout = tile_fanout();
        if (fanout > 1 && region.voxels >= (uint64_t)MIN_VOLUME * 2) {
            Region kids[MAX_TILE_FANOUT];
            const int n = split(region, kids, fanout);
            qsort(kids, n, sizeof(Region), cmp_kmin_desc);
            for (int i = 0; i < n; ++i)
                render8_tape(sub, ctx, kids[i], img, halt, callback);
        } else {
            Region A, B;

            bisect(region, &A, &B);

            render8_tape(sub, ctx, B, img, halt, callback);
            render8_tape(sub, ctx, A, img, halt, callback);
        }
    }

    tape_release(sub);
}

void render8(MathTree* tree, Region region,
             uint8_t** img, volatile int* halt,
             void (*callback)())
{
    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    render8_tape(deck_base(deck), ctx, region, img, halt, callback);
    tape_ctx_free(ctx);
    deck_free(deck);
}

////////////////////////////////////////////////////////////////////////////////

static
void region8(const Tape* tape, TapeCtx* ctx, Region region, uint8_t** img)
{
    float *X = malloc(region.voxels*sizeof(float)),
          *Y = malloc(region.voxels*sizeof(float)),
          *Z = malloc(region.voxels*sizeof(float));

    // Copy the X, Y, Z vectors into a flattened matrix form.
    int q = 0;
    for (int k = region.nk - 1; k >= 0; --k) {
        for (unsigned j = 0; j < region.nj; ++j) {
            for (unsigned i = 0; i < region.ni; ++i) {
                X[q] = region.X[i];
                Y[q] = region.Y[j];
                Z[q] = region.Z[k];
                q++;
            }
        }
    }
    region.X = X;
    region.Y = Y;
    region.Z = Z;

    const float* result = tape_eval_r(tape, ctx, region);

    // Free the allocated matrices
    free(X);
    free(Y);
    free(Z);

    for (int k = region.nk - 1; k >= 0; --k) {
        uint8_t L = region.L[k+1] >> 8;

        for (unsigned j = 0; j < region.nj; ++j) {
            int row = j + region.jmin;

            for (unsigned i = 0; i < region.ni; ++i) {
                int col = i + region.imin;

               if (*(result++) < 0 && img[row][col] < L) {
                    img[row][col] = L;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

static
void get_normals8(const Tape* tape, TapeCtx* ctx,
                  float* restrict X, float* restrict Y, float* restrict Z,
                  unsigned count, float (*normals)[3])
{
    Region dummy;
    dummy.X = X;
    dummy.Y = Y;
    dummy.Z = Z;
    dummy.voxels = count;

    const derivative* result = tape_eval_g(tape, ctx, dummy);

    // Calculate normals and copy over.
    for (unsigned i=0; i < count; ++i)
    {
        const float x = result[i].dx;
        const float y = result[i].dy;
        const float z = result[i].dz;

        const float dist = sqrt(pow(x,2) + pow(y, 2) + pow(z,2));
        normals[i][0] = dist ? x/dist : 0;
        normals[i][1] = dist ? y/dist : 0;
        normals[i][2] = dist ? z/dist : 0;
    }
}

static
void shade_pixels8(unsigned count, float (*normals)[3],
                   unsigned* is, unsigned* js, uint8_t (**out)[3])
{
    for (unsigned a=0; a < count; ++a)
    {
        for (int b=0; b < 3; ++b)
        {
            out[js[a]][is[a]][b] = normals[a][b]*127 + 128;
        }
    }
}

/*  Shades one xy tile of the depth buffer.  The tape has already
 *  been specialized against the tile's bounding box, so every
 *  gradient batch walks a pruned tape instead of the full deck.  */
static
void shaded8_block(Tape* tape, TapeCtx* ctx, const Region* region,
                   unsigned i0, unsigned i1, unsigned j0, unsigned j1,
                   uint16_t** depth, uint8_t (**out)[3],
                   volatile int* halt,
                   float* X, float* Y, float* Z,
                   unsigned* is, unsigned* js, float (*normals)[3])
{
    unsigned count = 0;
    for (unsigned j = j0; j < j1 && !*halt; ++j)
    {
        for (unsigned i = i0; i < i1 && !*halt; ++i)
        {
            // Load this pixel into the set of pixels to render.
            // (Image buffers are indexed absolutely, matching
            // render16, so sub-regions shade the right pixels.)
            const unsigned row = region->jmin + j;
            const unsigned col = region->imin + i;
            if (depth[row][col])
            {
                X[count] = region->X[i];
                Y[count] = region->Y[j];
                Z[count] = region->Z[0] + depth[row][col] / 65535.0f *
                            (region->Z[region->nk] - region->Z[0]);

                is[count] = col;
                js[count] = row;
                count++;
            }

            if (count == MIN_VOLUME/4 ||
                    (count && j == j1 - 1 && i == i1 - 1))
            {
                get_normals8(tape, ctx, X, Y, Z, count, normals);
                shade_pixels8(count, normals, is, js, out);
                count = 0;
            }
        }
    }
}

/*  Shading tile edge in pixels: big enough that the per-tile
 *  interval eval + push amortizes over many gradient batches,
 *  small enough that the pruned tape is genuinely short.  */
#define SHADE_TILE 64

void shaded8_tape(Tape* tape, TapeCtx* ctx, Region region,
                  uint16_t** depth, uint8_t (**out)[3],
                  volatile int* halt, void (*callback)())
{
    float *X = malloc(MIN_VOLUME*sizeof(float)),
          *Y = malloc(MIN_VOLUME*sizeof(float)),
          *Z = malloc(MIN_VOLUME*sizeof(float));

    unsigned *is = malloc(MIN_VOLUME*sizeof(unsigned));
    unsigned *js = malloc(MIN_VOLUME*sizeof(unsigned));

    float (*normals)[3] = malloc(MIN_VOLUME*sizeof(float[3]));

    for (unsigned j0 = 0; j0 < region.nj && !*halt; j0 += SHADE_TILE)
    {
        if (callback)   (*callback)();
        const unsigned j1 = j0 + SHADE_TILE < region.nj
                ? j0 + SHADE_TILE : region.nj;

        for (unsigned i0 = 0; i0 < region.ni && !*halt; i0 += SHADE_TILE)
        {
            const unsigned i1 = i0 + SHADE_TILE < region.ni
                    ? i0 + SHADE_TILE : region.ni;

            // The tile's depth range, from the depth buffer itself
            uint16_t dmin = UINT16_MAX, dmax = 0;
            for (unsigned j = j0; j < j1; ++j)
                for (unsigned i = i0; i < i1; ++i)
                {
                    const uint16_t d = depth[region.jmin + j]
                                            [region.imin + i];
                    if (!d) continue;
                    if (d < dmin) dmin = d;
                    if (d > dmax) dmax = d;
                }
            if (!dmax)      // nothing to shade in this tile
                continue;

            // Specialize the tape against the tile's bounding box
            // (STANDARD push: gradients need exact distances).  A
            // little slack on z since sample depths are quantized.
            const float zspan = region.Z[region.nk] - region.Z[0];
            const Interval Xi = { region.X[i0], region.X[i1] },
                           Yi = { region.Y[j0], region.Y[j1] },
                           Zi = { region.Z[0] +
                                    (dmin / 65535.0f) * zspan - 1e-3f,
                                  region.Z[0] +
                                    (dmax / 65535.0f) * zspan + 1e-3f };
            tape_eval_i(tape, ctx, Xi, Yi, Zi);
            Tape* sub = tape_push(tape, ctx, Xi, Yi, Zi,
                                  TAPE_PUSH_STANDARD);

            shaded8_block(sub, ctx, &region, i0, i1, j0, j1,
                          depth, out, halt, X, Y, Z, is, js, normals);

            tape_release(sub);
        }
    }

    free(X);
    free(Y);
    free(Z);

    free(is);
    free(js);

    free(normals);
}

void shaded8(MathTree* tree, Region region, uint16_t** depth,
             uint8_t (**out)[3], volatile int* halt,
             void (*callback)())
{
    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    shaded8_tape(deck_base(deck), ctx, region, depth, out, halt, callback);
    tape_ctx_free(ctx);
    deck_free(deck);
}


////////////////////////////////////////////////////////////////////////////////
void render16_tape(Tape* tape, TapeCtx* ctx, Region region,
                   uint16_t** img, volatile int* halt,
                   void (*callback)())
{
    // Special interrupt system, set asynchronously by on high
    if (*halt)  return;

    // Render pixel-by-pixel if we're below a certain size.
    if (region.voxels > 0 && region.voxels < MIN_VOLUME) {
        if (callback)
            (*callback)();
        region16(tape, ctx, region, img);
        return;
    }

    // Pre-emptively halt evaluation if all the points in this
    // region are already light.
    uint16_t L = region.L[region.nk];
    bool cull = true;
    for (unsigned row = region.jmin; cull && row < region.jmin + region.nj; ++row) {
        for (unsigned col = region.imin; cull && col < region.imin + region.ni; ++col) {
            if (L > img[row][col]) {
                cull = false;
                break;
            }
        }
    }
    if (cull) return;

    Interval X = {region.X[0], region.X[region.ni]},
             Y = {region.Y[0], region.Y[region.nj]},
             Z = {region.Z[0], region.Z[region.nk]};

    Interval result = tape_eval_i(tape, ctx, X, Y, Z);

    // If we're inside the object, fill with color.
    if (result.upper < 0) {
        for (unsigned row = region.jmin; row < region.jmin + region.nj; ++row) {
            for (unsigned col = region.imin; col < region.imin + region.ni; ++col) {
                if (L > img[row][col])  img[row][col] = L;
            }
        }
    }

    // In unambiguous cases, return immediately
    if (result.upper < 0 || result.lower >= 0)  return;

#if PRUNE
    Tape* sub = tape_push(tape, ctx, X, Y, Z, TAPE_PUSH_BINARY);
#else
    Tape* sub = tape_retain(tape);
#endif

    // Subdivide and recurse if we're not at voxel size.
    if (region.ni*region.nj*region.nk > 1) {
        const int fanout = tile_fanout();
        if (fanout > 1 && region.voxels >= (uint64_t)MIN_VOLUME * 2) {
            Region kids[MAX_TILE_FANOUT];
            const int n = split(region, kids, fanout);
            qsort(kids, n, sizeof(Region), cmp_kmin_desc);
            for (int i = 0; i < n; ++i)
                render16_tape(sub, ctx, kids[i], img, halt, callback);
        } else {
            Region A, B;
            bisect(region, &A, &B);

            render16_tape(sub, ctx, B, img, halt, callback);
            render16_tape(sub, ctx, A, img, halt, callback);
        }
    }

    tape_release(sub);
}

void render16(MathTree* tree, Region region,
              uint16_t** img, volatile int* halt,
              void (*callback)())
{
    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    render16_tape(deck_base(deck), ctx, region, img, halt, callback);
    tape_ctx_free(ctx);
    deck_free(deck);
}

////////////////////////////////////////////////////////////////////////////////

static
void region16(const Tape* tape, TapeCtx* ctx, Region region, uint16_t** img)
{
    float *X = malloc(region.voxels*sizeof(float)),
          *Y = malloc(region.voxels*sizeof(float)),
          *Z = malloc(region.voxels*sizeof(float));

    // Copy the X, Y, Z vectors into a flattened matrix form.
    int q = 0;
    for (int k = region.nk - 1; k >= 0; --k) {
        for (unsigned j = 0; j < region.nj; ++j) {
            for (unsigned i = 0; i < region.ni; ++i) {
                X[q] = region.X[i];
                Y[q] = region.Y[j];
                Z[q] = region.Z[k];
                q++;
            }
        }
    }
    region.X = X;
    region.Y = Y;
    region.Z = Z;

    const float* result = tape_eval_r(tape, ctx, region);

    // Free the allocated matrices
    free(X);
    free(Y);
    free(Z);

    for (int k = region.nk - 1; k >= 0; --k) {
        uint16_t L = region.L[k+1];

        for (unsigned j = 0; j < region.nj; ++j) {
            int row = j + region.jmin;

            for (unsigned i = 0; i < region.ni; ++i) {
                int col = i + region.imin;

               if (*(result++) < 0 && img[row][col] < L) {
                    img[row][col] = L;
                }
            }
        }
    }
}

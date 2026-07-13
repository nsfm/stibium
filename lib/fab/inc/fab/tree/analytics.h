#ifndef ANALYTICS_H
#define ANALYTICS_H

#include <cstdint>

struct MathTree_;

/*
 *  Field statistics from dense grid integration (cell centers).
 *  For flat (2D) fields, volume holds the area and z fields are 0.
 */
struct FieldStats
{
    double volume;          /*  model units^3 (units^2 when flat)  */
    double com[3];          /*  center of mass, uniform density  */
    float tight[6];         /*  tight bounds of the filled region:
                                xmin ymin zmin xmax ymax zmax
                                (accurate to one cell)  */
    uint64_t samples;       /*  total samples taken  */
    uint64_t inside;        /*  samples with f < 0  */
    double cell;            /*  sample spacing (model units)  */
};

/*
 *  Integrates the field over the given box.  resolution is samples
 *  per model unit; <= 0 picks a default targeting ~4M samples.
 *  flat means sample the z=0 plane only (2D shapes).
 *  Returns false if halted.
 */
bool analyze_field(struct MathTree_* tree,
                   float xmin, float ymin, float zmin,
                   float xmax, float ymax, float zmax,
                   float resolution, bool flat,
                   int threads, volatile int* halt,
                   FieldStats* out);

#endif

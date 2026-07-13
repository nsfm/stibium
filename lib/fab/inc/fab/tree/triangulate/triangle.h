#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Eigen/Dense"

#include <cstdint>

typedef Eigen::Vector3d Vec3f;

/*
 *  Indexed triangle: three indices into the Mesher's vertex table.
 *  (An index of UINT32_MAX in 'a' marks an erased triangle.)
 */
struct Tri {
    uint32_t a, b, c;
};

#endif // TRIANGLE_H

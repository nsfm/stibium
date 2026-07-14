#ifndef TREE_H
#define TREE_H

#include <stdint.h>

#include "fab/tree/node/opcodes.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @struct MathTree_
    @brief A structure containing nodes organized by rank and opcode.
*/
typedef struct MathTree_ {
    /** @var nodes
    Node pointers, stored in rows indexed by level */
    struct Node_*** nodes;

    /** @var active
    Number of nodes, indexed by level */
    unsigned* active;

    /** @var num_levels
    Number of levels in this tree */
    unsigned num_levels;

    /** @var constants
    Array of constant nodes */
    struct Node_ **constants;

    /** @var num_constants
    Size of contants array */
    unsigned       num_constants;

    /** @var head
    Root of this tree */
    struct Node_* head;
} MathTree;


/** @brief Creates a new tree, allocating the appropriate arrays.
    @param num_levels Number of levels in tree
    @param num_constants Number of constants in tree
 */
MathTree* new_tree(unsigned num_levels, unsigned num_constants);


/** @brief Frees a tree and all of its nodes
    @param tree Target tree
*/
void free_tree(MathTree* tree);


/** @brief Prints a math tree to stdout. */
void print_tree(MathTree* tree);


/** @brief Prints a math tree to a given file descriptor. */
void fprint_tree(MathTree* tree, int fd);

/*  Spatial specialization (formerly disable_nodes / enable_nodes /
 *  clone_tree) now lives in the immutable shortened-tape machinery:
 *  see fab/tree/tape.h and doc/TAPE-DESIGN.md.  */

/** @brief Returns a bit mask containing active axes in a tree.
    @details
    The bit mask is of the form (x_active << 2) | (y_active << 1) | (z_active)
*/
uint8_t active_axes(const MathTree* const tree);

#ifdef __cplusplus
}
#endif

#endif

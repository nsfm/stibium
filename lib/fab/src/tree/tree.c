#include <stdlib.h>

#include "fab/tree/tree.h"

#include "fab/tree/node/node.h"
#include "fab/tree/node/printers.h"
#include "fab/util/switches.h"

MathTree* new_tree(unsigned num_levels, unsigned num_constants)
{
    MathTree* tree = malloc(sizeof(MathTree));

    *tree = (MathTree){
        .nodes      = num_levels ?
                        calloc(num_levels, sizeof(Node**)) : NULL,
        .active     = num_levels ?
                        calloc(num_levels, sizeof(unsigned)) : NULL,
        .constants  = num_constants ?
                        malloc(sizeof(Node*)*num_constants) : NULL,

        .num_constants = num_constants,
        .head = NULL,
        .num_levels = num_levels,
    };

    return tree;
}


void free_tree(MathTree* tree)
{
    if (tree == NULL)   return;

    for (unsigned level=0; level < tree->num_levels; ++level) {
        for (unsigned n=0; n < tree->active[level]; ++n) {
            free_node(tree->nodes[level][n]);
        }
        free(tree->nodes[level]);
    }

    for (unsigned c=0; c < tree->num_constants; ++c) {
        free_node(tree->constants[c]);
    }

    free(tree->nodes);
    free(tree->active);
    free(tree->constants);

    free(tree);
}


void print_tree(MathTree* tree)
{
    print_node(tree->head);
}


void fdprint_tree(MathTree* tree, int fd)
{
    fdprint_node(tree->head, fd);
}

unsigned count_nodes(MathTree* tree)
{
    unsigned count = tree->num_constants;
    for (unsigned level=0; level < tree->num_levels; ++level) {
        count += tree->active[level];
    }
    return count;
}

////////////////////////////////////////////////////////////////////////////////

uint8_t active_axes(const MathTree* const tree)
{
    if (!tree->num_levels)  return 0;

    uint8_t active = 0;
    if (tree->num_levels) {
        for (unsigned a=0; a < tree->active[0]; ++a) {
            switch (tree->nodes[0][a]->opcode) {
                case OP_X:  active |= (1 << 2); break;
                case OP_Y:  active |= (1 << 1); break;
                case OP_Z:  active |= (1 << 0); break;
                default: ;
            }
        }
    }

    return active;
}

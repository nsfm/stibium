#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fab/tree/node/node.h"
#include "fab/tree/node/printers.h"
#include "fab/tree/math/math_f.h"
#include "fab/tree/grid.h"

#include "fab/util/switches.h"

// Non-recursively clone a node.
Node* clone_node(Node* n)
{
    // Allocate memory and copy everything over
    Node* clone = malloc(sizeof(Node));
    memcpy(clone, n, sizeof(Node));

    // Update children clone pointers
    if (n->lhs) clone->lhs = n->lhs->clone_address;
    if (n->rhs) clone->rhs = n->rhs->clone_address;
    if (n->mhs) clone->mhs = n->mhs->clone_address;

    // The clone shares the original's payload (grids are immutable
    // and refcounted, so clones across threads read the same data)
    if (n->opcode == OP_GRID)
        grid_retain((MeshGrid*)clone->payload);

    // Record the address of the new clone, so that clones of
    // its parents can be adjusted to point in the right place
    n->clone_address = clone;

    return clone;
}

void free_node(Node* n)
{
    if (n == NULL)  return;
    if (n->opcode == OP_GRID)
        grid_release((MeshGrid*)n->payload);
    free(n);
}

////////////////////////////////////////////////////////////////////////////////

Node* binary_n(Node* lhs, Node* rhs, float (*f)(float, float), Opcode op)
{
    Node* n = malloc(sizeof(Node));

    _Bool constant = (lhs->flags & NODE_CONSTANT) &&
                     (rhs->flags & NODE_CONSTANT);

    *n = (Node) {
        .opcode     = constant ? OP_CONST : op,
        .rank       = constant ? 0 : 1 + (lhs->rank > rhs->rank ?
                                          lhs->rank : rhs->rank),
        .flags      = constant ? NODE_CONSTANT : 0,
        .lhs        = constant ? NULL : lhs,
        .rhs        = constant ? NULL : rhs,
        .clone_address = NULL,
    };

    if (constant) {
        fill_results(n, (*f)(lhs->results.f, rhs->results.f));
    }

    return n;
}

Node* add_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, add_f, OP_ADD); }
Node* sub_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, sub_f, OP_SUB); }
Node* mul_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, mul_f, OP_MUL); }
Node* div_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, div_f, OP_DIV); }

Node* min_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, min_f, OP_MIN); }
Node* mod_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, mod_f, OP_MOD); }
Node* max_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, max_f, OP_MAX); }
Node* pow_n(Node* lhs, Node* rhs) { return binary_n(lhs, rhs, pow_f, OP_POW); }

Node* atan2_n(Node* a, Node* b) {return binary_n(a, b, atan2_f, OP_ATAN2); }

////////////////////////////////////////////////////////////////////////////////

Node* unary_n(Node* arg, float (*f)(float), Opcode op)
{
    Node* n = malloc(sizeof(Node));

    _Bool constant = arg->flags & NODE_CONSTANT;

    *n = (Node) {
        .opcode     = constant ? OP_CONST : op,
        .rank       = constant ? 0 : 1 + arg->rank,
        .flags      = constant ? NODE_CONSTANT : 0,
        .lhs        = constant ? NULL : arg,
        .rhs        = NULL,
        .clone_address = NULL,
    };

    if (constant) {
        fill_results(n, (*f)(arg->results.f));
    }
    return n;
}

Node* abs_n(Node* child) { return unary_n(child, abs_f, OP_ABS); }
Node* square_n(Node* child) { return unary_n(child, square_f, OP_SQUARE); }
Node* sqrt_n(Node* child) { return unary_n(child, sqrt_f, OP_SQRT); }
Node* floor_n(Node* child) { return unary_n(child, floor_f, OP_FLOOR); }
Node* log_n(Node* child) { return unary_n(child, log_f, OP_LOG); }
Node* sin_n(Node* child) { return unary_n(child, sin_f, OP_SIN); }
Node* cos_n(Node* child) { return unary_n(child, cos_f, OP_COS); }
Node* tan_n(Node* child) { return unary_n(child, tan_f, OP_TAN); }
Node* asin_n(Node* child) { return unary_n(child, asin_f, OP_ASIN); }
Node* acos_n(Node* child) { return unary_n(child, acos_f, OP_ACOS); }
Node* atan_n(Node* child) { return unary_n(child, atan_f, OP_ATAN); }
Node* neg_n(Node* child) { return unary_n(child, neg_f, OP_NEG); }
Node* exp_n(Node* child) { return unary_n(child, exp_f, OP_EXP); }

////////////////////////////////////////////////////////////////////////////////

Node* nonary_n(Opcode op)
{
    Node* n = malloc(sizeof(Node));

    *n = (Node) {
        .opcode     = op,
        .rank       = 0,
        .flags      = 0,
        .lhs        = NULL,
        .rhs        = NULL,
        .clone_address = NULL,
    };

    return n;
}

Node* constant_n(float value)
{
    Node* n = nonary_n(OP_CONST);
    n->flags = NODE_CONSTANT;
    fill_results(n, value);
    return n;
}

Node* grid_n(struct MeshGrid_* grid, Node* x, Node* y, Node* z)
{
    Node* n = malloc(sizeof(Node));

    int rank = x->rank > y->rank ? x->rank : y->rank;
    if (z->rank > rank)   rank = z->rank;

    /*  Never constant-folded, even with constant children: the
     *  payload lookup is cheap and folding would complicate the
     *  refcount story for no real win. */
    *n = (Node) {
        .opcode     = OP_GRID,
        .rank       = 1 + rank,
        .flags      = 0,
        .lhs        = x,
        .rhs        = y,
        .mhs        = z,
        .payload    = grid,
        .clone_address = NULL,
    };
    grid_retain(grid);

    return n;
}

Node* X_n()
{
    return nonary_n(OP_X);
}

Node* Y_n()
{
    return nonary_n(OP_Y);
}

Node* Z_n()
{
    return nonary_n(OP_Z);
}

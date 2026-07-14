#ifndef NODE_H
#define NODE_H

#include <stdbool.h>

#include "fab/tree/node/results.h"
#include "fab/tree/node/opcodes.h"
#include "fab/util/interval.h"
#include "fab/util/region.h"

#define NODE_CONSTANT   1
#define NODE_IN_TREE    8

#ifdef __cplusplus
extern "C" {
#endif

/** @struct Node_
    @brief Recursive data structure defining a node in a math tree.
*/
typedef struct Node_ {
    /** @var opcode
    Node operation */
    Opcode opcode;

    /** @var results
    Saved results from most recent evaluation */
    Results  results;

    /** @var rank
    Rank of the node in the tree. */
    int rank;

    /** @var flags
    Flags (combination of NODE_CONSTANT and NODE_IN_TREE) */
    uint8_t flags;

    /** @var lhs
    Left-hand child node (or NULL)
    */
    struct Node_* lhs;

    /** @var rhs
    Right-hand child node (or NULL)
    */
    struct Node_* rhs;

    /** @var mhs
    Third child node, used only by OP_GRID (the z sample coordinate;
    lhs and rhs carry x and y).  NULL everywhere else. */
    struct Node_* mhs;

    /** @var payload
    Refcounted heap payload, used only by OP_GRID (a MeshGrid*).
    deck_from_tree retains it; free_node releases it. */
    void* payload;
} Node;


/** @brief Frees a single node, releasing its payload if it has one.
    @details All node deallocation must go through here (not free)
    so OP_GRID payload refcounts stay balanced. */
void free_node(Node* n);

////////////////////////////////////////////////////////////////////////////////
// Node constructors
////////////////////////////////////////////////////////////////////////////////

// Binary operations
Node* add_n(Node* left, Node* right);
Node* sub_n(Node* left, Node* right);
Node* mul_n(Node* left, Node* right);
Node* div_n(Node* left, Node* right);

Node* min_n(Node* left, Node* right);
Node* max_n(Node* left, Node* right);
Node* pow_n(Node* left, Node* right);

Node* atan2_n(Node* a, Node* b);
Node* mod_n(Node* a, Node* b);
Node* floor_n(Node* child);
Node* log_n(Node* child);

// Unary arithmetic operators
Node* abs_n(Node* n);
Node* square_n(Node* n);
Node* sqrt_n(Node* n);
Node* sin_n(Node* n);
Node* cos_n(Node* n);
Node* tan_n(Node* n);
Node* asin_n(Node* n);
Node* acos_n(Node* n);
Node* atan_n(Node* n);
Node* neg_n(Node* n);
Node* exp_n(Node* n);

// Constants
Node* constant_n(float value);

// Variables
Node* X_n(void);
Node* Y_n(void);
Node* Z_n(void);

// Grid sample (ternary; retains the grid payload)
struct MeshGrid_;
Node* grid_n(struct MeshGrid_* grid, Node* x, Node* y, Node* z);

#ifdef __cplusplus
}
#endif

#endif

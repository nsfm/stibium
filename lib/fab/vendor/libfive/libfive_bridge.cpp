/*  Stibium <-> libfive bridge implementation (Stibium's own work,
 *  AGPL-3.0-or-later).  Compiled against libfive's Eigen in an isolated
 *  shared library; see libfive_bridge.h for the boundary contract.
 */
#include "libfive_bridge.h"

#include <unordered_map>

// Stibium kernel headers (plain C, Eigen-free).
extern "C" {
#include "fab/tree/node/node.h"
#include "fab/tree/node/opcodes.h"
}

// libfive.
#include "libfive/tree/tree.hpp"
#include "libfive/tree/opcode.hpp"
#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/region.hpp"
#include "libfive/render/brep/settings.hpp"

namespace stibium_libfive {

namespace {

// Recursive Stibium-Node -> libfive::Tree translation, memoized on the Node
// pointer (the Stibium tree is a DAG - subexpressions are shared - so a naive
// recursion would blow up exponentially).
struct Translator {
    std::unordered_map<const Node*, libfive::Tree> memo;
    bool ok = true;
    std::string err;

    libfive::Tree go(const Node* n) {
        if (!ok) {
            return libfive::Tree::invalid();
        }
        auto it = memo.find(n);
        if (it != memo.end()) {
            return it->second;
        }
        libfive::Tree t = build(n);
        memo.emplace(n, t);
        return t;
    }

    void fail(const std::string& m) {
        if (ok) {
            ok = false;
            err = m;
        }
    }

    libfive::Tree build(const Node* n) {
        using namespace libfive;
        switch (n->opcode) {
            // Nullary
            case OP_X: return Tree::X();
            case OP_Y: return Tree::Y();
            case OP_Z: return Tree::Z();
            case OP_CONST: return Tree(double(n->results.f));

            // Binary (lhs, rhs)
            case OP_ADD:   return Tree::binary(Opcode::OP_ADD,   go(n->lhs), go(n->rhs));
            case OP_SUB:   return Tree::binary(Opcode::OP_SUB,   go(n->lhs), go(n->rhs));
            case OP_MUL:   return Tree::binary(Opcode::OP_MUL,   go(n->lhs), go(n->rhs));
            case OP_DIV:   return Tree::binary(Opcode::OP_DIV,   go(n->lhs), go(n->rhs));
            case OP_MIN:   return Tree::binary(Opcode::OP_MIN,   go(n->lhs), go(n->rhs));
            case OP_MAX:   return Tree::binary(Opcode::OP_MAX,   go(n->lhs), go(n->rhs));
            case OP_POW:   return Tree::binary(Opcode::OP_POW,   go(n->lhs), go(n->rhs));
            case OP_MOD:   return Tree::binary(Opcode::OP_MOD,   go(n->lhs), go(n->rhs));
            case OP_ATAN2: return Tree::binary(Opcode::OP_ATAN2, go(n->lhs), go(n->rhs));

            // Unary (child in lhs)
            case OP_ABS:    return Tree::unary(Opcode::OP_ABS,    go(n->lhs));
            case OP_SQUARE: return Tree::unary(Opcode::OP_SQUARE, go(n->lhs));
            case OP_SQRT:   return Tree::unary(Opcode::OP_SQRT,   go(n->lhs));
            case OP_SIN:    return Tree::unary(Opcode::OP_SIN,    go(n->lhs));
            case OP_COS:    return Tree::unary(Opcode::OP_COS,    go(n->lhs));
            case OP_TAN:    return Tree::unary(Opcode::OP_TAN,    go(n->lhs));
            case OP_ASIN:   return Tree::unary(Opcode::OP_ASIN,   go(n->lhs));
            case OP_ACOS:   return Tree::unary(Opcode::OP_ACOS,   go(n->lhs));
            case OP_ATAN:   return Tree::unary(Opcode::OP_ATAN,   go(n->lhs));
            case OP_NEG:    return Tree::unary(Opcode::OP_NEG,    go(n->lhs));
            case OP_EXP:    return Tree::unary(Opcode::OP_EXP,    go(n->lhs));
            case OP_LOG:    return Tree::unary(Opcode::OP_LOG,    go(n->lhs));

            // No libfive equivalent.
            case OP_FLOOR:
                fail("OP_FLOOR has no libfive opcode (tilings/patterns "
                     "unsupported by the libfive port)");
                return Tree::invalid();
            case OP_GRID:
                fail("OP_GRID (imported-mesh / heightmap trilinear sample) "
                     "has no libfive opcode; would need an OracleClause");
                return Tree::invalid();
            default:
                fail("unsupported Stibium opcode " +
                     std::to_string(int(n->opcode)));
                return Tree::invalid();
        }
    }
};

}   // namespace

STIBIUM_LF_EXPORT
int mesh_shape(const Node_* root,
               double x0, double y0, double z0,
               double x1, double y1, double z1,
               double min_feature, double max_err,
               unsigned workers,
               std::vector<float>& verts,
               std::vector<uint32_t>& tris,
               std::string& err) {
    if (!root) {
        err = "null root node";
        return 1;
    }

    Translator tr;
    libfive::Tree shape = tr.go(reinterpret_cast<const Node*>(root));
    if (!tr.ok) {
        err = tr.err;
        return 2;
    }

    libfive::Region<3> region(
        {x0, y0, z0}, {x1, y1, z1});

    libfive::BRepSettings settings;
    settings.min_feature = min_feature;
    settings.max_err = max_err;
    settings.workers = workers ? workers : 1;
    settings.alg = libfive::DUAL_CONTOURING;

    std::unique_ptr<libfive::Mesh> mesh;
    try {
        mesh = libfive::Mesh::render(shape, region, settings);
    } catch (const std::exception& e) {
        err = std::string("libfive render threw: ") + e.what();
        return 3;
    }
    if (!mesh) {
        err = "libfive render returned null (cancelled or empty)";
        return 4;
    }

    verts.clear();
    tris.clear();
    verts.reserve(mesh->verts.size() * 3);
    tris.reserve(mesh->branes.size() * 3);
    // verts[0] is libfive's reserved marker vertex; keep it so brane indices
    // stay valid (it is never referenced by a triangle).
    for (const auto& v : mesh->verts) {
        verts.push_back(v.x());
        verts.push_back(v.y());
        verts.push_back(v.z());
    }
    for (const auto& b : mesh->branes) {
        tris.push_back(b[0]);
        tris.push_back(b[1]);
        tris.push_back(b[2]);
    }
    return 0;
}

}   // namespace stibium_libfive

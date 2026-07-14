#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QRegularExpression>
#include <QTextStream>

#include <cstdio>
#include <cstring>
#include <functional>

#include "app/vm_import.h"

#include "fab/tree/tree.h"
#include "fab/tree/parser.h"
#include "fab/tree/eval.h"
#include "fab/util/interval.h"

namespace {

/*
 *  A .vm file is a flat SSA tape: one op per line, arguments are the
 *  NAMES of previously-defined lines, the last line is the output.
 *  Translation emits each node's V1 prefix substring (memoized, so a
 *  shared node is built once but inlines textually at each use) and
 *  tracks expression depth, which is what bounds the parser stack.
 */
struct VmNode {
    QByteArray math;
    int depth;
};

/*  The parser stack is 4096 entries; stay well clear. */
const int MAX_DEPTH = 3500;

const char* SCRIPT_TEMPLATE_NOTE =
    "Imported from a Fidget .vm tape (frozen math; not parametric)";

/*  ops with a direct single-character V1 spelling */
QByteArray v1_unary(const QString& op)
{
    if (op == "abs")    return "b";
    if (op == "neg")    return "n";
    if (op == "sqrt")   return "r";
    if (op == "square") return "q";
    if (op == "floor")  return "F";
    if (op == "sin")    return "s";
    if (op == "cos")    return "c";
    if (op == "tan")    return "t";
    if (op == "asin")   return "S";
    if (op == "acos")   return "C";
    if (op == "atan")   return "T";
    if (op == "ln")     return "l";
    if (op == "exp")    return "x";
    return "";
}

QByteArray v1_binary(const QString& op)
{
    if (op == "add")    return "+";
    if (op == "sub")    return "-";
    if (op == "mul")    return "*";
    if (op == "div")    return "/";
    if (op == "min")    return "i";
    if (op == "max")    return "a";
    if (op == "mod")    return "M";   // both sides floor-mod (euclidean)
    return "";
}

}  // namespace

int importVmHeadless(const QString& vm_path, const QString& sb_path)
{
    QFile in(vm_path);
    if (!in.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        fprintf(stderr, "import-vm: cannot read '%s'\n",
                vm_path.toLocal8Bit().constData());
        return 1;
    }

    QHash<QString, VmNode> nodes;
    QString last;
    bool uses_z = false;
    int line_no = 0;

    QTextStream ts(&in);
    while (!ts.atEnd())
    {
        const QString line = ts.readLine();
        ++line_no;
        const QString trimmed = line.trimmed();
        if (trimmed.isEmpty() || trimmed.startsWith('#'))
            continue;

        const QStringList tok = trimmed.split(
                QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        if (tok.size() < 2)
        {
            fprintf(stderr, "import-vm: line %d: expected "
                    "'<name> <op> [args...]'\n", line_no);
            return 1;
        }
        const QString& name = tok[0];
        const QString& op = tok[1];

        // Resolves an argument name to its node, or fails loudly
        auto arg = [&](int idx) -> const VmNode* {
            if (idx >= tok.size() || !nodes.contains(tok[idx]))
            {
                fprintf(stderr, "import-vm: line %d: unknown or "
                        "missing operand for '%s'\n", line_no,
                        op.toLocal8Bit().constData());
                return nullptr;
            }
            return &nodes[tok[idx]];
        };

        VmNode out;
        if (op == "var-x")      out = {"X", 1};
        else if (op == "var-y") out = {"Y", 1};
        else if (op == "var-z") { out = {"Z", 1}; uses_z = true; }
        else if (op == "const")
        {
            if (tok.size() < 3)
            {
                fprintf(stderr, "import-vm: line %d: const needs a "
                        "value\n", line_no);
                return 1;
            }
            QString v = tok[2];
            if (v.startsWith('+'))
                v.remove(0, 1);
            out = {"f" + v.toUtf8(), 1};
        }
        else if (!v1_unary(op).isEmpty())
        {
            const VmNode* a = arg(2);
            if (!a)     return 1;
            out = {v1_unary(op) + a->math, a->depth + 1};
        }
        else if (op == "ceil")
        {
            // ceil(a) = -floor(-a)
            const VmNode* a = arg(2);
            if (!a)     return 1;
            out = {"nFn" + a->math, a->depth + 3};
        }
        else if (op == "round")
        {
            // round(a) = floor(a + 0.5)
            const VmNode* a = arg(2);
            if (!a)     return 1;
            out = {"F+" + a->math + "f0.5", a->depth + 2};
        }
        else if (op == "recip")
        {
            const VmNode* a = arg(2);
            if (!a)     return 1;
            out = {"/f1" + a->math, a->depth + 1};
        }
        else if (!v1_binary(op).isEmpty())
        {
            const VmNode* a = arg(2);
            const VmNode* b = arg(3);
            if (!a || !b)   return 1;
            out = {v1_binary(op) + a->math + b->math,
                   1 + qMax(a->depth, b->depth)};
        }
        else if (op == "atan2")
        {
            // No V1 spelling; embed the V2 form (grammar: '=expr;'
            // is a v1_expr, and '[...]' re-enters V1 inside V2)
            const VmNode* a = arg(2);
            const VmNode* b = arg(3);
            if (!a || !b)   return 1;
            out = {"=atan2([" + a->math + "],[" + b->math + "]);",
                   2 + qMax(a->depth, b->depth)};
        }
        else if (op == "compare" || op == "and" || op == "or" ||
                 op == "not")
        {
            fprintf(stderr, "import-vm: line %d: opcode '%s' has no "
                    "Stibium equivalent (discontinuous/boolean); "
                    "cannot translate this model\n", line_no,
                    op.toLocal8Bit().constData());
            return 1;
        }
        else
        {
            fprintf(stderr, "import-vm: line %d: unknown opcode "
                    "'%s'\n", line_no, op.toLocal8Bit().constData());
            return 1;
        }

        if (out.depth > MAX_DEPTH)
        {
            fprintf(stderr, "import-vm: expression depth %d exceeds "
                    "the parser budget (%d); this tape is too deeply "
                    "nested to freeze into a math string\n",
                    out.depth, MAX_DEPTH);
            return 1;
        }
        nodes[name] = out;
        last = name;
    }

    if (last.isEmpty())
    {
        fprintf(stderr, "import-vm: no ops found in '%s'\n",
                vm_path.toLocal8Bit().constData());
        return 1;
    }

    const QByteArray math = nodes[last].math;

    // Prove the string parses before writing anything
    MathTree* tree = parse(math.constData());
    if (!tree)
    {
        fprintf(stderr, "import-vm: internal error: generated math "
                "string failed to parse (please report)\n");
        return 1;
    }

    /*
     *  Auto-bounds: recursive interval descent from a big seed box.
     *  Boxes proven empty (lower > 0) are dropped; boxes proven or
     *  possibly solid grow the bounding box, ambiguous ones refine.
     *  2D tapes (no var-z) use a flat quadtree at z = 0.
     */
    float bb[6] = {1e30f, 1e30f, 1e30f, -1e30f, -1e30f, -1e30f};
    bool found = false;
    const std::function<void(float, float, float, float, float,
                             float, int)> descend =
        [&](float x0, float y0, float z0,
            float x1, float y1, float z1, int depth) {
        Interval X = {x0, x1}, Y = {y0, y1}, Z = {z0, z1};
        const Interval out = eval_i(tree, X, Y, Z);
        if (out.lower > 0)
            return;                     // proven empty
        if (out.upper < 0 || depth == 0)
        {
            // solid, or ambiguous at the resolution floor: keep
            bb[0] = qMin(bb[0], x0);  bb[3] = qMax(bb[3], x1);
            bb[1] = qMin(bb[1], y0);  bb[4] = qMax(bb[4], y1);
            bb[2] = qMin(bb[2], z0);  bb[5] = qMax(bb[5], z1);
            found = true;
            return;
        }
        const float xm = (x0 + x1)/2, ym = (y0 + y1)/2,
                    zm = (z0 + z1)/2;
        for (int i = 0; i < (uses_z ? 8 : 4); ++i)
        {
            descend(i & 1 ? xm : x0, i & 2 ? ym : y0,
                    uses_z ? (i & 4 ? zm : z0) : z0,
                    i & 1 ? x1 : xm, i & 2 ? y1 : ym,
                    uses_z ? (i & 4 ? z1 : zm) : z1,
                    depth - 1);
        }
    };
    /*  Two passes: a coarse sweep of the big seed box finds the rough
     *  extent (floor: 1 unit), then a second descent seeded on that
     *  box tightens to ~0.2% of the model's size — scale-free. */
    const float SEED = 512;
    descend(-SEED, -SEED, uses_z ? -SEED : 0,
            SEED, SEED, uses_z ? SEED : 0, 10);
    if (found)
    {
        float rough[6];
        memcpy(rough, bb, sizeof(bb));
        bb[0] = bb[1] = bb[2] = 1e30f;
        bb[3] = bb[4] = bb[5] = -1e30f;
        found = false;
        float c[3], h = 0;
        for (int a = 0; a < 3; ++a)
        {
            c[a] = (rough[a] + rough[a + 3]) / 2;
            h = qMax(h, (rough[a + 3] - rough[a]) / 2 + 1);
        }
        descend(c[0] - h, c[1] - h, uses_z ? c[2] - h : 0,
                c[0] + h, c[1] + h, uses_z ? c[2] + h : 0, 10);
    }
    free_tree(tree);

    if (!found)
    {
        fprintf(stderr, "import-vm: warning: could not locate the "
                "surface within [-%g, %g]^3; defaulting bounds to "
                "[-1, 1]\n", SEED, SEED);
        bb[0] = bb[1] = bb[2] = -1;
        bb[3] = bb[4] = bb[5] = 1;
    }
    else
    {
        // Breathing room so the surface never touches the box
        for (int a = 0; a < 3; ++a)
        {
            const float pad = 0.05f * (bb[a + 3] - bb[a]) + 1e-3f;
            bb[a] -= pad;
            bb[a + 3] += pad;
        }
        if (bb[0] <= -SEED || bb[3] >= SEED)
            fprintf(stderr, "import-vm: warning: shape reaches the "
                    "bounds probe's seed box (unbounded lattice?); "
                    "bounds are clamped to [-%g, %g]\n", SEED, SEED);
    }

    // Assemble the .sb: protocol 7, one script node
    const QString stem = QFileInfo(vm_path).completeBaseName()
            .toLower().replace(QRegularExpression("[^a-z0-9_]"), "_");

    QString bounds_args;
    if (uses_z)
        bounds_args = QString("%1, %2, %3, %4, %5, %6")
                .arg(bb[0]).arg(bb[1]).arg(bb[2])
                .arg(bb[3]).arg(bb[4]).arg(bb[5]);
    else
        bounds_args = QString("%1, %2, %3, %4")
                .arg(bb[0]).arg(bb[1]).arg(bb[3]).arg(bb[4]);

    QJsonArray script;
    script.append("import fab.types");
    script.append("");
    script.append(QString("title('Fidget import (%1)')").arg(stem));
    script.append("");
    script.append(QString("# %1").arg(SCRIPT_TEMPLATE_NOTE));
    script.append(QString("output('shape', fab.types.Shape('%1', %2))")
                  .arg(QString::fromUtf8(math)).arg(bounds_args));

    QJsonObject datum;
    datum["uid"] = 0;
    datum["name"] = "shape";
    datum["expr"] = QString(QChar(0x12));   // SIGIL_OUTPUT stub
    datum["type"] = "_fabtypes.Shape";

    QJsonObject node;
    node["uid"] = 0;
    node["name"] = stem.isEmpty() ? QString("import0") : stem;
    node["inspector"] = QJsonArray({0.0, 0.0});
    node["script"] = script;
    node["datums"] = QJsonArray({datum});

    QJsonObject doc;
    doc["type"] = "sb";
    doc["protocol"] = 7;
    doc["nodes"] = QJsonArray({node});

    QFile out_file(sb_path);
    if (!out_file.open(QIODevice::WriteOnly))
    {
        fprintf(stderr, "import-vm: cannot write '%s'\n",
                sb_path.toLocal8Bit().constData());
        return 1;
    }
    out_file.write(QJsonDocument(doc).toJson());

    fprintf(stderr, "import-vm: %s: %lld ops -> %lld byte math "
            "string, depth %d, bounds [%g %g %g] - [%g %g %g]%s\n",
            stem.toLocal8Bit().constData(),
            (long long)nodes.size(), (long long)math.size(),
            nodes[last].depth,
            bb[0], bb[1], uses_z ? bb[2] : 0.f,
            bb[3], bb[4], uses_z ? bb[5] : 0.f,
            uses_z ? "" : " (2D)");
    return 0;
}

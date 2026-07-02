/* kflc — AST dump helper. Used by the kflc driver (--dump) and by
 * embedding tools that want to inspect a parsed form. */

#include "kflc.h"

#include <stdio.h>

static const char *kind_name(KflcNodeKind k)
{
    switch (k) {
    case KFLN_FORM:         return "FORM";
    case KFLN_WIDGET:       return "WIDGET";
    case KFLN_ARG:          return "ARG";
    case KFLN_FN:           return "FN";
    case KFLN_FN_ARG:       return "FN_ARG";
    case KFLN_STMT_LET:     return "STMT_LET";
    case KFLN_STMT_CONST:   return "STMT_CONST";
    case KFLN_STMT_ASSIGN:  return "STMT_ASSIGN";
    case KFLN_STMT_RETURN:  return "STMT_RETURN";
    case KFLN_STMT_EXPR:    return "STMT_EXPR";
    case KFLN_STMT_PRINT:   return "STMT_PRINT";
    case KFLN_STMT_IF:      return "STMT_IF";
    case KFLN_STMT_WHILE:   return "STMT_WHILE";
    case KFLN_FN_DATA:      return "FN_DATA";
    case KFLN_STMT_INDEX_ASSIGN: return "STMT_INDEX_ASSIGN";
    case KFLN_STMT_SERIES:  return "STMT_SERIES";
    case KFLN_FN_WORLD:       return "FN_WORLD";
    case KFLN_TICK:           return "TICK";
    case KFLN_FRAME:          return "FRAME";
    case KFLN_EPOCH_LITERAL:  return "EPOCH_LITERAL";
    case KFLN_STMT_ASTRO_BODY:return "STMT_ASTRO_BODY";
    case KFLN_STMT_STEP:      return "STMT_STEP";
    case KFLN_STMT_PROPAGATE: return "STMT_PROPAGATE";
    case KFLN_STMT_FOR_EACH:  return "STMT_FOR_EACH";
    case KFLN_STMT_OBSERVE:   return "STMT_OBSERVE";
    case KFLN_STMT_LVALUE_ASSIGN: return "STMT_LVALUE_ASSIGN";
    case KFLN_ARENA:          return "ARENA";
    case KFLN_ALLOCATOR_BIND: return "ALLOCATOR_BIND";
    }
    return "?";
}

static const char *widget_name(KflcWidgetKind w)
{
    switch (w) {
    case KFL_W_PLOT:         return "plot";
    }
    return "?";
}

static void print_indent(FILE *out, int n)
{
    for (int i = 0; i < n; i++) fputs("  ", out);
}

static void print_value(FILE *out, const KflcValue *v)
{
    switch (v->kind) {
    case KFLV_NONE:     fprintf(out, "(none)"); break;
    case KFLV_INT:      fprintf(out, "%ld", v->u.i); break;
    case KFLV_FLOAT:    fprintf(out, "%g", v->u.f); break;
    case KFLV_STR:      fprintf(out, "\"%s\"", v->u.s ? v->u.s : ""); break;
    case KFLV_IDENT:    fprintf(out, "%s", v->u.s ? v->u.s : ""); break;
    case KFLV_COLOR:    fprintf(out, "#%06x", v->u.rgb); break;
    case KFLV_SIZE:     fprintf(out, "%dx%d", v->u.size.w, v->u.size.h); break;
    case KFLV_SHORTCUT: fprintf(out, "shortcut[%s]", v->u.s ? v->u.s : ""); break;
    case KFLV_ALIGN:    fprintf(out, "align[0x%x]", v->u.align); break;
    }
}

void kflc_dump_node(FILE *out, const KflcNode *n, int indent)
{
    if (!n) return;
    print_indent(out, indent);
    fprintf(out, "%s", kind_name(n->kind));
    if (n->kind == KFLN_WIDGET)    fprintf(out, "/%s", widget_name(n->widget));
    if (n->name) fprintf(out, " name=%s", n->name);
    if (n->position.kind != KFLV_NONE) {
        fprintf(out, " pos=");
        print_value(out, &n->position);
    }
    fprintf(out, " @%d\n", n->line);

    for (KflcAttr *a = n->attrs; a; a = a->next) {
        print_indent(out, indent + 1);
        fprintf(out, ". %s = ", a->name);
        print_value(out, &a->value);
        fprintf(out, " @%d\n", a->line);
    }
    for (KflcNode *c = n->children; c; c = c->next) {
        kflc_dump_node(out, c, indent + 1);
    }
}

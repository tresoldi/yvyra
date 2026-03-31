/*
 *  yvyra - Bayesian Phylogenetics for Linguistic Data
 *  Based on MrBayes 3.2.7a by Ronquist et al.
 *
 *  yamlparser.c: Minimal YAML subset parser.
 *  Handles scalars, block/flow mappings/sequences, comments,
 *  quoted strings. No anchors, tags, or multi-document.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "yamlparser.h"

/* ---- Node allocation ---- */

static YamlNode *NewNode (YamlNodeType type)
{
    YamlNode *n = (YamlNode *) calloc(1, sizeof(YamlNode));
    if (!n) return NULL;
    n->type = type;
    n->maxChildren = 16;
    n->children = (YamlNode **) calloc(n->maxChildren, sizeof(YamlNode *));
    return n;
}

static int AddChild (YamlNode *parent, YamlNode *child)
{
    if (parent->nChildren >= parent->maxChildren)
        {
        parent->maxChildren *= 2;
        parent->children = (YamlNode **) realloc(parent->children,
            parent->maxChildren * sizeof(YamlNode *));
        if (!parent->children) return -1;
        }
    parent->children[parent->nChildren++] = child;
    return 0;
}

static char *DupStr (const char *s)
{
    if (!s) return NULL;
    char *d = (char *) malloc(strlen(s) + 1);
    if (d) strcpy(d, s);
    return d;
}

void YamlFreeNode (YamlNode *node)
{
    int i;
    if (!node) return;
    if (node->key) free(node->key);
    if (node->scalar) free(node->scalar);
    if (node->children)
        {
        for (i = 0; i < node->nChildren; i++)
            YamlFreeNode(node->children[i]);
        free(node->children);
        }
    free(node);
}

/* ---- Parsing state ---- */

typedef struct {
    const char  *text;      /* full input text */
    const char  *pos;       /* current position */
    int         line;       /* current line number (1-based) */
    char        errmsg[256];
} ParseState;

static void SkipSpaces (ParseState *ps)
{
    while (*ps->pos == ' ' || *ps->pos == '\t')
        ps->pos++;
}

static void SkipToEol (ParseState *ps)
{
    while (*ps->pos && *ps->pos != '\n')
        ps->pos++;
}

static void SkipComment (ParseState *ps)
{
    if (*ps->pos == '#')
        SkipToEol(ps);
}

static void SkipBlankLines (ParseState *ps)
{
    while (*ps->pos)
        {
        const char *start = ps->pos;
        SkipSpaces(ps);
        SkipComment(ps);
        if (*ps->pos == '\n')
            { ps->pos++; ps->line++; }
        else if (*ps->pos == '\r')
            { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
        else
            { ps->pos = start; return; }
        }
}

static int GetIndent (ParseState *ps)
{
    const char *p = ps->pos;
    int indent = 0;
    while (*p == ' ') { indent++; p++; }
    return indent;
}

/* ---- Scalar parsing ---- */

static char *ParseQuotedScalar (ParseState *ps, char quote)
{
    char buf[YAML_MAX_SCALAR_LEN];
    int len = 0;

    ps->pos++;  /* skip opening quote */
    while (*ps->pos && *ps->pos != quote && len < YAML_MAX_SCALAR_LEN - 1)
        {
        if (*ps->pos == '\\' && ps->pos[1])
            {
            ps->pos++;
            if (*ps->pos == 'n') buf[len++] = '\n';
            else if (*ps->pos == 't') buf[len++] = '\t';
            else buf[len++] = *ps->pos;
            }
        else
            buf[len++] = *ps->pos;
        ps->pos++;
        }
    if (*ps->pos == quote) ps->pos++;
    buf[len] = '\0';
    return DupStr(buf);
}

static char *ParseBareScalar (ParseState *ps, const char *terminators)
{
    char buf[YAML_MAX_SCALAR_LEN];
    int len = 0;
    const char *start = ps->pos;

    while (*ps->pos && *ps->pos != '\n' && *ps->pos != '\r' &&
           !strchr(terminators, *ps->pos) && len < YAML_MAX_SCALAR_LEN - 1)
        {
        /* Stop at inline comment (space + #) */
        if (*ps->pos == '#' && ps->pos > start && *(ps->pos - 1) == ' ')
            { len--; break; }  /* remove trailing space before # */
        buf[len++] = *ps->pos;
        ps->pos++;
        }
    /* Trim trailing whitespace */
    while (len > 0 && (buf[len-1] == ' ' || buf[len-1] == '\t'))
        len--;
    buf[len] = '\0';
    return DupStr(buf);
}

/* ---- Flow parsing (inline { } and [ ]) ---- */

static YamlNode *ParseFlowValue (ParseState *ps);

static YamlNode *ParseFlowMapping (ParseState *ps)
{
    YamlNode *map = NewNode(YAML_MAPPING);
    ps->pos++;  /* skip { */

    while (*ps->pos && *ps->pos != '}')
        {
        SkipSpaces(ps);
        if (*ps->pos == '}') break;
        if (*ps->pos == ',') { ps->pos++; continue; }

        /* Parse key */
        char *key;
        if (*ps->pos == '"' || *ps->pos == '\'')
            key = ParseQuotedScalar(ps, *ps->pos);
        else
            key = ParseBareScalar(ps, ":,}");

        SkipSpaces(ps);
        if (*ps->pos == ':') ps->pos++;
        SkipSpaces(ps);

        /* Parse value */
        YamlNode *val = ParseFlowValue(ps);
        if (val)
            {
            val->key = key;
            AddChild(map, val);
            }
        else
            free(key);

        SkipSpaces(ps);
        if (*ps->pos == ',') ps->pos++;
        }
    if (*ps->pos == '}') ps->pos++;
    return map;
}

static YamlNode *ParseFlowSequence (ParseState *ps)
{
    YamlNode *seq = NewNode(YAML_SEQUENCE);
    ps->pos++;  /* skip [ */

    while (*ps->pos && *ps->pos != ']')
        {
        SkipSpaces(ps);
        if (*ps->pos == ']') break;
        if (*ps->pos == ',') { ps->pos++; continue; }

        YamlNode *val = ParseFlowValue(ps);
        if (val) AddChild(seq, val);

        SkipSpaces(ps);
        if (*ps->pos == ',') ps->pos++;
        }
    if (*ps->pos == ']') ps->pos++;
    return seq;
}

static YamlNode *ParseFlowValue (ParseState *ps)
{
    SkipSpaces(ps);
    if (*ps->pos == '{')
        return ParseFlowMapping(ps);
    if (*ps->pos == '[')
        return ParseFlowSequence(ps);
    if (*ps->pos == '"' || *ps->pos == '\'')
        {
        YamlNode *n = NewNode(YAML_SCALAR);
        n->scalar = ParseQuotedScalar(ps, *ps->pos);
        return n;
        }
    /* Bare scalar */
    YamlNode *n = NewNode(YAML_SCALAR);
    n->scalar = ParseBareScalar(ps, ",}]");
    return n;
}

/* ---- Block parsing (indentation-based) ---- */

static YamlNode *ParseBlockValue (ParseState *ps, int minIndent);

static YamlNode *ParseBlockMapping (ParseState *ps, int indent)
{
    YamlNode *map = NewNode(YAML_MAPPING);

    while (*ps->pos)
        {
        SkipBlankLines(ps);
        if (!*ps->pos) break;

        int curIndent = GetIndent(ps);
        if (curIndent < indent) break;
        if (curIndent > indent)
            {
            /* Unexpected extra indentation */
            break;
            }

        const char *lineStart = ps->pos;
        ps->pos += indent;
        SkipSpaces(ps);

        /* Check for sequence item */
        if (*ps->pos == '-' && (ps->pos[1] == ' ' || ps->pos[1] == '\n'))
            {
            ps->pos = lineStart;
            break;  /* This is a sequence, not a mapping */
            }

        /* Parse key */
        char *key;
        if (*ps->pos == '"' || *ps->pos == '\'')
            key = ParseQuotedScalar(ps, *ps->pos);
        else
            key = ParseBareScalar(ps, ":");

        if (!key || !*key)
            { if (key) free(key); ps->pos = lineStart; break; }

        SkipSpaces(ps);
        if (*ps->pos != ':')
            {
            /* Not a mapping entry */
            free(key);
            ps->pos = lineStart;
            break;
            }
        ps->pos++;  /* skip : */
        SkipSpaces(ps);

        /* Parse value: inline or next-line block */
        YamlNode *val = NULL;
        if (*ps->pos == '\n' || *ps->pos == '\r' || *ps->pos == '#' || !*ps->pos)
            {
            /* Value on next line(s) */
            SkipComment(ps);
            if (*ps->pos == '\n') { ps->pos++; ps->line++; }
            else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
            val = ParseBlockValue(ps, indent + 1);
            }
        else if (*ps->pos == '{')
            val = ParseFlowMapping(ps);
        else if (*ps->pos == '[')
            val = ParseFlowSequence(ps);
        else
            {
            /* Inline scalar */
            val = NewNode(YAML_SCALAR);
            if (*ps->pos == '"' || *ps->pos == '\'')
                val->scalar = ParseQuotedScalar(ps, *ps->pos);
            else
                val->scalar = ParseBareScalar(ps, "");
            SkipComment(ps);
            if (*ps->pos == '\n') { ps->pos++; ps->line++; }
            else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
            }

        if (val)
            {
            val->key = key;
            AddChild(map, val);
            }
        else
            free(key);
        }

    return map;
}

static YamlNode *ParseBlockSequence (ParseState *ps, int indent)
{
    YamlNode *seq = NewNode(YAML_SEQUENCE);

    while (*ps->pos)
        {
        SkipBlankLines(ps);
        if (!*ps->pos) break;

        int curIndent = GetIndent(ps);
        if (curIndent < indent) break;

        const char *lineStart = ps->pos;
        ps->pos += curIndent;

        if (*ps->pos != '-' || (ps->pos[1] != ' ' && ps->pos[1] != '\n'))
            {
            ps->pos = lineStart;
            break;
            }

        ps->pos++;  /* skip - */
        if (*ps->pos == ' ') ps->pos++;  /* skip space after - */

        /* Parse item value */
        YamlNode *val = NULL;
        SkipSpaces(ps);

        if (*ps->pos == '\n' || *ps->pos == '\r' || *ps->pos == '#' || !*ps->pos)
            {
            /* Block value on next line */
            SkipComment(ps);
            if (*ps->pos == '\n') { ps->pos++; ps->line++; }
            else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
            val = ParseBlockValue(ps, curIndent + 2);
            }
        else if (*ps->pos == '{')
            val = ParseFlowMapping(ps);
        else if (*ps->pos == '[')
            val = ParseFlowSequence(ps);
        else
            {
            /* Check if this is an inline mapping (key: value on same line as -) */
            const char *colon = NULL;
            const char *scan = ps->pos;
            int inQuote = 0;
            while (*scan && *scan != '\n' && *scan != '#')
                {
                if (*scan == '"' || *scan == '\'') inQuote = !inQuote;
                if (*scan == ':' && !inQuote && (scan[1] == ' ' || scan[1] == '\n'))
                    { colon = scan; break; }
                scan++;
                }

            if (colon)
                {
                /* Inline mapping starting on - line.
                   Parse first key:value here, then continue with
                   ParseBlockMapping for subsequent lines at same indent. */
                int itemIndent = (int)(ps->pos - lineStart);
                val = NewNode(YAML_MAPPING);

                /* Parse the first key:value on this line */
                char *firstKey;
                if (*ps->pos == '"' || *ps->pos == '\'')
                    firstKey = ParseQuotedScalar(ps, *ps->pos);
                else
                    firstKey = ParseBareScalar(ps, ":");
                SkipSpaces(ps);
                if (*ps->pos == ':') ps->pos++;
                SkipSpaces(ps);

                YamlNode *firstVal = NULL;
                if (*ps->pos == '{')
                    firstVal = ParseFlowMapping(ps);
                else if (*ps->pos == '[')
                    firstVal = ParseFlowSequence(ps);
                else if (*ps->pos == '\n' || *ps->pos == '\r' || *ps->pos == '#' || !*ps->pos)
                    {
                    SkipComment(ps);
                    if (*ps->pos == '\n') { ps->pos++; ps->line++; }
                    else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
                    firstVal = ParseBlockValue(ps, itemIndent + 1);
                    }
                else
                    {
                    firstVal = NewNode(YAML_SCALAR);
                    if (*ps->pos == '"' || *ps->pos == '\'')
                        firstVal->scalar = ParseQuotedScalar(ps, *ps->pos);
                    else
                        firstVal->scalar = ParseBareScalar(ps, "");
                    SkipComment(ps);
                    if (*ps->pos == '\n') { ps->pos++; ps->line++; }
                    else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
                    }
                if (firstVal)
                    {
                    firstVal->key = firstKey;
                    AddChild(val, firstVal);
                    }
                else if (firstKey)
                    free(firstKey);

                /* Continue reading subsequent keys at the same indent */
                while (*ps->pos)
                    {
                    SkipBlankLines(ps);
                    if (!*ps->pos) break;
                    int nextIndent = GetIndent(ps);
                    if (nextIndent < itemIndent) break;
                    if (nextIndent != itemIndent) break;

                    const char *nl = ps->pos;
                    ps->pos += nextIndent;

                    /* Check for sequence item (next - at parent indent) */
                    if (*ps->pos == '-' && (ps->pos[1] == ' ' || ps->pos[1] == '\n'))
                        { ps->pos = nl; break; }

                    char *nextKey;
                    if (*ps->pos == '"' || *ps->pos == '\'')
                        nextKey = ParseQuotedScalar(ps, *ps->pos);
                    else
                        nextKey = ParseBareScalar(ps, ":");

                    if (!nextKey || !*nextKey)
                        { if (nextKey) free(nextKey); ps->pos = nl; break; }

                    SkipSpaces(ps);
                    if (*ps->pos != ':')
                        { free(nextKey); ps->pos = nl; break; }
                    ps->pos++;
                    SkipSpaces(ps);

                    YamlNode *nextVal = NULL;
                    if (*ps->pos == '{')
                        nextVal = ParseFlowMapping(ps);
                    else if (*ps->pos == '[')
                        nextVal = ParseFlowSequence(ps);
                    else if (*ps->pos == '\n' || *ps->pos == '\r' || *ps->pos == '#' || !*ps->pos)
                        {
                        SkipComment(ps);
                        if (*ps->pos == '\n') { ps->pos++; ps->line++; }
                        else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
                        nextVal = ParseBlockValue(ps, itemIndent + 1);
                        }
                    else
                        {
                        nextVal = NewNode(YAML_SCALAR);
                        if (*ps->pos == '"' || *ps->pos == '\'')
                            nextVal->scalar = ParseQuotedScalar(ps, *ps->pos);
                        else
                            nextVal->scalar = ParseBareScalar(ps, "");
                        SkipComment(ps);
                        if (*ps->pos == '\n') { ps->pos++; ps->line++; }
                        else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
                        }
                    if (nextVal)
                        {
                        nextVal->key = nextKey;
                        AddChild(val, nextVal);
                        }
                    else if (nextKey)
                        free(nextKey);
                    }
                }
            else
                {
                /* Plain scalar */
                val = NewNode(YAML_SCALAR);
                val->scalar = ParseBareScalar(ps, "");
                SkipComment(ps);
                if (*ps->pos == '\n') { ps->pos++; ps->line++; }
                else if (*ps->pos == '\r') { ps->pos++; if (*ps->pos == '\n') ps->pos++; ps->line++; }
                }
            }

        if (val) AddChild(seq, val);
        }

    return seq;
}

static YamlNode *ParseBlockValue (ParseState *ps, int minIndent)
{
    SkipBlankLines(ps);
    if (!*ps->pos) return NULL;

    int indent = GetIndent(ps);
    if (indent < minIndent) return NULL;

    const char *lineStart = ps->pos;
    const char *p = ps->pos + indent;

    /* Detect if this is a sequence (starts with -) or mapping (has key:) */
    if (*p == '-' && (p[1] == ' ' || p[1] == '\n'))
        return ParseBlockSequence(ps, indent);

    /* Check for mapping (look for key: pattern) */
    int inQuote = 0;
    const char *scan = p;
    while (*scan && *scan != '\n' && *scan != '#')
        {
        if (*scan == '"' || *scan == '\'') inQuote = !inQuote;
        if (*scan == ':' && !inQuote && (scan[1] == ' ' || scan[1] == '\n' || scan[1] == '\0'))
            return ParseBlockMapping(ps, indent);
        scan++;
        }

    /* Fallback: scalar */
    ps->pos = p;
    YamlNode *n = NewNode(YAML_SCALAR);
    n->scalar = ParseBareScalar(ps, "");
    SkipComment(ps);
    if (*ps->pos == '\n') { ps->pos++; ps->line++; }
    return n;
}

/* ---- Public API ---- */

YamlNode *YamlParseString (const char *text)
{
    ParseState ps;
    ps.text = text;
    ps.pos = text;
    ps.line = 1;
    ps.errmsg[0] = '\0';

    YamlNode *root = ParseBlockValue(&ps, 0);
    return root;
}

YamlNode *YamlParseFile (const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) return NULL;

    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    char *text = (char *) malloc(size + 1);
    if (!text) { fclose(fp); return NULL; }

    size_t nread = fread(text, 1, size, fp);
    text[nread] = '\0';
    fclose(fp);

    YamlNode *root = YamlParseString(text);
    free(text);
    return root;
}

/* ---- Lookup helpers ---- */

YamlNode *YamlGetChild (YamlNode *node, const char *key)
{
    int i;
    if (!node || node->type != YAML_MAPPING) return NULL;
    for (i = 0; i < node->nChildren; i++)
        if (node->children[i]->key && strcmp(node->children[i]->key, key) == 0)
            return node->children[i];
    return NULL;
}

const char *YamlGetScalar (YamlNode *node, const char *key)
{
    YamlNode *child = YamlGetChild(node, key);
    if (child && child->type == YAML_SCALAR)
        return child->scalar;
    return NULL;
}

int YamlGetInt (YamlNode *node, const char *key, int defaultVal)
{
    const char *s = YamlGetScalar(node, key);
    if (s) return atoi(s);
    return defaultVal;
}

double YamlGetDouble (YamlNode *node, const char *key, double defaultVal)
{
    const char *s = YamlGetScalar(node, key);
    if (s) return atof(s);
    return defaultVal;
}

int YamlIsMissing (const char *scalar)
{
    if (!scalar) return 1;
    if (strcmp(scalar, "~") == 0) return 1;
    if (strcmp(scalar, "?") == 0) return 1;
    if (strcmp(scalar, "null") == 0) return 1;
    if (strcmp(scalar, "missing") == 0) return 1;
    return 0;
}

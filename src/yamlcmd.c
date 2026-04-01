/*
 *  yvyra - Bayesian Phylogenetics for Linguistic Data
 *  Based on MrBayes 3.2.7a by Ronquist et al.
 *
 *  yamlcmd.c: Convert a parsed YAML configuration into a NEXUS string
 *  that the existing yvyra parser can execute. Handles taxa, characters
 *  with named states, cost matrices, weights, tipweights, and analysis
 *  settings.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "yamlparser.h"
#include "yamlcmd.h"

#define BUF_GROW    65536
#define MAX_STATES  24

/* ---- Built-in model templates ---- */

#define TPL_MAX_STATES  6

typedef struct {
    const char  *name;
    const char  *description;
    int         nStates;
    const char  *stateNames[TPL_MAX_STATES];
    double      costs[TPL_MAX_STATES][TPL_MAX_STATES];
} ModelTemplate;

static const ModelTemplate TEMPLATES[] = {
    /* --- General --- */
    {"symmetric",
     "Equal costs both directions",
     2, {"state_0", "state_1"},
     {{0, 1}, {1, 0}}},

    {"irreversible",
     "Forward only, reversal near-impossible",
     2, {"absent", "present"},
     {{0, 1}, {100, 0}}},

    {"biased",
     "Forward common, reverse rare (3:1)",
     2, {"state_0", "state_1"},
     {{0, 1}, {3, 0}}},

    {"strong_bias",
     "Forward common, reverse very rare (10:1)",
     2, {"state_0", "state_1"},
     {{0, 1}, {10, 0}}},

    {"ordered_3",
     "3-state ordered: adjacent cheap, jumps expensive",
     3, {"state_0", "state_1", "state_2"},
     {{0, 1, 10}, {1, 0, 1}, {10, 1, 0}}},

    {"dollo",
     "Trait gained once (rare), lost many times (Dollo parsimony)",
     2, {"absent", "present"},
     {{0, 10}, {1, 0}}},

    {"camin_sokal",
     "Trait gained many times, never lost",
     2, {"absent", "present"},
     {{0, 1}, {100, 0}}},

    /* --- Sound change --- */
    {"loss",
     "Irreversible phonological loss",
     2, {"retained", "lost"},
     {{0, 1}, {100, 0}}},

    {"lenition",
     "Consonant weakening common, fortition rare",
     2, {"fortis", "lenis"},
     {{0, 1}, {2, 0}}},

    {"lenition_3",
     "3-stage lenition: fortis > lenis > lost",
     3, {"fortis", "lenis", "lost"},
     {{0, 2, 4}, {8, 0, 1}, {100, 100, 0}}},

    {"fortition",
     "Consonant strengthening (reverse of lenition)",
     2, {"weak", "strong"},
     {{0, 2}, {1, 0}}},

    {"chain_shift",
     "Cyclic sound changes (forward cycle cheap, reverse expensive)",
     3, {"stage_a", "stage_b", "stage_c"},
     {{0, 1, 10}, {10, 0, 1}, {1, 10, 0}}},

    {"merger",
     "Phonological merger: distinct > merged (near-irreversible)",
     2, {"distinct", "merged"},
     {{0, 1}, {100, 0}}},

    /* --- Morphology --- */
    {"simplification",
     "Morphological simplification (4 stages)",
     4, {"complex", "intermediate", "simple", "absent"},
     {{0, 2, 20, 100}, {15, 0, 1, 3}, {100, 20, 0, 1}, {100, 100, 100, 0}}},

    {"grammaticalization",
     "Lexical > grammatical > clitic > affix (5 stages)",
     5, {"lexical", "semi_grammatical", "grammatical", "clitic", "affix"},
     {{0, 2, 10, 20, 40}, {25, 0, 1, 10, 20}, {100, 25, 0, 2, 10}, {100, 100, 25, 0, 1}, {100, 100, 100, 20, 0}}},

    {"analogy",
     "Regularization by analogy: irregular > regular",
     3, {"irregular", "semi_regular", "regular"},
     {{0, 2, 10}, {10, 0, 1}, {100, 8, 0}}},

    {"morphological_loss",
     "Loss of productive morphology (4 stages)",
     4, {"productive", "residual", "fossilized", "absent"},
     {{0, 3, 10, 100}, {100, 0, 2, 10}, {100, 100, 0, 1}, {100, 100, 100, 0}}},

    /* --- Syntax --- */
    {"word_order",
     "6-way word order (SOV, SVO, VSO, VOS, OVS, OSV)",
     6, {"SOV", "SVO", "VSO", "VOS", "OVS", "OSV"},
     {{0, 2, 5, 8, 10, 12}, {3, 0, 3, 8, 10, 12}, {5, 2, 0, 5, 8, 10},
      {8, 5, 5, 0, 5, 8}, {10, 8, 8, 5, 0, 5}, {12, 10, 10, 8, 5, 0}}},

    {"configurationality",
     "Word order freedom: free > semi_free > semi_fixed > fixed",
     4, {"free", "semi_free", "semi_fixed", "fixed"},
     {{0, 2, 8, 20}, {8, 0, 2, 8}, {20, 8, 0, 1}, {100, 20, 6, 0}}},

    /* --- Manuscript / Stemmatology --- */
    {"transmission",
     "Manuscript copying errors (4 stages of corruption)",
     4, {"original", "minor_error", "major_error", "corrupted"},
     {{0, 2, 8, 20}, {10, 0, 2, 8}, {100, 12, 0, 1}, {100, 100, 15, 0}}},

    {"contamination",
     "Cross-manuscript contamination (4 stages)",
     4, {"pure", "minor_contamination", "mixed", "heavily_contaminated"},
     {{0, 3, 10, 20}, {15, 0, 2, 10}, {100, 20, 0, 1}, {100, 100, 100, 0}}},

    {"lacuna",
     "Text loss / lacunae (4 stages)",
     4, {"complete", "minor_gaps", "major_gaps", "fragmentary"},
     {{0, 2, 8, 20}, {8, 0, 2, 8}, {100, 10, 0, 1}, {100, 100, 12, 0}}},

    {"scribal_error",
     "Scribal copying errors by type",
     4, {"correct", "orthographic", "lexical", "syntactic"},
     {{0, 1, 2, 3}, {3, 0, 2, 3}, {6, 6, 0, 3}, {10, 10, 10, 0}}},

    /* sentinel */
    {NULL, NULL, 0, {NULL}, {{0}}}
};

static const ModelTemplate *FindTemplate (const char *name)
{
    int i;
    for (i = 0; TEMPLATES[i].name != NULL; i++)
        if (strcmp(TEMPLATES[i].name, name) == 0)
            return &TEMPLATES[i];
    return NULL;
}

/* ---- Extern globals from command.c for state name storage ---- */
extern char ***charStateNames;
extern int  *charNStateNames;

/* ---- Globals for estimated weight system ---- */
int     *isWeightEstimated = NULL;  /* per-original-char: 1 if estimated */
double  *initialWeights = NULL;     /* per-original-char: initial weight */

void SetCharStateNames (int charIdx, int nStates, const char **stateNames)
{
    int i;
    if (!charStateNames || !charNStateNames) return;
    if (charIdx < 0) return;

    /* Free old names if any */
    if (charStateNames[charIdx])
        {
        for (i = 0; i < charNStateNames[charIdx]; i++)
            if (charStateNames[charIdx][i])
                free(charStateNames[charIdx][i]);
        free(charStateNames[charIdx]);
        }

    charStateNames[charIdx] = (char **) calloc(nStates, sizeof(char *));
    charNStateNames[charIdx] = nStates;
    for (i = 0; i < nStates; i++)
        {
        charStateNames[charIdx][i] = (char *) malloc(strlen(stateNames[i]) + 1);
        strcpy(charStateNames[charIdx][i], stateNames[i]);
        }
}

#define MAX_TAXA    200
#define MAX_CHARS   500

/* ---- Dynamic string buffer ---- */

typedef struct {
    char    *data;
    int     len;
    int     cap;
} Buf;

static void BufInit (Buf *b)
{
    b->cap = BUF_GROW;
    b->data = (char *) malloc(b->cap);
    b->data[0] = '\0';
    b->len = 0;
}

static void BufAppend (Buf *b, const char *fmt, ...)
{
    va_list ap;
    int needed;

    va_start(ap, fmt);
    needed = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);

    while (b->len + needed + 1 > b->cap)
        {
        b->cap += BUF_GROW;
        b->data = (char *) realloc(b->data, b->cap);
        }

    va_start(ap, fmt);
    vsnprintf(b->data + b->len, b->cap - b->len, fmt, ap);
    va_end(ap);
    b->len += needed;
}

/* ---- State name to index mapping ---- */

typedef struct {
    char    *names[MAX_STATES];
    int     n;
} StateMap;

static int StateIndex (StateMap *sm, const char *name)
{
    int i;
    for (i = 0; i < sm->n; i++)
        if (strcmp(sm->names[i], name) == 0)
            return i;
    return -1;
}

static int StateAdd (StateMap *sm, const char *name)
{
    int idx = StateIndex(sm, name);
    if (idx >= 0) return idx;
    if (sm->n >= MAX_STATES) return -1;
    sm->names[sm->n] = (char *) malloc(strlen(name) + 1);
    strcpy(sm->names[sm->n], name);
    return sm->n++;
}

static void StateFree (StateMap *sm)
{
    int i;
    for (i = 0; i < sm->n; i++)
        free(sm->names[i]);
    sm->n = 0;
}

/* ---- Build state map from model costs or character data ---- */

static void BuildStateMapFromCosts (YamlNode *costsNode, StateMap *sm)
{
    int i;
    sm->n = 0;
    if (!costsNode || costsNode->type != YAML_MAPPING) return;
    for (i = 0; i < costsNode->nChildren; i++)
        if (costsNode->children[i]->key)
            StateAdd(sm, costsNode->children[i]->key);
}

static void BuildStateMapFromData (YamlNode *dataNode, StateMap *sm)
{
    int i;
    sm->n = 0;
    if (!dataNode || dataNode->type != YAML_MAPPING) return;
    for (i = 0; i < dataNode->nChildren; i++)
        {
        YamlNode *val = dataNode->children[i];
        if (val->type == YAML_SCALAR && !YamlIsMissing(val->scalar))
            StateAdd(sm, val->scalar);
        else if (val->type == YAML_MAPPING)
            {
            /* Polymorphic: {state: prob, ...} */
            int j;
            for (j = 0; j < val->nChildren; j++)
                if (val->children[j]->key)
                    StateAdd(sm, val->children[j]->key);
            }
        }
}

static void BuildStateMapFromExplicit (YamlNode *statesNode, StateMap *sm)
{
    int i;
    sm->n = 0;
    if (!statesNode || statesNode->type != YAML_SEQUENCE) return;
    for (i = 0; i < statesNode->nChildren; i++)
        if (statesNode->children[i]->type == YAML_SCALAR)
            StateAdd(sm, statesNode->children[i]->scalar);
}

/* ---- Validation ---- */

static int ValidateYaml (YamlNode *root)
{
    YamlNode *taxaNode, *charsNode, *modelsNode, *analysisNode;
    int i, j, k, nTaxa, nChars, nErrors = 0;

    taxaNode = YamlGetChild(root, "taxa");
    charsNode = YamlGetChild(root, "characters");
    modelsNode = YamlGetChild(root, "models");
    analysisNode = YamlGetChild(root, "analysis");

    /* Required sections */
    if (!taxaNode || taxaNode->type != YAML_SEQUENCE)
        { fprintf(stderr, "Error: 'taxa' section missing or not a list\n"); nErrors++; }
    if (!charsNode || charsNode->type != YAML_SEQUENCE)
        { fprintf(stderr, "Error: 'characters' section missing or not a list\n"); nErrors++; }
    if (nErrors > 0) return nErrors;

    nTaxa = taxaNode->nChildren;
    nChars = charsNode->nChildren;

    /* Minimum taxa */
    if (nTaxa < 4)
        { fprintf(stderr, "Error: need at least 4 taxa (found %d)\n", nTaxa); nErrors++; }
    if (nChars < 1)
        { fprintf(stderr, "Error: need at least 1 character\n"); nErrors++; }
    if (nErrors > 0) return nErrors;

    /* Collect taxon IDs for cross-referencing */
    char *taxaIds[MAX_TAXA];
    for (i = 0; i < nTaxa && i < MAX_TAXA; i++)
        {
        YamlNode *t = taxaNode->children[i];
        const char *id = NULL;
        if (t->type == YAML_MAPPING)
            id = YamlGetScalar(t, "id");
        else if (t->type == YAML_SCALAR)
            id = t->scalar;
        taxaIds[i] = (char *)(id ? id : "");

        /* Check duplicate taxon IDs */
        for (j = 0; j < i; j++)
            if (strcmp(taxaIds[i], taxaIds[j]) == 0)
                { fprintf(stderr, "Error: duplicate taxon id '%s'\n", taxaIds[i]); nErrors++; }
        }

    /* Validate each character */
    for (i = 0; i < nChars && i < MAX_CHARS; i++)
        {
        YamlNode *ch = charsNode->children[i];
        if (ch->type != YAML_MAPPING)
            { fprintf(stderr, "Error: character %d is not a mapping\n", i+1); nErrors++; continue; }

        const char *charId = YamlGetScalar(ch, "id");
        if (!charId || !*charId)
            { fprintf(stderr, "Error: character %d missing 'id'\n", i+1); nErrors++; continue; }

        /* Duplicate character ID check */
        for (j = 0; j < i; j++)
            {
            const char *prevId = YamlGetScalar(charsNode->children[j], "id");
            if (prevId && strcmp(charId, prevId) == 0)
                { fprintf(stderr, "Error: duplicate character id '%s'\n", charId); nErrors++; }
            }

        /* Model reference check */
        const char *model = YamlGetScalar(ch, "model");
        if (model)
            {
            int found = 0;
            if (modelsNode && YamlGetChild(modelsNode, model)) found = 1;
            if (!found && FindTemplate(model)) found = 1;
            if (!found)
                { fprintf(stderr, "Error: character '%s' references unknown model '%s'\n", charId, model); nErrors++; }
            }

        /* Weight check */
        YamlNode *weightNode = YamlGetChild(ch, "weight");
        if (weightNode && weightNode->type == YAML_SCALAR)
            {
            double w = atof(weightNode->scalar);
            if (w < 0.0)
                { fprintf(stderr, "Error: character '%s' has negative weight %g\n", charId, w); nErrors++; }
            }

        /* Data section check */
        YamlNode *dataNode = YamlGetChild(ch, "data");
        if (!dataNode || dataNode->type != YAML_MAPPING)
            { fprintf(stderr, "Error: character '%s' missing 'data' section\n", charId); nErrors++; continue; }

        /* Build state map for this character */
        StateMap sm;
        sm.n = 0;
        YamlNode *statesNode = YamlGetChild(ch, "states");
        if (statesNode && statesNode->type == YAML_SEQUENCE)
            BuildStateMapFromExplicit(statesNode, &sm);
        else if (model)
            {
            if (modelsNode)
                {
                YamlNode *modelNode = YamlGetChild(modelsNode, model);
                if (modelNode)
                    {
                    YamlNode *costsNode = YamlGetChild(modelNode, "costs");
                    if (costsNode) BuildStateMapFromCosts(costsNode, &sm);
                    }
                }
            if (sm.n == 0)
                {
                const ModelTemplate *tpl = FindTemplate(model);
                if (tpl)
                    for (k = 0; k < tpl->nStates; k++)
                        StateAdd(&sm, tpl->stateNames[k]);
                }
            }
        if (sm.n == 0)
            BuildStateMapFromData(dataNode, &sm);

        /* Check every taxon in data exists in taxa list */
        for (j = 0; j < dataNode->nChildren; j++)
            {
            const char *taxon = dataNode->children[j]->key;
            if (!taxon) continue;
            int found = 0;
            for (k = 0; k < nTaxa; k++)
                if (strcmp(taxaIds[k], taxon) == 0) { found = 1; break; }
            if (!found)
                { fprintf(stderr, "Error: character '%s' data references unknown taxon '%s'\n", charId, taxon); nErrors++; }
            }

        /* Check every taxon appears in data (or warn) */
        for (j = 0; j < nTaxa; j++)
            {
            if (!YamlGetChild(dataNode, taxaIds[j]))
                { fprintf(stderr, "Error: character '%s' has no data for taxon '%s' (use ~ for missing)\n", charId, taxaIds[j]); nErrors++; }
            }

        /* Check state values in data match declared states */
        if (sm.n > 0)
            {
            for (j = 0; j < dataNode->nChildren; j++)
                {
                YamlNode *val = dataNode->children[j];
                if (val->type == YAML_SCALAR)
                    {
                    if (!YamlIsMissing(val->scalar) && StateIndex(&sm, val->scalar) < 0)
                        { fprintf(stderr, "Error: character '%s', taxon '%s': unknown state '%s' (valid: ",
                                  charId, val->key ? val->key : "?", val->scalar);
                          for (k = 0; k < sm.n; k++) fprintf(stderr, "%s%s", k ? ", " : "", sm.names[k]);
                          fprintf(stderr, ")\n"); nErrors++; }
                    }
                else if (val->type == YAML_MAPPING)
                    {
                    /* Polymorphic: check keys are valid states */
                    double probSum = 0.0;
                    for (k = 0; k < val->nChildren; k++)
                        {
                        const char *sn = val->children[k]->key;
                        if (sn && StateIndex(&sm, sn) < 0)
                            { fprintf(stderr, "Error: character '%s', taxon '%s': unknown state '%s' in probability\n",
                                      charId, val->key ? val->key : "?", sn); nErrors++; }
                        if (val->children[k]->type == YAML_SCALAR)
                            {
                            double p = atof(val->children[k]->scalar);
                            if (p < 0.0 || p > 1.0)
                                { fprintf(stderr, "Error: character '%s', taxon '%s': probability %g out of [0, 1]\n",
                                          charId, val->key ? val->key : "?", p); nErrors++; }
                            probSum += p;
                            }
                        }
                    if (probSum > 1.0 + 1e-10)
                        { fprintf(stderr, "Error: character '%s', taxon '%s': probabilities sum to %g (> 1.0)\n",
                                  charId, val->key ? val->key : "?", probSum); nErrors++; }
                    }
                }
            }
        StateFree(&sm);
        }

    /* Validate root_prior */
    if (analysisNode)
        {
        const char *rp = YamlGetScalar(analysisNode, "root_prior");
        if (rp && strcmp(rp, "uniform") != 0 && strcmp(rp, "first_state") != 0)
            {
            /* Must be a valid state name in at least one character */
            int found = 0;
            for (j = 0; j < nChars && !found; j++)
                {
                YamlNode *ch2 = charsNode->children[j];
                const char *model2 = YamlGetScalar(ch2, "model");
                StateMap sm2;
                sm2.n = 0;
                if (model2 && modelsNode)
                    {
                    YamlNode *mn = YamlGetChild(modelsNode, model2);
                    if (mn) { YamlNode *cn = YamlGetChild(mn, "costs"); if (cn) BuildStateMapFromCosts(cn, &sm2); }
                    }
                if (sm2.n == 0 && model2)
                    { const ModelTemplate *tpl = FindTemplate(model2); if (tpl) for (k=0;k<tpl->nStates;k++) StateAdd(&sm2, tpl->stateNames[k]); }
                if (sm2.n == 0)
                    { YamlNode *dn = YamlGetChild(ch2, "data"); if (dn) BuildStateMapFromData(dn, &sm2); }
                if (StateIndex(&sm2, rp) >= 0) found = 1;
                StateFree(&sm2);
                }
            if (!found)
                { fprintf(stderr, "Error: root_prior '%s' is not a valid state name in any character\n", rp); nErrors++; }
            }
        }

    /* Validate tree_model */
    {
        const char *tm = analysisNode ? YamlGetScalar(analysisNode, "tree_model") : NULL;
        if (tm && strcmp(tm, "clock") != 0 && strcmp(tm, "unrooted") != 0
               && strcmp(tm, "unconstrained") != 0)
            { fprintf(stderr, "Error: tree_model must be 'clock' or 'unrooted' (got '%s')\n", tm); nErrors++; }

        /* Validate clock model */
        const char *cm = analysisNode ? YamlGetScalar(analysisNode, "clock") : NULL;
        if (cm && strcmp(cm,"strict")!=0 && strcmp(cm,"igr")!=0 && strcmp(cm,"iln")!=0
               && strcmp(cm,"tk02")!=0 && strcmp(cm,"cpp")!=0 && strcmp(cm,"wn")!=0)
            { fprintf(stderr, "Error: clock must be one of: strict, igr, iln, tk02, cpp, wn (got '%s')\n", cm); nErrors++; }

        /* Determine if clock model is active (for cross-validation below) */
        int isClock = 1;
        if (tm && (strcmp(tm, "unrooted") == 0 || strcmp(tm, "unconstrained") == 0))
            isClock = 0;

        /* Validate outgroup */
        const char *og = analysisNode ? YamlGetScalar(analysisNode, "outgroup") : NULL;
        if (og)
            {
            if (isClock)
                { fprintf(stderr, "Error: 'outgroup' cannot be used with tree_model 'clock' (clock trees are rooted by the model)\n"); nErrors++; }
            else
                {
                int found = 0;
                for (j = 0; j < nTaxa; j++)
                    if (strcmp(taxaIds[j], og) == 0) { found = 1; break; }
                if (!found)
                    { fprintf(stderr, "Error: outgroup '%s' is not a valid taxon id\n", og); nErrors++; }
                }
            }
        else if (!isClock)
            {
            /* No explicit outgroup in unrooted mode — check reserved name is not used */
            for (j = 0; j < nTaxa; j++)
                if (strcmp(taxaIds[j], "_outgroup_") == 0)
                    { fprintf(stderr, "Error: taxon id '_outgroup_' is reserved (use 'outgroup:' in analysis section to set a custom outgroup)\n"); nErrors++; }
            }
    }

    /* Validate models section */
    if (modelsNode && modelsNode->type == YAML_MAPPING)
        {
        for (i = 0; i < modelsNode->nChildren; i++)
            {
            YamlNode *model = modelsNode->children[i];
            if (!model->key || model->type != YAML_MAPPING) continue;

            YamlNode *costsNode = YamlGetChild(model, "costs");
            if (!costsNode || costsNode->type != YAML_MAPPING)
                { fprintf(stderr, "Error: model '%s' missing 'costs' section\n", model->key); nErrors++; continue; }

            int nStates = costsNode->nChildren;

            /* Check each row */
            for (j = 0; j < costsNode->nChildren; j++)
                {
                YamlNode *row = costsNode->children[j];
                int nCols = 0;

                if (row->type == YAML_SEQUENCE)
                    nCols = row->nChildren;
                else if (row->type == YAML_SCALAR)
                    {
                    /* Count space-separated values */
                    const char *p = row->scalar;
                    while (*p)
                        {
                        while (*p == ' ' || *p == '\t') p++;
                        if (!*p) break;
                        nCols++;
                        while (*p && *p != ' ' && *p != '\t') p++;
                        }
                    }

                if (nCols != nStates)
                    { fprintf(stderr, "Error: model '%s', row '%s': expected %d values, got %d\n",
                              model->key, row->key ? row->key : "?", nStates, nCols); nErrors++; }

                /* Check diagonal is zero */
                if (row->type == YAML_SEQUENCE && j < row->nChildren)
                    {
                    double diag = atof(row->children[j]->scalar ? row->children[j]->scalar : "0");
                    if (fabs(diag) > 1e-10)
                        { fprintf(stderr, "Error: model '%s', row '%s': diagonal must be 0 (got %g)\n",
                                  model->key, row->key ? row->key : "?", diag); nErrors++; }
                    }
                }
            }
        }

    return nErrors;
}

/* ---- Main conversion ---- */

char *YamlToNexus (const char *yamlFilename)
{
    YamlNode *root, *taxaNode, *charsNode, *modelsNode, *analysisNode, *diagNode;
    int i, j, k, nTaxa, nChars;
    Buf buf;

    /* Parse YAML file */
    root = YamlParseFile(yamlFilename);
    if (!root)
        {
        fprintf(stderr, "Error: could not parse YAML file '%s'\n", yamlFilename);
        return NULL;
        }

    /* Validate before conversion — all errors are fatal */
    {
        int nErrors = ValidateYaml(root);
        if (nErrors > 0)
            {
            fprintf(stderr, "\n%d error(s) found in '%s'. Fix all errors before running.\n", nErrors, yamlFilename);
            YamlFreeNode(root);
            return NULL;
            }
    }

    taxaNode = YamlGetChild(root, "taxa");
    charsNode = YamlGetChild(root, "characters");
    modelsNode = YamlGetChild(root, "models");
    analysisNode = YamlGetChild(root, "analysis");
    diagNode = YamlGetChild(root, "diagnostics");

    /* Tree model: clock (default) or unrooted */
    const char *treeModel = analysisNode ? YamlGetScalar(analysisNode, "tree_model") : NULL;
    int useClock = 1;  /* default: clock */
    if (treeModel && (strcmp(treeModel, "unrooted") == 0 || strcmp(treeModel, "unconstrained") == 0))
        useClock = 0;

    /* Detect polymorphic/confidence data (tipweights) — incompatible with clock trees */
    if (useClock && !treeModel)
        {
        int hasTipweights = 0;
        for (i = 0; i < charsNode->nChildren && !hasTipweights; i++)
            {
            YamlNode *ch = charsNode->children[i];
            YamlNode *dataNode = ch->type == YAML_MAPPING ? YamlGetChild(ch, "data") : NULL;
            if (!dataNode) continue;
            for (j = 0; j < dataNode->nChildren && !hasTipweights; j++)
                if (dataNode->children[j]->type == YAML_MAPPING)
                    hasTipweights = 1;
            }
        if (hasTipweights)
            {
            fprintf(stderr, "Note: confidence-weighted data detected; using unrooted tree model\n");
            fprintf(stderr, "      (clock trees do not yet support tipweights; set tree_model: clock to override)\n");
            useClock = 0;
            }
        }

    /* Clock model: igr (default), strict, iln, tk02, cpp, wn */
    const char *clockModel = analysisNode ? YamlGetScalar(analysisNode, "clock") : NULL;
    if (!clockModel) clockModel = "igr";

    /* Outgroup handling: dummy only for unrooted trees without explicit outgroup */
    const char *userOutgroup = analysisNode ? YamlGetScalar(analysisNode, "outgroup") : NULL;
    int useDummyOutgroup = (!useClock && userOutgroup == NULL) ? 1 : 0;

    /* Validation passed — safe to proceed */
    nTaxa = taxaNode->nChildren;
    nChars = charsNode->nChildren;

    /* Collect taxa IDs and names */
    char *taxaIds[MAX_TAXA];
    char *taxaNames[MAX_TAXA];
    for (i = 0; i < nTaxa && i < MAX_TAXA; i++)
        {
        YamlNode *t = taxaNode->children[i];
        const char *id = NULL;
        const char *name = NULL;
        if (t->type == YAML_MAPPING)
            {
            id = YamlGetScalar(t, "id");
            name = YamlGetScalar(t, "name");
            }
        else if (t->type == YAML_SCALAR)
            id = t->scalar;
        taxaIds[i] = (char *)(id ? id : "unknown");
        taxaNames[i] = (char *)(name ? name : id ? id : "unknown");
        }

    /* Build per-character state maps and model references */
    StateMap charStates[MAX_CHARS];
    char *charModelRef[MAX_CHARS];
    char *charIds[MAX_CHARS];
    double charWeights[MAX_CHARS];

    for (i = 0; i < nChars && i < MAX_CHARS; i++)
        {
        YamlNode *ch = charsNode->children[i];
        charStates[i].n = 0;
        charModelRef[i] = NULL;
        charIds[i] = NULL;
        charWeights[i] = 1.0;

        if (ch->type != YAML_MAPPING) continue;

        const char *id = YamlGetScalar(ch, "id");
        charIds[i] = (char *)(id ? id : "?");

        const char *model = YamlGetScalar(ch, "model");
        charModelRef[i] = (char *)model;

        /* Parse weight: number or "estimated" */
        {
            const char *wStr = YamlGetScalar(ch, "weight");
            if (wStr && strcmp(wStr, "estimated") == 0)
                {
                charWeights[i] = 1.0;  /* initial value for MCMC */
                /* Mark for estimation (allocated later) */
                if (!isWeightEstimated)
                    {
                    isWeightEstimated = (int *) calloc(MAX_CHARS, sizeof(int));
                    initialWeights = (double *) calloc(MAX_CHARS, sizeof(double));
                    }
                if (isWeightEstimated && i < MAX_CHARS)
                    {
                    isWeightEstimated[i] = 1;
                    initialWeights[i] = 1.0;
                    }
                }
            else
                charWeights[i] = YamlGetDouble(ch, "weight", 1.0);
        }

        /* Build state map: explicit states > model costs > template > data */
        YamlNode *statesNode = YamlGetChild(ch, "states");
        if (statesNode && statesNode->type == YAML_SEQUENCE)
            BuildStateMapFromExplicit(statesNode, &charStates[i]);
        else if (model && modelsNode)
            {
            YamlNode *modelNode = YamlGetChild(modelsNode, model);
            if (modelNode)
                {
                YamlNode *costsNode = YamlGetChild(modelNode, "costs");
                if (costsNode)
                    BuildStateMapFromCosts(costsNode, &charStates[i]);
                }
            }

        /* Template fallback: if model name matches a built-in template */
        if (charStates[i].n == 0 && model)
            {
            const ModelTemplate *tpl = FindTemplate(model);
            if (tpl)
                {
                int s;
                for (s = 0; s < tpl->nStates; s++)
                    StateAdd(&charStates[i], tpl->stateNames[s]);
                }
            }

        if (charStates[i].n == 0)
            {
            YamlNode *dataNode = YamlGetChild(ch, "data");
            if (dataNode)
                BuildStateMapFromData(dataNode, &charStates[i]);
            }
        }

    /* Determine max number of states across all characters */
    int maxStates = 0;
    for (i = 0; i < nChars; i++)
        if (charStates[i].n > maxStates)
            maxStates = charStates[i].n;

    /* Build symbols string */
    char symbols[MAX_STATES + 1];
    for (i = 0; i < maxStates && i < 10; i++)
        symbols[i] = '0' + i;
    symbols[maxStates < 10 ? maxStates : 10] = '\0';

    /* ---- Generate NEXUS ---- */
    BufInit(&buf);

    BufAppend(&buf, "#NEXUS\n\n");
    BufAppend(&buf, "Begin data;\n");
    BufAppend(&buf, "    Dimensions ntax=%d nchar=%d;\n",
              useDummyOutgroup ? nTaxa + 1 : nTaxa, nChars);
    BufAppend(&buf, "    Format datatype=standard symbols=\"%s\" missing=? gap=-;\n", symbols);
    BufAppend(&buf, "    Matrix\n");

    /* Dummy outgroup taxon with all-missing data (absorbs calculation root) */
    if (useDummyOutgroup)
        {
        BufAppend(&buf, "%-12s ", "_outgroup_");
        for (j = 0; j < nChars; j++)
            BufAppend(&buf, "?");
        BufAppend(&buf, "\n");
        }

    /* Build matrix row for each taxon */
    for (i = 0; i < nTaxa; i++)
        {
        BufAppend(&buf, "%-12s ", taxaNames[i]);
        for (j = 0; j < nChars; j++)
            {
            YamlNode *ch = charsNode->children[j];
            YamlNode *dataNode = YamlGetChild(ch, "data");
            YamlNode *taxonData = dataNode ? YamlGetChild(dataNode, taxaIds[i]) : NULL;

            if (!taxonData || (taxonData->type == YAML_SCALAR && YamlIsMissing(taxonData->scalar)))
                {
                BufAppend(&buf, "?");
                }
            else if (taxonData->type == YAML_SCALAR)
                {
                int idx = StateIndex(&charStates[j], taxonData->scalar);
                if (idx >= 0)
                    BufAppend(&buf, "%d", idx);
                else
                    BufAppend(&buf, "?");
                }
            else if (taxonData->type == YAML_MAPPING)
                {
                /* Polymorphic — encode as ? (handled via tipweights) */
                BufAppend(&buf, "?");
                }
            }
        BufAppend(&buf, "\n");
        }

    BufAppend(&buf, "    ;\nEnd;\n\n");

    /* ---- yvyra block ---- */
    BufAppend(&buf, "Begin yvyra;\n");
    BufAppend(&buf, "    set autoclose=yes nowarnings=yes quitonerror=no;\n");

    /* Seed */
    if (analysisNode)
        {
        int seed = YamlGetInt(analysisNode, "seed", 0);
        if (seed > 0)
            BufAppend(&buf, "    set seed=%d;\n", seed);
        }

    /* Outgroup */
    if (useDummyOutgroup)
        BufAppend(&buf, "    outgroup _outgroup_;\n");
    else if (userOutgroup)
        {
        /* Map taxon id to the display name used in the matrix */
        const char *ogDisplayName = userOutgroup;
        for (i = 0; i < nTaxa; i++)
            if (strcmp(taxaIds[i], userOutgroup) == 0)
                { ogDisplayName = taxaNames[i]; break; }
        BufAppend(&buf, "    outgroup %s;\n", ogDisplayName);
        }

    /* Character labels */
    BufAppend(&buf, "    charlabels");
    for (j = 0; j < nChars; j++)
        BufAppend(&buf, " %d=%s", j + 1, charIds[j]);
    BufAppend(&buf, ";\n\n");

    /* Models (usercost commands) */
    if (modelsNode && modelsNode->type == YAML_MAPPING)
        {
        for (i = 0; i < modelsNode->nChildren; i++)
            {
            YamlNode *model = modelsNode->children[i];
            if (!model->key || model->type != YAML_MAPPING) continue;

            YamlNode *costsNode = YamlGetChild(model, "costs");
            if (!costsNode || costsNode->type != YAML_MAPPING) continue;

            StateMap modelStates;
            modelStates.n = 0;
            BuildStateMapFromCosts(costsNode, &modelStates);

            BufAppend(&buf, "    usercost %s (Standard) = %d\n", model->key, modelStates.n);

            /* Write cost matrix rows */
            for (j = 0; j < costsNode->nChildren; j++)
                {
                YamlNode *row = costsNode->children[j];
                BufAppend(&buf, "        ");

                if (row->type == YAML_SEQUENCE)
                    {
                    for (k = 0; k < row->nChildren; k++)
                        {
                        if (k > 0) BufAppend(&buf, " ");
                        BufAppend(&buf, "%s", row->children[k]->scalar ? row->children[k]->scalar : "0");
                        }
                    }
                else if (row->type == YAML_SCALAR)
                    {
                    /* Row as space-separated values in a single scalar */
                    BufAppend(&buf, "%s", row->scalar);
                    }
                BufAppend(&buf, "\n");
                }
            BufAppend(&buf, "    ;\n");
            StateFree(&modelStates);
            }
        }

    /* Emit usercost for template models not defined in the models: section */
    {
        /* Collect which template names are referenced by characters */
        int t;
        for (j = 0; j < nChars; j++)
            {
            if (!charModelRef[j]) continue;
            /* Skip if already defined in models: section */
            if (modelsNode && YamlGetChild(modelsNode, charModelRef[j])) continue;

            const ModelTemplate *tpl = FindTemplate(charModelRef[j]);
            if (!tpl) continue;

            /* Check if we already emitted this template */
            int alreadyEmitted = 0;
            for (k = 0; k < j; k++)
                if (charModelRef[k] && strcmp(charModelRef[k], charModelRef[j]) == 0)
                    { alreadyEmitted = 1; break; }
            if (alreadyEmitted) continue;

            BufAppend(&buf, "    usercost %s (Standard) = %d\n", tpl->name, tpl->nStates);
            for (i = 0; i < tpl->nStates; i++)
                {
                BufAppend(&buf, "        ");
                for (k = 0; k < tpl->nStates; k++)
                    {
                    if (k > 0) BufAppend(&buf, " ");
                    BufAppend(&buf, "%g", tpl->costs[i][k]);
                    }
                BufAppend(&buf, "\n");
                }
            BufAppend(&buf, "    ;\n");
            }
    }

    /* Ctype assignments */
    if (modelsNode)
        {
        for (i = 0; i < modelsNode->nChildren; i++)
            {
            if (!modelsNode->children[i]->key) continue;
            const char *modelName = modelsNode->children[i]->key;

            BufAppend(&buf, "    ctype %s:", modelName);
            int count = 0;
            for (j = 0; j < nChars; j++)
                {
                if (charModelRef[j] && strcmp(charModelRef[j], modelName) == 0)
                    {
                    BufAppend(&buf, " %d", j + 1);
                    count++;
                    }
                }
            if (count > 0)
                BufAppend(&buf, ";\n");
            else
                {
                /* Remove the partial line */
                buf.len -= (int)strlen(modelName) + 12 + (count == 0 ? 0 : 0);
                /* Actually just skip */
                buf.data[buf.len] = '\0';
                }
            }
        }

    /* Ctype assignments for template models */
    {
        for (j = 0; j < nChars; j++)
            {
            if (!charModelRef[j]) continue;
            if (modelsNode && YamlGetChild(modelsNode, charModelRef[j])) continue;
            if (!FindTemplate(charModelRef[j])) continue;

            /* Check if already emitted ctype for this template */
            int done = 0;
            for (k = 0; k < j; k++)
                if (charModelRef[k] && strcmp(charModelRef[k], charModelRef[j]) == 0)
                    { done = 1; break; }
            if (done) continue;

            BufAppend(&buf, "    ctype %s:", charModelRef[j]);
            for (k = 0; k < nChars; k++)
                if (charModelRef[k] && strcmp(charModelRef[k], charModelRef[j]) == 0)
                    BufAppend(&buf, " %d", k + 1);
            BufAppend(&buf, ";\n");
            }
    }

    /* Character weights */
    {
        int hasWeights = 0;
        for (j = 0; j < nChars; j++)
            if (fabs(charWeights[j] - 1.0) > 1e-10)
                hasWeights = 1;

        if (hasWeights)
            {
            BufAppend(&buf, "    wtset * =");
            for (j = 0; j < nChars; j++)
                if (fabs(charWeights[j] - 1.0) > 1e-10)
                    BufAppend(&buf, " %d:%g", j + 1, charWeights[j]);
            BufAppend(&buf, ";\n");
            }
    }

    /* Tipweights from polymorphic data */
    for (j = 0; j < nChars; j++)
        {
        YamlNode *ch = charsNode->children[j];
        YamlNode *dataNode = YamlGetChild(ch, "data");
        if (!dataNode) continue;

        for (i = 0; i < nTaxa; i++)
            {
            YamlNode *td = YamlGetChild(dataNode, taxaIds[i]);
            if (!td || td->type != YAML_MAPPING) continue;

            /* This is a polymorphic or confidence-weighted entry */
            double probs[MAX_STATES];
            int stateSet[MAX_STATES];
            int nSet = 0;
            double sumProb = 0.0;

            memset(probs, 0, sizeof(probs));

            for (k = 0; k < td->nChildren; k++)
                {
                if (!td->children[k]->key) continue;
                int sIdx = StateIndex(&charStates[j], td->children[k]->key);
                if (sIdx >= 0 && td->children[k]->type == YAML_SCALAR)
                    {
                    double p = atof(td->children[k]->scalar);
                    probs[sIdx] = p;
                    stateSet[nSet++] = sIdx;
                    sumProb += p;
                    }
                }

            /* Distribute remaining probability to unspecified states */
            if (sumProb < 1.0 - 1e-10)
                {
                double remain = 1.0 - sumProb;
                int nUnspec = 0;
                for (k = 0; k < charStates[j].n; k++)
                    {
                    int found = 0;
                    int m;
                    for (m = 0; m < td->nChildren; m++)
                        if (td->children[m]->key && StateIndex(&charStates[j], td->children[m]->key) == k)
                            { found = 1; break; }
                    if (!found) nUnspec++;
                    }
                if (nUnspec > 0)
                    {
                    double each = remain / nUnspec;
                    for (k = 0; k < charStates[j].n; k++)
                        {
                        int found = 0;
                        int m;
                        for (m = 0; m < td->nChildren; m++)
                            if (td->children[m]->key && StateIndex(&charStates[j], td->children[m]->key) == k)
                                { found = 1; break; }
                        if (!found) probs[k] = each;
                        }
                    }
                }

            /* Emit tipweights command */
            BufAppend(&buf, "    tipweights %s %d =", taxaNames[i], j + 1);
            for (k = 0; k < charStates[j].n; k++)
                if (probs[k] > 0.0)
                    BufAppend(&buf, " %d:%g", k, probs[k]);
            BufAppend(&buf, ";\n");
            }
        }

    BufAppend(&buf, "\n");

    /* Analysis settings */
    {
        const char *coding = analysisNode ? YamlGetScalar(analysisNode, "coding") : NULL;
        if (coding)
            BufAppend(&buf, "    lset coding=%s;\n", coding);
        else
            BufAppend(&buf, "    lset coding=variable;\n");

        /* Clock tree priors */
        if (useClock)
            {
            BufAppend(&buf, "    prset brlenspr=clock:birthdeath;\n");
            BufAppend(&buf, "    prset clockratepr=fixed(1.0);\n");
            BufAppend(&buf, "    prset treeagepr=gamma(1,1);\n");
            BufAppend(&buf, "    prset nodeagepr=unconstrained;\n");

            if (strcmp(clockModel, "strict") == 0)
                BufAppend(&buf, "    prset clockvarpr=strict;\n");
            else if (strcmp(clockModel, "igr") == 0)
                {
                BufAppend(&buf, "    prset clockvarpr=igr;\n");
                BufAppend(&buf, "    prset igrvarpr=exponential(1.0);\n");
                }
            else if (strcmp(clockModel, "iln") == 0)
                {
                BufAppend(&buf, "    prset clockvarpr=iln;\n");
                BufAppend(&buf, "    prset ilnvarpr=exponential(1.0);\n");
                }
            else if (strcmp(clockModel, "tk02") == 0)
                {
                BufAppend(&buf, "    prset clockvarpr=tk02;\n");
                BufAppend(&buf, "    prset tk02varpr=exponential(1.0);\n");
                }
            else if (strcmp(clockModel, "cpp") == 0)
                BufAppend(&buf, "    prset clockvarpr=cpp;\n");
            else if (strcmp(clockModel, "wn") == 0)
                {
                BufAppend(&buf, "    prset clockvarpr=wn;\n");
                BufAppend(&buf, "    prset wnvarpr=exponential(1.0);\n");
                }
            }

        /* Root state prior.
           Upstream MrBayes's statefreqpr for variable-state Mk only supports symmetric
           Dirichlet or Fixed(Equal). Asymmetric per-state priors require
           per-character partitions which are not currently supported.

           For "first_state": use a symmetric Dirichlet with low concentration
           (alpha=0.5) which allows the data to pull frequencies away from uniform,
           combined with a note that true first-state priors need future work.
           For "uniform": use the default Fixed(Equal). */
        {
            const char *rootPrior = analysisNode ? YamlGetScalar(analysisNode, "root_prior") : NULL;

            if (rootPrior)
                {
                if (strcmp(rootPrior, "uniform") == 0 || strcmp(rootPrior, "equal") == 0)
                    BufAppend(&buf, "    prset statefreqpr=Fixed(Equal);\n");
                else if (strcmp(rootPrior, "first_state") == 0)
                    {
                    /* Low Dirichlet alpha allows data-driven frequency estimation,
                       which tends to favor the more common (ancestral) state */
                    BufAppend(&buf, "    prset statefreqpr=Dirichlet(1);\n");
                    BufAppend(&buf, "    [ root_prior=first_state: Dirichlet prior on state frequencies ]\n");
                    }
                else
                    {
                    /* Named state — same treatment as first_state for now */
                    BufAppend(&buf, "    prset statefreqpr=Dirichlet(1);\n");
                    BufAppend(&buf, "    [ root_prior=%s: Dirichlet prior on state frequencies ]\n", rootPrior);
                    }
                }

            /* Per-character overrides noted but deferred */
            YamlNode *rootOverrides = analysisNode ? YamlGetChild(analysisNode, "root_prior_overrides") : NULL;
            if (rootOverrides && rootOverrides->type == YAML_MAPPING && rootOverrides->nChildren > 0)
                BufAppend(&buf, "    [ Note: per-character root_prior_overrides require multi-partition support ]\n");
        }

        const char *sitelikes = analysisNode ? YamlGetScalar(analysisNode, "sitelikes") : NULL;
        if (sitelikes && (strcmp(sitelikes, "yes") == 0 || strcmp(sitelikes, "Yes") == 0))
            BufAppend(&buf, "    report sitelikes=yes;\n");

        int ngen = useClock ? 100000 : 50000;
        int samplefr = 50, nrun = 1, nchain = 1;
        int printfr = 5000, diagnfr = 5000;

        if (analysisNode)
            {
            ngen = YamlGetInt(analysisNode, "iterations",
                   YamlGetInt(analysisNode, "ngen", ngen));
            samplefr = YamlGetInt(analysisNode, "sample_frequency",
                       YamlGetInt(analysisNode, "samplefr", samplefr));
            nrun = YamlGetInt(analysisNode, "runs",
                   YamlGetInt(analysisNode, "nrun", nrun));
            nchain = YamlGetInt(analysisNode, "chains",
                     YamlGetInt(analysisNode, "nchain", nchain));
            printfr = YamlGetInt(analysisNode, "print_frequency",
                      YamlGetInt(analysisNode, "printfr", printfr));
            diagnfr = YamlGetInt(analysisNode, "diagnostic_frequency",
                      YamlGetInt(analysisNode, "diagnfr", diagnfr));
            }

        BufAppend(&buf, "    mcmcp nrun=%d nchain=%d ngen=%d samplefr=%d;\n",
                  nrun, nchain, ngen, samplefr);
        BufAppend(&buf, "    mcmcp printfr=%d diagnfr=%d;\n", printfr, diagnfr);
        BufAppend(&buf, "    mcmc;\n");
        if (useDummyOutgroup)
            BufAppend(&buf, "    delete _outgroup_;\n");
        BufAppend(&buf, "    sump;\n");
        BufAppend(&buf, "    sumt;\n");
        BufAppend(&buf, "    convergence;\n");
    }

    /* Diagnostics */
    if (diagNode && diagNode->type == YAML_SEQUENCE)
        {
        for (i = 0; i < diagNode->nChildren; i++)
            {
            YamlNode *item = diagNode->children[i];
            if (item->type == YAML_SCALAR)
                {
                /* bare command like "asr_entropy" */
                if (strcmp(item->scalar, "asr_entropy") == 0)
                    BufAppend(&buf, "    asrentropy;\n");
                }
            else if (item->type == YAML_MAPPING && item->nChildren > 0)
                {
                const char *cmd = item->children[0]->key;
                YamlNode *cladeNode = item->children[0];

                if (cladeNode->type == YAML_SEQUENCE)
                    {
                    BufAppend(&buf, "    %s clade={", cmd);
                    for (k = 0; k < cladeNode->nChildren; k++)
                        {
                        if (k > 0) BufAppend(&buf, ",");
                        if (cladeNode->children[k]->type == YAML_SCALAR)
                            BufAppend(&buf, "%s", cladeNode->children[k]->scalar);
                        }
                    BufAppend(&buf, "};\n");
                    }
                }
            }
        }

    BufAppend(&buf, "End;\n");

    /* Store per-character state names for diagnostic output.
       Allocate globals if not already done (they'll be used after
       the NEXUS parser reads the matrix and sets numChar). */
    {
        /* Allocate globals — use nChars from YAML (matches numChar after parse) */
        if (!charStateNames)
            charStateNames = (char ***) calloc(nChars, sizeof(char **));
        if (!charNStateNames)
            charNStateNames = (int *) calloc(nChars, sizeof(int));

        if (charStateNames && charNStateNames)
            {
            for (j = 0; j < nChars; j++)
                {
                if (charStates[j].n > 0)
                    SetCharStateNames(j, charStates[j].n,
                                      (const char **)charStates[j].names);
                }
            }
    }

    /* Clean up */
    for (i = 0; i < nChars; i++)
        StateFree(&charStates[i]);
    YamlFreeNode(root);

    return buf.data;
}

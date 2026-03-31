/*
 *  yvyra - Bayesian Phylogenetics for Linguistic Data
 *  Based on MrBayes 3.2.7a by Ronquist et al.
 *
 *  Copyright (C) 2026 Tiago Tresoldi
 *
 *  diagnostics.c: Character diagnostics, sensitivity analysis, and
 *  consistency index. Post-processing commands that operate on
 *  MCMC output files (.p and .t).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bayes.h"
#include "charweight.h"
#include "command.h"
#include "mcmc.h"
#include "utils.h"
#include "diagnostics.h"

/* Maximum line length for reading .p/.t files */
#define DIAG_MAX_LINE   1000000
#define DIAG_MAX_COLS   1000

/* Clade specification for chardiag/sensitivity */
static char     *diagCladeTaxa[100];
static int      diagNCladeTaxa = 0;
static char     diagCladeStr[1000];


/*--------------------------------------------------------------
|   ReadSiteLnL: Read per-site log-likelihoods from .p file
|   Applies relative burnin (discards first relBurnin fraction)
---------------------------------------------------------------*/
int ReadSiteLnL (char *filename, YFlt relBurnin, SiteLnLData *data)
{
    FILE    *fp;
    char    *line, *tok;
    int     i, j, nCols, nLnLCols, lnLColIdx[DIAG_MAX_COLS];
    int     totalLnLCol = -1;
    int     nTotalLines, nBurnin, lineCount;
    YFlt  *matrix;
    char    slkFilename[310];
    char    *actualFile;

    data->nSamples = 0;
    data->nChars = 0;
    data->siteLnL = NULL;
    data->totalLnL = NULL;
    data->charNames = NULL;

    /* Try .slk file first (new format), fall back to .p file */
    {
        int len = strlen(filename);
        if (len > 2 && strcmp(filename + len - 2, ".p") == 0)
            {
            strncpy(slkFilename, filename, len - 2);
            slkFilename[len - 2] = '\0';
            strcat(slkFilename, ".slk");
            }
        else
            slkFilename[0] = '\0';
    }

    fp = NULL;
    actualFile = filename;
    if (slkFilename[0] != '\0')
        {
        fp = fopen(slkFilename, "r");
        if (fp)
            actualFile = slkFilename;
        }
    if (!fp)
        {
        fp = fopen(filename, "r");
        actualFile = filename;
        }
    if (!fp)
        {
        YvyraPrint ("%s   Could not open file '%s'\n", spacer, filename);
        return (ERROR);
        }

    line = (char *) SafeMalloc(DIAG_MAX_LINE);
    if (!line) { fclose(fp); return (ERROR); }

    /* Skip comment lines (start with [ ) */
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        if (line[0] != '[' && line[0] != '\n' && line[0] != '\r')
            break;
        }

    /* Parse header line */
    nCols = 0;
    nLnLCols = 0;
    tok = strtok(line, "\t\n\r");
    while (tok)
        {
        if (strncmp(tok, "lnL(", 4) == 0)
            {
            lnLColIdx[nLnLCols++] = nCols;
            }
        else if (strcmp(tok, "lnLike") == 0)
            {
            totalLnLCol = nCols;
            }
        nCols++;
        tok = strtok(NULL, "\t\n\r");
        }

    if (nLnLCols == 0)
        {
        YvyraPrint ("%s   No per-site lnL columns found in '%s'\n", spacer, filename);
        YvyraPrint ("%s   Run the analysis with 'report sitelikes=yes' first\n", spacer);
        free(line); fclose(fp);
        return (ERROR);
        }

    /* Store column names */
    data->nChars = nLnLCols;
    data->charNames = (char **) SafeCalloc(nLnLCols, sizeof(char *));

    /* Re-read header for names */
    rewind(fp);
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        if (line[0] != '[' && line[0] != '\n' && line[0] != '\r')
            break;
        }
    nCols = 0;
    j = 0;
    tok = strtok(line, "\t\n\r");
    while (tok)
        {
        if (strncmp(tok, "lnL(", 4) == 0)
            {
            data->charNames[j] = (char *) SafeMalloc(strlen(tok) + 1);
            strcpy(data->charNames[j], tok);
            j++;
            }
        nCols++;
        tok = strtok(NULL, "\t\n\r");
        }

    /* Count data lines */
    nTotalLines = 0;
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        if (line[0] != '\n' && line[0] != '\r' && line[0] != '[')
            nTotalLines++;
        }

    nBurnin = (int)(relBurnin * nTotalLines);
    data->nSamples = nTotalLines - nBurnin;

    if (data->nSamples <= 0)
        {
        YvyraPrint ("%s   No samples after burnin in '%s'\n", spacer, filename);
        free(line); fclose(fp);
        return (ERROR);
        }

    /* Allocate */
    data->siteLnL = (YFlt *) SafeCalloc(data->nSamples * nLnLCols, sizeof(YFlt));
    data->totalLnL = (YFlt *) SafeCalloc(data->nSamples, sizeof(YFlt));
    if (!data->siteLnL || !data->totalLnL)
        {
        free(line); fclose(fp);
        return (ERROR);
        }

    /* Re-read and extract data */
    rewind(fp);
    /* Skip comment + header */
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        if (line[0] != '[' && line[0] != '\n' && line[0] != '\r')
            break;  /* header */
        }

    lineCount = 0;
    i = 0;
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '[')
            continue;
        lineCount++;
        if (lineCount <= nBurnin)
            continue;

        /* Parse columns */
        nCols = 0;
        j = 0;
        tok = strtok(line, "\t\n\r");
        while (tok)
            {
            if (nCols == totalLnLCol)
                data->totalLnL[i] = atof(tok);
            if (j < nLnLCols && nCols == lnLColIdx[j])
                {
                data->siteLnL[i * nLnLCols + j] = atof(tok);
                j++;
                }
            nCols++;
            tok = strtok(NULL, "\t\n\r");
            }
        i++;
        if (i >= data->nSamples)
            break;
        }

    free(line);
    fclose(fp);

    YvyraPrint ("%s   Read %d post-burnin samples with %d per-site lnL columns from '%s'\n",
                  spacer, data->nSamples, data->nChars, actualFile);

    /* If we read from .slk, totalLnL is missing — read it from the .p file */
    if (actualFile != filename && totalLnLCol < 0)
        {
        FILE *fpP = fopen(filename, "r");
        if (fpP)
            {
            char *pline = (char *) SafeMalloc(DIAG_MAX_LINE);
            if (pline)
                {
                int pTotalCol = -1, pnCols;

                /* Skip comments, find header */
                while (fgets(pline, DIAG_MAX_LINE, fpP))
                    if (pline[0] != '[' && pline[0] != '\n' && pline[0] != '\r')
                        break;

                /* Find lnLike column */
                pnCols = 0;
                tok = strtok(pline, "\t\n\r");
                while (tok)
                    {
                    if (strcmp(tok, "lnLike") == 0)
                        pTotalCol = pnCols;
                    pnCols++;
                    tok = strtok(NULL, "\t\n\r");
                    }

                if (pTotalCol >= 0)
                    {
                    lineCount = 0;
                    i = 0;
                    nBurnin = (int)(relBurnin * nTotalLines);
                    while (fgets(pline, DIAG_MAX_LINE, fpP))
                        {
                        if (pline[0] == '\n' || pline[0] == '\r' || pline[0] == '[')
                            continue;
                        lineCount++;
                        if (lineCount <= nBurnin)
                            continue;
                        pnCols = 0;
                        tok = strtok(pline, "\t\n\r");
                        while (tok)
                            {
                            if (pnCols == pTotalCol && i < data->nSamples)
                                data->totalLnL[i] = atof(tok);
                            pnCols++;
                            tok = strtok(NULL, "\t\n\r");
                            }
                        i++;
                        if (i >= data->nSamples)
                            break;
                        }
                    YvyraPrint ("%s   Read total lnL from '%s'\n", spacer, filename);
                    }
                free(pline);
                }
            fclose(fpP);
            }
        }
    return (NO_ERROR);
}


void FreeSiteLnL (SiteLnLData *data)
{
    int i;
    if (data->siteLnL) free(data->siteLnL);
    if (data->totalLnL) free(data->totalLnL);
    if (data->charNames)
        {
        for (i = 0; i < data->nChars; i++)
            if (data->charNames[i]) free(data->charNames[i]);
        free(data->charNames);
        }
    data->siteLnL = NULL;
    data->totalLnL = NULL;
    data->charNames = NULL;
    data->nSamples = 0;
    data->nChars = 0;
}


/*--------------------------------------------------------------
|   ReadTreeSamples: Read Newick trees from .t file
---------------------------------------------------------------*/
int ReadTreeSamples (char *filename, YFlt relBurnin, TreeData *data)
{
    FILE    *fp;
    char    *line, *p;
    int     i, nTotalTrees, nBurnin, treeCount;
    int     inTranslate = NO;

    data->nSamples = 0;
    data->nTaxa = 0;
    data->taxaNames = NULL;
    data->newickTrees = NULL;

    fp = fopen(filename, "r");
    if (!fp)
        {
        YvyraPrint ("%s   Could not open file '%s'\n", spacer, filename);
        return (ERROR);
        }

    line = (char *) SafeMalloc(DIAG_MAX_LINE);
    if (!line) { fclose(fp); return (ERROR); }

    /* Parse translate block for taxon names */
    data->taxaNames = (char **) SafeCalloc(100, sizeof(char *));
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        /* Trim whitespace */
        p = line;
        while (*p == ' ' || *p == '\t') p++;

        if (strncmp(p, "translate", 9) == 0)
            { inTranslate = YES; continue; }

        if (inTranslate == YES)
            {
            if (*p == ';')
                { inTranslate = NO; continue; }
            /* Parse "number name[,;]" */
            int idx;
            char name[200];
            if (sscanf(p, "%d %199s", &idx, name) == 2)
                {
                /* Remove trailing comma or semicolon */
                int len = strlen(name);
                if (len > 0 && (name[len-1] == ',' || name[len-1] == ';'))
                    name[len-1] = '\0';
                if (idx > 0 && idx <= 100)
                    {
                    data->taxaNames[idx-1] = (char *) SafeMalloc(strlen(name) + 1);
                    strcpy(data->taxaNames[idx-1], name);
                    if (idx > data->nTaxa)
                        data->nTaxa = idx;
                    }
                }
            }
        }

    /* Count trees */
    rewind(fp);
    nTotalTrees = 0;
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (strncmp(p, "tree ", 5) == 0 || strncmp(p, "tree\t", 5) == 0)
            nTotalTrees++;
        }

    nBurnin = (int)(relBurnin * nTotalTrees);
    data->nSamples = nTotalTrees - nBurnin;

    if (data->nSamples <= 0)
        {
        YvyraPrint ("%s   No trees after burnin in '%s'\n", spacer, filename);
        free(line); fclose(fp);
        return (ERROR);
        }

    /* Allocate tree strings */
    data->newickTrees = (char **) SafeCalloc(data->nSamples, sizeof(char *));

    /* Re-read trees */
    rewind(fp);
    treeCount = 0;
    i = 0;
    while (fgets(line, DIAG_MAX_LINE, fp))
        {
        p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (strncmp(p, "tree ", 5) != 0 && strncmp(p, "tree\t", 5) != 0)
            continue;
        treeCount++;
        if (treeCount <= nBurnin)
            continue;

        /* Extract Newick string (after "= [&...] " or "= ") */
        char *eq = strchr(p, '=');
        if (eq)
            {
            eq++;
            while (*eq == ' ') eq++;
            if (*eq == '[')
                {
                eq = strchr(eq, ']');
                if (eq) eq++;
                while (eq && *eq == ' ') eq++;
                }
            /* Remove trailing whitespace/newline */
            int len = strlen(eq);
            while (len > 0 && (eq[len-1] == '\n' || eq[len-1] == '\r' || eq[len-1] == ' '))
                len--;
            data->newickTrees[i] = (char *) SafeMalloc(len + 1);
            strncpy(data->newickTrees[i], eq, len);
            data->newickTrees[i][len] = '\0';
            i++;
            }
        if (i >= data->nSamples)
            break;
        }

    free(line);
    fclose(fp);

    YvyraPrint ("%s   Read %d post-burnin trees (%d taxa) from '%s'\n",
                  spacer, data->nSamples, data->nTaxa, filename);
    return (NO_ERROR);
}


void FreeTreeSamples (TreeData *data)
{
    int i;
    if (data->taxaNames)
        {
        for (i = 0; i < data->nTaxa; i++)
            if (data->taxaNames[i]) free(data->taxaNames[i]);
        free(data->taxaNames);
        }
    if (data->newickTrees)
        {
        for (i = 0; i < data->nSamples; i++)
            if (data->newickTrees[i]) free(data->newickTrees[i]);
        free(data->newickTrees);
        }
    data->taxaNames = NULL;
    data->newickTrees = NULL;
    data->nSamples = 0;
    data->nTaxa = 0;
}


/*--------------------------------------------------------------
|   TreeHasClade: Check if a Newick string contains a clade
|   Uses taxon numbers from the translate block.
|   Returns YES if all cladeTaxa form a monophyletic group.
---------------------------------------------------------------*/
int TreeHasClade (char *newick, char **cladeTaxa, int nCladeTaxa,
                  char **allTaxa, int nTaxa)
{
    /* Simple approach: extract all bipartitions from the Newick string
       and check if any bipartition matches the clade.

       For small trees (5-50 taxa), a recursive descent suffices.
       We parse the Newick, and for each internal node, collect the
       set of descendant taxa. If any descendant set exactly equals
       the clade set, the clade is monophyletic. */

    int     i, j;
    BitsLong cladeBits;         /* bit set for clade taxa (up to 64 taxa) */
    BitsLong stack[200];        /* stack of taxon bit sets */
    int     stackTop;
    char    *p;

    if (nCladeTaxa <= 1 || nCladeTaxa >= nTaxa)
        return YES;  /* trivial clades */
    if (nTaxa > 64)
        return NO;   /* bit set overflow guard */

    /* Build clade bit set */
    cladeBits = 0;
    for (i = 0; i < nCladeTaxa; i++)
        {
        for (j = 0; j < nTaxa; j++)
            {
            if (allTaxa[j] && strcmp(allTaxa[j], cladeTaxa[i]) == 0)
                {
                cladeBits |= ((BitsLong)1 << j);
                break;
                }
            }
        }

    /* Parse Newick: track descendant taxa using bit sets */
    stackTop = 0;
    memset(stack, 0, sizeof(stack));

    p = newick;
    while (*p)
        {
        if (*p == '(')
            {
            stackTop++;
            if (stackTop >= 200) return NO;
            stack[stackTop] = 0;
            }
        else if (*p == ')')
            {
            if (stackTop > 0)
                {
                /* Check if current descendant set matches clade */
                if (stack[stackTop] == cladeBits)
                    return YES;

                /* Merge into parent */
                stack[stackTop-1] |= stack[stackTop];
                stackTop--;
                }
            }
        else if (*p >= '0' && *p <= '9')
            {
            if (p > newick && *(p-1) == ':')
                {
                while (*p && *p != ',' && *p != ')' && *p != '(')
                    p++;
                p--;
                }
            else
                {
                int taxNum = 0;
                while (*p >= '0' && *p <= '9')
                    {
                    taxNum = taxNum * 10 + (*p - '0');
                    p++;
                    }
                p--;
                if (taxNum >= 1 && taxNum <= nTaxa)
                    stack[stackTop] |= ((BitsLong)1 << (taxNum - 1));
                }
            }
        p++;
        }

    return NO;
}


/*--------------------------------------------------------------
|   FitchParsimonyScore: Compute Fitch parsimony score for one
|   character on a Newick tree with numbered tips.
|
|   tipStates[i] is a bit set of observed states for taxon i+1
|   (0-based array, 1-based taxon numbers in Newick).
|   Missing data (tipStates[i]==0) means all states possible.
|   Returns score (number of changes) via *score, NO_ERROR/ERROR.
---------------------------------------------------------------*/
int FitchParsimonyScore (char *newick, int nTaxa, int *tipStates, int *score)
{
    /* Stack-based Fitch on a Newick string.
       '(' pushes a sentinel (-1).
       Taxon numbers push their state bit sets.
       ')' pops back to the sentinel, merges children via Fitch, pushes result. */
    int     stack[200];
    int     top = -1;
    int     changes = 0;
    char    *p;
    int     allStatesMask;
    int     SENTINEL = -1;

    /* Build mask for "all states" (used for missing data) */
    allStatesMask = 0;
    {
        int i;
        for (i = 0; i < nTaxa; i++)
            allStatesMask |= tipStates[i];
        if (allStatesMask == 0)
            allStatesMask = 0x3; /* at least binary */
    }

    p = newick;
    while (*p)
        {
        if (*p == '(')
            {
            top++;
            if (top >= 200) return (ERROR);
            stack[top] = SENTINEL;
            }
        else if (*p == ')')
            {
            /* Pop children back to sentinel, merge via Fitch */
            int merged = 0;
            int nChildren = 0;
            while (top >= 0 && stack[top] != SENTINEL)
                {
                int child = stack[top--];
                if (nChildren == 0)
                    merged = child;
                else
                    {
                    int intersect = merged & child;
                    if (intersect != 0)
                        merged = intersect;
                    else
                        {
                        merged = merged | child;
                        changes++;
                        }
                    }
                nChildren++;
                }
            /* Pop the sentinel itself */
            if (top >= 0 && stack[top] == SENTINEL)
                top--;
            /* Push the merged result */
            top++;
            stack[top] = merged;
            }
        else if (*p >= '0' && *p <= '9')
            {
            /* Skip if this is a branch length (preceded by ':') */
            if (p > newick && *(p-1) == ':')
                {
                while (*p && *p != ',' && *p != ')' && *p != '(')
                    p++;
                p--;
                }
            else
                {
                /* Parse taxon number */
                int taxNum = 0;
                while (*p >= '0' && *p <= '9')
                    {
                    taxNum = taxNum * 10 + (*p - '0');
                    p++;
                    }
                p--;
                /* Push taxon state set */
                top++;
                if (top >= 200) return (ERROR);
                if (taxNum >= 1 && taxNum <= nTaxa && tipStates[taxNum-1] != 0)
                    stack[top] = tipStates[taxNum-1];
                else
                    stack[top] = allStatesMask;  /* missing data */
                }
            }
        /* commas, colons, semicolons, letters — skip */
        p++;
        }

    /* Handle unrooted case: merge remaining entries */
    while (top >= 1)
        {
        int child2 = stack[top--];
        int child1 = stack[top];
        int intersect = child1 & child2;
        if (intersect != 0)
            stack[top] = intersect;
        else
            {
            stack[top] = child1 | child2;
            changes++;
            }
        }

    *score = changes;
    return (NO_ERROR);
}


/*--------------------------------------------------------------
|   ShannonEntropy: Compute entropy in bits
---------------------------------------------------------------*/
YFlt ShannonEntropy (YFlt *probs, int n)
{
    int     i;
    YFlt  h = 0.0;

    for (i = 0; i < n; i++)
        {
        if (probs[i] > 0.0)
            h -= probs[i] * log2(probs[i]);
        }
    return h;
}


/*--------------------------------------------------------------
|   ParseNewickTree: Parse a Newick string into a SimpleTree.
|   Tips are identified by their taxon numbers (1-based in Newick).
|   Returns NO_ERROR on success.
---------------------------------------------------------------*/
int ParseNewickTree (char *newick, int nTaxa, SimpleTree *tree)
{
    int     nodeCount, tipCount, intCount;
    int     stackIdx[200], stackTop;
    int     curParent;
    char    *p;

    tree->nTips = nTaxa;
    tree->nIntNodes = nTaxa;    /* one extra for resolving trifurcations */
    tree->nNodes = 2 * nTaxa;
    tree->left    = (int *)     SafeCalloc(tree->nNodes, sizeof(int));
    tree->right   = (int *)     SafeCalloc(tree->nNodes, sizeof(int));
    tree->parent  = (int *)     SafeCalloc(tree->nNodes, sizeof(int));
    tree->taxon   = (int *)     SafeCalloc(tree->nNodes, sizeof(int));
    tree->descTips = (BitsLong *) SafeCalloc(tree->nNodes, sizeof(BitsLong));

    if (!tree->left || !tree->right || !tree->parent || !tree->taxon || !tree->descTips)
        return (ERROR);

    {
        int i;
        for (i = 0; i < tree->nNodes; i++)
            {
            tree->left[i] = tree->right[i] = tree->parent[i] = -1;
            tree->taxon[i] = -1;
            tree->descTips[i] = 0;
            }
    }

    /* Nodes 0..nTaxa-1 are reserved for tips (filled when encountered).
       Nodes nTaxa..2*nTaxa-2 are internal nodes (allocated on '(').
       For unrooted trees with root trifurcation, we resolve by
       creating a virtual internal node for the first two root children. */
    tipCount = 0;
    intCount = 0;
    stackTop = -1;
    nodeCount = 0;  /* suppress warning */

    /* Start with a virtual root internal node */
    curParent = nTaxa + intCount++;

    p = newick;
    while (*p)
        {
        if (*p == '(')
            {
            /* Push current parent, create new internal node */
            stackTop++;
            stackIdx[stackTop] = curParent;
            if (intCount < tree->nIntNodes)
                {
                int newNode = nTaxa + intCount++;
                /* Attach to parent */
                if (tree->left[curParent] == -1)
                    tree->left[curParent] = newNode;
                else if (tree->right[curParent] == -1)
                    tree->right[curParent] = newNode;
                else
                    {
                    /* Trifurcation: create extra internal node to hold
                       existing left+right, then attach new child to parent.right */
                    if (intCount < tree->nIntNodes)
                        {
                        int extraNode = nTaxa + intCount++;
                        tree->left[extraNode] = tree->left[curParent];
                        tree->right[extraNode] = tree->right[curParent];
                        tree->parent[extraNode] = curParent;
                        tree->parent[tree->left[curParent]] = extraNode;
                        tree->parent[tree->right[curParent]] = extraNode;
                        tree->left[curParent] = extraNode;
                        tree->right[curParent] = newNode;
                        }
                    }
                tree->parent[newNode] = curParent;
                curParent = newNode;
                }
            }
        else if (*p == ')')
            {
            /* Return to parent */
            if (stackTop >= 0)
                curParent = stackIdx[stackTop--];
            }
        else if (*p >= '0' && *p <= '9')
            {
            /* Skip branch lengths */
            if (p > newick && *(p-1) == ':')
                {
                while (*p && *p != ',' && *p != ')' && *p != '(' && *p != ';')
                    p++;
                p--;
                }
            else
                {
                /* Parse taxon number */
                int taxNum = 0;
                while (*p >= '0' && *p <= '9')
                    { taxNum = taxNum * 10 + (*p - '0'); p++; }
                p--;

                if (taxNum >= 1 && taxNum <= nTaxa && tipCount < nTaxa)
                    {
                    int tipNode = tipCount++;
                    tree->taxon[tipNode] = taxNum;
                    tree->descTips[tipNode] = (BitsLong)1 << (taxNum - 1);
                    /* Attach to parent */
                    if (tree->left[curParent] == -1)
                        tree->left[curParent] = tipNode;
                    else if (tree->right[curParent] == -1)
                        tree->right[curParent] = tipNode;
                    else
                        {
                        /* Third child at this node (root trifurcation with tip):
                           create extra node for existing children */
                        if (intCount < tree->nIntNodes)
                            {
                            int extraNode = nTaxa + intCount++;
                            tree->left[extraNode] = tree->left[curParent];
                            tree->right[extraNode] = tree->right[curParent];
                            tree->parent[extraNode] = curParent;
                            tree->parent[tree->left[curParent]] = extraNode;
                            tree->parent[tree->right[curParent]] = extraNode;
                            tree->left[curParent] = extraNode;
                            tree->right[curParent] = tipNode;
                            }
                        }
                    tree->parent[tipNode] = curParent;
                    }
                }
            }
        p++;
        }

    /* Update actual internal node count */
    tree->nIntNodes = intCount;
    tree->nNodes = nTaxa + intCount;

    /* Compute descTips for internal nodes (bottom-up) */
    {
        int i, changed;
        do {
            changed = 0;
            for (i = nTaxa; i < tree->nNodes; i++)
                {
                BitsLong bits = 0;
                if (tree->left[i] >= 0)  bits |= tree->descTips[tree->left[i]];
                if (tree->right[i] >= 0) bits |= tree->descTips[tree->right[i]];
                if (bits != tree->descTips[i])
                    { tree->descTips[i] = bits; changed = 1; }
                }
        } while (changed);
    }

    return (NO_ERROR);
}


void FreeSimpleTree (SimpleTree *tree)
{
    if (tree->left)     free(tree->left);
    if (tree->right)    free(tree->right);
    if (tree->parent)   free(tree->parent);
    if (tree->taxon)    free(tree->taxon);
    if (tree->descTips) free(tree->descTips);
    tree->left = tree->right = tree->parent = tree->taxon = NULL;
    tree->descTips = NULL;
}


/*--------------------------------------------------------------
|   FitchDownpassSets: Run Fitch downpass on a SimpleTree and
|   store the state bit set at each internal node.
|
|   tipStates[i]: bit set for tip with taxon number i+1 (0-based).
|   nodeStateSets[i]: output bit set for node i (internal nodes
|   are at indices nTaxa..nNodes-1).
---------------------------------------------------------------*/
int FitchDownpassSets (SimpleTree *tree, int *tipStates, int *nodeStateSets)
{
    int     i, allMask;

    /* Compute allMask for missing data */
    allMask = 0;
    for (i = 0; i < tree->nTips; i++)
        allMask |= tipStates[i];
    if (allMask == 0)
        allMask = 0x3;

    /* Initialize tips */
    for (i = 0; i < tree->nTips; i++)
        {
        if (tree->taxon[i] >= 1 && tree->taxon[i] <= tree->nTips)
            nodeStateSets[i] = tipStates[tree->taxon[i] - 1];
        else
            nodeStateSets[i] = allMask;
        if (nodeStateSets[i] == 0)
            nodeStateSets[i] = allMask;
        }

    /* Bottom-up pass for internal nodes */
    for (i = tree->nTips; i < tree->nNodes; i++)
        {
        int lft = tree->left[i];
        int rht = tree->right[i];

        if (lft < 0 || rht < 0)
            {
            /* Shouldn't happen in a well-formed tree */
            nodeStateSets[i] = allMask;
            continue;
            }

        int intersect = nodeStateSets[lft] & nodeStateSets[rht];
        if (intersect != 0)
            nodeStateSets[i] = intersect;
        else
            nodeStateSets[i] = nodeStateSets[lft] | nodeStateSets[rht];
        }

    return (NO_ERROR);
}


/*--------------------------------------------------------------
|   CompareYFlt: qsort comparator for YFlt
---------------------------------------------------------------*/
static int CompareYFlt (const void *a, const void *b)
{
    YFlt va = *(const YFlt *)a;
    YFlt vb = *(const YFlt *)b;
    if (va < vb) return -1;
    if (va > vb) return  1;
    return 0;
}


/*--------------------------------------------------------------
|   BuildNodeLabel: Create a human-readable label for an internal
|   node based on its descendant tip set.
---------------------------------------------------------------*/
static void BuildNodeLabel (BitsLong descBits, int nTaxa, char **taxaNames,
                            char *buf, int bufSize)
{
    int     i, count;
    BitsLong allBits;

    /* Check if this is the root (all taxa) */
    allBits = 0;
    for (i = 0; i < nTaxa; i++)
        allBits |= ((BitsLong)1 << i);

    if (descBits == allBits)
        {
        snprintf(buf, bufSize, "(root)");
        return;
        }

    /* Count descendants */
    count = 0;
    for (i = 0; i < nTaxa; i++)
        if (descBits & ((BitsLong)1 << i))
            count++;

    /* If small enough, list all names */
    if (count <= 3)
        {
        int pos = 0;
        int first = YES;
        pos += snprintf(buf + pos, bufSize - pos, "(");
        for (i = 0; i < nTaxa && pos < bufSize - 2; i++)
            {
            if (descBits & ((BitsLong)1 << i))
                {
                if (first == NO)
                    pos += snprintf(buf + pos, bufSize - pos, ",");
                pos += snprintf(buf + pos, bufSize - pos, "%s",
                                taxaNames[i] ? taxaNames[i] : "?");
                first = NO;
                }
            }
        snprintf(buf + pos, bufSize - pos, ")");
        }
    else
        {
        /* Show first two + count */
        int pos = 0, shown = 0;
        pos += snprintf(buf + pos, bufSize - pos, "(");
        for (i = 0; i < nTaxa && shown < 2 && pos < bufSize - 2; i++)
            {
            if (descBits & ((BitsLong)1 << i))
                {
                if (shown > 0)
                    pos += snprintf(buf + pos, bufSize - pos, ",");
                pos += snprintf(buf + pos, bufSize - pos, "%s",
                                taxaNames[i] ? taxaNames[i] : "?");
                shown++;
                }
            }
        snprintf(buf + pos, bufSize - pos, "+%d)", count - 2);
        }
}


/*--------------------------------------------------------------
|   DoAsrentropy: Ancestral state entropy from posterior trees.
|
|   Parsimony-based ASR: for each posterior tree, runs Fitch
|   downpass to get state sets at internal nodes. Aggregates
|   across samples with equal weight within state sets.
|   Reports entropy, confidence, and credible intervals.
---------------------------------------------------------------*/
int DoAsrentropy (void)
{
    TreeData        trData;
    SiteLnLData     slData = {0, 0, NULL, NULL, NULL};
    SimpleTree      refTree, sampleTree;
    char            tFilename[300], pFilename[300];
    YFlt          relBurnin = 0.25;
    int             s, j, i, k, n;
    int             nIntNodes, maxStates;
    int             *tipStates, *nodeStateSets;
    int             mapIdx;

    /* --- Read tree and likelihood data --- */
    if (defMatrix == NO || matrix == NULL)
        {
        YvyraPrint ("%s   A matrix must be defined before running asrentropy\n", spacer);
        return (ERROR);
        }

    if (chainParams.numRuns == 1)
        sprintf(tFilename, "%s.t", chainParams.chainFileName);
    else
        sprintf(tFilename, "%s.run1.t", chainParams.chainFileName);
    if (chainParams.numRuns == 1)
        sprintf(pFilename, "%s.p", chainParams.chainFileName);
    else
        sprintf(pFilename, "%s.run1.p", chainParams.chainFileName);

    if (ReadTreeSamples(tFilename, relBurnin, &trData) == ERROR)
        return (ERROR);

    if (trData.nSamples < 2)
        {
        YvyraPrint ("%s   Need at least 2 post-burnin tree samples\n", spacer);
        FreeTreeSamples(&trData);
        return (ERROR);
        }

    /* Try to read site lnL for MAP tree identification; if unavailable,
       read just totalLnL from .p file; if that also fails, use last sample */
    {
        int haveLnL = NO;
        if (ReadSiteLnL(pFilename, relBurnin, &slData) == NO_ERROR)
            haveLnL = YES;

        mapIdx = trData.nSamples - 1;  /* default: last sample */
        if (haveLnL == YES && slData.totalLnL != NULL)
            {
            mapIdx = 0;
            int nComp = (trData.nSamples < slData.nSamples) ? trData.nSamples : slData.nSamples;
            for (s = 1; s < nComp; s++)
                if (slData.totalLnL[s] > slData.totalLnL[mapIdx])
                    mapIdx = s;
            }
        else
            {
            /* No site lnL available — try reading totalLnL from .p file directly */
            FILE *fpP = fopen(pFilename, "r");
            if (fpP)
                {
                char *pline = (char *) SafeMalloc(DIAG_MAX_LINE);
                if (pline)
                    {
                    int pTotalCol = -1, nCols, lineCount = 0, idx = 0;
                    int nTotal = 0, nBurnin;
                    YFlt *totalLnL;
                    char *tok;

                    /* Find header */
                    while (fgets(pline, DIAG_MAX_LINE, fpP))
                        if (pline[0] != '[' && pline[0] != '\n' && pline[0] != '\r')
                            break;
                    nCols = 0;
                    tok = strtok(pline, "\t\n\r");
                    while (tok)
                        {
                        if (strcmp(tok, "lnLike") == 0)
                            pTotalCol = nCols;
                        nCols++;
                        tok = strtok(NULL, "\t\n\r");
                        }

                    /* Count data lines */
                    while (fgets(pline, DIAG_MAX_LINE, fpP))
                        if (pline[0] != '\n' && pline[0] != '\r' && pline[0] != '[')
                            nTotal++;
                    nBurnin = (int)(relBurnin * nTotal);

                    if (pTotalCol >= 0 && nTotal - nBurnin > 0)
                        {
                        totalLnL = (YFlt *) SafeCalloc(nTotal - nBurnin, sizeof(YFlt));
                        rewind(fpP);
                        while (fgets(pline, DIAG_MAX_LINE, fpP))
                            if (pline[0] != '[' && pline[0] != '\n' && pline[0] != '\r')
                                break;
                        lineCount = 0; idx = 0;
                        while (fgets(pline, DIAG_MAX_LINE, fpP))
                            {
                            if (pline[0] == '\n' || pline[0] == '\r' || pline[0] == '[')
                                continue;
                            lineCount++;
                            if (lineCount <= nBurnin) continue;
                            nCols = 0;
                            tok = strtok(pline, "\t\n\r");
                            while (tok)
                                {
                                if (nCols == pTotalCol && idx < nTotal - nBurnin)
                                    totalLnL[idx] = atof(tok);
                                nCols++;
                                tok = strtok(NULL, "\t\n\r");
                                }
                            idx++;
                            if (idx >= nTotal - nBurnin) break;
                            }
                        mapIdx = 0;
                        for (s = 1; s < idx && s < trData.nSamples; s++)
                            if (totalLnL[s] > totalLnL[mapIdx])
                                mapIdx = s;
                        free(totalLnL);
                        }
                    free(pline);
                    }
                fclose(fpP);
                }
            }
    }

    /* Parse MAP tree as reference */
    if (ParseNewickTree(trData.newickTrees[mapIdx], trData.nTaxa, &refTree) == ERROR)
        {
        YvyraPrint ("%s   Could not parse MAP tree\n", spacer);
        FreeSiteLnL(&slData); FreeTreeSamples(&trData);
        return (ERROR);
        }

    nIntNodes = refTree.nIntNodes;
    maxStates = MAX_STD_STATES;  /* upper bound */

    YvyraPrint ("\n%s   ANCESTRAL STATE ENTROPY (parsimony-based)\n", spacer);
    YvyraPrint ("%s   %d post-burnin trees, %d taxa, %d internal nodes\n",
                  spacer, trData.nSamples, trData.nTaxa, nIntNodes);
    YvyraPrint ("%s   Reference: MAP tree (sample %d, lnL = %.4f)\n",
                  spacer, mapIdx + 1, slData.totalLnL[mapIdx]);

    /* --- Allocate per-sample storage for credible intervals ---
       stateFreq[node * numChar * maxStates + char * maxStates + state]
       is the accumulated weighted count.
       perSampleProbs[node * numChar * maxStates * nSamples + char * maxStates * nSamples + state * nSamples + sample]
       stores per-sample probability for credible intervals. */

    int nSamples = (trData.nSamples < slData.nSamples) ? trData.nSamples : slData.nSamples;
    size_t freqSize = (size_t)nIntNodes * numChar * maxStates;
    YFlt *stateFreq = (YFlt *) SafeCalloc(freqSize, sizeof(YFlt));
    YFlt *perSampleProbs = (YFlt *) SafeCalloc(freqSize * nSamples, sizeof(YFlt));
    int *sampleCount = (int *) SafeCalloc((size_t)nIntNodes, sizeof(int));

    if (!stateFreq || !perSampleProbs || !sampleCount)
        {
        YvyraPrint ("%s   Problem allocating ASR storage\n", spacer);
        FreeSimpleTree(&refTree);
        FreeSiteLnL(&slData); FreeTreeSamples(&trData);
        if (stateFreq) free(stateFreq);
        if (perSampleProbs) free(perSampleProbs);
        if (sampleCount) free(sampleCount);
        return (ERROR);
        }

    tipStates = (int *) SafeCalloc(trData.nTaxa, sizeof(int));
    nodeStateSets = (int *) SafeCalloc(refTree.nNodes, sizeof(int));

    /* --- Process each posterior tree sample --- */
    for (s = 0; s < nSamples; s++)
        {
        memset(&sampleTree, 0, sizeof(SimpleTree));
        if (ParseNewickTree(trData.newickTrees[s], trData.nTaxa, &sampleTree) == ERROR)
            continue;

        /* For each original character */
        for (j = 0; j < numChar; j++)
            {
            int nStates = charInfo[j].numStates;
            int validMask = (1 << nStates) - 1;

            if (charInfo[j].isExcluded == YES || nStates < 2)
                continue;

            /* Build tip states from matrix */
            for (i = 0; i < trData.nTaxa; i++)
                {
                /* Find matrix taxon index for .t file taxon i */
                int matIdx = -1;
                int k2;
                for (k2 = 0; k2 < numTaxa; k2++)
                    {
                    if (trData.taxaNames[i] && taxaNames[k2] &&
                        strcmp(trData.taxaNames[i], taxaNames[k2]) == 0)
                        { matIdx = k2; break; }
                    }
                if (matIdx >= 0)
                    tipStates[i] = matrix[pos(matIdx, j, numChar)] & validMask;
                else
                    tipStates[i] = validMask;
                if (tipStates[i] == 0)
                    tipStates[i] = validMask;
                }

            /* Run Fitch downpass */
            FitchDownpassSets(&sampleTree, tipStates, nodeStateSets);

            /* Map sample internal nodes to reference tree by descTips */
            for (n = sampleTree.nTips; n < sampleTree.nNodes; n++)
                {
                BitsLong sampleDesc = sampleTree.descTips[n];
                int stateSet = nodeStateSets[n];
                int setSize = 0;
                int nStatesInSet;

                /* Count bits in stateSet */
                {
                    int tmp = stateSet;
                    while (tmp) { setSize += tmp & 1; tmp >>= 1; }
                }
                if (setSize == 0) continue;
                nStatesInSet = setSize;

                /* Find matching reference node */
                for (k = refTree.nTips; k < refTree.nNodes; k++)
                    {
                    if (refTree.descTips[k] == sampleDesc)
                        {
                        int refIntIdx = k - refTree.nTips;
                        size_t baseIdx = (size_t)refIntIdx * numChar * maxStates + (size_t)j * maxStates;

                        /* Add weighted state counts */
                        YFlt weight = 1.0 / nStatesInSet;
                        int st;
                        for (st = 0; st < nStates; st++)
                            {
                            if (stateSet & (1 << st))
                                {
                                stateFreq[baseIdx + st] += weight;
                                perSampleProbs[((size_t)baseIdx + st) * nSamples + s] = weight;
                                }
                            }
                        if (refIntIdx < nIntNodes)
                            sampleCount[refIntIdx]++;
                        break;
                        }
                    }
                }
            }

        FreeSimpleTree(&sampleTree);
        }

    /* --- Compute and print results per character --- */
    for (j = 0; j < numChar; j++)
        {
        int nStates = charInfo[j].numStates;
        char charName[64];

        if (charInfo[j].isExcluded == YES || nStates < 2)
            continue;

        if (charLabels != NULL && charLabels[j] != NULL)
            snprintf(charName, sizeof(charName), "%s", charLabels[j]);
        else
            snprintf(charName, sizeof(charName), "char_%d", j + 1);

        YvyraPrint ("\n%s   Character: %s (%d states)\n", spacer, charName, nStates);
        YvyraPrint ("%s   %-24s", spacer, "Node");

        /* Print state headers (named if available) */
        for (k = 0; k < nStates; k++)
            {
            if (charStateNames && charNStateNames && j < numChar &&
                charNStateNames[j] > k && charStateNames[j] && charStateNames[j][k])
                {
                char hdr[64];
                snprintf(hdr, sizeof(hdr), "p(%s)", charStateNames[j][k]);
                YvyraPrint ("  %-7s", hdr);
                }
            else
                YvyraPrint ("  p(%d)  ", k);
            }
        YvyraPrint (" Entropy  Conf   Best  95%%CI_lo  95%%CI_hi\n");

        YvyraPrint ("%s   %-24s", spacer, "----");
        for (k = 0; k < nStates; k++)
            YvyraPrint ("  -----  ");
        YvyraPrint (" -------  ----   ----  --------  --------\n");

        /* For each reference internal node */
        for (n = 0; n < nIntNodes; n++)
            {
            int refNodeIdx = refTree.nTips + n;
            size_t baseIdx = (size_t)n * numChar * maxStates + (size_t)j * maxStates;
            YFlt probs[MAX_STD_STATES];
            YFlt totalWeight = 0.0;
            YFlt entropy, confidence, maxEntropy;
            int bestState = 0;
            YFlt bestProb = 0.0;
            char nodeLabel[128];

            /* Normalize frequencies to probabilities */
            for (k = 0; k < nStates; k++)
                totalWeight += stateFreq[baseIdx + k];

            if (totalWeight <= 0.0)
                continue;

            for (k = 0; k < nStates; k++)
                probs[k] = stateFreq[baseIdx + k] / totalWeight;

            /* Find best state */
            for (k = 0; k < nStates; k++)
                if (probs[k] > bestProb)
                    { bestProb = probs[k]; bestState = k; }

            /* Entropy and confidence */
            entropy = ShannonEntropy(probs, nStates);
            maxEntropy = log2((double)nStates);
            confidence = (maxEntropy > 0.0) ? 1.0 - entropy / maxEntropy : 1.0;

            /* Credible intervals for best state */
            YFlt ciLo = 0.0, ciHi = 1.0;
            {
                size_t probBaseIdx = ((size_t)baseIdx + bestState) * nSamples;
                YFlt *sortBuf = (YFlt *) SafeCalloc(nSamples, sizeof(YFlt));
                if (sortBuf)
                    {
                    int nValid = 0;
                    for (s = 0; s < nSamples; s++)
                        {
                        /* For samples where this node exists, use recorded prob;
                           for samples where node doesn't exist, use 0 */
                        sortBuf[nValid++] = perSampleProbs[probBaseIdx + s];
                        }
                    qsort(sortBuf, nValid, sizeof(YFlt), CompareYFlt);
                    int loIdx = (int)(0.025 * nValid);
                    int hiIdx = (int)(0.975 * nValid);
                    if (hiIdx >= nValid) hiIdx = nValid - 1;
                    ciLo = sortBuf[loIdx];
                    ciHi = sortBuf[hiIdx];
                    free(sortBuf);
                    }
            }

            /* Build node label */
            BuildNodeLabel(refTree.descTips[refNodeIdx], trData.nTaxa,
                           trData.taxaNames, nodeLabel, sizeof(nodeLabel));

            YvyraPrint ("%s   %-24s", spacer, nodeLabel);
            for (k = 0; k < nStates; k++)
                YvyraPrint ("  %5.3f  ", probs[k]);

            /* Show best state as name if available */
            if (charStateNames && charNStateNames && j < numChar &&
                charNStateNames[j] > bestState && charStateNames[j] && charStateNames[j][bestState])
                YvyraPrint (" %7.4f  %4.2f   %-12s  %8.4f  %8.4f\n",
                              entropy, confidence, charStateNames[j][bestState], ciLo, ciHi);
            else
                YvyraPrint (" %7.4f  %4.2f   %4d  %8.4f  %8.4f\n",
                              entropy, confidence, bestState, ciLo, ciHi);
            }
        }

    /* Clean up */
    free(stateFreq);
    free(perSampleProbs);
    free(sampleCount);
    free(tipStates);
    free(nodeStateSets);
    FreeSimpleTree(&refTree);
    FreeSiteLnL(&slData);
    FreeTreeSamples(&trData);

    return (NO_ERROR);
}


/*--------------------------------------------------------------
|   DoChardiag: Character contribution analysis
---------------------------------------------------------------*/
int DoChardiag (void)
{
    SiteLnLData     slData;
    TreeData        trData;
    int             i, j, s;
    YFlt          *meanLL, *proportion, *support;
    YFlt          totalAbs, sum;
    int             *hasClade;
    int             nWith, nWithout;
    char            pFilename[300], tFilename[300];
    YFlt          relBurnin = 0.25;

    if (diagNCladeTaxa < 2)
        {
        YvyraPrint ("%s   Please specify a clade with at least 2 taxa\n", spacer);
        YvyraPrint ("%s   Example: chardiag clade={Luwian,Lycian}\n", spacer);
        return (ERROR);
        }

    /* Build filenames from chainParams */
    if (chainParams.numRuns == 1)
        sprintf(pFilename, "%s.p", chainParams.chainFileName);
    else
        sprintf(pFilename, "%s.run1.p", chainParams.chainFileName);
    if (chainParams.numRuns == 1)
        sprintf(tFilename, "%s.t", chainParams.chainFileName);
    else
        sprintf(tFilename, "%s.run1.t", chainParams.chainFileName);

    /* Read data */
    if (ReadSiteLnL(pFilename, relBurnin, &slData) == ERROR)
        return (ERROR);
    if (ReadTreeSamples(tFilename, relBurnin, &trData) == ERROR)
        { FreeSiteLnL(&slData); return (ERROR); }

    /* Replace .p/.slk column names with charLabels if available.
       Column lnL(d.c) has c as the 1-based index among non-dummy chars
       in division d. The absolute compressed char index is
       compCharStart + numDummyChars + (c-1). */
    if (charLabels != NULL && origChar != NULL)
        {
        for (j = 0; j < slData.nChars; j++)
            {
            int compIdx = -1, divIdx = -1;
            char *dot = strrchr(slData.charNames[j], '.');
            char *paren = strchr(slData.charNames[j], '(');
            if (dot)
                sscanf(dot + 1, "%d", &compIdx);
            if (paren)
                sscanf(paren + 1, "%d", &divIdx);

            /* Convert from non-dummy 1-based index to absolute compressed index */
            int absCompIdx = -1;
            if (compIdx >= 1 && divIdx >= 1 && divIdx <= numCurrentDivisions)
                {
                ModelInfo *mi = &modelSettings[divIdx - 1];
                absCompIdx = mi->compCharStart + mi->numDummyChars + (compIdx - 1);
                }
            else if (compIdx >= 1)
                {
                /* Fallback: assume single division, offset by first division's dummies */
                absCompIdx = modelSettings[0].compCharStart + modelSettings[0].numDummyChars + (compIdx - 1);
                }

            if (absCompIdx >= 0 && absCompIdx < numCompressedChars)
                {
                int origIdx = origChar[absCompIdx];
                if (origIdx >= 0 && origIdx < numChar && charLabels[origIdx] != NULL)
                    {
                    free(slData.charNames[j]);
                    slData.charNames[j] = (char *) SafeMalloc(strlen(charLabels[origIdx]) + 1);
                    if (slData.charNames[j])
                        strcpy(slData.charNames[j], charLabels[origIdx]);
                    }
                else if (origIdx < 0)
                    {
                    /* Ascertainment correction dummy pattern */
                    char buf[32];
                    sprintf(buf, "(asc_corr_%d)", compIdx);
                    free(slData.charNames[j]);
                    slData.charNames[j] = (char *) SafeMalloc(strlen(buf) + 1);
                    if (slData.charNames[j])
                        strcpy(slData.charNames[j], buf);
                    }
                }
            }
        }

    /* Ensure sample counts match */
    if (slData.nSamples != trData.nSamples)
        {
        YvyraPrint ("%s   WARNING: .p has %d samples but .t has %d; using minimum\n",
                      spacer, slData.nSamples, trData.nSamples);
        if (trData.nSamples < slData.nSamples)
            slData.nSamples = trData.nSamples;
        else
            trData.nSamples = slData.nSamples;
        }

    /* --- 1. Likelihood contributions --- */
    meanLL = (YFlt *) SafeCalloc(slData.nChars, sizeof(YFlt));
    proportion = (YFlt *) SafeCalloc(slData.nChars, sizeof(YFlt));

    for (j = 0; j < slData.nChars; j++)
        {
        sum = 0.0;
        for (s = 0; s < slData.nSamples; s++)
            sum += slData.siteLnL[s * slData.nChars + j];
        meanLL[j] = sum / slData.nSamples;
        }

    totalAbs = 0.0;
    for (j = 0; j < slData.nChars; j++)
        totalAbs += fabs(meanLL[j]);

    for (j = 0; j < slData.nChars; j++)
        proportion[j] = (totalAbs > 0.0) ? fabs(meanLL[j]) / totalAbs : 0.0;

    /* Build sort index (descending by proportion) */
    {
        int *sortIdx = (int *) SafeCalloc(slData.nChars, sizeof(int));
        for (j = 0; j < slData.nChars; j++) sortIdx[j] = j;
        for (i = 0; i < slData.nChars - 1; i++)
            for (j = i + 1; j < slData.nChars; j++)
                if (proportion[sortIdx[j]] > proportion[sortIdx[i]])
                    { int tmp = sortIdx[i]; sortIdx[i] = sortIdx[j]; sortIdx[j] = tmp; }

        YvyraPrint ("\n%s   CHARACTER LIKELIHOOD CONTRIBUTIONS\n", spacer);
        YvyraPrint ("%s   %-20s %12s %10s\n", spacer, "Character", "Mean lnL", "Proportion");
        YvyraPrint ("%s   %-20s %12s %10s\n", spacer, "---------", "--------", "----------");
        for (i = 0; i < slData.nChars; i++)
            {
            j = sortIdx[i];
            YvyraPrint ("%s   %-20s %12.4f %9.1f%%\n", spacer,
                          slData.charNames[j], meanLL[j], proportion[j] * 100.0);
            }
        free(sortIdx);
    }

    /* --- 2. Clade-character support --- */
    {
        int maxSamples = (slData.nSamples > trData.nSamples) ? slData.nSamples : trData.nSamples;
        hasClade = (int *) SafeCalloc(maxSamples, sizeof(int));
    }
    support = (YFlt *) SafeCalloc(slData.nChars, sizeof(YFlt));

    nWith = nWithout = 0;
    for (s = 0; s < trData.nSamples; s++)
        {
        hasClade[s] = TreeHasClade(trData.newickTrees[s], diagCladeTaxa,
                                    diagNCladeTaxa, trData.taxaNames, trData.nTaxa);
        if (hasClade[s] == YES)
            nWith++;
        else
            nWithout++;
        }

    YFlt cladePosterior = (YFlt)nWith / trData.nSamples;
    YvyraPrint ("\n%s   CLADE-CHARACTER SUPPORT for {%s}\n", spacer, diagCladeStr);
    YvyraPrint ("%s   Clade posterior probability: %.4f (%d/%d samples)\n",
                  spacer, cladePosterior, nWith, trData.nSamples);

    if (nWith > 0 && nWithout > 0)
        {
        for (j = 0; j < slData.nChars; j++)
            {
            YFlt sumWith = 0.0, sumWithout = 0.0;
            for (s = 0; s < slData.nSamples; s++)
                {
                if (hasClade[s] == YES)
                    sumWith += slData.siteLnL[s * slData.nChars + j];
                else
                    sumWithout += slData.siteLnL[s * slData.nChars + j];
                }
            support[j] = (sumWith / nWith) - (sumWithout / nWithout);
            }

        /* Sort by absolute support (descending) */
        {
            int *sortIdx = (int *) SafeCalloc(slData.nChars, sizeof(int));
            for (j = 0; j < slData.nChars; j++) sortIdx[j] = j;
            for (i = 0; i < slData.nChars - 1; i++)
                for (j = i + 1; j < slData.nChars; j++)
                    if (fabs(support[sortIdx[j]]) > fabs(support[sortIdx[i]]))
                        { int tmp = sortIdx[i]; sortIdx[i] = sortIdx[j]; sortIdx[j] = tmp; }

            YvyraPrint ("%s   %-20s %12s %10s\n", spacer, "Character", "Support", "Direction");
            YvyraPrint ("%s   %-20s %12s %10s\n", spacer, "---------", "-------", "---------");
            for (i = 0; i < slData.nChars; i++)
                {
                j = sortIdx[i];
                YvyraPrint ("%s   %-20s %12.4f %10s\n", spacer,
                              slData.charNames[j], support[j],
                              support[j] > 0.01 ? "SUPPORTS" : (support[j] < -0.01 ? "CONFLICTS" : "neutral"));
                }
            free(sortIdx);
        }
        }
    else
        {
        YvyraPrint ("%s   Cannot compute support: clade is %s in all samples\n",
                      spacer, nWith > 0 ? "present" : "absent");
        }

    /* --- 3. Weight sensitivity (importance sampling) --- */
    YvyraPrint ("\n%s   WEIGHT SENSITIVITY for {%s}\n", spacer, diagCladeStr);
    YvyraPrint ("%s   Baseline clade probability: %.4f\n", spacer, cladePosterior);

    if (nWith > 0 && nWithout > 0)
        {
        YFlt *origTotalLL = (YFlt *) SafeCalloc(slData.nSamples, sizeof(YFlt));
        YFlt *logRatios = (YFlt *) SafeCalloc(slData.nSamples, sizeof(YFlt));
        YFlt *impWeights = (YFlt *) SafeCalloc(slData.nSamples, sizeof(YFlt));
        YFlt *sensProb = (YFlt *) SafeCalloc(slData.nChars, sizeof(YFlt));
        YFlt *sensDelta = (YFlt *) SafeCalloc(slData.nChars, sizeof(YFlt));

        for (s = 0; s < slData.nSamples; s++)
            origTotalLL[s] = slData.totalLnL[s];

        for (j = 0; j < slData.nChars; j++)
            {
            YFlt maxRatio = -1e300;
            for (s = 0; s < slData.nSamples; s++)
                {
                logRatios[s] = -slData.siteLnL[s * slData.nChars + j];
                if (logRatios[s] > maxRatio)
                    maxRatio = logRatios[s];
                }
            YFlt sumW = 0.0;
            for (s = 0; s < slData.nSamples; s++)
                {
                impWeights[s] = exp(logRatios[s] - maxRatio);
                sumW += impWeights[s];
                }
            for (s = 0; s < slData.nSamples; s++)
                impWeights[s] /= sumW;
            YFlt rwProb = 0.0;
            for (s = 0; s < slData.nSamples; s++)
                rwProb += impWeights[s] * (hasClade[s] == YES ? 1.0 : 0.0);
            sensProb[j] = rwProb;
            sensDelta[j] = cladePosterior - rwProb;
            }

        /* Sort by absolute delta (descending) */
        {
            int *sortIdx = (int *) SafeCalloc(slData.nChars, sizeof(int));
            for (j = 0; j < slData.nChars; j++) sortIdx[j] = j;
            for (i = 0; i < slData.nChars - 1; i++)
                for (j = i + 1; j < slData.nChars; j++)
                    if (fabs(sensDelta[sortIdx[j]]) > fabs(sensDelta[sortIdx[i]]))
                        { int tmp = sortIdx[i]; sortIdx[i] = sortIdx[j]; sortIdx[j] = tmp; }

            YvyraPrint ("%s   %-20s %10s %12s\n", spacer, "Character", "Prob w/o", "Delta");
            YvyraPrint ("%s   %-20s %10s %12s\n", spacer, "---------", "--------", "-----");
            for (i = 0; i < slData.nChars; i++)
                {
                j = sortIdx[i];
                YvyraPrint ("%s   %-20s %9.4f %+11.4f\n", spacer,
                              slData.charNames[j], sensProb[j], sensDelta[j]);
                }
            free(sortIdx);
        }

        free(origTotalLL);
        free(logRatios);
        free(impWeights);
        free(sensProb);
        free(sensDelta);
        }

    /* --- 4. Consistency Index (Fitch parsimony on MAP tree) --- */
    if (defMatrix == YES && matrix != NULL)
        {
        /* Find MAP tree (sample with highest total lnL) */
        int mapIdx = 0;
        for (s = 1; s < slData.nSamples; s++)
            if (slData.totalLnL[s] > slData.totalLnL[mapIdx])
                mapIdx = s;

        YvyraPrint ("\n%s   CONSISTENCY INDEX (Fitch parsimony on MAP tree)\n", spacer);
        YvyraPrint ("%s   MAP tree from sample %d (lnL = %.4f)\n",
                      spacer, mapIdx + 1, slData.totalLnL[mapIdx]);

        /* Build taxon number mapping: .t file taxon i → matrix taxon index.
           The .t translate block has taxaNames in numbered order;
           we need to match these to the global taxaNames array. */
        int *tFileToMatrix = (int *) SafeCalloc(trData.nTaxa, sizeof(int));
        for (i = 0; i < trData.nTaxa; i++)
            {
            tFileToMatrix[i] = -1;
            for (j = 0; j < numTaxa; j++)
                {
                if (trData.taxaNames[i] && taxaNames[j] &&
                    strcmp(trData.taxaNames[i], taxaNames[j]) == 0)
                    {
                    tFileToMatrix[i] = j;
                    break;
                    }
                }
            }

        /* Buffer CI results for sorting */
        typedef struct { char name[64]; int steps; int minSteps; YFlt ci; } CIResult;
        CIResult *ciResults = (CIResult *) SafeCalloc(numChar, sizeof(CIResult));
        int nCI = 0;

        /* For each original character, compute Fitch score */
        for (j = 0; j < numChar; j++)
            {
            int *tipStates;
            int fitchScore, nObsStates, minSteps;
            YFlt ci;
            int statesSeen;

            if (charInfo[j].isExcluded == YES)
                continue;

            /* Build tip state bit sets for this character.
               matrix[pos(taxon, char, numChar)] is a bit set of observed states. */
            tipStates = (int *) SafeCalloc(trData.nTaxa, sizeof(int));
            statesSeen = 0;
            for (i = 0; i < trData.nTaxa; i++)
                {
                if (tFileToMatrix[i] >= 0)
                    {
                    tipStates[i] = matrix[pos(tFileToMatrix[i], j, numChar)];
                    statesSeen |= tipStates[i];
                    }
                }

            /* Mask tip states to only valid states for this character */
            {
                int validMask = (1 << charInfo[j].numStates) - 1;
                for (i = 0; i < trData.nTaxa; i++)
                    if (tipStates[i] != 0)
                        tipStates[i] &= validMask;
            }

            /* Count actually observed states (only unambiguous tips) */
            {
                int obsBits = 0;
                for (i = 0; i < trData.nTaxa; i++)
                    {
                    int ts = tipStates[i];
                    /* Only count if exactly one bit set (unambiguous) */
                    if (ts > 0 && (ts & (ts - 1)) == 0)
                        obsBits |= ts;
                    }
                nObsStates = 0;
                while (obsBits) { nObsStates += obsBits & 1; obsBits >>= 1; }
            }

            if (nObsStates <= 1)
                {
                /* Invariant character — CI = 1.0 by definition */
                free(tipStates);
                continue;
                }

            minSteps = nObsStates - 1;

            if (FitchParsimonyScore(trData.newickTrees[mapIdx], trData.nTaxa,
                                    tipStates, &fitchScore) == ERROR)
                {
                free(tipStates);
                continue;
                }

            ci = (fitchScore > 0) ? (YFlt)minSteps / fitchScore : 1.0;

            /* Store result */
            if (charLabels != NULL && charLabels[j] != NULL)
                snprintf(ciResults[nCI].name, sizeof(ciResults[nCI].name), "%s", charLabels[j]);
            else
                snprintf(ciResults[nCI].name, sizeof(ciResults[nCI].name), "char_%d", j + 1);
            ciResults[nCI].steps = fitchScore;
            ciResults[nCI].minSteps = minSteps;
            ciResults[nCI].ci = ci;
            nCI++;

            free(tipStates);
            }

        /* Sort by steps descending (most homoplasy first) */
        for (i = 0; i < nCI - 1; i++)
            for (j = i + 1; j < nCI; j++)
                if (ciResults[j].steps > ciResults[i].steps)
                    { CIResult tmp = ciResults[i]; ciResults[i] = ciResults[j]; ciResults[j] = tmp; }

        YvyraPrint ("%s   %-20s %6s %8s %6s\n", spacer, "Character", "Steps", "MinSteps", "CI");
        YvyraPrint ("%s   %-20s %6s %8s %6s\n", spacer, "---------", "-----", "--------", "----");
        for (i = 0; i < nCI; i++)
            YvyraPrint ("%s   %-20s %6d %8d %6.3f\n",
                          spacer, ciResults[i].name, ciResults[i].steps, ciResults[i].minSteps, ciResults[i].ci);

        free(ciResults);
        free(tFileToMatrix);
        }

    /* Clean up */
    free(meanLL);
    free(proportion);
    free(support);
    free(hasClade);
    FreeSiteLnL(&slData);
    FreeTreeSamples(&trData);

    return (NO_ERROR);
}


/*--------------------------------------------------------------
|   DoChardiagParm: Parse chardiag command parameters
|   Syntax: chardiag clade={Luwian,Lycian};
---------------------------------------------------------------*/
int DoChardiagParm (char *parmName, char *tkn)
{
    char    *p, *tok;

    if (expecting == Expecting(PARAMETER))
        {
        /* Parse the whole token as "clade={taxon1,taxon2,...}" */
        if (strncmp(tkn, "clade", 5) == 0 || strncmp(tkn, "Clade", 5) == 0)
            {
            expecting = Expecting(EQUALSIGN);
            }
        else
            {
            YvyraPrint ("%s   Unknown parameter '%s' for chardiag\n", spacer, tkn);
            return (ERROR);
            }
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(LEFTCURL);
        }
    else if (expecting == Expecting(LEFTCURL))
        {
        /* Start collecting clade taxa */
        diagNCladeTaxa = 0;
        diagCladeStr[0] = '\0';
        expecting = Expecting(ALPHA) | Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA) || expecting == Expecting(NUMBER))
        {
        /* Add taxon to clade */
        if (diagNCladeTaxa < 100)
            {
            diagCladeTaxa[diagNCladeTaxa] = (char *) SafeMalloc(strlen(tkn) + 1);
            strcpy(diagCladeTaxa[diagNCladeTaxa], tkn);
            if (diagNCladeTaxa > 0)
                strcat(diagCladeStr, ",");
            strcat(diagCladeStr, tkn);
            diagNCladeTaxa++;
            }
        expecting = Expecting(COMMA) | Expecting(RIGHTCURL);
        }
    else if (expecting == Expecting(COMMA))
        {
        expecting = Expecting(ALPHA) | Expecting(NUMBER);
        }
    else if (expecting == Expecting(RIGHTCURL))
        {
        expecting = Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


/*--------------------------------------------------------------
|   DoSensitivity: Alias for chardiag (same implementation)
---------------------------------------------------------------*/
int DoSensitivity (void)
{
    return DoChardiag();
}

int DoSensitivityParm (char *parmName, char *tkn)
{
    return DoChardiagParm(parmName, tkn);
}


/*--------------------------------------------------------------
|   ComputeESS: Compute Effective Sample Size from an array
|   of MCMC samples using the initial positive sequence estimator.
---------------------------------------------------------------*/
static YFlt ComputeESS (YFlt *samples, int n)
{
    int     i, maxLag;
    YFlt  mean, var, sum, rho, prevRho;

    if (n < 4) return (YFlt)n;

    /* Compute mean */
    mean = 0.0;
    for (i = 0; i < n; i++) mean += samples[i];
    mean /= n;

    /* Compute variance */
    var = 0.0;
    for (i = 0; i < n; i++) var += (samples[i] - mean) * (samples[i] - mean);
    var /= (n - 1);
    if (var < 1e-300) return (YFlt)n;

    /* Initial positive sequence estimator (Geyer 1992) */
    maxLag = n / 2;
    sum = 0.0;
    prevRho = 1.0;
    for (i = 1; i < maxLag; i++)
        {
        rho = 0.0;
        int j;
        for (j = 0; j < n - i; j++)
            rho += (samples[j] - mean) * (samples[j + i] - mean);
        rho /= (n - i) * var;

        /* Stop when autocorrelation goes negative */
        if (rho < 0.0) break;
        sum += rho;
        prevRho = rho;
        }

    return (YFlt)n / (1.0 + 2.0 * sum);
}


/*--------------------------------------------------------------
|   DoConvergence: Comprehensive convergence diagnostics.
|   Reads .pstat, .tstat, and .w files to report ESS, PSRF,
|   split frequency diagnostics, and weight posteriors.
---------------------------------------------------------------*/
int DoConvergence (void)
{
    char    filename[300];
    FILE    *fp;
    char    *line;
    int     nIssues = 0;
    YFlt  minESS = 1e30;
    char    minESSParam[100] = "";

    line = (char *) SafeMalloc(DIAG_MAX_LINE);
    if (!line) return (ERROR);

    YvyraPrint ("\n%s   ====================================================================\n", spacer);
    YvyraPrint ("%s   CONVERGENCE DIAGNOSTICS\n", spacer);
    YvyraPrint ("%s   ====================================================================\n\n", spacer);

    /* --- 1. Read .pstat for parameter ESS and PSRF --- */
    if (chainParams.numRuns == 1)
        sprintf(filename, "%s.pstat", chainParams.chainFileName);
    else
        sprintf(filename, "%s.pstat", chainParams.chainFileName);

    fp = fopen(filename, "r");
    if (fp)
        {
        YvyraPrint ("%s   Parameter convergence (from %s):\n", spacer, filename);
        YvyraPrint ("%s   %-20s %10s %10s %10s\n", spacer, "Parameter", "min ESS", "avg ESS", "PSRF");
        YvyraPrint ("%s   %-20s %10s %10s %10s\n", spacer, "---------", "-------", "-------", "----");

        while (fgets(line, DIAG_MAX_LINE, fp))
            {
            if (line[0] == '[' || line[0] == '\n' || line[0] == '\r') continue;
            if (strncmp(line, "Parameter", 9) == 0) continue;  /* header */

            char param[100];
            YFlt mean, variance, lower, upper, median, ess1, ess2, psrf;
            int nread = sscanf(line, "%99s %lf %lf %lf %lf %lf %lf %lf %lf",
                               param, &mean, &variance, &lower, &upper, &median,
                               &ess1, &ess2, &psrf);

            /* Handle both single-run (7 cols: ESS) and multi-run (9 cols: minESS avgESS PSRF) */
            if (nread >= 7)
                {
                YFlt essVal = ess1;  /* single-run ESS or multi-run minESS */
                char status[20] = "";

                if (essVal < 100)
                    { strcpy(status, "LOW"); nIssues++; }
                else if (essVal < 200)
                    strcpy(status, "marginal");
                else
                    strcpy(status, "ok");

                if (nread >= 9 && psrf > 1.05)
                    { strcpy(status, "NOT CONVERGED"); nIssues++; }

                YvyraPrint ("%s   %-20s %10.1f", spacer, param, essVal);
                if (nread >= 8)
                    YvyraPrint (" %10.1f", ess2);
                else
                    YvyraPrint (" %10s", "N/A");
                if (nread >= 9)
                    YvyraPrint (" %10.3f", psrf);
                else
                    YvyraPrint (" %10s", "N/A");
                YvyraPrint ("   %s\n", status);

                if (essVal < minESS)
                    { minESS = essVal; strncpy(minESSParam, param, 99); }
                }
            }
        fclose(fp);
        YvyraPrint ("\n");
        }

    /* --- 2. Read .tstat for split frequency diagnostics --- */
    sprintf(filename, "%s.tstat", chainParams.chainFileName);

    fp = fopen(filename, "r");
    if (fp)
        {
        /* Detect format from header: multi-run has "Stddev" column */
        int multiRun = NO;
        int nSplits = 0;
        YFlt maxStdDev = 0.0;

        while (fgets(line, DIAG_MAX_LINE, fp))
            {
            if (line[0] == '[' || line[0] == '\n' || line[0] == '\r') continue;
            if (strncmp(line, "ID", 2) == 0)
                {
                if (strstr(line, "Stddev") || strstr(line, "stddev") || strstr(line, "StdDev"))
                    multiRun = YES;
                break;
                }
            }

        YvyraPrint ("%s   Split frequency diagnostics (from %s):\n", spacer, filename);
        if (multiRun)
            {
            YvyraPrint ("%s   %-10s %10s %10s %10s\n", spacer, "Split ID", "Prob", "Std Dev", "Status");
            YvyraPrint ("%s   %-10s %10s %10s %10s\n", spacer, "--------", "----", "-------", "------");
            }
        else
            {
            YvyraPrint ("%s   %-10s %10s\n", spacer, "Split ID", "Prob");
            YvyraPrint ("%s   %-10s %10s\n", spacer, "--------", "----");
            }

        while (fgets(line, DIAG_MAX_LINE, fp))
            {
            if (line[0] == '[' || line[0] == '\n' || line[0] == '\r') continue;
            if (strncmp(line, "ID", 2) == 0) continue;

            int id, nobs;
            YFlt prob, stddev, minp, maxp;
            int nread = sscanf(line, "%d %d %lf %lf %lf %lf",
                               &id, &nobs, &prob, &stddev, &minp, &maxp);

            if (multiRun && nread >= 4)
                {
                char status[20] = "ok";
                if (stddev > 0.05)
                    { strcpy(status, "HIGH"); nIssues++; }
                else if (stddev > 0.01)
                    strcpy(status, "marginal");

                YvyraPrint ("%s   %-10d %10.4f %10.4f   %s\n", spacer, id, prob, stddev, status);
                if (stddev > maxStdDev) maxStdDev = stddev;
                nSplits++;
                }
            else if (nread >= 3)
                {
                YvyraPrint ("%s   %-10d %10.4f\n", spacer, id, prob);
                nSplits++;
                }
            }
        fclose(fp);
        if (nSplits == 0)
            YvyraPrint ("%s   (no informative splits found)\n", spacer);
        if (nSplits > 0 && multiRun)
            YvyraPrint ("%s   Max split std dev: %.6f\n", spacer, maxStdDev);
        YvyraPrint ("\n");
        }

    /* --- 3. Read .w file for weight ESS --- */
    {
        char wFilename[300];
        if (chainParams.numRuns == 1)
            sprintf(wFilename, "%s.w", chainParams.chainFileName);
        else
            sprintf(wFilename, "%s.run1.w", chainParams.chainFileName);

        fp = fopen(wFilename, "r");
        if (fp)
            {
            /* Count columns and rows */
            char *header = NULL;
            int nCols = 0, nRows = 0;
            char *colNames[MAX_ESTIMATED_WEIGHTS + 1];

            /* Read header */
            while (fgets(line, DIAG_MAX_LINE, fp))
                {
                if (line[0] == '[') continue;
                header = line;
                break;
                }

            if (header)
                {
                char *tok = strtok(header, "\t\n\r");
                while (tok)
                    {
                    if (nCols > 0 && nCols <= MAX_ESTIMATED_WEIGHTS)
                        {
                        colNames[nCols - 1] = (char *) SafeMalloc(strlen(tok) + 1);
                        if (colNames[nCols - 1]) strcpy(colNames[nCols - 1], tok);
                        }
                    nCols++;
                    tok = strtok(NULL, "\t\n\r");
                    }
                nCols--;  /* subtract Gen column */
                }

            /* Count data rows */
            long dataStart = ftell(fp);
            while (fgets(line, DIAG_MAX_LINE, fp))
                if (line[0] != '[' && line[0] != '\n' && line[0] != '\r')
                    nRows++;

            /* Apply burnin (25%) */
            int burnin = nRows / 4;
            int nPost = nRows - burnin;

            if (nCols > 0 && nPost > 2)
                {
                YFlt *samples = (YFlt *) SafeCalloc(nPost, sizeof(YFlt));
                int col;

                YvyraPrint ("%s   Estimated weight convergence (from %s):\n", spacer, wFilename);
                YvyraPrint ("%s   %-20s %10s %10s %10s %10s\n", spacer, "Weight", "Mean", "Std Dev", "ESS", "Status");
                YvyraPrint ("%s   %-20s %10s %10s %10s %10s\n", spacer, "------", "----", "-------", "---", "------");

                for (col = 0; col < nCols; col++)
                    {
                    fseek(fp, dataStart, SEEK_SET);
                    int row = 0, postRow = 0;
                    while (fgets(line, DIAG_MAX_LINE, fp))
                        {
                        if (line[0] == '[' || line[0] == '\n' || line[0] == '\r') continue;
                        row++;
                        if (row <= burnin) continue;
                        char *tok2 = strtok(line, "\t\n\r");
                        int c = 0;
                        while (tok2)
                            {
                            if (c == col + 1 && postRow < nPost)  /* +1 to skip Gen */
                                samples[postRow++] = atof(tok2);
                            c++;
                            tok2 = strtok(NULL, "\t\n\r");
                            }
                        }

                    /* Compute stats */
                    YFlt mean = 0, var = 0, ess;
                    int s;
                    for (s = 0; s < postRow; s++) mean += samples[s];
                    mean /= postRow;
                    for (s = 0; s < postRow; s++) var += (samples[s] - mean) * (samples[s] - mean);
                    var = (postRow > 1) ? sqrt(var / (postRow - 1)) : 0;
                    ess = ComputeESS(samples, postRow);

                    char status[20] = "ok";
                    if (ess < 100)
                        { strcpy(status, "LOW"); nIssues++; }
                    else if (ess < 200)
                        strcpy(status, "marginal");

                    YvyraPrint ("%s   %-20s %10.4f %10.4f %10.1f   %s\n",
                                spacer, col < nCols ? colNames[col] : "?", mean, var, ess, status);

                    if (ess < minESS)
                        { minESS = ess; strncpy(minESSParam, col < nCols ? colNames[col] : "w(?)", 99); }
                    }

                free(samples);
                for (col = 0; col < nCols; col++)
                    if (colNames[col]) free(colNames[col]);
                }
            fclose(fp);
            YvyraPrint ("\n");
            }
    }

    /* --- 4. Overall assessment --- */
    YvyraPrint ("%s   ====================================================================\n", spacer);
    YvyraPrint ("%s   CONVERGENCE ASSESSMENT\n", spacer);
    YvyraPrint ("%s   ====================================================================\n", spacer);

    if (minESS < 1e29)
        YvyraPrint ("%s   Minimum ESS: %.1f (%s)\n", spacer, minESS, minESSParam);

    if (nIssues == 0)
        {
        YvyraPrint ("%s   Status: CONVERGED\n", spacer);
        YvyraPrint ("%s   All parameters have ESS >= 200 and PSRF <= 1.05\n", spacer);
        }
    else
        {
        YvyraPrint ("%s   Status: %d ISSUE(S) DETECTED\n", spacer, nIssues);
        if (minESS < 200)
            YvyraPrint ("%s   - ESS below 200 for some parameters (run longer chain)\n", spacer);
        YvyraPrint ("%s   Recommendations:\n", spacer);
        YvyraPrint ("%s   - Increase iterations (current chain may be too short)\n", spacer);
        YvyraPrint ("%s   - Check for multimodality (multiple peaks in posterior)\n", spacer);
        YvyraPrint ("%s   - Consider using multiple runs (nrun=2) for PSRF\n", spacer);
        }

    YvyraPrint ("%s   ====================================================================\n\n", spacer);

    free(line);
    return (NO_ERROR);
}



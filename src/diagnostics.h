#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include "bayes.h"

/* Per-site log-likelihood samples read from .p file */
typedef struct {
    int     nSamples;       /* number of post-burnin samples              */
    int     nChars;         /* number of per-site lnL columns             */
    YFlt  *siteLnL;       /* [nSamples * nChars] row-major matrix       */
    YFlt  *totalLnL;      /* [nSamples] total log-likelihood per sample */
    char    **charNames;    /* [nChars] column header names               */
} SiteLnLData;

/* Tree samples read from .t file */
typedef struct {
    int     nSamples;       /* number of post-burnin tree samples         */
    int     nTaxa;          /* number of taxa                             */
    char    **taxaNames;    /* [nTaxa] taxon names from translate block   */
    char    **newickTrees;  /* [nSamples] Newick strings                  */
} TreeData;

/* Read per-site log-likelihoods from .p file(s) */
int  ReadSiteLnL (char *filename, YFlt relBurnin, SiteLnLData *data);
void FreeSiteLnL (SiteLnLData *data);

/* Read tree samples from .t file(s) */
int  ReadTreeSamples (char *filename, YFlt relBurnin, TreeData *data);
void FreeTreeSamples (TreeData *data);

/* Check if a Newick tree contains a clade (monophyletic group) */
int  TreeHasClade (char *newick, char **cladeTaxa, int nCladeTaxa,
                   char **allTaxa, int nTaxa);

/* Simple parsed tree for post-processing */
typedef struct {
    int     nTips;          /* number of tip nodes (= nTaxa)                */
    int     nIntNodes;      /* number of internal nodes (nTips - 1)         */
    int     nNodes;         /* total nodes (2 * nTips - 1)                  */
    int     *left;          /* [nNodes] left child index (-1 if tip)        */
    int     *right;         /* [nNodes] right child index (-1 if tip)       */
    int     *parent;        /* [nNodes] parent index (-1 if root)           */
    int     *taxon;         /* [nNodes] taxon number for tips (1-based), -1 for internal */
    BitsLong *descTips;     /* [nNodes] bit set of descendant tip indices   */
} SimpleTree;

int  ParseNewickTree (char *newick, int nTaxa, SimpleTree *tree);
void FreeSimpleTree  (SimpleTree *tree);

/* Diagnostic computations */
int  DoChardiag (void);
int  DoChardiagParm (char *parmName, char *tkn);
int  DoSensitivity (void);
int  DoSensitivityParm (char *parmName, char *tkn);
int  DoAsrentropy (void);
int  DoConvergence (void);

/* Fitch parsimony on a Newick tree */
int  FitchParsimonyScore (char *newick, int nTaxa, int *tipStates, int *score);

/* Fitch downpass with state set extraction at internal nodes */
int  FitchDownpassSets (SimpleTree *tree, int *tipStates, int *nodeStateSets);

/* Math utilities */
YFlt ShannonEntropy (YFlt *probs, int n);

#endif /* DIAGNOSTICS_H_ */

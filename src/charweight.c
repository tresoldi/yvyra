/*
 *  yvyra - Bayesian Phylogenetics for Linguistic Data
 *  Based on MrBayes 3.2.7a by Ronquist et al.
 *
 *  charweight.c: Estimated character weight MCMC sampling.
 *  Uses log-normal multiplier proposals with Gamma(2,2) priors.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bayes.h"
#include "charweight.h"
#include "utils.h"

/* ---- State ---- */

int     hasEstimatedWeights = NO;

/* Per-original-character data */
static int      nEstimated = 0;             /* number of estimated-weight chars */
static int      estimatedIdx[MAX_ESTIMATED_WEIGHTS]; /* original char indices */
static YFlt   weights[MAX_ESTIMATED_WEIGHTS]; /* current weight values */
static YFlt   priorAlpha = 2.0;            /* Gamma prior shape */
static YFlt   priorBeta = 2.0;             /* Gamma prior rate */

/* Proposal state (for accept/reject) */
static int      lastProposedIdx = -1;       /* index into estimatedIdx[] */
static YFlt   lastOldWeight = 1.0;
static YFlt   tuning = 0.5;                /* multiplier tuning parameter */

/* Per-compressed-character weight lookup */
static YFlt   *compCharWeight = NULL;       /* [numCompressedChars] */


/* ---- Gamma log-density ---- */

static YFlt LnGammaDensity (YFlt x, YFlt alpha, YFlt beta)
{
    /* log[ beta^alpha / Gamma(alpha) * x^(alpha-1) * exp(-beta*x) ] */
    if (x <= 0.0) return -1e300;
    return alpha * log(beta) - LnGamma(alpha) + (alpha - 1.0) * log(x) - beta * x;
}


/* ---- Public API ---- */

int InitCharWeights (void)
{
    int i;
    extern int      *isWeightEstimated;   /* from yamlcmd.c globals */
    extern YFlt   *initialWeights;        /* from yamlcmd.c globals */
    extern int      numChar;
    extern int      numCompressedChars;
    extern int      *origChar;

    nEstimated = 0;
    hasEstimatedWeights = NO;

    if (!isWeightEstimated) return NO_ERROR;

    /* Collect estimated-weight characters */
    for (i = 0; i < numChar && nEstimated < MAX_ESTIMATED_WEIGHTS; i++)
        {
        if (isWeightEstimated[i])
            {
            estimatedIdx[nEstimated] = i;
            weights[nEstimated] = initialWeights ? initialWeights[i] : 1.0;
            nEstimated++;
            }
        }

    if (nEstimated == 0) return NO_ERROR;
    hasEstimatedWeights = YES;

    /* Build per-compressed-character weight lookup */
    compCharWeight = (YFlt *) calloc(numCompressedChars, sizeof(YFlt));
    if (!compCharWeight) return ERROR;

    /* Initialize all to 1.0 */
    for (i = 0; i < numCompressedChars; i++)
        compCharWeight[i] = 1.0;

    /* Map estimated weights to compressed characters */
    for (i = 0; i < numCompressedChars; i++)
        {
        if (origChar[i] < 0) continue;  /* dummy char */
        int oc = origChar[i];
        int j;
        for (j = 0; j < nEstimated; j++)
            {
            if (estimatedIdx[j] == oc)
                {
                compCharWeight[i] = weights[j];
                break;
                }
            }
        }

    return NO_ERROR;
}


void FreeCharWeights (void)
{
    if (compCharWeight)
        {
        free(compCharWeight);
        compCharWeight = NULL;
        }
    nEstimated = 0;
    hasEstimatedWeights = NO;
}


YFlt ProposeCharWeight (int chain, RandLong *seed)
{
    YFlt oldW, newW, lnPriorRatio, lnProposalRatio;
    int idx, origIdx, i;
    extern int numCompressedChars;
    extern int *origChar;

    if (nEstimated == 0) return 0.0;

    /* Pick a random estimated-weight character */
    idx = (int)(RandomNumber(seed) * nEstimated);
    origIdx = estimatedIdx[idx];
    oldW = weights[idx];

    /* Log-normal multiplier proposal */
    newW = oldW * exp(tuning * (RandomNumber(seed) - 0.5));
    if (newW < 0.001) newW = 0.001;  /* lower bound */
    if (newW > 1000.0) newW = 1000.0; /* upper bound */

    /* Store for accept/reject */
    lastProposedIdx = idx;
    lastOldWeight = oldW;

    /* Update weight */
    weights[idx] = newW;

    /* Update compressed character lookup */
    for (i = 0; i < numCompressedChars; i++)
        if (origChar[i] >= 0 && origChar[i] == origIdx)
            compCharWeight[i] = newW;

    /* Compute log Hastings ratio (proposal ratio for multiplier) */
    lnProposalRatio = log(newW / oldW);

    /* Compute log prior ratio: Gamma(alpha, beta) */
    lnPriorRatio = LnGammaDensity(newW, priorAlpha, priorBeta)
                 - LnGammaDensity(oldW, priorAlpha, priorBeta);

    return lnPriorRatio + lnProposalRatio;
}


void AcceptCharWeight (void)
{
    /* Weight already updated in ProposeCharWeight */
    lastProposedIdx = -1;
}


void RejectCharWeight (void)
{
    int i, origIdx;
    extern int numCompressedChars;
    extern int *origChar;

    if (lastProposedIdx < 0) return;

    origIdx = estimatedIdx[lastProposedIdx];
    weights[lastProposedIdx] = lastOldWeight;

    /* Revert compressed character lookup */
    for (i = 0; i < numCompressedChars; i++)
        if (origChar[i] >= 0 && origChar[i] == origIdx)
            compCharWeight[i] = lastOldWeight;

    lastProposedIdx = -1;
}


YFlt GetCharWeightMultiplier (int compCharIdx)
{
    if (!compCharWeight) return 1.0;
    return compCharWeight[compCharIdx];
}


void PrintCharWeightHeader (FILE *fp)
{
    int i;
    extern char **charLabels;
    extern int numChar;

    for (i = 0; i < nEstimated; i++)
        {
        int oc = estimatedIdx[i];
        if (charLabels && oc < numChar && charLabels[oc])
            fprintf(fp, "\tw(%s)", charLabels[oc]);
        else
            fprintf(fp, "\tw(%d)", oc + 1);
        }
}


void PrintCharWeightValues (FILE *fp)
{
    int i;
    for (i = 0; i < nEstimated; i++)
        fprintf(fp, "\t%f", weights[i]);
}

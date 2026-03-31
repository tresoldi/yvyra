/*
 *  yvyra - Bayesian Phylogenetics for Linguistic Data
 *  Based on MrBayes 3.2.7a by Ronquist et al.
 *
 *  Copyright (C) 2026 Tiago Tresoldi
 *  Copyright (C) 2002-2023
 *
 *  John P. Huelsenbeck
 *  Dept. Integrative Biology
 *  University of California, Berkeley
 *  Berkeley, CA 94720-3140
 *  johnh@berkeley.edu
 *
 *  Fredrik Ronquist
 *  Swedish Museum of Natural History
 *  Box 50007
 *  SE-10405 Stockholm, SWEDEN
 *  fredrik.ronquist@nrm.se
 *
 *  With important contributions by
 *
 *  Paul van der Mark (paulvdm@sc.fsu.edu)
 *  Maxim Teslenko (maxkth@gmail.com)
 *  Chi Zhang (zhangchicool@gmail.com)
 *
 *  and by many users (run 'acknowledgments' to see more info)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 */

#include "bayes.h"
#include "charweight.h"
#include "likelihood.h"
#include "model.h"
#include "utils.h"

#define LIKE_EPSILON                1.0e-300

/* global variables declared here */
CLFlt     *preLikeL;                  /* precalculated cond likes for left descendant */
CLFlt     *preLikeR;                  /* precalculated cond likes for right descendant*/
CLFlt     *preLikeA;                  /* precalculated cond likes for ancestor        */

/* global variables used here but declared elsewhere */
extern int      *chainId;
extern int      numLocalChains;
extern int      rateProbRowSize;            /* size of rate probs for one chain one state   */
extern YFlt   **rateProbs;                /* pointers to rate probs used by adgamma model */

/* local prototypes */
YFlt    GetRate (int division, int chain);
int       RemoveNodeScalers(TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       RemoveNodeScalers_SSE(TreeNode *p, int division, int chain);
#endif
#if defined (AVX_ENABLED)
int       RemoveNodeScalers_AVX(TreeNode *p, int division, int chain);
#endif
int       SetStdQMatrix (YFlt **a, int nStates, YFlt *bs, int cType);
int       UpDateCijk (int whichPart, int whichChain);


/*----------------------------------------------------------------
|
|   CondLikeDown_Bin: binary model with or without rate variation
|
-----------------------------------------------------------------*/
int CondLikeDown_Std (TreeNode *p, int division, int chain)
{
    int             a, c, h, i, j, k, nStates, nCats, tmp;
    CLFlt           *clL, *clR, *clP, *pL, *pR, *tiPL, *tiPR, likeL, likeR;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* Flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* Conditional likelihood space is assumed to be arranged in numGammaCats blocks of data. Each block contains all data for one gamma category.
    Each gamma cat block consist of numChars sequences of data, each of this sequences corresponds to a character of data matrix. 
    A sequence consists of nStates for all non-binary data, otherwise length of sequence is nStates*numBetaCats (i.e. 2*numBetaCats) */

    /* calculate ancestral probabilities */
    for (k=h=0; k<m->numRateCats; k++)
        {
        /* calculate ancestral probabilities */
        for (c=0; c<m->numChars; c++)
            {
            nStates = m->nStates[c];
        
            /* the following lines ensure that nCats is 1 unless */
            /* the character is binary and beta categories are used  */
            if (nStates == 2)
                nCats = m->numBetaCats;
            else
                nCats = 1;

            tmp = k*nStates*nStates; /* tmp contains offset to skip rate cats that already processed*/
            tiPL = pL + m->tiIndex[c] + tmp;
            tiPR = pR + m->tiIndex[c] + tmp;
            tmp = (m->numRateCats-1)*2*2; /* tmp contains size of block of tpi matrices across all rate cats (minus one) for single beta category. Further used only if character is binary to jump to next beta category */
                
            for (j=0; j<nCats;j++)
                {
                for (a=0; a<nStates; a++)
                    {
                    likeL = likeR = 0.0;
                    for (i=0; i<nStates; i++)
                        {
                        likeL += *(tiPL++) * clL[i];
                        likeR += *(tiPR++) * clR[i];
                        }
                    clP[h++] = likeL * likeR;
                    }
                clL += nStates;
                clR += nStates;
        
                tiPL += tmp;
                tiPR += tmp;
                }
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeRoot_Bin: binary model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_Std (TreeNode *p, int division, int chain)
{
    int             a, c, h, i, j, k, nStates=0, nCats=0, tmp;
    CLFlt           *clL, *clR, *clP, *clA, *pL, *pR, *pA, *tiPL, *tiPR, *tiPA,
                    likeL, likeR, likeA;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* calculate ancestral probabilities */
    for (k=h=0; k<m->numRateCats; k++)
        {
        /* calculate ancestral probabilities */
        for (c=0; c<m->numChars; c++)
            {
            nStates = m->nStates[c];
        
            /* the following lines ensure that nCats is 1 unless */
            /* the character is binary and beta categories are used  */
            if (nStates == 2)
                nCats = m->numBetaCats;
            else
                nCats = 1;

            tmp = k*nStates*nStates; /* tmp contains offset to skip gamma cats that already processed*/
            tiPL = pL + m->tiIndex[c] + tmp;
            tiPR = pR + m->tiIndex[c] + tmp;
            tiPA = pA + m->tiIndex[c] + tmp;
            tmp = (m->numRateCats-1)*2*2; /* tmp contains size of block of tpi matrices across all rate cats (minus one) for single beta category. Further used only if character is binary to jump to next beta category */
                
            for (j=0; j<nCats;j++)
                {
                for (a=0; a<nStates; a++)
                    {
                    likeL = likeR = likeA = 0.0;
                    for (i=0; i<nStates; i++)
                        {
                        likeL += *(tiPL++) * clL[i];
                        likeR += *(tiPR++) * clR[i];
                        likeA += *(tiPA++) * clA[i];
                        }
                    clP[h++] = likeL * likeR * likeA;
                    }
                clL += nStates;
                clR += nStates;
                clA += nStates;
        
                tiPL += tmp;
                tiPR += tmp;
                tiPA += tmp;
                }
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeUp_Bin: pull likelihoods up and calculate scaled
|       finals, binary model with or without rate variation
|
-----------------------------------------------------------------*/
int     CondLikeUp_Std (TreeNode *p, int division, int chain)
{
    int             a, c, i, j, k, t, nStates, nCats, coppySize,tmp;
    CLFlt           *clFA, *clFP, *clDP, *pA, *tiP, condLikeUp[MAX_STD_STATES], sum;
    ModelInfo       *m;
    
    /* find model settings for this division */
    m = &modelSettings[division];
    
    /* calculate final states */
    if (p->anc->anc == NULL)
        {
        /* this is the root node */
        /* find conditional likelihood pointers = down cond likes */
        /* use conditional likelihood scratch space for final cond likes */
        clDP = m->condLikes[m->condLikeIndex[chain][p->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index]];
        
        coppySize=0;
        /* final cond likes = downpass cond likes */
        for (c=0; c<m->numChars; c++)
            {
            /* calculate nStates and nCats */
            nStates = m->nStates[c];
            
            /* the following lines ensure that nCats is 1 unless */
            /* the character is binary and beta categories are used  */
            if (nStates == 2)
                nCats = m->numBetaCats;
            else
                nCats = 1;

            coppySize+=nCats*nStates;
            }

        /* finally multiply with the rate cats */
        coppySize *= m->numRateCats;

        /* copy cond likes */ 
        for (k=0; k<coppySize; k++)
            *(clFP++) = *(clDP++);
        }
    else
        {
        /* find conditional likelihood pointers */
        /* use conditional likelihood scratch space for final cond likes */
        clFA = m->condLikes[m->condLikeScratchIndex[p->anc->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index     ]];
        clDP = m->condLikes[m->condLikeIndex[chain][p->index     ]];

        /* find transition probabilities */
        pA = m->tiProbs[m->tiProbsIndex[chain][p->index]];
        
        for (k=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {

                /* calculate nStates and nCats */
                nStates = m->nStates[c];
                
                /* the following lines ensure that nCats is 1 unless */
                /* the character is binary and beta categories are used  */
                if (nStates == 2)
                    nCats = m->numBetaCats;
                else
                    nCats = 1;

                tmp = k*nStates*nStates; /* tmp contains offset to skip rate cats that already processed*/
                tiP = pA + m->tiIndex[c] + tmp;
                tmp = (m->numRateCats-1)*2*2; /* tmp contains size of block of tpi matrices across all rate cats (minus one) for single beta category. Further used only if character is binary to jump to next beta category */

                /* now calculate the final cond likes */
                for (t=0; t<nCats; t++)
                    {
                    for (a=j=0; a<nStates; a++)
                        {
                        sum = 0.0;
                        for (i=0; i<nStates; i++)
                            sum += tiP[j++]*clDP[i];
                        if (sum == 0.0)
                            condLikeUp[a] = 0.0;    /* we lost the conditional likelihood in the downpass (can occur in gamma model) */
                        else
                            condLikeUp[a] = clFA[a] / sum;
                        }
                        
                    for (a=j=0; a<nStates; a++)
                        {
                        sum = 0.0;
                        for (i=0; i<nStates; i++)
                            {
                            sum += condLikeUp[i] * tiP[j++];
                            }
                        clFP[a] = sum * clDP[a];
                        }

                    clFP += nStates;
                    clFA += nStates;
                    clDP += nStates;
                    tiP += tmp;
                    }
                }
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeScaler_Gen: general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_Std (TreeNode *p, int division, int chain)
{
    int             c, n, k, nStates, numReps;
    CLFlt           scaler, *clPtr, **clP, *scP, *lnScaler;
    ModelInfo       *m;

    m = &modelSettings[division];

    numReps=0;
    for (c=0; c<m->numChars; c++)
        {
        if (m->nStates[c] == 2)
            numReps += m->numBetaCats * 2;
        else
            numReps += m->nStates[c];
        }

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += numReps;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numChars; c++)
        {
        scaler = 0.0;
        nStates = m->nStates[c];
        if (nStates == 2)
            nStates = m->numBetaCats * 2;

        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                if (clP[k][n] > scaler)
                    scaler = clP[k][n];
                }
            }

        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                clP[k][n] /= scaler;
            clP[k] += nStates;
            }

        scP[c]       = (CLFlt) log (scaler);    /* store node scaler */
        lnScaler[c] += scP[c];                  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;
        
    return NO_ERROR;
}


/* FlipCijkSpace: Flip space for cijks with scratch area */
void FlipCijkSpace (ModelInfo* m, int chain)
{
    int         temp;

    temp                = m->cijkIndex[chain];
    m->cijkIndex[chain] = m->cijkScratchIndex;
    m->cijkScratchIndex = temp;
}


/* FlipCondLikeSpace: Flip space for conditional likelihoods with scratch area */
void FlipCondLikeSpace (ModelInfo* m, int chain, int nodeIndex)
{
    int         temp;

    temp                               = m->condLikeIndex[chain][nodeIndex];
    m->condLikeIndex[chain][nodeIndex] = m->condLikeScratchIndex[nodeIndex];
    m->condLikeScratchIndex[nodeIndex] = temp;
}


/* FlipNodeScalerSpace: Flip space for node scalers and scaler flag with scratch area */
void FlipNodeScalerSpace (ModelInfo* m, int chain, int nodeIndex)
{
    int         temp;

    temp                                 = m->nodeScalerIndex[chain][nodeIndex];
    m->nodeScalerIndex[chain][nodeIndex] = m->nodeScalerScratchIndex[nodeIndex];
    m->nodeScalerScratchIndex[nodeIndex] = temp;

    temp                                 = m->unscaledNodes[chain][nodeIndex];
    m->unscaledNodes[chain][nodeIndex]   = m->unscaledNodesScratch[nodeIndex];
    m->unscaledNodesScratch[nodeIndex]   = temp;
}


/* FlipSiteScalerSpace: Flip space for ln site scalers */
void FlipSiteScalerSpace (ModelInfo *m, int chain)
{
    int  temp;


    temp = m->siteScalerIndex[chain];
    m->siteScalerIndex[chain] = m->siteScalerScratchIndex;
    m->siteScalerScratchIndex = temp;

}


/* FlipTiProbsSpace: Flip space for ti probs with scratch area */
void FlipTiProbsSpace (ModelInfo* m, int chain, int nodeIndex)
{
    int         temp;

    temp                              = m->tiProbsIndex[chain][nodeIndex];
    m->tiProbsIndex[chain][nodeIndex] = m->tiProbsScratchIndex[nodeIndex];
    m->tiProbsScratchIndex[nodeIndex] = temp;
}


/*------------------------------------------------------------------
|
|   Likelihood_Adgamma: all n-state models with autocorrelated
|        discrete gamma rate variation, NOT morph, restriction,
|        codon or doublet models; just fill in rateProbs
|
-------------------------------------------------------------------*/
/* Find the correction group index for a character based on its model properties */
static int FindCorrGroup (ModelInfo *m, int c)
{
    int g;
    for (g = 0; g < m->numCorrGroups; g++)
        {
        if (m->nStates[c] == m->corrGroupNStates[g] &&
            m->cType[c] == m->corrGroupCType[g])
            {
            if (m->cType[c] == USERTYPE)
                {
                if (m->userTypeIdx[c] == m->corrGroupUserTypeIdx[g])
                    return g;
                }
            else
                return g;
            }
        }
    return 0;  /* fallback to first group */
}


int Likelihood_Std (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats)
{
    int             b, c, g, j, k, nBetaCats, nRateCats, nStates, numReps;
    YFlt          catLike, catFreq, rateFreq, like, *bs, *bsBase,
                    pUnobserved, pObserved;
    YFlt          groupPUnobs[50], groupPObs[50];
    CLFlt           *clPtr, **clP, *lnScaler, *nSitesOfPat;
    ModelInfo       *m;
    YFlt          *siteLnL;

    m = &modelSettings[division];
    siteLnL = m->siteLnL;  /* may be NULL if not requested */

    numReps=0;
    for (c=0; c<m->numChars; c++)
        {
        if (m->nStates[c] == 2)
            numReps += m->numBetaCats * 2;
        else
            numReps += m->nStates[c];
        }
    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += numReps;
        }
    
    /* find base frequencies */
    bsBase = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);

    /* find rate category number and frequencies */
    nRateCats = m->numRateCats;
    rateFreq = 1.0 / nRateCats;

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    *lnL = 0.0; /* reset lnL */

    /* Initialize per-group pUnobserved accumulators */
    for (g = 0; g < m->numCorrGroups; g++)
        groupPUnobs[g] = 0.0;

    /* --- Phase 1: Compute dummy character likelihoods (per correction group) --- */
    for (c=j=0; c<m->numDummyChars; c++)
        {
        like = 0.0;
        nStates = m->nStates[c];
        bs = bsBase + m->bsIndex[c];
        if (nStates == 2 && m->numBetaCats > 1)
            {
            nBetaCats = m->numBetaCats;
            catFreq = rateFreq / nBetaCats;
            }
        else
            {
            nBetaCats = 1;
            catFreq = rateFreq;
            }
        for (b=0; b<nBetaCats; b++)
            {
            for (k=0; k<nRateCats; k++)
                {
                catLike = 0.0;
                for (j=0; j<nStates; j++)
                    catLike += clP[k][j] * bs[j];
                like += catLike * catFreq;
                clP[k] += nStates;
                }
            bs += nStates;
            }

        /* Find which correction group this dummy belongs to */
        g = FindCorrGroup(m, c);
        groupPUnobs[g] += like * exp(lnScaler[c]);
        }

    /* Compute per-group pObserved */
    for (g = 0; g < m->numCorrGroups; g++)
        {
        groupPObs[g] = 1.0 - groupPUnobs[g];
        if (groupPObs[g] < LIKE_EPSILON)
            groupPObs[g] = LIKE_EPSILON;
        }

    /* Fallback for partitions with no correction groups (coding=all) */
    if (m->numCorrGroups == 0)
        {
        pObserved = 1.0;
        }

    /* --- Phase 2: Compute real character likelihoods with per-group correction --- */
    for (c=m->numDummyChars; c<m->numChars; c++)
        {
        like = 0.0;
        nStates = m->nStates[c];
        bs = bsBase + m->bsIndex[c];
        if (nStates == 2 && m->numBetaCats > 1)
            {
            nBetaCats = m->numBetaCats;
            catFreq = rateFreq / nBetaCats;
            }
        else
            {
            nBetaCats = 1;
            catFreq = rateFreq;
            }
        for (b=0; b<nBetaCats; b++)
            {
            for (k=0; k<nRateCats; k++)
                {
                catLike = 0.0;
                for (j=0; j<nStates; j++)
                    catLike += clP[k][j] * bs[j];
                like += catLike * catFreq;
                clP[k] += nStates;
                }
            bs += nStates;
            }
        /* check against LIKE_EPSILON (values close to zero are problematic) */
        if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
            YvyraPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
            (*lnL) = YFLT_NEG_MAX;
            abortMove = YES;
            return ERROR;
            }
        else
            {
            YFlt wt = hasEstimatedWeights ? GetCharWeightMultiplier(c + m->compCharStart) : 1.0;
            (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c] * wt;
            if (siteLnL != NULL && (c - m->numDummyChars) >= 0 && (c - m->numDummyChars) < (m->numChars - m->numDummyChars))
                siteLnL[c - m->numDummyChars] = lnScaler[c] + log(like);
            }

        /* Apply per-group ascertainment correction */
        if (m->numCorrGroups > 0)
            {
            YFlt wt = hasEstimatedWeights ? GetCharWeightMultiplier(c + m->compCharStart) : 1.0;
            g = FindCorrGroup(m, c);
            (*lnL) -= log(groupPObs[g]) * nSitesOfPat[c] * wt;
            }
        }

    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   Likelihood_Cont: likelihood for continuous traits
|
|   This function calculates the restricted maximum likelihood (REML)
|      using phylogenetic independent contrasts (PICs) (Felsenstein 1985).
|   These standardized contrasts are, under a BM model, both independent and identically distributed.
|   The “restricted” part of REML refers to the fact that it calculates likelihood based on a transformed
|      set of data where the effect of nuisance parameters (the root state in this case) has been removed.
|
-------------------------------------------------------------------*/
int Likelihood_Cont (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats)
{
    /* start calculating phylogenetic independent contrasts (PICs) and REML */

    /* reset log likelihood */
    (*lnL) = 0.0;

	//chi TODO
	
    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   Likelihood_Pars: likelihood under the Tuffley and Steel (1997)
|       model for characters with constant number of states. The idea
|       is described in:
|
|       Tuffley, C., and M. Steel. 1997. Links between maximum likelihood
|          and maximum parsimony under a simple model of site substitution.
|          Bull. Math. Bio. 59:581-607.
|
|       The likelihood under the Tuffley and Steel (1997) model is:
|       
|       L = k^[-(T + n)]
|      
|       where L is the likelihood
|             k is the number of character states
|             T is the parsimony tree length
|             n is the number of characters 
|
|   The parsimony calculator does not use character packing; this is
|       to enable reweighting of characters 
|
|   Note that this is an empirical Bayes approach in that it uses the
|       maximum likelihood branch length.
|
-------------------------------------------------------------------*/
int Likelihood_Pars (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats)
{
    int             c, i, nStates;
    BitsLong        *pL, *pR, *pP, *pA, x;
    CLFlt           nParsChars, treeLength;
    CLFlt           length, *nSitesOfPat, *newNodeLengthPtr, oldNodeLength;
    Tree            *t;
    ModelInfo       *m;

    /* Find model settings */
    m = &modelSettings[division];

    /* Get tree */
    t = GetTree(m->brlens,chain,state[chain]);
    
    /* Get parsimony tree length */
    treeLength = (CLFlt) m->parsTreeLength[2 * chain + state[chain]];
    
    /* Get number of states */
    nStates = m->numStates;

    /* Get number of sites of pat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;

    /* Now make downpass node by node */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];

        /* continue if no work needs to be done */
        if (p->upDateCl == NO)
            continue;

        /* flip space */
        FlipCondLikeSpace(m, chain, p->index);
        
        /* find parsimony sets for the node and its environment */
        pL    = m->parsSets[m->condLikeIndex[chain][p->left->index ]];
        pR    = m->parsSets[m->condLikeIndex[chain][p->right->index]];
        pP    = m->parsSets[m->condLikeIndex[chain][p->index       ]];
        // oldpP = m->parsSets[m->condLikeScratchIndex[p->index    ]];

        /* find old and new node lengths */
        oldNodeLength    =  m->parsNodeLens[m->condLikeScratchIndex[p->index]];
        newNodeLengthPtr = &m->parsNodeLens[m->condLikeIndex[chain][p->index]];
        
        if (t->isRooted == NO && p->anc->anc == NULL)
            {
            pA = m->parsSets[m->condLikeIndex[chain][p->anc->index]];
            length = 0.0;
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    length += nSitesOfPat[c];
                    }
                if ((x & pA[c]) == 0)
                    length += nSitesOfPat[c];
                pP[c] = x;
                }
            treeLength += (length - oldNodeLength);
            *newNodeLengthPtr = length;
            }
        else
            {
            length = 0.0;  // done = 0;
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    length += nSitesOfPat[c];
                    }
                pP[c] = x;
                }
            treeLength += (length - oldNodeLength);
            *newNodeLengthPtr = length;
            }
        }

    /* Count number of characters in the partition. It is calculated
       on the fly because this number is going to differ for
       different chains if character reweighting is used. */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    nParsChars = 0.0;
    for (c=0; c<m->numChars; c++)
        nParsChars += nSitesOfPat[c];

    /* Calculate likelihood from parsimony tree length */
    *lnL = - ((treeLength + nParsChars) *  log (nStates));

    /* Store current parsimony tree length */
    m->parsTreeLength[2 * chain + state[chain]] = treeLength;

    return (NO_ERROR);
}


/*------------------------------------------------------------------
|
|   Likelihood_Pars: likelihood under the Tuffley and Steel (1997)
|       model for characters with constant number of states.
|
|   This variant of Likelihood_Pars assumes that the number of states
|       is variable. It does not take state order into account.
|
-------------------------------------------------------------------*/
int Likelihood_ParsStd (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats)
{
    int             c, i, *nStates;
    BitsLong        *pL, *pR, *pP, *pA, x;
    CLFlt           *treeLength;
    CLFlt           *nSitesOfPat;
    Tree            *t;
    ModelInfo       *m;

    /* Find model settings */
    m = &modelSettings[division];

    /* Get tree */
    t = GetTree(m->brlens,chain,state[chain]);
    
    /* Allocate space for parsimony tree length */
    treeLength = (CLFlt *) SafeCalloc (m->numChars, sizeof (CLFlt));
    
    /* Get number of states */
    nStates = m->nStates;

    /* Get number of sites of pat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;

    /* Make downpass node by node; do not skip any nodes */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];

        /* flip space */
        FlipCondLikeSpace(m, chain, p->index);
        
        /* find parsimony sets for the node and its environment */
        pL    = m->parsSets[m->condLikeIndex[chain][p->left->index ]];
        pR    = m->parsSets[m->condLikeIndex[chain][p->right->index]];
        pP    = m->parsSets[m->condLikeIndex[chain][p->index       ]];

        if (t->isRooted == NO && p->anc->anc == NULL)
            {
            pA = m->parsSets[m->condLikeIndex[chain][p->anc->index]];
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    treeLength[c] += nSitesOfPat[c];
                    }
                if ((x & pA[c]) == 0)
                    treeLength[c] += nSitesOfPat[c];
                pP[c] = x;
                }
            }
        else
            {
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    treeLength[c] += nSitesOfPat[c];
                    }
                pP[c] = x;
                }
            }
        }

    /* Calculate the likelihood one character at a time */
    *lnL = 0.0;
    for (c=0; c<m->numChars; c++)
        {
        *lnL -= ((treeLength[c] + nSitesOfPat[c]) * log (nStates[c]));
        }

    /* Free space for parsimony character states */
    free (treeLength);

    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   LaunchLogLikeForDivision: calculate the log likelihood of the
|       new state of the chain for a single division
|
-----------------------------------------------------------------*/
void LaunchLogLikeForDivision(int chain, int d, YFlt* lnL)
{
    int i;
    TreeNode        *p;
    ModelInfo       *m;
    Tree            *tree;
#   if defined (TIMING_ANALIZ)
    clock_t         CPUTimeStart;
#   endif
    
    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);
    
    if (m->upDateCijk == YES)
        {
        if (UpDateCijk(d, chain)== ERROR)
            {
            (*lnL) = YFLT_NEG_MAX; /* effectively abort the move */
            return;
            }
        m->upDateAll = YES;
        }
    
        
    if (m->parsModelId == NO && m->dataType != CONTINUOUS)
        {
        /* get site scalers ready */
        FlipSiteScalerSpace(m, chain);
        if (m->upDateAll == YES)
            ResetSiteScalers(m, chain);
        else
            CopySiteScalers(m, chain);

        /* pass over tree */
        for (i=0; i<tree->nIntNodes; i++)
            {
            p = tree->intDownPass[i];
            
            if (p->left->upDateTi == YES)
                {
                /* shift state of ti probs for node */
                FlipTiProbsSpace (m, chain, p->left->index);
                m->TiProbs (p->left, d, chain);
                }
            
            if (p->right->upDateTi == YES)
                {
                /* shift state of ti probs for node */
                FlipTiProbsSpace (m, chain, p->right->index);
                m->TiProbs (p->right, d, chain);
                }
            
            if (tree->isRooted == NO)
                {
                if (p->anc->anc == NULL /* && p->upDateTi == YES */)
                    {
                    /* shift state of ti probs for node */
                    FlipTiProbsSpace (m, chain, p->index);
                    m->TiProbs (p, d, chain);
                    }
                }
            
            if (p->upDateCl == YES)
                {
                if (tree->isRooted == NO)
                    {
                    if (p->anc->anc == NULL)
                        {
                        TIME(m->CondLikeRoot (p, d, chain),CPUCondLikeRoot);
                        }
                    else
                        {
                        TIME(m->CondLikeDown (p, d, chain),CPUCondLikeDown);                        
                        }
                    }
                else
                    {
                    TIME(m->CondLikeDown (p, d, chain),CPUCondLikeDown);
                    }

                if (m->unscaledNodes[chain][p->index] == 0 && m->upDateAll == NO)
                    {
#if defined (SSE_ENABLED)
                    if (m->useVec == VEC_SSE)
                        {
                        TIME(RemoveNodeScalers_SSE (p, d, chain),CPUScalersRemove);
                        }
#if defined (AVX_ENABLED)
                    else if (m->useVec == VEC_AVX)
                        {
                        TIME(RemoveNodeScalers_AVX (p, d, chain),CPUScalersRemove);
                        }
#endif
                    else
                        {
                        TIME(RemoveNodeScalers (p, d, chain),CPUScalersRemove);
                        }
#   else
                    TIME(RemoveNodeScalers (p, d, chain),CPUScalersRemove);
#   endif
                    }
                FlipNodeScalerSpace (m, chain, p->index);
                m->unscaledNodes[chain][p->index] = 1 + m->unscaledNodes[chain][p->left->index] + m->unscaledNodes[chain][p->right->index];
                
                if (m->unscaledNodes[chain][p->index] >= m->rescaleFreq[chain] && p->anc->anc != NULL)
                    {
                    TIME(m->CondLikeScaler (p, d, chain),CPUScalers);
                    }
                }
            }
        }
    
    /* call likelihood function to summarize result */
    TIME(m->Likelihood (tree->root->left, d, chain, lnL, (chainId[chain] % chainParams.numChains)),CPULilklihood);
    return;
}


/*----------------------------------------------------------------
|
|   RemoveNodeScalers: Remove node scalers
|
-----------------------------------------------------------------*/
int RemoveNodeScalers (TreeNode *p, int division, int chain)
{
    int             c;
    CLFlt           *scP, *lnScaler;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    assert (m->unscaledNodes[chain][p->index] == 0);

    /* find scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* remove scalers */
    for (c=0; c<m->numChars; c++)
        lnScaler[c] -= scP[c];

    return NO_ERROR;
}


#if defined (AVX_ENABLED)
/*----------------------------------------------------------------
 |
 |   RemoveNodeScalers_AVX: Remove node scalers, AVX code
 |
 -----------------------------------------------------------------*/
int RemoveNodeScalers_AVX (TreeNode *p, int division, int chain)
{
    int             c;
    __m256          *scP_AVX, *lnScaler_AVX;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    assert (m->unscaledNodes[chain][p->index] == 0);
    
    /* find scalers */
    scP_AVX = (__m256*)(m->scalers[m->nodeScalerIndex[chain][p->index]]);
    
    /* find site scalers */
    lnScaler_AVX = (__m256*)(m->scalers[m->siteScalerIndex[chain]]);
    
    /* remove scalers */
    for (c=0; c<m->numVecChars; c++)
    {
        lnScaler_AVX[c] = _mm256_sub_ps(lnScaler_AVX[c], scP_AVX[c]);
    }
    
    return NO_ERROR;
    
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   RemoveNodeScalers_SSE: Remove node scalers, SSE code
|
-----------------------------------------------------------------*/
int RemoveNodeScalers_SSE (TreeNode *p, int division, int chain)
{
    int             c;
    __m128          *scP_SSE, *lnScaler_SSE;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    assert (m->unscaledNodes[chain][p->index] == 0);

    /* find scalers */
    scP_SSE = (__m128*)(m->scalers[m->nodeScalerIndex[chain][p->index]]);

    /* find site scalers */
    lnScaler_SSE = (__m128*)(m->scalers[m->siteScalerIndex[chain]]);

    /* remove scalers */
    for (c=0; c<m->numVecChars; c++)
        {
        lnScaler_SSE[c] = _mm_sub_ps(lnScaler_SSE[c], scP_SSE[c]);
        }

    return NO_ERROR;
    
}
#endif


/* CopySiteScalers: Copy site scalers from scratch space into current space */
void CopySiteScalers (ModelInfo *m, int chain)
{
    CLFlt       *from, *to;

    from = m->scalers[m->siteScalerScratchIndex];
    to   = m->scalers[m->siteScalerIndex[chain]];
    memcpy ((void*) to, (void*) from, (size_t)(m->numChars) * sizeof(CLFlt));
}


/*----------------------------------------------------------------------
 |
 |   ResetSiteScalers: Set log site scalers to 0.0.
 |
 ------------------------------------------------------------------------*/
void ResetSiteScalers (ModelInfo *m, int chain)
{
    int     c;
    CLFlt   *lnScaler;

    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    for (c=0; c<m->numChars; c++)
        lnScaler[c] = 0.0;
}
int SetStdQMatrix (YFlt **a, int nStates, YFlt *bs, int cType)
{
    register int    i, j;
    YFlt          scaler;

    /* This function sets up rate matrices for standard characters.
       For UNORD: symmetric rates proportional to stationary frequencies.
       For ORD: only adjacent states have non-zero rates.
       For USERTYPE: rates are taken from the user-defined rate matrix
       (passed via bs pointer, which is repurposed for this case). */

    /* set Q matrix to 0 */
    for (i=0; i<nStates; i++)
        for (j=0; j<nStates; j++)
            a[i][j] = 0.0;

    /* initialize Q matrix */
    scaler = 0.0;
    if (cType == UNORD)
        {
        /* unordered characters */
        for (i=0; i<nStates; i++)
            {
            for (j=0; j<nStates; j++)
                {
                if (i != j)
                    {
                    a[i][i] -= (a[i][j] = bs[j]);
                    scaler += bs[i] * a[i][j];
                    }
                }
            }
        }
    else if (cType == ORD)
        {
        /* ordered characters */
        for (i=0; i<nStates; i++)
            {
            for (j=0; j<nStates; j++)
                {
                if (abs(i - j) == 1)
                    {
                    a[i][i] -= (a[i][j] = bs[j]);
                    scaler += bs[i] * a[i][j];
                    }
                }
            }
        }
    else
        {
        /* This shouldn't be reached for USERTYPE (handled in TiProbs_Std) */
        YvyraPrint ("%s   ERROR: Unknown cType %d in SetStdQMatrix\n", spacer, cType);
        return (ERROR);
        }

    /* rescale Q matrix */
    for (i=0; i<nStates; i++)
        for (j=0; j<nStates; j++)
            a[i][j] /= scaler;

#   if defined DEBUG_SETSTDQMATRIX
    for (i=0; i<nStates; i++)
        {
        for (j=0; j<nStates; j++)
            printf ("%0.5lf ", a[i][j]);
        printf ("\n");
        }
#   endif

    return (NO_ERROR);
}
int TiProbs_Std (TreeNode *p, int division, int chain)
{
    int         b, c, i, j, k, n, s, nStates, index=0, index2;
    YFlt      v, eV1, eV2, eV3, eV4, eV5, *catRate,
                baseRate, theRate, pi, f1, f2, f3, f4, f5, f6, f7, root,
                *eigenValues, *cijk, sum, *bs, mu, length;
    CLFlt       pNoChange, pChange, *tiP;
    ModelInfo   *m;
#   if defined (DEBUG_TIPROBS_STD)
    int         index3;
#   endif

    m = &modelSettings[division];

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
    
    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;
    
    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->wnBranchRates != NULL)
        {
        length = GetParamSubVals (m->wnBranchRates, chain, state[chain])[p->index];
        }
    else if (m->ilnBranchRates != NULL)
        {
        length = GetParamSubVals (m->ilnBranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths, which might
       occur in relaxed clock models; an elegant solution would be to substitute the stationary
       probs and initial probs but for now we truncate lengths at small or large values TODO */
    if (length > BRLENS_MAX)
        length = BRLENS_MAX;
    else if (length < BRLENS_MIN)
        length = BRLENS_MIN;

    /* find base frequencies */
    bs = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);
    
    /* fill in values; this has to be done differently if state freqs are not equal
       or if usertype characters are present (which need eigendecomposition) */
    if (m->stateFreq->paramId == SYMPI_EQUAL)
        {
        /* equal state frequencies */
        /* fill in values for unordered characters */
        index = 0;
#   if defined (DEBUG_TIPROBS_STD)
        index3 = 0;
#   endif
        for (nStates=2; nStates<=MAX_STD_STATES; nStates++)
            {
            if (m->isTiNeeded[nStates-2] == NO)
                continue;
            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * baseRate * catRate[k];
                eV1 =  exp(-(nStates / (nStates -  1.0)) * v);
                pChange   = (CLFlt) ((1.0 / nStates) - ((1.0 / nStates) * eV1));
                pNoChange = (CLFlt) ((1.0 / nStates) + ((nStates - 1.0) / nStates) * eV1);
                if (pChange<0.0)
                    pChange = (CLFlt) 0.0;
                for (i=0; i<nStates; i++)
                    {
                    for (j=0; j<nStates; j++)
                        {
                        if (i == j)
                            tiP[index++] = pNoChange;
                        else
                            tiP[index++] = pChange;
                        }
                    }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* TODO: need a general algorithm for ordered characters */
        /* 3-state ordered character */
        if (m->isTiNeeded[MAX_STD_STATES-1] == YES)
            {
            nStates = 3;
            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                eV1 =  exp (-(3.0 / 4.0) * v);
                eV2 =  exp (-(9.0 / 4.0) * v);
                
                /* pij(0,0) */
                tiP[index] = (CLFlt) ((1.0 / 3.0) + (eV1 / 2.0) + (eV2 / 6.0));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+3] = (CLFlt) ((1.0 / 3.0) - (eV2 / 3.0));
                /* pij(0,2) */
                tiP[index+2] = (CLFlt) ((1.0 / 3.0) - (eV1 / 2.0) + (eV2 / 6.0));
                /* pij(1,1) */
                tiP[index+4] = (CLFlt) ((1.0 / 3.0) + (2.0 * eV2 / 3.0));
                
                /* fill in mirror part of matrix */
                index += 5;
                index2 = index - 2;
                for (i=0; i<4; i++)
                    tiP[index++] = tiP[index2--];

                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }

#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* 4-state ordered character */
        if (m->isTiNeeded[MAX_STD_STATES] == YES)
            {
            nStates = 4;
            pi = 1.0 / 4.0;
            root =  sqrt (2.0);
            f1 = root +  1.0;
            f2 = root -  1.0;

            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                eV1 =  1.0 / (exp ((4.0 * v) / 3.0));
                eV2 =  exp ((2.0 * (root - 2.0) * v) / 3.0) / root;
                eV3 =  1.0 / (root *  exp ((2.0 * (root + 2.0) * v) / 3.0));
                
                /* pij(0,0) */
                tiP[index] = (CLFlt) (pi * (1.0 + eV1 + (f1*eV2) + (f2*eV3)));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+4] = (CLFlt) (pi * (1.0 - eV1 + eV2 - eV3));
                /* pij(0,2) = tiP(1,3) */
                tiP[index+2] = tiP[index+7] = (CLFlt) (pi * (1.0 - eV1 - eV2 + eV3));
                /* pij(0,3) */
                tiP[index+3] = (CLFlt) (pi * (1.0 + eV1 - (f1*eV2) - (f2*eV3)));
                /* pij(1,1) */
                tiP[index+5] = (CLFlt) (pi * (1.0 + eV1 + (f2*eV2) + (f1*eV3)));
                /* pij(1,2) */
                tiP[index+6] = (CLFlt) (pi * (1.0 + eV1 - (f2*eV2) - (f1*eV3)));

                /* fill in mirror part of matrix */
                index += 8;
                index2 = index - 1;
                for (i=0; i<8; i++)
                    tiP[index++] = tiP[index2--];
        
                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* 5-state ordered character */
        if (m->isTiNeeded[MAX_STD_STATES+1] == YES)
            {
            nStates = 5;
            pi = 1.0 / 5.0;
            root =  sqrt (5.0);

            f5 = root /  4.0;
            f1 =  0.75 + f5;;
            f2 =  1.25 + f5;
            f3 =  1.25 - f5;
            f4 =  0.75 - f5;
            f5 = f5 *  2.0;
            f6 = f5 +  0.5;
            f7 = f5 -  0.5;

            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                v *=  5.0 /  16.0;

                eV1 =  exp ((root -  3.0) * v);
                eV2 =  exp (-(root +  3.0) * v);
                eV3 =  exp ((root -  5.0) * v);
                eV4 =  exp (-(root +  5.0) * v);

                /* pij(0,0) */
                tiP[index] = (CLFlt) (pi* (1.0 + (f1*eV3) + (f2*eV1) + (f3*eV2) + (f4*eV4)));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+5] = (CLFlt) (pi*(1.0 - (eV3/2.0) + (f5*eV1) - (f5*eV2) - (eV4/2.0)));
                /* pij(0,2) = pij(2,0) */
                tiP[index+2] = tiP[index+10] = (CLFlt) (pi*(1.0 - (f6*eV3) + (f7*eV4)));
                /* pij(0,3) = pij(1,4) */
                tiP[index+3] = tiP[index+9] = (CLFlt) (pi*(1.0 - (eV3/2.0) - (f5*eV1) + (f5*eV2) - (eV4/2.0)));
                /* pij(0,4) */
                tiP[index+4] = (CLFlt) (pi*(1.0 + (f1*eV3) - (f2*eV1) - (f3*eV2) + (f4*eV4)));
                /* pij(1,1) */
                tiP[index+6] = (CLFlt) (pi*(1.0 + (f4*eV3) + (f3*eV1) + (f2*eV2) + (f1*eV4)));
                /* pij(1,2) = pij(2,1) */
                tiP[index+7] = tiP[index+11] = (CLFlt) (pi*(1.0 + (f7*eV3) - (f6*eV4)));
                /* pij(1,3) */
                tiP[index+8] = (CLFlt) (pi*(1.0 + (f4*eV3) - (f3*eV1) - (f2*eV2) + (f1*eV4)));
                /* pij(2,2) */
                tiP[index+12] = (CLFlt) (pi*(1.0 + (2.0*eV3) + (2.0*eV4)));

                /* fill in mirror part of matrix */
                index += 13;
                index2 = index - 2;
                for (i=0; i<12; i++)
                    tiP[index++] = tiP[index2--];

                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* 6-state ordered character */
        if (m->isTiNeeded[MAX_STD_STATES+2] == YES)
            {
            nStates = 6;
            pi =  1.0 / 6.0;
            root =  sqrt (3.0);

            f4 = (3.0 / (2.0 * root));
            f1 =  1.0 + f4;
            f2 =  1.0 - f4;
            f3 =  0.5 + f4;
            f4 =  0.5 - f4;

            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                v /=  5.0;

                eV1 =  exp (-9 * v);
                eV2 =  exp (-6 * v);
                eV3 =  exp (-3 * v);
                eV4 =  exp (3.0 * (root - 2.0) * v);
                eV5 =  exp (-3.0 * (root + 2.0) * v);

                /* pij(0,0) */
                tiP[index] = (CLFlt) (pi* (1.0 + (0.5*eV1) + eV2 + (1.5*eV3) + (f1*eV4) + (f2*eV5)));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+6] = (CLFlt) (pi*(1.0 - eV1 - eV2 + (f3*eV4) + (f4*eV5)));
                /* pij(0,2) = pij(2,0) */
                tiP[index+2] = tiP[index+12] = (CLFlt) (pi*(1.0 + (0.5*eV1) - eV2 - (1.5*eV3) + (0.5*eV4) + (0.5*eV5)));
                /* pij(0,3) = pij(2,5) */
                tiP[index+3] = tiP[index+17] = (CLFlt) (pi*(1.0 + (0.5*eV1) + eV2 - (1.5*eV3) - (0.5*eV4) - (0.5*eV5)));
                /* pij(0,4) = pij(1,5) */
                tiP[index+4] = tiP[index+11] = (CLFlt) (pi*(1.0 - eV1 + eV2 - (f3*eV4) - (f4*eV5)));
                /* pij(0,5) */
                tiP[index+5] = (CLFlt) (pi*(1.0 + (0.5*eV1) - eV2 + (1.5*eV3) - (f1*eV4) - (f2*eV5)));
                /* pij(1,1) */
                tiP[index+7] = (CLFlt) (pi*(1.0 + (2.0*eV1) + eV2 + eV4 + eV5));
                /* pij(1,2) = pij(2,1) */
                tiP[index+8] = tiP[index+13] = (CLFlt) (pi*(1.0 - eV1 + eV2 - (f4*eV4) - (f3*eV5)));
                /* pij(1,3) = pij(2,4) */
                tiP[index+9] = tiP[index+16] = (CLFlt) (pi*(1.0 - eV1 - eV2 + (f4*eV4) + (f3*eV5)));
                /* pij(1,4) */
                tiP[index+10] = (CLFlt) (pi*(1.0 + (2.0*eV1) - eV2 - eV4 - eV5));
                /* pij(2,2) */
                tiP[index+14] = (CLFlt) (pi*(1.0 + (0.5*eV1) + eV2 + (1.5*eV3) + (f2*eV4) + (f1*eV5)));
                /* pij(2,3) */
                tiP[index+15] = (CLFlt) (pi*(1.0 + (0.5*eV1) - eV2 + (1.5*eV3) - (f2*eV4) - (f1*eV5)));

                /* fill in mirror part of matrix */
                index += 18;
                index2 = index - 1;
                for (i=0; i<18; i++)
                    tiP[index++] = tiP[index2--];

                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }
        }
    else
        {
        /* unequal state frequencies */
        index = 0;

        /* first fill in for binary characters using beta categories if needed */
        if (m->isTiNeeded[0] == YES)
            {
            /* cycle through beta and gamma cats */
            for (b=0; b<m->numBetaCats; b++)
                {
                mu =  1.0 / (2.0 * bs[0] * bs[1]);
                for (k=0; k<m->numRateCats; k++)
                    {
                    /* calculate probabilities */
                    v =  length * baseRate * catRate[k];
                    eV1 =  exp(- mu * v);
                    tiP[index++] = (CLFlt) (bs[0] + (bs[1] * eV1));
                    tiP[index++] = (CLFlt) (bs[1] - (bs[1] * eV1));
                    tiP[index++] = (CLFlt) (bs[0] - (bs[0] * eV1));
                    tiP[index++] = (CLFlt) (bs[1] + (bs[0] * eV1));
                    }
                /* update stationary state frequency pointer */
                bs += 2;
                }
            }

        /* now use general algorithm for the other cases */
        if (m->cijkLength > 0)
            {
            /* first update cijk if necessary */
            if (m->cijkLength > 0 && m->upDateCijk == YES)
                {
                if (UpDateCijk (division, chain) == ERROR)
                    return (ERROR);
                }

            /* then get first set of eigenvalues */
            eigenValues = m->cijks[m->cijkIndex[chain]];

            /* and cycle through the relevant characters */
            for (c=0; c<m->stateFreq->nSympi; c++)
                {
                n = m->stateFreq->sympinStates[c];

                /* skip USERTYPE chars (handled by direct matrix exp above) */
                if (m->stateFreq->sympiCType[c] == USERTYPE)
                    {
                    eigenValues += (n * n * n) + (2 * n);
                    continue;
                    }

                /* fill in values */
                for (k=0; k<m->numRateCats; k++)
                    {
                    v =  length * baseRate * catRate[k];
                    cijk = eigenValues + (2 * n);

                    for (i=0; i<n; i++)
                        {
                        for (j=0; j<n; j++)
                            {
                            sum = 0.0;
                            for (s=0; s<n; s++)
                                sum += (*cijk++) * exp(eigenValues[s] * v);
                            tiP[index++] = (CLFlt) ((sum <  0.0) ?  0.0 : sum);
                            }
                        }
                    }

                /* update eigenValues pointer */
                eigenValues += (n * n * n) + (2 * n);
                }
            }
        }

    /* handle USERTYPE characters via direct matrix exponentiation
       (runs after both the equal-freq and non-equal-freq paths above) */
    if (m->hasUserType == YES)
        {
        for (c=0; c<m->numChars; c++)
            {
            if (m->cType[c] != USERTYPE)
                continue;

            n = m->nStates[c];
            int utIdx = m->userTypeIdx[c];
            if (utIdx < 0 || utIdx >= numUserTypes)
                return (ERROR);
            UserType *ut = &userTypes[utIdx];

            /* build Q matrix from usertype rates */
            YFlt Q[MAX_STD_STATES][MAX_STD_STATES];
            YFlt *bsC = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]) + m->bsIndex[c];
            YFlt scaler = 0.0;
            for (i=0; i<n; i++)
                {
                Q[i][i] = 0.0;
                for (j=0; j<n; j++)
                    if (i != j)
                        {
                        Q[i][j] = ut->rates[i*n + j];
                        Q[i][i] -= Q[i][j];
                        scaler += bsC[i] * Q[i][j];
                        }
                }
            if (scaler > 0.0)
                for (i=0; i<n; i++)
                    for (j=0; j<n; j++)
                        Q[i][j] /= scaler;

            /* Write transition probs at tiIndex[c] for each rate cat */
            index = m->tiIndex[c];
            for (k=0; k<m->numRateCats; k++)
                {
                v = length * baseRate * catRate[k];

                /* exp(Q*v) via scaling and squaring + Taylor series */
                YFlt maxRate = 0.0;
                for (i=0; i<n; i++)
                    if (-Q[i][i] > maxRate) maxRate = -Q[i][i];
                int nSquare = 0;
                YFlt scaledV = v;
                while (maxRate * scaledV > 0.5)
                    { scaledV /= 2.0; nSquare++; }

                YFlt P[MAX_STD_STATES][MAX_STD_STATES];
                YFlt Tk[MAX_STD_STATES][MAX_STD_STATES];
                for (i=0; i<n; i++)
                    for (j=0; j<n; j++)
                        { P[i][j] = (i==j) ? 1.0 : 0.0; Tk[i][j] = P[i][j]; }
                for (int term=1; term<=20; term++)
                    {
                    YFlt Tnew[MAX_STD_STATES][MAX_STD_STATES];
                    for (i=0; i<n; i++)
                        for (j=0; j<n; j++)
                            { sum=0; for (s=0; s<n; s++) sum += Tk[i][s]*Q[s][j]*scaledV; Tnew[i][j]=sum/term; }
                    YFlt maxT = 0.0;
                    for (i=0; i<n; i++)
                        for (j=0; j<n; j++)
                            { Tk[i][j]=Tnew[i][j]; P[i][j]+=Tk[i][j]; if(fabs(Tk[i][j])>maxT) maxT=fabs(Tk[i][j]); }
                    if (maxT < 1e-15) break;
                    }
                for (int sq=0; sq<nSquare; sq++)
                    {
                    YFlt P2[MAX_STD_STATES][MAX_STD_STATES];
                    for (i=0; i<n; i++)
                        for (j=0; j<n; j++)
                            { sum=0; for (s=0; s<n; s++) sum += P[i][s]*P[s][j]; P2[i][j]=sum; }
                    for (i=0; i<n; i++)
                        for (j=0; j<n; j++) P[i][j] = P2[i][j];
                    }
                for (i=0; i<n; i++)
                    for (j=0; j<n; j++)
                        tiP[index++] = (CLFlt)((P[i][j] < 0.0) ? 0.0 : P[i][j]);
                }
            }
        }

    return NO_ERROR;
}


int UpDateCijk (int whichPart, int whichChain)
{
    int         c, i, j, k, n, n3, isComplex, sizeOfSingleCijk, cType, numQAllocated;
    YFlt      **q[100], **eigvecs, **inverseEigvecs;
    YFlt      *eigenValues, *eigvalsImag, *cijk;
    YFlt      *bs, *bsBase, *rateOmegaValues=NULL, rA=0.0, rS=0.0, posScaler, *omegaCatFreq=NULL;
    YComplex     **Ceigvecs, **CinverseEigvecs;
    ModelInfo   *m;
    Param       *p;

    /* get a pointer to the model settings for this partition */
    m = &modelSettings[whichPart];
    assert (m->upDateCijk == YES);
    
    /* we should only go through here if we have cijk information available for the partition */
    if (m->cijkLength > 0) 
        {
        /* flip cijk space */
        FlipCijkSpace(m, whichChain);
        
        /* Only STANDARD datatype supported */
        if (m->dataType == STANDARD)
            {
            /* set pointers and other stuff needed */
            numQAllocated = 1;
            p = m->stateFreq;
            eigenValues = m->cijks[m->cijkIndex[whichChain]];
            q[0] = AllocateSquareDoubleMatrix (MAX_STD_STATES);
            eigvecs = AllocateSquareDoubleMatrix (MAX_STD_STATES);
            inverseEigvecs = AllocateSquareDoubleMatrix (MAX_STD_STATES);
            Ceigvecs = AllocateSquareComplexMatrix (MAX_STD_STATES);
            CinverseEigvecs = AllocateSquareComplexMatrix (MAX_STD_STATES);
            bsBase = GetParamStdStateFreqs (m->stateFreq, whichChain, state[whichChain]);
            
            /* cycle over characters needing cijks */
            for (c=0; c<p->nSympi; c++)
                {
                n = p->sympinStates[c];
                bs = bsBase + p->sympiBsIndex[c];
                cType = p->sympiCType[c];
                n3 = n * n * n;
                eigvalsImag = eigenValues + n;
                cijk = eigenValues + (2 * n);
                if (cType == USERTYPE)
                    {
                    /* USERTYPE chars use direct matrix exp in TiProbs_Std;
                       skip eigendecomposition here (may have complex eigenvalues) */
                    eigenValues += (n3 + (2 * n));
                    continue;
                    }
                else if (SetStdQMatrix (q[0], n, bs, cType) == ERROR)
                    return (ERROR);
                isComplex = GetEigens (n, q[0], eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
                if (isComplex == NO)
                    {
                    CalcCijk (n, cijk, eigvecs, inverseEigvecs);
                    }
                else
                    {
                    if (isComplex == YES)
                        YvyraPrint ("%s   ERROR: Complex eigenvalues found!\n", spacer);
                    else
                        YvyraPrint ("%s   ERROR: Computing eigenvalues problem!\n", spacer);
                    goto errorExit;
                    }
                eigenValues += (n3 + (2 * n));
                }
            }
        for (k=0; k<numQAllocated; k++)
            FreeSquareDoubleMatrix (q[k]);
        FreeSquareDoubleMatrix (eigvecs);
        FreeSquareDoubleMatrix (inverseEigvecs);
        FreeSquareComplexMatrix (Ceigvecs);
        FreeSquareComplexMatrix (CinverseEigvecs);
        }
        
    return (NO_ERROR);

    errorExit:      
        for (k=0; k<numQAllocated; k++)
            FreeSquareDoubleMatrix (q[k]);
        FreeSquareDoubleMatrix (eigvecs);
        FreeSquareDoubleMatrix (inverseEigvecs);
        FreeSquareComplexMatrix (Ceigvecs);
        FreeSquareComplexMatrix (CinverseEigvecs);

        return ERROR;
}


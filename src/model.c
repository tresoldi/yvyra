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
/* best.h removed: species tree (BEST) not supported in yvyra */
#include "command.h"
#include "model.h"
#include "mcmc.h"
#include "proposal.h"
#include "sumpt.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

#undef  DEBUG_ADDDUMMYCHARS
#undef  DEBUG_CONSTRAINTS
#undef  DEBUG_COMPRESSDATA
#undef  DEBUG_CPP

/* local prototypes */
int       AddDummyChars (void);
void      AllocateCppEvents (Param *p);
MCMCMove *AllocateMove (MoveType *moveType, Param *param);
int       AllocateNormalParams (void);
int       AllocateTreeParams (void);
void      CheckCharCodingType (Matrix *m, CharInfo *ci);
int       CompressData (void);
int       InitializeChainTrees (Param *p, int from, int to, int isRooted);
int       FillBrlensSubParams (Param *param, int chn, int state);
void      FreeCppEvents (Param *p);
void      FreeMove (MCMCMove *mv);
void      GetPossibleRestrictionSites (int resSiteCode, int *sites);
int       GetUserTreeFromName (int *index, char *treeName);
void      InitializeMcmcTrees (Param *p);
int       IsApplicable (Param *param);
int       IsApplicable_FiveTaxaOrMore (Param *param);
int       IsApplicable_FourTaxaOrMore (Param *param);
int       IsApplicable_ThreeTaxaOrMore (Param *param);
int       IsApplicable_TreeAgeMove (Param *param);
int       IsModelSame (int whichParam, int part1, int part2, int *isApplic1, int *isApplic2);
int       LargestMovableSubtree(Param *treeParam);
int       NumActiveParts (void);
int       NumInformativeHardConstraints (ModelParams *mp);
int       NumNonExcludedChar (void);
int       NumStates (int part);
int       PrintCompMatrix (void);
int       PrintMatrix (void);
int       ProcessStdChars (RandLong *seed);
int       SetModelParams (void);
int       SetPopSizeParam (Param *param, int chn, int state, PolyTree *pt);
int       SetRelaxedClockParam (Param *param, int chn, int state, PolyTree *pt);
int       SetUpLinkTable (void);
int       ShowMoves (int used);
int       ShowParameters (int showStartVals, int showMoves, int showAllAvailable);
int       UpdateCppEvolLength (int *nEvents, YFlt **pos, YFlt **rateMult, YFlt *evolLength, TreeNode *p, YFlt baseRate);

/* globals */
int             *activeParams[NUM_LINKED];              /* a table holding the parameter link status        */
int             localOutGroup;                          /* outgroup for non-excluded taxa                   */
Calibration     **localTaxonCalibration = NULL;         /* stores local taxon calibrations (ages)           */
char            **localTaxonNames = NULL;               /* points to names of non-excluded taxa             */
Model           *modelParams;                           /* holds model params                               */
ModelInfo       *modelSettings;                         /* stores important info on model params            */
MCMCMove        **moves;                                /* vector of pointers to applicable moves           */
int             numApplicableMoves;                     /* number of moves applicable to model parameters   */
int             numCurrentDivisions;                    /* number of partitions of data                     */
int             numGlobalChains;                        /* number of global chains                          */
int             numLocalTaxa;                           /* number of non-excluded taxa                      */
int             numLocalChar;                           /* number of non-excluded characters                */
int             numParams;                              /* number of parameters in model                    */
int             numTopologies;                          /* number of topologies for one chain and state     */
int             numTrees;                               /* number of trees for one chain and state          */
Param           *params;                                /* vector of parameters                             */
Param           *printParams;                           /* vector of subst model parameters to print        */
ShowmovesParams showmovesParams;                        /* holds parameters for Showmoves command           */
Param           *treePrintparams;                       /* vector of tree parameters to print               */
int             setUpAnalysisSuccess;                   /* Set to YES if analysis is set without error      */

/* globals used to describe and change the current model; allocated in AllocCharacters and SetPartition  */
int         *numVars;                                   /* number of variables in setting arrays         */
int         *activeParts;                               /* partitions changes should apply to            */
int         *linkTable[NUM_LINKED];                     /* how parameters are linked across parts        */
int         *tempLinkUnlink[NUM_LINKED];                /* for changing parameter linkage                */
int         *tempLinkUnlinkVec;                         /* for changing parameter linkage                */
YFlt      *tempNum;                                   /* vector of numbers used for setting arrays     */

/* parser flags and variables only used in this file */
int         foundComma, foundBeta, foundAaSetting, modelIsFixed, linkNum, foundLeftPar, tempNumStates,
            foundBSNum[999], foundDSNum[999], foundFSNum[999], foundBSTime[999], foundDSTime[999], foundFSTime[999];
YFlt      tempStateFreqs[200], tempAaModelPrs[10];
char        colonPr[100], clockPr[30];

/* parser flags and variables shared with command.c */
extern int  fromI, toJ, foundDash, foundExp, foundEqual, isNegative;

/* other local variables (this file) */
YFlt          empiricalFreqs[200];         /* empirical base frequencies for partition                 */
int             intValsRowSize = 0;          /* row size of intValues matrix                            */
int             *intValues = NULL;           /* stores int values of chain parameters                   */
Tree            **mcmcTree;                  /* pointers to trees for mcmc                              */
int             paramValsRowSize = 0;        /* row size of paramValues matrix                          */
YFlt          *paramValues = NULL;         /* stores actual values and subvalues of chain parameters  */
int             *relevantParts = NULL;       /* partitions that are affected by this move               */
Param           **subParamPtrs;              /* pointer to subparams for topology params                */
int             *stateSize;                  /* # states for each compressed char                       */
// int          foundCurly;
// char         *plotTokenP;                 /* plotToken[CMD_STRING_LENGTH];*/


/*-----------------------------------------------------------------------
|
|   AddDummyChars: Add dummy characters to relevant partitions
|
------------------------------------------------------------------------*/
int AddDummyChars (void)
{
    int         i, j, k, d, numIncompatible, numDeleted, oldRowSize,
                newRowSize, numDummyChars, newColumn, newChar, oldColumn, oldChar, 
                isCompat, *tempChar, numIncompatibleChars;
    BitsLong    *tempMatrix, bitsLongOne = 1;
    CLFlt       *tempSitesOfPat;
    ModelInfo   *m;
    ModelParams *mp;
    CharInfo    cinfo;
    Matrix      matrix;

    /* set pointers to NULL */
    tempMatrix = NULL;
    tempSitesOfPat = NULL;
    tempChar = NULL;

    /* check how many dummy characters needed in total */
    numDummyChars = 0;  // numStdChars = 0;
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        mp = &modelParams[d];

        m->numDummyChars = 0;

        if (0) /* RESTRICTION removed */
            {
            if (mp->coding & NOABSENCESITES)
                m->numDummyChars++;
            if (mp->coding & NOPRESENCESITES)
                m->numDummyChars++;
            if (mp->coding & NOSINGLETONABSENCE)
                m->numDummyChars += numLocalTaxa;
            if (mp->coding & NOSINGLETONPRESENCE)
                m->numDummyChars += numLocalTaxa;
            }

        if (mp->dataType == STANDARD && !strcmp(mp->parsModel,"No"))
            {
            if (mp->coding & VARIABLE)
                {
                /* Determine distinct model groups from real characters.
                   Each group defined by (nStates, ctype, userTypeIndex).
                   Create nStates dummy chars per group for correct Mkv correction. */
                int g, nGroups = 0, oldCol;
                m->numCorrGroups = 0;

                for (oldCol = m->compMatrixStart; oldCol < m->compMatrixStop; oldCol++)
                    {
                    int oc = origChar[oldCol];
                    if (oc < 0) continue;  /* skip existing dummies */
                    int ns = charInfo[oc].numStates;
                    int ct = charInfo[oc].ctype;
                    int ut = (ct == USERTYPE) ? charInfo[oc].userTypeIndex : -1;

                    /* Check if this group already exists */
                    int found = NO;
                    for (g = 0; g < nGroups; g++)
                        {
                        if (m->corrGroupNStates[g] == ns &&
                            m->corrGroupCType[g] == ct &&
                            m->corrGroupUserTypeIdx[g] == ut)
                            { found = YES; break; }
                        }
                    if (found == NO && nGroups < 50)
                        {
                        m->corrGroupNStates[nGroups] = ns;
                        m->corrGroupCType[nGroups] = ct;
                        m->corrGroupUserTypeIdx[nGroups] = ut;
                        nGroups++;
                        }
                    }

                m->numCorrGroups = nGroups;
                for (g = 0; g < nGroups; g++)
                    m->numDummyChars += m->corrGroupNStates[g];
                }
            if (mp->coding & NOSINGLETONS)
                m->numDummyChars += 2*numLocalTaxa;
            }

        numDummyChars += m->numDummyChars;
        m->numChars += m->numDummyChars;
        }

    /* exit if dummy characters not needed */
    if (numDummyChars == 0)
        return NO_ERROR;

    /* print original compressed matrix */
#   if  defined (DEBUG_ADDDUMMYCHARS)
    YvyraPrint ("Compressed matrix before adding dummy characters...\n");
    PrintCompMatrix();
#   endif

    /* set row sizes for old and new matrices */
    oldRowSize = compMatrixRowSize;
    compMatrixRowSize += numDummyChars;
    newRowSize = compMatrixRowSize;
    numCompressedChars += numDummyChars;

    /* allocate space for new data */
    tempMatrix = (BitsLong *) SafeCalloc (numLocalTaxa * newRowSize, sizeof(BitsLong));
    tempSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    tempChar = (int *) SafeCalloc (compMatrixRowSize, sizeof(int));
    if (!tempMatrix || !tempSitesOfPat || !tempChar)
        {
        YvyraPrint ("%s   Problem allocating temporary variables in AddDummyChars\n", spacer);
        goto errorExit;
        }

    /* initialize indices */
    oldChar = newChar = newColumn = numDeleted = 0;

    /* set up matrix struct */
    matrix.origin = compMatrix;
    matrix.nRows = numLocalTaxa;
    matrix.rowSize = oldRowSize;

    /* loop over divisions */
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        mp = &modelParams[d];

        /* insert the dummy characters first for each division */
        if (m->numDummyChars > 0)
            {
            for (k=0; k<2; k++)
                {
                if (((mp->coding & NOSINGLETONPRESENCE) && k == 0) || ((mp->coding & NOSINGLETONABSENCE) && k == 1))
                    {
                    for (i=0; i< numLocalTaxa; i++)
                        {
                        for (j=0; j<numLocalTaxa; j++)
                            {
                            if (j == i)
                                tempMatrix[pos(j,newColumn,newRowSize)] = (bitsLongOne << k) ^ 3;
                            else
                                tempMatrix[pos(j,newColumn,newRowSize)] = bitsLongOne << k;
                            }
                        tempSitesOfPat[newChar] = 0;
                        tempChar[newColumn] = -1;
                        newChar++;
                        newColumn++;
                        }
                    }
                }

            /* Create per-model-group dummy characters for coding=variable.
               Each group gets nStates dummy patterns (one per constant state).
               Dummy char for state s has all taxa = (1 << s). */
            if (mp->coding & VARIABLE)
                {
                int g, s;
                for (g = 0; g < m->numCorrGroups; g++)
                    {
                    m->corrGroupDummyStart[g] = newChar - (newChar - newColumn + newColumn);
                    /* corrGroupDummyStart is relative to division start;
                       we record newChar offset for later use in ProcessStdChars */
                    m->corrGroupDummyStart[g] = newChar;
                    for (s = 0; s < m->corrGroupNStates[g]; s++)
                        {
                        for (i = 0; i < numLocalTaxa; i++)
                            tempMatrix[pos(i, newColumn, newRowSize)] = (BitsLong)1 << s;
                        tempSitesOfPat[newChar] = 0;
                        tempChar[newColumn] = -1;
                        newChar++;
                        newColumn++;
                        }
                    }
                }
            }

        /* add the normal characters */
        numIncompatible = numIncompatibleChars = 0;
        for (oldColumn=m->compMatrixStart; oldColumn<m->compMatrixStop; oldColumn++)
            {
            isCompat = YES;
            /* first check if the character is supposed to be present */
            if (m->numDummyChars > 0)
                {
                /* set up matrix struct */
                matrix.column = oldColumn;
                /* set up charinfo struct */
                cinfo.dType = mp->dataType;
                cinfo.cType = charInfo[origChar[oldChar]].ctype;
                cinfo.nStates = charInfo[origChar[oldChar]].numStates;
                CheckCharCodingType(&matrix, &cinfo);

                if (mp->coding & VARIABLE)
                    {
                    if((mp->coding & VARIABLE) == VARIABLE && cinfo.variable == NO)
                        {
                        isCompat = NO;
                        }
                    else if((mp->coding & NOABSENCESITES) && cinfo.constant[0] == YES)
                        {
                        isCompat = NO;
                        }
                    else if((mp->coding & NOPRESENCESITES) && cinfo.constant[1] == YES)
                        {
                        isCompat = NO;
                        }
                    }
                if (mp->coding & NOSINGLETONS)
                    {
                    if((mp->coding & NOSINGLETONS) == NOSINGLETONS && cinfo.informative == NO && cinfo.variable == YES)
                        {
                        isCompat = NO;
                        }
                    else if((mp->coding & NOSINGLETONABSENCE) && cinfo.singleton[0] == YES && cinfo.informative == NO)
                        {
                        isCompat = NO;
                        }
                    else if((mp->coding & NOSINGLETONPRESENCE) && cinfo.singleton[1] == YES && cinfo.informative == NO)
                        {
                        isCompat = NO;
                        }
                    }
                }

            if (isCompat == NO)
                {
                numIncompatible++;  // printf("%d\n",  origChar[oldChar]+1);
                numIncompatibleChars += (int) numSitesOfPat[oldChar];
                oldChar++;
                }
            else
                {
                /* add character */
                for (i=0; i<numLocalTaxa; i++)
                    tempMatrix[pos(i,newColumn,newRowSize)] = compMatrix[pos(i,oldColumn,oldRowSize)];
                /* set indices */
                compCharPos[origChar[oldColumn]] = newChar;
                compColPos[origChar[oldColumn]] = newColumn;
                tempSitesOfPat[newChar] = numSitesOfPat[oldChar];
                tempChar[newColumn] = origChar[oldColumn];
                newColumn++;
                if ((oldColumn-m->compMatrixStart+1) % m->nCharsPerSite == 0)
                    {
                    newChar++;
                    oldChar++;
                    }
                }
            }

        /* print a warning if there are incompatible characters */
        if (numIncompatible > 0)
            {
            m->numChars -= numIncompatible;
            m->numUncompressedChars -= numIncompatibleChars;
            numDeleted += numIncompatible;
            if (numIncompatibleChars > 1)
                {
                YvyraPrint ("%s   WARNING: There are %d characters incompatible with the specified\n", spacer, numIncompatibleChars);
                YvyraPrint ("%s            coding bias. These characters will be excluded.\n", spacer);
                }
            else
                {
                YvyraPrint ("%s   WARNING: There is one character incompatible with the specified\n", spacer);
                YvyraPrint ("%s            coding bias. This character will be excluded.\n", spacer);
                }
            }

        /* update division comp matrix and comp char pointers */
        m->compCharStop = newChar;
        m->compMatrixStop = newColumn;
        m->compCharStart = newChar - m->numChars;
        m->compMatrixStart = newColumn - m->nCharsPerSite * m->numChars;

        }   /* next division */

    /* compress matrix if necessary */
    if (numDeleted > 0)
        {
        for (i=k=0; i<numLocalTaxa; i++)
            {
            for (j=0; j<newRowSize-numDeleted; j++)
                {
                tempMatrix[k++] = tempMatrix[j+i*newRowSize];
                }
            }
        numCompressedChars -= numDeleted;
        compMatrixRowSize -= numDeleted;
        }

    /* free old data, set pointers to new data */
    free (compMatrix);
    free (numSitesOfPat);
    free (origChar);
    
    compMatrix = tempMatrix;
    numSitesOfPat = tempSitesOfPat;
    origChar = tempChar;
    
    tempMatrix = NULL;
    tempSitesOfPat = NULL;
    tempChar = NULL;
    
    /* print new compressed matrix */
#   if  defined (DEBUG_ADDDUMMYCHARS)
    YvyraPrint ("After adding dummy characters...\n");
    PrintCompMatrix();
#   endif

    return NO_ERROR;

    errorExit:
        if (tempMatrix)
            free (tempMatrix);
        if (tempSitesOfPat)
            free (tempSitesOfPat);
        if (tempChar)
            free (tempChar);

        return ERROR;   
}


/* Allocate space for cpp events */
void AllocateCppEvents (Param *p)
{
    int     i;

    p->nEvents = (int **) SafeCalloc (2*numGlobalChains, sizeof (int *));
    p->nEvents[0] = (int *) SafeCalloc (2*numGlobalChains*(2*numLocalTaxa), sizeof (int));
    for (i=1; i<2*numGlobalChains; i++)
        p->nEvents[i] = p->nEvents[i-1] + (2*numLocalTaxa);
    p->position = (YFlt ***) SafeCalloc (2*numGlobalChains, sizeof (YFlt **));
    p->position[0] = (YFlt **) SafeCalloc (2*numGlobalChains*(2*numLocalTaxa), sizeof (YFlt *));
    for (i=1; i<2*numGlobalChains; i++)
        p->position[i] = p->position[i-1] + (2*numLocalTaxa);
    p->rateMult = (YFlt ***) SafeCalloc (2*numGlobalChains, sizeof (YFlt **));
    p->rateMult[0] = (YFlt **) SafeCalloc (2*numGlobalChains*(2*numLocalTaxa), sizeof (YFlt *));
    for (i=1; i<2*numGlobalChains; i++)
        p->rateMult[i] = p->rateMult[i-1] + (2*numLocalTaxa);

}


/*----------------------------------------------------------------------------
|
|   AllocateMove: Allocate space for and initialize one applicable move
|
-----------------------------------------------------------------------------*/
MCMCMove *AllocateMove (MoveType *moveType, Param *param)
{
    int         i, j, nameLength;
    char        *partitionDescriptor = "";
    MCMCMove    *temp;

    if ((temp = (MCMCMove *) SafeCalloc (1, sizeof (MCMCMove))) == NULL)
        return (NULL);

    /* calculate length */
    if (strcmp (moveType->paramName, "") == 0)
        nameLength = (int) (strlen (moveType->shortName) + strlen (param->name)) + 10;
    else
        {
        partitionDescriptor = param->name;
        while (*partitionDescriptor != '\0')
            {
            if (*partitionDescriptor == '{')
                break;
            partitionDescriptor++;
            }
        nameLength = (int) (strlen (moveType->shortName) + strlen (moveType->paramName) + strlen (partitionDescriptor)) + 10;
        }
    /* add length of names of subparams */
    if (moveType->subParams == YES)
        {
        for (i=0; i<param->nSubParams; i++)
            nameLength += (int)(strlen(param->subParams[i]->name)) + 1;
        }

    if ((temp->name = (char *) SafeCalloc (nameLength, sizeof (char))) == NULL)
        {
        free (temp);
        return NULL;
        }

    if ((temp->nAccepted = (int *) SafeCalloc (5*numGlobalChains, sizeof (int))) == NULL)
        {
        free (temp->name);
        free (temp);
        return NULL;
        }
    temp->nTried       = temp->nAccepted + numGlobalChains;
    temp->nBatches     = temp->nAccepted + 2*numGlobalChains;
    temp->nTotAccepted = temp->nAccepted + 3*numGlobalChains;
    temp->nTotTried    = temp->nAccepted + 4*numGlobalChains; 
    
    if ((temp->relProposalProb = (YFlt *) SafeCalloc (4*numGlobalChains, sizeof (YFlt))) == NULL)
        {
        free (temp->nAccepted);
        free (temp->name);
        free (temp);
        return NULL;
        }
    temp->cumProposalProb = temp->relProposalProb + numGlobalChains;
    temp->targetRate = temp->relProposalProb + 2*numGlobalChains;
    temp->lastAcceptanceRate = temp->relProposalProb + 3*numGlobalChains;

    if ((temp->tuningParam = (YFlt **) SafeCalloc (numGlobalChains, sizeof (YFlt *))) == NULL)
        {
        free (temp->relProposalProb);
        free (temp->nAccepted);
        free (temp->name);
        free (temp);
        return NULL;
        }
    if (moveType->numTuningParams > 0 && (temp->tuningParam[0] = (YFlt *) SafeCalloc (moveType->numTuningParams*numGlobalChains, sizeof (YFlt))) == NULL)
        {
        free (temp->tuningParam);
        free (temp->relProposalProb);
        free (temp->nAccepted);
        free (temp->name);
        free (temp);
        return NULL;
        }
    for (i=1; i<numGlobalChains; i++)
        temp->tuningParam[i] = temp->tuningParam[0] + (i * moveType->numTuningParams);

    /* set default values */
    if (strcmp(moveType->paramName, "") != 0)
        sprintf (temp->name, "%s(%s%s)", moveType->shortName, moveType->paramName, partitionDescriptor);
    else
        {
        sprintf (temp->name, "%s(%s", moveType->shortName, param->name);
        if (moveType->subParams == YES)
            {
            for (i=0; i<param->nSubParams; i++)
                {
                strcat(temp->name,",");
                strcat(temp->name,param->subParams[i]->name);
                }
            }
        strcat (temp->name,")");
        }
        
    temp->moveType = moveType;
    temp->moveFxn = moveType->moveFxn;
    for (i=0; i<numGlobalChains; i++)
        {
        temp->relProposalProb[i] = moveType->relProposalProb;
        temp->cumProposalProb[i] = 0.0;
        temp->nAccepted[i] = 0;
        temp->nTried[i] = 0;
        temp->nBatches[i] = 0;
        temp->nTotAccepted[i] = 0;
        temp->nTotTried[i] = 0;
        temp->targetRate[i] = moveType->targetRate;
        temp->lastAcceptanceRate[i] = 0.0;
        for (j=0; j<moveType->numTuningParams; j++)
            temp->tuningParam[i][j] = moveType->tuningParam[j];
        }
    return (temp);
}


/*----------------------------------------------------------------------
|
|   AllocateNormalParams: Allocate space for normal parameters
|
-----------------------------------------------------------------------*/
int AllocateNormalParams (void)
{
    int         i, k, nOfParams, nOfIntParams;
    Param       *p;
    
    /* Count the number of param values and subvalues */
    nOfParams = 0;
    nOfIntParams = 0;
    for (k=0; k<numParams; k++)
        {
        nOfParams += params[k].nValues;
        nOfParams += params[k].nSubValues;
        nOfIntParams += params[k].nIntValues;
        }

    /* Set row size and find total number of values */
    paramValsRowSize = nOfParams;
    intValsRowSize = nOfIntParams;
    nOfParams *= (2 * numGlobalChains);
    nOfIntParams *= (2 * numGlobalChains);

    if (memAllocs[ALLOC_PARAMVALUES] == YES)
        {
        if (nOfParams > 1)
            {
            paramValues = (YFlt *) SafeRealloc ((void *) paramValues, nOfParams * sizeof (YFlt));
            for (i=0; i<nOfParams; i++)
                paramValues[i] = 0.0;
            }
        if (nOfIntParams > 0)
            intValues = (int *) SafeRealloc ((void *) intValues, nOfIntParams * sizeof(int));
        }
    else
        {
        if (nOfParams > 1)
            paramValues = (YFlt *) SafeCalloc (nOfParams, sizeof(YFlt));
        else
            paramValues = NULL;
        if (nOfIntParams > 0)
            intValues = (int *) SafeCalloc (nOfIntParams, sizeof(int));
        else
            intValues = NULL;
        }
    if ((nOfParams > 0 && !paramValues) || (nOfIntParams > 0 && !intValues))
        {
        YvyraPrint ("%s   Problem allocating paramValues\n", spacer);
        if (paramValues)
            free (paramValues);
        if (intValues)
            free (intValues);
        return ERROR;
        }
    else
        memAllocs[ALLOC_PARAMVALUES] = YES;

    /* set pointers to values for chain 1 state 0            */
    /* this scheme keeps the chain and state values together */
    nOfParams = 0;
    nOfIntParams = 0;
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->nValues > 0)
            p->values = paramValues + nOfParams;
        else
            p->values = NULL;
        nOfParams += p->nValues;
        if (p->nSubValues > 0)
            p->subValues = paramValues + nOfParams;
        else
            p->subValues = NULL;
        nOfParams += p->nSubValues;
        if (p->nIntValues > 0)
            p->intValues = intValues + nOfIntParams;
        else
            p->intValues = NULL;
        nOfIntParams += p->nIntValues;
        }
    
    /* allocate space for cpp events */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_CPPEVENTS)
           AllocateCppEvents(p);
        }

    return (NO_ERROR);
}


/*----------------------------------------------------------------------
|
|   AllocateTreeParams: Allocate space for tree parameters
|
-----------------------------------------------------------------------*/
int AllocateTreeParams (void)
{
    int         i, j, k, n, nOfParams, nOfTrees, isRooted, numSubParamPtrs, allModelsStationary;
    Param       *p, *q;

    /* Count the number of trees and dated trees */
    /* based on branch length parameters */
    /*  One tree is needed for each brlen parameter.
        A topology may apply to several trees; a topology parameter
        contains pointers to all trees it applies to */
    numTrees = 0;
    numTopologies = 0;
    for (k=0; k<numParams; k++)
        {
        if (params[k].paramType == P_BRLENS)
            numTrees++;
        else if (params[k].paramType == P_TOPOLOGY)
            numTopologies++;
        }

    /* We need to add the trees that do not have any branch lengths */
    /* that is, the pure parsimony model trees, and the species trees */
    for (k=0; k<numParams; k++)
        {
        if (params[k].paramType == P_TOPOLOGY)
            {
            if (!strcmp(modelParams[params[k].relParts[0]].parsModel, "Yes"))
                numTrees++;
            }
        else if (params[k].paramType == P_SPECIESTREE)
            {
            numTopologies++;
            numTrees++;
            }
        }

    /* Finally add subparam pointers for relaxed clock parameters and species trees */
    numSubParamPtrs = 0;
    for (k=0; k<numParams; k++)
        {
        if (params[k].paramType == P_TOPOLOGY && params[k].paramId == TOPOLOGY_SPECIESTREE)
            numSubParamPtrs += 1;
        else if (params[k].paramType == P_TOPOLOGY && !strcmp(modelParams[params[k].relParts[0]].parsModel, "Yes"))
            numSubParamPtrs += 1;
        else if (params[k].paramType == P_BRLENS)
            numSubParamPtrs += 1;
        else if (params[k].paramType == P_CPPEVENTS)
            numSubParamPtrs += 3;
        else if (params[k].paramType == P_TK02BRANCHRATES)
            numSubParamPtrs += 2;
        else if (params[k].paramType == P_WNBRANCHRATES)
            numSubParamPtrs += 2;
        else if (params[k].paramType == P_IGRBRANCHRATES)
            numSubParamPtrs += 2;
        else if (params[k].paramType == P_ILNBRANCHRATES)
            numSubParamPtrs += 2;
        else if (params[k].paramType == P_MIXEDBRCHRATES)
            numSubParamPtrs += 2;
        }
        
    /* Allocate space for trees and subparam pointers */
    if (memAllocs[ALLOC_MCMCTREES] == YES)
        {
        free (subParamPtrs);
        free (mcmcTree);
        subParamPtrs = NULL;
        mcmcTree = NULL;
        memAllocs[ALLOC_MCMCTREES] = NO;
        }
    if (numSubParamPtrs > 0)
        subParamPtrs = (Param **) SafeCalloc (numSubParamPtrs, sizeof (Param *));
    mcmcTree = (Tree **) SafeCalloc (numTrees * 2 * numGlobalChains, sizeof (Tree *));
    if ((numSubParamPtrs>0 && !subParamPtrs) || !mcmcTree)
        {
        if (subParamPtrs) free (subParamPtrs);
        if (mcmcTree) free (mcmcTree);
        subParamPtrs = NULL;
        mcmcTree = NULL;
        YvyraPrint ("%s   Problem allocating MCMC trees\n", spacer);
        printf("subparams: %d -- trees: %d \n", numSubParamPtrs, numTrees);
        return (ERROR);
        }
    else
        memAllocs[ALLOC_MCMCTREES] = YES;

    /* Initialize number of subparams, just in case */
    for (k=0; k<numParams; k++)
        params[k].nSubParams = 0;
    
    /* Count number of trees (brlens) for each topology or species tree */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (params[k].paramType == P_BRLENS)
            {
            q = modelSettings[p->relParts[0]].topology;
            q->nSubParams++;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        else if (params[k].paramType == P_TOPOLOGY)
            {
            q = modelSettings[p->relParts[0]].speciesTree;
            if (q != NULL)
                q->nSubParams++;
            }
        }

    /* Make sure there is also one subparam for a parsimony tree */
    for (k=0; k<numParams; k++)
        if (params[k].paramType == P_TOPOLOGY)
            {
            p = &params[k];
            if (p->nSubParams == 0)
                p->nSubParams = 1;
            }

    /* Count subparams for relaxed clock parameters */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_CPPEVENTS)
            {
            q = modelSettings[p->relParts[0]].cppRate;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].cppMultDev;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].brlens;
            q->nSubParams++;
            }
        else if (p->paramType == P_TK02BRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].tk02var;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].brlens;
            q->nSubParams++;
            }
        else if (p->paramType == P_WNBRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].wnvar;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].brlens;
            q->nSubParams++;
            }
        else if (p->paramType == P_IGRBRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].igrvar;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].brlens;
            q->nSubParams++;
            }
        else if (p->paramType == P_ILNBRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].ilnvar;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].brlens;
            q->nSubParams++;
            }
        else if (p->paramType == P_MIXEDBRCHRATES)
            {
            q = modelSettings[p->relParts[0]].mixedvar;
            q->nSubParams++;
            q = modelSettings[p->relParts[0]].brlens;
            q->nSubParams++;
            }
        }

    /* set pointers to subparams */
    nOfParams = 0;
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->nSubParams > 0)
            {
            p->subParams = subParamPtrs + nOfParams;
            nOfParams += p->nSubParams;
            }
        }
    assert (nOfParams == numSubParamPtrs);

    /* Set brlens param pointers and tree values */
    /* the scheme below keeps trees for the same state and chain together */
    nOfTrees = 0;
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if ((p->paramType == P_BRLENS) ||
            (p->paramType == P_TOPOLOGY && p->paramId == TOPOLOGY_PARSIMONY_UNIFORM) ||
            (p->paramType == P_TOPOLOGY && p->paramId == TOPOLOGY_PARSIMONY_CONSTRAINED) ||
            (p->paramType == P_SPECIESTREE))
            {
            /* allocate space for trees and initialize trees */
            p->treeIndex = nOfTrees;
            p->tree = mcmcTree + nOfTrees;
            nOfTrees++;
            }
        }

    /* Set topology params and associated brlen subparams */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_TOPOLOGY)
            {
            if (p->paramId == TOPOLOGY_PARSIMONY_UNIFORM ||
                p->paramId == TOPOLOGY_PARSIMONY_CONSTRAINED)
                /* pure parsimony topology case */
                {
                /* there is no brlen subparam */
                /* so let subparam point to the param itself */
                q = p->subParams[0] = p;
                /* p->tree and p->treeIndex have been set above */
                }
            else
                {
                /* first set brlens pointers for any parsimony partitions */
                for (i=j=0; i<p->nRelParts; i++)
                    {
                    if (modelSettings[p->relParts[i]].parsModelId == YES)
                        {
                        modelSettings[p->relParts[i]].brlens = p;
                        }
                    }

                /* now proceed with pointer assignment */
                q = modelSettings[p->relParts[0]].brlens;
                n = 0;  /* number of stored subParams */
                i = 0;  /* relevant partition number  */
                while (i < p->nRelParts)
                    {
                    for (j=0; j<n; j++)
                        if (q == p->subParams[j])
                            break;
                    
                    if (j == n && q != p)   /* a new tree (brlens) for this topology */
                        {
                        p->subParams[n++] = q;
                        }
                    q = modelSettings[p->relParts[++i]].brlens;
                    }
                
                p->tree = p->subParams[0]->tree;
                p->treeIndex = p->subParams[0]->treeIndex;
                }
            }
        else if (p->paramType == P_SPECIESTREE)
            {
            /* now proceed with pointer assignment */
            q = modelSettings[p->relParts[0]].topology;
            n = 0;  /* number of stored subParams */
            i = 0;  /* relevant partition number  */
            while (i < p->nRelParts)
                {
                for (j=0; j<n; j++)
                    if (q == p->subParams[j])
                        break;

                if (j == n && q != p)   /* a new topology for this species tree */
                    {
                    p->subParams[n++] = q;
                    }
                q = modelSettings[p->relParts[++i]].topology;
                }
            }
        }

    /* Check for constraints */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_TOPOLOGY)
            {
            if (!strcmp(modelParams[p->relParts[0]].topologyPr, "Constraints"))
                {
                for (i=0; i<p->nSubParams; i++)
                    {
                    q = p->subParams[i];
                    q->checkConstraints = YES;
                    }
                }
            else
                {
                for (i=0; i<p->nSubParams; i++)
                    {
                    q = p->subParams[i];
                    q->checkConstraints = NO;
                    }
                }
            }
        }

    /* update paramId */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_TOPOLOGY)
            {
            if (p->nSubParams > 1)
                {
                if (p->paramId == TOPOLOGY_NCL_UNIFORM_HOMO)
                    p->paramId = TOPOLOGY_NCL_UNIFORM_HETERO;
                else if (p->paramId == TOPOLOGY_NCL_CONSTRAINED_HOMO)
                    p->paramId = TOPOLOGY_NCL_CONSTRAINED_HETERO;
                else if (p->paramId == TOPOLOGY_NCL_FIXED_HOMO)
                    p->paramId = TOPOLOGY_NCL_FIXED_HETERO;
                else
                    {
                    YvyraPrint ("%s   A clock tree cannot have more than one set of branch lengths\n", spacer);
                    printf ("nparam:%d paramid:%d",p->nSubParams,p->paramId);
                    return (ERROR);
                    }
                }
            }
        }

    /* finally initialize trees */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_BRLENS)
            {
            /* find type of tree */
            allModelsStationary = YES;
            for (i=0; i<p->nRelParts; i++)
                {
                if (strcmp(modelParams[p->relParts[i]].statefreqModel, "Stationary"      ))
                    {
                     allModelsStationary = NO;
                     break;
                    }
                }
            
            if (!strcmp(modelParams[p->relParts[0]].brlensPr,"Clock") || allModelsStationary == NO)
                isRooted = YES;
            else
                isRooted = NO;

            if (InitializeChainTrees (p, 0, numGlobalChains, isRooted) == ERROR)
                return (ERROR);
            }
        else if (p->paramType == P_SPECIESTREE)
            {
            if (InitializeChainTrees (p, 0, numGlobalChains, YES) == ERROR)
                return (ERROR);
            }
        else if (p->paramType == P_TOPOLOGY && p->subParams[0]==p)
            {
            if (InitializeChainTrees (p, 0, numGlobalChains, NO) == ERROR)
                return (ERROR);
            }
        }

    /* now initialize subparam pointers for relaxed clock models */
    /* use nSubParams to point to the next available subParam by first
       resetting all nSubParams to 0 */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_CPPRATE ||
            p->paramType == P_CPPMULTDEV ||
            p->paramType == P_BRLENS ||
            p->paramType == P_TK02VAR ||
            p->paramType == P_WNVAR ||
            p->paramType == P_IGRVAR ||
            p->paramType == P_ILNVAR ||
            p->paramType == P_MIXEDVAR)
            p->nSubParams = 0;
        }
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType == P_CPPEVENTS)
            {
            q = modelSettings[p->relParts[0]].cppRate;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].cppMultDev;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].brlens;
            q->subParams[q->nSubParams++] = p;
            p->treeIndex = q->treeIndex;
            p->tree = q->tree;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        else if (p->paramType == P_TK02BRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].tk02var;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].brlens;
            q->subParams[q->nSubParams++] = p;
            p->treeIndex = q->treeIndex;
            p->tree = q->tree;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        else if (p->paramType == P_WNBRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].wnvar;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].brlens;
            q->subParams[q->nSubParams++] = p;
            p->treeIndex = q->treeIndex;
            p->tree = q->tree;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        else if (p->paramType == P_IGRBRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].igrvar;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].brlens;
            q->subParams[q->nSubParams++] = p;
            p->treeIndex = q->treeIndex;
            p->tree = q->tree;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        else if (p->paramType == P_ILNBRANCHRATES)
            {
            q = modelSettings[p->relParts[0]].ilnvar;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].brlens;
            q->subParams[q->nSubParams++] = p;
            p->treeIndex = q->treeIndex;
            p->tree = q->tree;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        else if (p->paramType == P_MIXEDBRCHRATES)
            {
            q = modelSettings[p->relParts[0]].mixedvar;
            q->subParams[q->nSubParams++] = p;
            q = modelSettings[p->relParts[0]].brlens;
            q->subParams[q->nSubParams++] = p;
            p->treeIndex = q->treeIndex;
            p->tree = q->tree;
            if (p->printParam == YES)
                q->nPrintSubParams++;
            }
        }

    return (NO_ERROR);
}


int AreDoublesEqual (YFlt x, YFlt y, YFlt tol)
{
    if ((x - y) < -tol || (x - y) > tol)
        return (NO);
    else
        return (YES);
}


int ChangeNumChains (int from, int to)
{
    int         i, i1, j, k, m, st, nRuns, fromIndex, toIndex, run, chn, *tempIntVals, nCppEventParams, *toEvents, *fromEvents;
    MCMCMove    **tempMoves, *fromMove, *toMove;
    Tree        **tempTrees;
    YFlt      *tempVals, **toRateMult, **toPosition, **fromRateMult, **fromPosition, *stdStateFreqsOld;
    Param       *p, *q, *cppEventParams = NULL;
    Tree        **oldMcmcTree, *tree;

    if (from == to)
        return (NO_ERROR);

    /* set new number of chains */
    chainParams.numChains = to;
    nRuns = chainParams.numRuns;
    numGlobalChains = chainParams.numRuns * chainParams.numChains;

    /* Do the normal parameters */  
    /* first save old values */
    tempVals = paramValues;
    paramValues = NULL;
    tempIntVals = intValues;
    intValues = NULL;
    memAllocs[ALLOC_PARAMS] = NO;
    /* .. and old cpp events parameters */
    nCppEventParams = 0;
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_CPPEVENTS)
            nCppEventParams++;
        }
    if (nCppEventParams > 0)
        {
        cppEventParams = (Param *) SafeCalloc ((size_t)(nCppEventParams), sizeof(Param));
        for (i=0; i<nCppEventParams; i++)
            {
            cppEventParams[i].paramType = P_CPPEVENTS;
            AllocateCppEvents (&cppEventParams[i]);
            }
        for (i=j=0; i<numParams; i++)
            {
            p = &params[i];
            if (p->paramType == P_CPPEVENTS)
                {
                cppEventParams[j].nEvents = p->nEvents;
                p->nEvents = NULL;
                cppEventParams[j].position = p->position;
                p->position = NULL;
                cppEventParams[j].rateMult = p->rateMult;
                p->rateMult = NULL;
                j++;
                }
            }
        }
    if (AllocateNormalParams () == ERROR)
        return (ERROR);

    /* then fill all params */
    FillNormalParams (&globalSeed, 0, numGlobalChains);
    
    /* finally overwrite with old values if present */
    for (run=0; run<nRuns; run++)
        {
        for (chn=0; chn<from; chn++)
            {
            if (chn < to)
                {
                fromIndex = (run*from + chn)*2*paramValsRowSize;
                toIndex = (run*to + chn)*2*paramValsRowSize;
                for (i=0; i<2*paramValsRowSize; i++)
                    paramValues[toIndex++] = tempVals[fromIndex++];
                fromIndex = (run*from + chn)*2*intValsRowSize;
                toIndex = (run*to + chn)*2*intValsRowSize;
                for (i=0; i<2*intValsRowSize; i++)
                    intValues[toIndex++] = tempIntVals[fromIndex++];
                for (i=i1=0; i<numParams; i++)
                    {
                    p = &params[i];
                    if (p->paramType == P_CPPEVENTS)
                        {
                        fromIndex = 2*(run*from + chn);
                        toIndex = 2*(run*to + chn);
                        fromEvents = cppEventParams[i1].nEvents[fromIndex];
                        toEvents = p->nEvents[toIndex];
                        fromPosition = cppEventParams[i1].position[fromIndex];
                        toPosition = p->position[toIndex];
                        fromRateMult = cppEventParams[i1].rateMult[fromIndex];
                        toRateMult = p->rateMult[toIndex];
                        for (j=0; j<2*numLocalTaxa; j++)
                            {
                            toEvents[j] = fromEvents[j];
                            toPosition[j] = (YFlt *) SafeRealloc ((void *)toPosition[j], toEvents[j]*sizeof(YFlt));
                            toRateMult[j] = (YFlt *) SafeRealloc ((void *)toRateMult[j], toEvents[j]*sizeof(YFlt));
                            for (k=0; k<toEvents[j]; k++)
                                {
                                toPosition[j][k] = fromPosition[j][k];
                                toRateMult[j][k] = fromRateMult[j][k];
                                }
                            }
                        i1++;
                        }
                    }
                assert (nCppEventParams==i1);
                }
            }
        }
    
    /* and free up space */
    free (tempVals);
    if (intValsRowSize > 0)
        free (tempIntVals);
    for (i=0; i<nCppEventParams; i++)
        {
        numGlobalChains = chainParams.numRuns * from; /* Revert to the old value to clean old Cpp events in FreeCppEvents() */
        FreeCppEvents(&cppEventParams[i]);
        numGlobalChains = chainParams.numRuns * chainParams.numChains; /*Set to proper value again*/
        }
    if (nCppEventParams > 0)
        free (cppEventParams);

    /* then do the trees (they cannot be done before the parameters because otherwise FillTreeParams will overwrite
       relaxed clock parameters that need to be saved) */

    /* reallocate trees */
    tempTrees = (Tree **) SafeCalloc (2*nRuns*from*numTrees, sizeof(Tree *));
    for (i=0; i<2*nRuns*from*numTrees; i++)
        tempTrees[i] = mcmcTree[i];
    oldMcmcTree = mcmcTree;
    mcmcTree = (Tree **) SafeRealloc ((void *)(mcmcTree), 2 * (size_t)numGlobalChains * (size_t)numTrees * sizeof(Tree*));
    for (i=0; i<2*nRuns*to*numTrees; i++)
        mcmcTree[i] = NULL;

    /* move the old trees over */
    for (run=0; run<nRuns; run++)
        {
        for (chn=0; chn<from; chn++)
            {
            if (chn >= to)
                continue;
            /*Here we move only one tree per chain/state?! Should not we move numTrees??*/
            fromIndex = 2*(run*from + chn)  * numTrees;
            toIndex   = 2*(run*to   + chn)  * numTrees;
            for (k=0;k<2*numTrees;k++)
                {
                mcmcTree[toIndex+k]    = tempTrees[fromIndex+k];
                tempTrees[fromIndex+k] = NULL;
                }
            }
        }

    /* remove any remaining old trees */
    for (i=0; i<2*nRuns*from*numTrees; i++)
        if (tempTrees[i] != NULL)
            FreeTree (tempTrees[i]);
    free (tempTrees);

    /* now fill in the tree parameters */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_TOPOLOGY)
            {
            p->tree += (mcmcTree - oldMcmcTree);    /* calculate new address */
            for (j=0; j<p->nSubParams; j++)
                {
                q = p->subParams[j];
                assert (q->paramType==P_BRLENS);
                q->tree += (mcmcTree - oldMcmcTree);    /* calculate new address */
                if (to > from)
                    for (run=0; run<nRuns; run++)
                        {
                        /*rename old trees because each run has more chains*/
                        for (m=0; m<from; m++)
                            {
                            for (st=0; st<2; st++)
                                {
                                tree = GetTree (q,run*to + m, st);
                                if (numTrees > 1)
                                    sprintf (tree->name, "mcmc.tree%d_%d", p->treeIndex+1, run*to + m +1);
                                else /* if (numTrees == 1) */
                                    sprintf (tree->name, "mcmc.tree_%d", run*to + m +1);
                                }
                            }
                        InitializeChainTrees (q, run*to + from, run*to + to , GetTree (q, 0, 0)->isRooted);
                        }
                }
            }
        else if (p->paramType == P_CPPEVENTS || p->paramType == P_TK02BRANCHRATES || p->paramType == P_IGRBRANCHRATES ||
                 p->paramType == P_ILNBRANCHRATES || p->paramType == P_MIXEDBRCHRATES || p->paramType == P_WNBRANCHRATES)
            p->tree += (mcmcTree - oldMcmcTree);
        else
            assert (p->paramType==P_BRLENS || p->tree==NULL);
        }

    /* fill new tree parameters */
    if (to > from)
        {
        for (run=0; run<nRuns; run++)
            {
            for (chn=from; chn<to; chn++)
                {
                toIndex = run*to + chn;
                FillTreeParams (&globalSeed, toIndex, toIndex+1);
                }
            }
        }

    /* fix stationary frequencies for standard data */
    if (stdStateFreqsRowSize > 0)
        {
        assert (memAllocs[ALLOC_STDSTATEFREQS] == YES);
        stdStateFreqsOld=stdStateFreqs;
        stdStateFreqs = (YFlt *) SafeMalloc (2 * (size_t)stdStateFreqsRowSize * (size_t)numGlobalChains * sizeof (YFlt));
        if (!stdStateFreqs)
            {
            YvyraPrint ("%s   Problem reallocating stdStateFreqs\n", spacer);
            return (ERROR);
            }

        /* set pointers */
        for (k=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->paramType != P_PI || modelParams[p->relParts[0]].dataType != STANDARD)
                continue;
            p->stdStateFreqs += stdStateFreqs-stdStateFreqsOld;
            }
        
        for (run=0; run<nRuns; run++)
            {
            /* copy old chains values*/
            for (chn=0; chn<from; chn++)
                {
                if (chn >= to)
                    break;

                fromIndex = 2*(run*from + chn)*stdStateFreqsRowSize;
                toIndex = 2*(run*to + chn)*stdStateFreqsRowSize;
                for (k=0;k<2*stdStateFreqsRowSize;k++)
                    {
                    stdStateFreqs[toIndex+k]=stdStateFreqsOld[fromIndex+k];
                    }
                }
            /* set new chains */
            FillStdStateFreqs (run*to+from, run*to+to, &globalSeed);
            }
        free(stdStateFreqsOld);
    }

    /* Do the moves */
    /* first allocate space and set up default moves */
    tempMoves = moves;
    moves = NULL;
    memAllocs[ALLOC_MOVES] = NO;
    SetMoves ();
    
    /* then overwrite with old values if present */
    for (i=0; i<numApplicableMoves; i++)
        {
        toMove = moves[i];
        fromMove = tempMoves[i];
        for (run=0; run<nRuns; run++)
            {
            for (chn=0; chn<from; chn++)
                {
                if (chn < to)
                    {
                    fromIndex = run*from + chn;
                    toIndex = run*to + chn;
                    toMove->relProposalProb[toIndex] = fromMove->relProposalProb[fromIndex];
                    for (j=0; j<toMove->moveType->numTuningParams; j++)
                        {
                        toMove->tuningParam[toIndex][j] = fromMove->tuningParam[fromIndex][j];
                        }
                    }
                }
            }
        }
    
    /* and free up space */
    for (i=0; i<numApplicableMoves; i++)
        FreeMove (tempMoves[i]);
    free (tempMoves);
    
    return (NO_ERROR);
}


int ChangeNumRuns (int from, int to)
{
    int         i, i1, j, k, n, nChains;
    Param       *p, *q;
    MoveType    *mvt;
    Tree        **oldMcmcTree;
    YFlt      *oldParamValues;
    YFlt      *stdStateFreqsOld;
    int         *oldintValues;

    if (from == to)
        return (NO_ERROR);

#if 0
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_CPPEVENTS)
            {
            printf ("Trees before changing number of runs\n");
            for (j=0; j<numGlobalChains; j++)
                {
                printf ("Event tree for chain %d\n", j+1);
                for (k=0; k<2*numLocalTaxa; k++)
                    {
                    printf ("%d -- %d:", k, p->nEvents[2*j][k]);
                    for (i1=0; i1<p->nEvents[2*j][k]; i1++)
                        {
                        if (i1 == 0)
                            printf (" (%lf %lf,", p->position[2*j][k], p->rateMult[2*j][k]);
                        else if (i1 == p->nEvents[2*j][k]-1)
                            printf (" %lf %lf)", p->position[2*j][k], p->rateMult[2*j][k]);
                        else
                            printf (" %lf %lf,", p->position[2*j][k], p->rateMult[2*j][k]);
                        }
                    printf ("\n");
                    }
                for (k=0; k<2*numLocalTaxa; k++)
                    {
                    printf ("%d -- %d:", k, p->nEvents[2*j][k]);
                    for (i1=0; i1<p->nEvents[2*j+1][k]; i1++)
                        {
                        if (i1 == 0)
                            printf (" (%lf %lf,", p->position[2*j+1][k], p->rateMult[2*j+1][k]);
                        else if (i1 == p->nEvents[2*j][k]-1)
                            printf (" %lf %lf)", p->position[2*j+1][k], p->rateMult[2*j+1][k]);
                        else
                            printf (" %lf %lf,", p->position[2*j+1][k], p->rateMult[2*j+1][k]);
                        }
                    printf ("\n");
                    }
                }
            }
        }
#endif

    /* set new number of runs */
    chainParams.numRuns = to;
    nChains = chainParams.numChains;
    numGlobalChains = chainParams.numRuns * chainParams.numChains;

    /* do the trees, free tree's memory if we reduce number of trees. */
    for (i=to*2*nChains*numTrees; i<from*2*nChains*numTrees; i++)
        {
        FreeTree (mcmcTree[i]);
        }
    oldMcmcTree = mcmcTree;
    mcmcTree = (Tree **) SafeRealloc ((void *) mcmcTree, numTrees * 2 * numGlobalChains * sizeof (Tree *));
    if (mcmcTree == NULL)
        {
        memAllocs[ALLOC_MCMCTREES] = NO;
        YvyraPrint ("%s   Problem reallocating mcmcTree\n", spacer);
        return (ERROR);
        }

    for (i=from*2*nChains*numTrees; i<to*2*nChains*numTrees; i++)
        {
        mcmcTree[i]=NULL;
        }  
    /* then the cppevents parameters */
    for (i1=0; i1<numParams; i1++)
        {
        p = &params[i1];
        if (p->paramType == P_CPPEVENTS)
            {
            p->nEvents = (int **) SafeRealloc ((void *)p->nEvents, 2*numGlobalChains*sizeof (int *));
            p->nEvents[0] = (int *) SafeRealloc ((void *)p->nEvents[0], 2*numGlobalChains*(2*numLocalTaxa)*sizeof (int));
            for (i=1; i<2*numGlobalChains; i++)
                p->nEvents[i] = p->nEvents[i-1] + (2*numLocalTaxa);
            if (from > to)
                {
                for (j=numGlobalChains; j<from*nChains; j++)
                    {
                    for (k=0; k<2*numLocalTaxa; k++)
                        {
                        free (p->position[2*j+0][k]);
                        p->position[2*j+0][k] = NULL;
                        free (p->rateMult[2*j+0][k]);
                        p->rateMult[2*j+0][k] = NULL;
                        }
                    }
                }
            p->position = (YFlt ***) SafeRealloc ((void *)p->position, 2*numGlobalChains*sizeof (YFlt **));
            p->position[0] = (YFlt **) SafeRealloc ((void *)p->position[0], 2*numGlobalChains*(2*numLocalTaxa)*sizeof (YFlt *));
            for (i=1; i<2*numGlobalChains; i++)
                p->position[i] = p->position[i-1] + (2*numLocalTaxa);
            p->rateMult = (YFlt ***) SafeRealloc ((void *)p->rateMult, 2*numGlobalChains*sizeof (YFlt **));
            p->rateMult[0] = (YFlt **) SafeRealloc ((void *)p->rateMult[0], 2*numGlobalChains*(2*numLocalTaxa)*sizeof (YFlt *));
            for (i=1; i<2*numGlobalChains; i++)
                p->rateMult[i] = p->rateMult[i-1] + (2*numLocalTaxa);
            if (to > from)
                {
                for (j=from*nChains; j<numGlobalChains; j++)
                    {
                    for (k=0; k<2*numLocalTaxa; k++)
                        {
                        p->nEvents[2*j+0][k] = 0;
                        p->position[2*j+0][k] = NULL;
                        p->rateMult[2*j+0][k] = NULL;
                        p->nEvents[2*j+1][k] = 0;
                        p->position[2*j+1][k] = NULL;
                        p->rateMult[2*j+1][k] = NULL;
                        }
                    }
                }
            }
        }
    /* and finally the normal parameters */
    oldParamValues = paramValues;
    if (paramValues != NULL)
        paramValues = (YFlt *) SafeRealloc ((void *) paramValues, (size_t)(paramValsRowSize * 2 * numGlobalChains * sizeof (YFlt)));
    oldintValues = intValues;
    if (intValues != NULL)
        intValues = (int *) SafeRealloc ((void *) intValues, (size_t)(intValsRowSize * 2 * numGlobalChains * sizeof (int)));
    if (paramValues == NULL)
        {
        YvyraPrint ("%s   Problem reallocating paramValues\n", spacer);
        return (ERROR);
        }
    if (intValues == NULL && intValsRowSize > 0)
        {
        YvyraPrint ("%s   Problem reallocating intValues\n", spacer);
        return (ERROR);
        }
    for (i=0; i<numParams; i++)
        {
        if (paramValues)
            {
            params[i].values += (paramValues - oldParamValues);
            params[i].subValues += (paramValues - oldParamValues);
            }
        if (intValues)
            params[i].intValues += (intValues - oldintValues);
        }

    /* fill new chains parameters with appropriate values */
    if (to > from)
        FillNormalParams (&globalSeed, from*nChains, to*nChains);

    /* now fill in the tree parameters */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_TOPOLOGY)
            {
            p->tree += (mcmcTree - oldMcmcTree);    /* calculate new address */
            for (j=0; j<p->nSubParams; j++)
                {
                q = p->subParams[j];
                assert (q->paramType==P_BRLENS);
                q->tree += (mcmcTree - oldMcmcTree);    /* calculate new address */
                InitializeChainTrees (q, from*nChains, to*nChains, GetTree (q, 0, 0)->isRooted);
                }
            }
        else if (p->paramType == P_CPPEVENTS || p->paramType == P_TK02BRANCHRATES || p->paramType == P_IGRBRANCHRATES ||
                 p->paramType == P_ILNBRANCHRATES || p->paramType == P_MIXEDBRCHRATES || p->paramType == P_WNBRANCHRATES)
            p->tree += (mcmcTree - oldMcmcTree);
        }

    FillTreeParams (&globalSeed, from*nChains, to*nChains);

    /* fix stationary frequencies for standard data */
    if (stdStateFreqsRowSize > 0)
        {
        assert (memAllocs[ALLOC_STDSTATEFREQS] == YES);
        stdStateFreqsOld=stdStateFreqs;
        stdStateFreqs = (YFlt *) SafeRealloc ((void *) stdStateFreqs, stdStateFreqsRowSize * 2 * numGlobalChains * sizeof (YFlt));
        if (!stdStateFreqs)
            {
            YvyraPrint ("%s   Problem reallocating stdStateFreqs\n", spacer);
            return (ERROR);
            }
        
        /* set pointers */
        for (k=n=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->paramType != P_PI || modelParams[p->relParts[0]].dataType != STANDARD)
                continue;
            p->stdStateFreqs += stdStateFreqs-stdStateFreqsOld;
            }
        
        FillStdStateFreqs (from*nChains, to*nChains, &globalSeed);
        }

    /* do the moves */
    for (i=0; i<numApplicableMoves; i++)
        {
        mvt = moves[i]->moveType;
        moves[i]->tuningParam = (YFlt **) SafeRealloc ((void *) moves[i]->tuningParam, (size_t)numGlobalChains * sizeof (YFlt *));
        if (mvt->numTuningParams > 0)
            moves[i]->tuningParam[0] = (YFlt *) SafeRealloc ((void *) moves[i]->tuningParam[0], (size_t)numGlobalChains * (size_t)(mvt->numTuningParams) * sizeof (YFlt));
        for (j=1; j<numGlobalChains; j++)
            moves[i]->tuningParam[j] = moves[i]->tuningParam[0] + j * mvt->numTuningParams;
        moves[i]->relProposalProb = (YFlt *) SafeRealloc ((void *) moves[i]->relProposalProb, 4 * (size_t)numGlobalChains * sizeof (YFlt));
        moves[i]->cumProposalProb = moves[i]->relProposalProb + numGlobalChains;
        moves[i]->targetRate = moves[i]->relProposalProb + 2*numGlobalChains;
        moves[i]->lastAcceptanceRate = moves[i]->relProposalProb + 3*numGlobalChains;
        moves[i]->nAccepted = (int *) SafeRealloc ((void *) moves[i]->nAccepted, 5* (size_t)numGlobalChains * sizeof (int));
        moves[i]->nTried = moves[i]->nAccepted + numGlobalChains;
        moves[i]->nBatches = moves[i]->nAccepted + 2*numGlobalChains;
        moves[i]->nTotAccepted = moves[i]->nAccepted + 3*numGlobalChains;
        moves[i]->nTotTried    = moves[i]->nAccepted + 4*numGlobalChains;
        /* initialize all values to default */
        for (j=0; j<numGlobalChains; j++)
            {
            moves[i]->nAccepted[j] = 0;
            moves[i]->nTried[j] = 0;
            moves[i]->nBatches[j] = 0;
            moves[i]->nTotAccepted[j] = 0;
            moves[i]->nTotTried[j] = 0;
            moves[i]->relProposalProb[j] = mvt->relProposalProb;
            moves[i]->cumProposalProb[j] = 0.0;
            moves[i]->lastAcceptanceRate[j] = 0.0;
            for (k=0; k<mvt->numTuningParams; k++)
                moves[i]->tuningParam[j][k] = mvt->tuningParam[k];
            moves[i]->targetRate[j] = mvt->targetRate;
            }
        }

#if 0
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_CPPEVENTS)
            {
            printf ("Trees after changing number of runs\n");
            for (j=0; j<numGlobalChains; j++)
                {
                printf ("Event tree for chain %d\n", j+1);
                for (k=0; k<2*numLocalTaxa; k++)
                    {
                    printf ("%d -- %d:", k, p->nEvents[2*j][k]);
                    assert (p->nEvents[2*j] >= 0);
                    for (i1=0; i1<p->nEvents[2*j][k]; i1++)
                        {
                        if (i1 == 0)
                            printf (" (%lf %lf,", p->position[2*j][k], p->rateMult[2*j][k]);
                        else if (i1 == p->nEvents[2*j][k]-1)
                            printf (" %lf %lf)", p->position[2*j][k], p->rateMult[2*j][k]);
                        else
                            printf (" %lf %lf,", p->position[2*j][k], p->rateMult[2*j][k]);
                        }
                    printf ("\n");
                    }
                for (k=0; k<2*numLocalTaxa; k++)
                    {
                    printf ("%d -- %d:", k, p->nEvents[2*j+1][k]);
                    assert (p->nEvents[2*j+1] >= 0);
                    for (i1=0; i1<p->nEvents[2*j+1][k]; i1++)
                        {
                        if (i1 == 0)
                            printf (" (%lf %lf,", p->position[2*j+1][k], p->rateMult[2*j+1][k]);
                        else if (i1 == p->nEvents[2*j][k]-1)
                            printf (" %lf %lf)", p->position[2*j+1][k], p->rateMult[2*j+1][k]);
                        else
                            printf (" %lf %lf,", p->position[2*j+1][k], p->rateMult[2*j+1][k]);
                        }
                    printf ("\n");
                    }
                }
            }
        }
#endif

    return (NO_ERROR);
}


/*-----------------------------------------------------------
|
|   CheckCharCodingType: check if character is parsimony-
|       informative, variable, or constant
|
-----------------------------------------------------------*/
void CheckCharCodingType (Matrix *m, CharInfo *ci)
{
    int         i, j, k, x, n1[MAX_STD_STATES], n2[MAX_STD_STATES], largest, smallest,
                numPartAmbig, numConsidered, numInformative, lastInformative=0, uniqueBits,
                newPoss, oldPoss;
    BitsLong    combinations[2048], *tempComb, *newComb, *oldComb, bitsLongOne=1;

    /* set up comb pointers */
    oldComb = combinations;
    newComb = oldComb + 1024;

    /* set counters to 0 */
    numPartAmbig = numConsidered = 0;

    /* set variable and informative to yes */
    ci->variable = ci->informative = YES;

    /* set constant to no and state counters to 0 for all states */
    for (i=0; i<MAX_STD_STATES; i++)
        {
        ci->constant[i] = ci->singleton[i] = NO;
        n1[i] = n2[i] = 0;
        }

    for (i=0; i<m->nRows; i++)
        {
        /* retrieve character */
        x = (int) m->origin[m->column + i*m->rowSize];

        /* add it to counters if not all ambiguous */
        if (NBits(x) < ci->nStates)
            {
            numConsidered++;
            if (NBits(x) > 1)
                numPartAmbig++;
            for (j=0; j<MAX_STD_STATES; j++)
                {
                if (((bitsLongOne<<j) & x) != 0)
                    {   
                    n1[j]++;
                    if (NBits(x) == 1)
                        n2[j]++;
                    }
                }
            }
        }

    /* if the ambig counter for any state is equal to the number of considered
       states, then set constant for that state and set variable and informative to no */
    for (i=0; i<MAX_STD_STATES; i++)
        {
        if (n1[i] == numConsidered)
            {
            ci->constant[i] = YES;
            ci->variable = ci->informative = NO;
            }
            else if (n1[i] == 1)
            {
            ci->singleton[i] = YES;
            }
        }

    /* return if variable is no */
    if (ci->variable == NO)
        return;

    /* the character is either (variable and uninformative) or informative */
    
    /* first consider unambiguous characters */
    /* find smallest and largest unambiguous state for this character */
    smallest = MAX_STD_STATES-1;
    largest = 0;
    for (i=0; i<MAX_STD_STATES; i++)
        {
        if (n2[i] > 0)
            {
            if (i < smallest)
                smallest = i;
            if (i > largest)
                largest = i;
            }
        }
        
    /* count the number of informative states in the unambiguous codings */
    for (i=numInformative=0; i<MAX_STD_STATES; i++)
        {
        if (ci->cType == ORD && n2[i] > 0 && i != smallest && i != largest)
            {
            numInformative++;
            lastInformative = i;
            }
        else if (n2[i] > 1)
            {
            numInformative++;
            lastInformative = i;
            }
        }

    /* set informative */
    if (numInformative > 1)
        ci->informative = YES;
    else
        ci->informative = NO;

    
    /* we can return now unless informative is no and numPartAmbig is not 0 */
    if (!(numPartAmbig > 0 && ci->informative == NO))
        return;

    /* check if partially ambiguous observations make this character informative
       after all */
    
    /* first set the bits for the taken states */
    x = 0;
    for (i=0; i<MAX_STD_STATES; i++)
        {
        if (n2[i] > 0 && i != lastInformative)
            x |= (bitsLongOne<<i);
        }
    oldPoss = 1;
    oldComb[0] = x;

    /* now go through all partambig chars and see if we can add them without
       making the character informative */
    for (i=0; i<m->nRows; i++)
        {
        x = (int) m->origin[m->column + i*m->rowSize];
        /* if partambig */ 
        if (NBits(x) > 1 && NBits(x) < ci->nStates)
            {
            /* remove lastInformative */
            x &= !(bitsLongOne<<lastInformative);
            /* reset newPoss */
            newPoss = 0;
            /* see if we can add it, store all possible combinations */
            for (j=0; j<oldPoss; j++)
                {
                uniqueBits = x & (!oldComb[j]);
                for (k=0; k<MAX_STD_STATES; k++)
                    {
                    if (((bitsLongOne<<k) & uniqueBits) != 0)
                        newComb[newPoss++] = oldComb[j] | (bitsLongOne<<k);
                    }
                }
            /* break out if we could not add it */
            if (newPoss == 0)
                break;
            
            /* prepare for next partAmbig */
            oldPoss = newPoss;
            tempComb = oldComb;
            oldComb = newComb;
            newComb = tempComb;
            }
        }

    if (i < m->nRows)
        ci->informative = YES;

    return;
}


/*-----------------------------------------------------------
|
|   CheckModel: check model and warn user if strange things
|      are discovered.
|
-------------------------------------------------------------*/
int CheckModel (void)
{
    int         ch, i, j, k, answer;
    Tree        *t = NULL;
    TreeNode    *p;
    YFlt      treeAge, clockRate;
 
    /* there should only be one calibrated tree */
    for (i=0; i<numTrees; i++)
        {
        t = GetTreeFromIndex(i,0,0);
        if (t->isCalibrated == YES)
            break;
        }
    
    if (i < numTrees)
        {
        if (!strcmp(modelParams[t->relParts[0]].clockRatePr, "Fixed") &&
            AreDoublesEqual (modelParams[t->relParts[0]].clockRateFix, 1.0, 1E-6) == YES)
            {
            YvyraPrint("%s   WARNING: You have calibrated the tree but the clock rate is fixed to 1.0.\n", spacer);
            YvyraPrint("%s      This means that time is measured in expected changes per time unit. If\n", spacer);
            YvyraPrint("%s      the calibrations use a different time scale, you need to modify the model\n", spacer);
            YvyraPrint("%s      by introducing a prior for the clock rate ('prset clockratepr').\n", spacer);

            if (noWarn == NO)
                {
                answer = WantTo("Do you want to continue with the run regardless");
                if (answer == YES)
                    {
                    YvyraPrint("%s   Continuing with the run...\n\n", spacer);
                    }
                else
                    {
                    YvyraPrint("%s   Stopping the run...\n\n", spacer);
                    return (ERROR);
                    }
                }
            }
        }

    /* 
     * Check that the clock rate is consistent with the tree age prior. We cannot check this earlier
     * because the clock rate and tree age are set separately, and we do not know in which order they are set.
     * We need to check all chains because start values are set separately for each chain.
     */
    for (i=0; i<numTrees; i++)
        {
        for (ch=0; ch<numGlobalChains; ch++)
            {
            t = GetTreeFromIndex(i,ch,0);
            if (t->isClock == YES && t->isCalibrated == YES)
                {
                clockRate = *GetParamVals(modelSettings[t->relParts[0]].clockRate, ch, 0);
                treeAge = t->root->left->nodeDepth / clockRate;
                if (!AreDoublesEqual(treeAge, t->root->left->age, 0.000001))
                    {
                    YvyraPrint("%s   ERROR: The tree age setting is inconsistent with the specified tree age prior.\n", spacer);
                    return (ERROR);
                    }
                if (modelParams[t->relParts[0]].treeAgePr.prior == fixed)
                    {
                    if (!AreDoublesEqual(treeAge, modelParams[t->relParts[0]].treeAgePr.priorParams[0], 0.000001))
                        {
                        YvyraPrint("%s   ERROR: The clock rate is inconsistent with the specified tree age prior.\n", spacer);
                        return (ERROR);
                        }
                    }
                else if (modelParams[t->relParts[0]].treeAgePr.prior == uniform)
                    {
                    if (treeAge < modelParams[t->relParts[0]].treeAgePr.priorParams[0] || treeAge > modelParams[t->relParts[0]].treeAgePr.priorParams[1])
                        {    
                        YvyraPrint("%s   ERROR: The clock rate is inconsistent with the specified tree age prior.\n", spacer);
                        return (ERROR);
                        }
                    }
                else if (modelParams[t->relParts[0]].treeAgePr.prior == offsetExponential ||
                         modelParams[t->relParts[0]].treeAgePr.prior == offsetGamma ||
                         modelParams[t->relParts[0]].treeAgePr.prior == truncatedNormal ||
                         modelParams[t->relParts[0]].treeAgePr.prior == offsetLogNormal)
                    {
                    if (treeAge < modelParams[t->relParts[0]].treeAgePr.priorParams[0])
                        {    
                        YvyraPrint("%s   ERROR: The clock rate is inconsistent with the specified tree age prior.\n", spacer);
                        return (ERROR);
                        }
                    }
                }
            }
        }

    /* check coalescence model */
    for (i=0; i<numTrees; i++)
        {
        t = GetTreeFromIndex(i, 0, 0);
        if ((!strcmp(modelParams[t->relParts[0]].clockPr,"Coalescence") ||
             !strcmp(modelParams[t->relParts[0]].clockPr,"Speciestreecoalescence"))
            && !strcmp(modelParams[t->relParts[0]].clockRatePr, "Fixed")
            && AreDoublesEqual (modelParams[t->relParts[0]].clockRateFix, 1.0, 1E-6) == YES)
            {
            if (i == 0) // We only warn for the first tree
                {
                YvyraPrint("%s   WARNING: You are using a coalescent model but the clock rate is fixed to 1.0.\n", spacer);
                YvyraPrint("%s      This is likely to be incorrect unless you have set the population size prior\n", spacer);
                YvyraPrint("%s      ('prset popsizepr') to reflect an appropriate prior on theta. Please check that \n", spacer);
                YvyraPrint("%s      the prior on theta is reasonable for your data.\n", spacer);

                if (noWarn == NO)
                    {
                    answer = WantTo("Do you want to continue with the run regardless");
                    if (answer == YES)
                        {
                        YvyraPrint("%s   Continuing with the run...\n\n", spacer);
                        }
                    else
                        {
                        YvyraPrint("%s   Stopping the run...\n\n", spacer);
                        return (ERROR);
                        }
                    }
                }
            }
        }

    /* Check consistency of best model. First we guarantee that if one topology has
       a species tree prior, then all topologies have the same prior. Then we make
       sure that all clock trees have a coalescence prior. */

    j = 0;
    for (i=0; i<numCurrentDivisions; i++)
        {
        if (!strcmp(modelParams[i].topologyPr, "Speciestree"))
            j++;
        }

    if (j > 0)
        {
        if (j != numCurrentDivisions)
            {
            YvyraPrint("%s   ERROR: If one gene tree has a speciestree prior then all\n", spacer);
            YvyraPrint("%s          gene trees must have the same prior.\n", spacer);
            return (ERROR);
            }
        for (i=0; i<numTrees-1; i++)
            {
            t = GetTreeFromIndex(i,0,0);
            if (strcmp(modelParams[t->relParts[0]].clockPr,"Speciestreecoalescence") != 0)
                {
                YvyraPrint("%s   ERROR: All gene trees must have a speciestreecoalescence prior\n", spacer);
                YvyraPrint("%s          if they fold into a species tree.\n", spacer);
                return (ERROR);
                }
            if (t->isCalibrated == YES)
                {
                for (k=0; k<t->nNodes-1; k++)
                    {
                    p = t->allDownPass[k];
                    if (p->calibration != NULL)
                        {
                        YvyraPrint("%s   ERROR: Gene trees cannot be individually calibrated\n", spacer);
                        YvyraPrint("%s          if they fold into a species tree.\n", spacer);
                        return (ERROR);
                        }
                    }
                }
            }
        }

    return NO_ERROR;
}


/*-----------------------------------------------------------
|
|   CheckExpandedModels: check data partitions that have
|   the codon or doublet model specified
|
-------------------------------------------------------------*/
void CodingToString(int coding, char* string)
{
    if(coding == ALL)
        strcpy(string, "All");
    else if(coding == INFORMATIVE)
        strcpy(string, "Informative");
    else if((coding & VARIABLE) == VARIABLE)
        {
        if (coding == VARIABLE)
            {
            strcpy(string, "Variable");
            }
        else if (coding & NOSINGLETONABSENCE)
            {
            strcpy(string, "Variable|Nosingletonabsence");
            }
        else if (coding & NOSINGLETONPRESENCE)
            {
            strcpy(string, "Variable|Nosingletonpresence");
            }
        }
    else if((coding & NOSINGLETONS) == NOSINGLETONS)
        {
        if (coding == NOSINGLETONS)
            {
            strcpy(string, "Nosingletons");
            }
        else if (coding & NOABSENCESITES)
            {
            strcpy(string, "Noabsencesites|Nosingletons");
            }
        else if (coding & NOPRESENCESITES)
            {
            strcpy(string, "Nopresencesites|Nosingletons");
            }
        }
    else if(coding == NOABSENCESITES)
        {
        strcpy(string, "Noabsencesites");
        }
    else if(coding == NOPRESENCESITES)
        {
        strcpy(string, "Nopresencesites");
        }
    else if(coding == NOSINGLETONABSENCE)
        {
        strcpy(string, "Nosingletonabsence");
        }
    else if(coding == NOSINGLETONPRESENCE)
        {
        strcpy(string, "Nosingletonpresence");
        }
    else if(coding == (NOABSENCESITES | NOSINGLETONABSENCE))
        {
        strcpy(string, "Noabsencesites|Nosingletonabsence");
        }
    else if(coding == (NOABSENCESITES | NOSINGLETONPRESENCE))
        {
        strcpy(string, "Noabsencesites|Nosingletonpresence");
        }
    else if(coding == (NOPRESENCESITES | NOSINGLETONABSENCE))
        {
        strcpy(string, "Nopresencesites|Nosingletonabsence");
        }
    else if(coding == (NOPRESENCESITES | NOSINGLETONPRESENCE))
        {
        strcpy(string, "Nopresencesites|Nosingletonpresence");
        }
}


/*-----------------------------------------------------------
|
|   CompressData: compress original data matrix
|
-------------------------------------------------------------*/
int CompressData (void)
{
    int             a, c, d, i, j, k, t, col[3], isSame, newRow, newColumn,
                    *isTaken, *tempSitesOfPat, *tempChar;
    BitsLong        *tempMatrix;
    ModelInfo       *m;
    ModelParams     *mp;

    /* set all pointers that will be allocated locally to NULL */
    isTaken = NULL;
    tempMatrix = NULL;
    tempSitesOfPat = NULL;
    tempChar = NULL;

#   if defined DEBUG_COMPRESSDATA
    if (PrintMatrix() == ERROR)
        goto errorExit;
    getchar();
#   endif
 
    /* allocate indices pointing from original to compressed matrix */
    if (memAllocs[ALLOC_COMPCOLPOS] == YES)
        {
        free (compColPos);
        compColPos = NULL;
        memAllocs[ALLOC_COMPCOLPOS] = NO;
        }
    compColPos = (int *)SafeMalloc((size_t)numChar * sizeof(int));
    if (!compColPos)
        {
        YvyraPrint ("%s   Problem allocating compColPos (%d)\n", spacer, numChar * sizeof(int));
        goto errorExit;
        }
    for (i=0; i<numChar; i++)
        compColPos[i] = 0;
    memAllocs[ALLOC_COMPCOLPOS] = YES;

    if (memAllocs[ALLOC_COMPCHARPOS] == YES)
        {
        free (compCharPos);
        compCharPos = NULL;
        memAllocs[ALLOC_COMPCHARPOS] = NO;
        }
    compCharPos = (int *)SafeMalloc((size_t)numChar * sizeof(int));
    if (!compCharPos)
        {
        YvyraPrint ("%s   Problem allocating compCharPos (%d)\n", spacer, numChar * sizeof(int));
        goto errorExit;
        }
    for (i=0; i<numChar; i++)
        compCharPos[i] = 0;
    memAllocs[ALLOC_COMPCHARPOS] = YES;

    /* allocate space for temporary matrix, tempSitesOfPat,             */
    /* vector keeping track of whether a character has been compressed, */
    /* and vector indexing first original char for each compressed char */
    tempMatrix = (BitsLong *) SafeCalloc (numLocalTaxa * numLocalChar, sizeof(BitsLong));
    tempSitesOfPat = (int *) SafeCalloc (numLocalChar, sizeof(int));
    isTaken = (int *) SafeCalloc (numChar, sizeof(int));
    tempChar = (int *) SafeCalloc (numLocalChar, sizeof(int));
    if (!tempMatrix || !tempSitesOfPat || !isTaken || !tempChar)
        {
        YvyraPrint ("%s   Problem allocating temporary variables in CompressData\n", spacer);
        goto errorExit;
        }

    /* initialize isTaken */
    for (c=0; c<numChar; c++)
        isTaken[c] = NO;

    /* set index to first empty column in temporary matrix */
    newColumn = 0;

    /* initialize number of compressed characters */
    numCompressedChars = 0;

    /* sort and compress data */
    for (d=0; d<numCurrentDivisions; d++)
        {
        /* set pointers to the model params and settings for this division */
        m = &modelSettings[d];
        mp = &modelParams[d];

        /* set column offset for this division in compressed matrix */
        m->compMatrixStart = newColumn;

        /* set compressed character offset for this division */
        m->compCharStart = numCompressedChars;

        /* set number of compressed characters to 0 for this division */
        m->numChars = 0;

        /* find the number of original characters per model site */
        m->nCharsPerSite = 1;
        if (mp->dataType == DNA || mp->dataType == RNA)
            {   
            if (!strcmp(mp->nucModel, "Doublet"))
                m->nCharsPerSite = 2;
            if (!strcmp(mp->nucModel, "Codon") || !strcmp(mp->nucModel, "Protein"))
                m->nCharsPerSite = 3;
            }
        
        /* sort and compress the characters for this division */
        for (c=0; c<numChar; c++)
            {
            if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != d+1 || isTaken[c] == YES)
                continue;

            col[0] = c;
            isTaken[c] = YES;
            
            /* find additional columns if more than one character per model site      */
            /* return error if the number of matching characters is smaller or larger */
            /* than the actual number of characters per model site                    */
            if (m->nCharsPerSite > 1)
                {
                j = 1;
                if (charInfo[c].charId == 0)
                    {
                    YvyraPrint ("%s   Character %d is not properly defined\n", spacer, c+1);
                    goto errorExit;
                    }
                for (i=c+1; i<numChar; i++)
                    {
                    if (charInfo[i].charId == charInfo[c].charId)
                        {
                        if (j >= m->nCharsPerSite)
                            {
                            YvyraPrint ("%s   Too many matches in charId (division %d char %d)\n", spacer, d, numCompressedChars);
                            goto errorExit;
                            }
                        else
                            {
                            col[j++] = i;
                            isTaken[i] = YES;
                            }
                        }
                    }
                if (j != m->nCharsPerSite)
                    {
                    YvyraPrint ("%s   Too few matches in charId (division %d char %d)\n", spacer, d, numCompressedChars);
                    goto errorExit;
                    }
                }
            
            /* add character to temporary matrix in column(s) at newColumn */
            for (t=newRow=0; t<numTaxa; t++)
                {
                if (taxaInfo[t].isDeleted == YES)
                    continue;

                for (k=0; k<m->nCharsPerSite; k++)
                    {
                    tempMatrix[pos(newRow,newColumn+k,numLocalChar)] = matrix[pos(t,col[k],numChar)];
                    }
                newRow++;
                }
            
            /* is it unique? */
            isSame = NO;
            if (mp->dataType != CONTINUOUS)
                {
                for (i=m->compMatrixStart; i<newColumn; i+=m->nCharsPerSite)
                    {
                    isSame = YES;
                    for (j=0; j<numLocalTaxa; j++)
                        {
                        for (k=0; k<m->nCharsPerSite; k++)
                            if (tempMatrix[pos(j,newColumn+k,numLocalChar)] != tempMatrix[pos(j,i+k,numLocalChar)])
                                {
                                isSame = NO;
                                break;
                                }
                        if (isSame == NO)
                            break;
                        }
                    if (isSame == YES)
                        break;
                    }
                }

            /* if subject to data augmentation, it is always unique */
            if (!strcmp(mp->augmentData, "Yes"))
                {
                for (k=0; k<m->nCharsPerSite; k++)
                    {
                    if (charInfo[col[k]].isMissAmbig == YES)
                        isSame = NO;
                    }
                }

            if (isSame == NO)
                {
                /* if it is unique then it should be added */
                tempSitesOfPat[numCompressedChars] = 1;
                for (k=0; k<m->nCharsPerSite; k++)
                    {
                    compColPos[col[k]] = newColumn + k;
                    compCharPos[col[k]] = numCompressedChars;
                    tempChar[newColumn + k] = col[k];
                    }
                newColumn+=m->nCharsPerSite;
                m->numChars++;
                numCompressedChars++;
                }
            else
                {
                /* if it is not unique then simply update tempSitesOfPat     */
                /* calculate compressed character position and put it into a */
                /* (i points to compressed column position)                  */
                a = m->compCharStart + ((i - m->compMatrixStart) / m->nCharsPerSite);
                tempSitesOfPat[a]++;
                for (k=0; k<m->nCharsPerSite; k++)
                    {
                    compColPos[col[k]] = i;
                    compCharPos[col[k]] = a;
                    /* tempChar (pointing from compressed to uncompressed) */
                    /* can only be set for first pattern */
                    }
                }
            }   /* next character */
            
        /* check that the partition has at least a single character */
        if (m->numChars <= 0)
            {
            YvyraPrint ("%s   You must have at least one site in a partition. Partition %d\n", spacer, d+1);
            YvyraPrint ("%s   has %d site patterns.\n", spacer, m->numChars);
            goto errorExit;
            }

        m->compCharStop = m->compCharStart + m->numChars;
        m->compMatrixStop = newColumn;

        } /* next division */

    compMatrixRowSize = newColumn;

    /* now we know the size, so we can allocate space for the compressed matrix ... */
    if (memAllocs[ALLOC_COMPMATRIX] == YES)
        {
        free (compMatrix);
        compMatrix = NULL;
        memAllocs[ALLOC_COMPMATRIX] = NO;
        }
    compMatrix = (BitsLong *) SafeCalloc (compMatrixRowSize * numLocalTaxa, sizeof(BitsLong));
    if (!compMatrix)
        {
        YvyraPrint ("%s   Problem allocating compMatrix (%d)\n", spacer, compMatrixRowSize * numLocalTaxa * sizeof(BitsLong));
        goto errorExit;
        }
    memAllocs[ALLOC_COMPMATRIX] = YES;
    
    if (memAllocs[ALLOC_NUMSITESOFPAT] == YES)
        {
        free (numSitesOfPat);
        numSitesOfPat = NULL;
        memAllocs[ALLOC_NUMSITESOFPAT] = NO;
        }
    numSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!numSitesOfPat)
        {
        YvyraPrint ("%s   Problem allocating numSitesOfPat (%d)\n", spacer, numCompressedChars * sizeof(YFlt));
        goto errorExit;
        }
    memAllocs[ALLOC_NUMSITESOFPAT] = YES;

    if (memAllocs[ALLOC_ORIGCHAR] == YES)
        {
        free (origChar);
        origChar = NULL;
        memAllocs[ALLOC_ORIGCHAR] = NO;
        }
    origChar = (int *)SafeMalloc((size_t)compMatrixRowSize * sizeof(int));
    if (!origChar)
        {
        YvyraPrint ("%s   Problem allocating origChar (%d)\n", spacer, numCompressedChars * sizeof(int));
        goto errorExit;
        }
    memAllocs[ALLOC_ORIGCHAR] = YES;

    /* ... and copy the data there */
    for (i=0; i<numLocalTaxa; i++)
        for (j=0; j<compMatrixRowSize; j++)
            compMatrix[pos(i,j,compMatrixRowSize)] = tempMatrix[pos(i,j,numLocalChar)];

    /* Apply character weights: sum weights of all original chars per pattern */
    {
    YFlt *weightAccum = (YFlt *) SafeCalloc (numCompressedChars, sizeof(YFlt));
    if (!weightAccum)
        goto errorExit;
    for (i=0; i<numChar; i++)
        {
        if (charInfo[i].isExcluded == YES)
            continue;
        if (compCharPos[i] >= 0 && compCharPos[i] < numCompressedChars)
            weightAccum[compCharPos[i]] += charInfo[i].weight;
        }
    for (i=0; i<numCompressedChars; i++)
        numSitesOfPat[i] = (CLFlt) weightAccum[i];
    free (weightAccum);
    }

    for (i=0; i<compMatrixRowSize; i++)
        origChar[i] = tempChar[i];

#   if defined (DEBUG_COMPRESSDATA)
    if (PrintCompMatrix() == ERROR)
        goto errorExit;
    getchar();
#   endif

    /* free the temporary variables */
    free (tempSitesOfPat);
    free (tempMatrix);
    free (isTaken);
    free (tempChar);

    return NO_ERROR;

    errorExit:
        if (tempMatrix)
            free (tempMatrix);
        if (tempSitesOfPat)
            free (tempSitesOfPat);
        if (isTaken)
            free (isTaken);
        if (tempChar)
            free (tempChar);

        return ERROR;
}


int DataType (int part)
{
    int     i;

    for (i=0; i<numChar; i++)
        {
        if (partitionId[i][partitionNum] == part + 1)
            break;
        }

    return (charInfo[i].charType);
}


int DoLink (void)
{
    int         i, j, newLine;
    
    YvyraPrint ("%s   Linking\n", spacer);
    
    /* update status of linkTable */
    for (j=0; j<NUM_LINKED; j++)
        {
        newLine = YES;
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (tempLinkUnlink[j][i] == YES)
                {
                if (newLine == YES)
                    {
                    linkNum++;
                    newLine = NO;
                    }
                linkTable[j][i] = linkNum;
                }
            }
        }
    
#   if 0
    for (j=0; j<NUM_LINKED; j++)
        {
        YvyraPrint ("%s   ", spacer);
        for (i=0; i<numCurrentDivisions; i++)
            YvyraPrint ("%d", linkTable[j][i]);
        YvyraPrint ("\n");
        }
#   endif

    /* reinitialize the temporary table */
    for (j=0; j<NUM_LINKED; j++)
        for (i=0; i<numCurrentDivisions; i++)
            tempLinkUnlink[j][i] = NO;

    /* set up parameters and moves */
    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);
    
    return (NO_ERROR);
}


int DoLinkParm (char *parmName, char *tkn)
{
    int         i, j, tempInt;

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
        return (ERROR);
        }
        
    if (inValidCommand == YES)
        {
        for (j=0; j<NUM_LINKED; j++)
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[j][i] = NO;
        inValidCommand = NO;
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(LEFTPAR);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        /* initialize tempLinkUnlinkVec to no */
        for (i=0; i<numCurrentDivisions; i++)
            tempLinkUnlinkVec[i] = NO;
        fromI = toJ = -1;
        foundDash = NO;
        expecting = Expecting(NUMBER) | Expecting(ALPHA);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        if (fromI != -1)
            tempLinkUnlinkVec[fromI-1] = YES;
        /* now copy tempLinkUnlinkVec to appropriate row of tempLinkUnlink */
        if (!strcmp(parmName, "Tratio"))
            {
            }
        else if (!strcmp(parmName, "Revmat"))
            {
            }
        else if (!strcmp(parmName, "Omega"))
            {
            }
        else if (!strcmp(parmName, "Statefreq"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_PI][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Mixturerates"))
        {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_MIXTURE_RATES][i] = tempLinkUnlinkVec[i];
        }
        else if (!strcmp(parmName, "Shape"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_SHAPE][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Pinvar"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_PINVAR][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Correlation"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_CORREL][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Ratemultiplier"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_RATEMULT][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Switchrates"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_SWITCH][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Topology"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_TOPOLOGY][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Brlens"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_BRLENS][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Speciationrate"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_SPECRATE][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Extinctionrate"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_EXTRATE][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Popsize"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_POPSIZE][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Growthrate"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_GROWTH][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Aamodel"))
            {
            }
        else if (!strcmp(parmName, "Cpprate"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_CPPRATE][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Cppmultdev"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_CPPMULTDEV][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Cppevents"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_CPPEVENTS][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "TK02var"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_TK02VAR][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "TK02branchrates"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_TK02BRANCHRATES][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "WNvar"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_WNVAR][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "WNbranchrates"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_WNBRANCHRATES][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "IGRvar"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_IGRVAR][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "IGRbranchrates"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_IGRBRANCHRATES][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "ILNvar"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_ILNVAR][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "ILNbranchrates"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_ILNBRANCHRATES][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Mixedvar"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_MIXEDVAR][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Mixedbrchrates"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_MIXEDBRCHRATES][i] = tempLinkUnlinkVec[i];
            }
        else if (!strcmp(parmName, "Browncorr"))
            {
            for (i=0; i<numCurrentDivisions; i++)
                tempLinkUnlink[P_BMCORR][i] = tempLinkUnlinkVec[i];
            }
        else
            {
            YvyraPrint ("%s   Couldn't find parameter %s to link/unlink\n", spacer, parmName);
            }
        
        expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(COMMA))
        {
        foundComma = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == DIFFERENT)
            {
            YvyraPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
            return (ERROR);
            }
        for (i=0; i<numCurrentDivisions; i++)
            tempLinkUnlinkVec[i] = YES;
        expecting  = Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%d", &tempInt);
        if (tempInt > numCurrentDivisions)
            {
            YvyraPrint ("%s   Partition delimiter is too large\n", spacer);
            return (ERROR);
            }
        if (fromI == -1)
            fromI = tempInt;
        else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
            {
            toJ = tempInt;
            for (i=fromI-1; i<toJ; i++)
                tempLinkUnlinkVec[i] = YES;
            fromI = toJ = -1;
            foundDash = NO;
            }
        else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
            {
            tempLinkUnlinkVec[fromI-1] = YES;
            fromI = tempInt;
            foundComma = NO;
            }
        expecting  = Expecting(COMMA);
        expecting |= Expecting(DASH);
        expecting |= Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoLset (void)
{
    int         i, nApplied, lastActive=0;
    
    nApplied = NumActiveParts ();
    for (i=numCurrentDivisions-1; i>=0; i--)
        {
        if (activeParts[i] == YES)
            {
            lastActive = i;
            break;
            }
        }
            
    /* YvyraPrint ("\n"); */
    if (numCurrentDivisions == 1)
        YvyraPrint ("%s   Successfully set likelihood model parameters\n", spacer);
    else 
        {
        if (nApplied == numCurrentDivisions || nApplied == 0)
            {
            YvyraPrint ("%s   Successfully set likelihood model parameters to all\n", spacer);
            YvyraPrint ("%s      applicable data partitions \n", spacer);
            }
        else
            {
            YvyraPrint ("%s   Successfully set likelihood model parameters to\n", spacer);
            if (nApplied == 1)
                YvyraPrint ("%s   partition", spacer);
            else
                YvyraPrint ("%s   partitions", spacer);
            for (i=0; i<numCurrentDivisions; i++)
                {
                if (activeParts[i] == YES)
                    {
                    if (i == lastActive && nApplied > 1)
                        YvyraPrint (" and %d", i+1);
                    else
                        YvyraPrint (" %d", i+1);
                    if (nApplied > 2 && i != lastActive)
                        YvyraPrint (",");
                    }
                }
            YvyraPrint (" (if applicable)\n");
            }
        }

    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);
    
    return (NO_ERROR);
}


int DoLsetParm (char *parmName, char *tkn)
{
    int         i, j, tempInt, nApplied;
    char        tempStr[100];

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
        return (ERROR);
        }
    if (inValidCommand == YES)
        {
        for (i=0; i<numCurrentDivisions; i++)
            activeParts[i] = NO;
        inValidCommand = NO;
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Applyto (Applyto) *************************************************************/
        if (!strcmp(parmName, "Applyto"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(LEFTPAR);
            else if (expecting == Expecting(LEFTPAR))
                {
                for (i=0; i<numCurrentDivisions; i++)
                    activeParts[i] = NO;
                fromI = toJ = -1;
                foundDash = NO;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                if (fromI != -1)
                    activeParts[fromI-1] = YES;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(COMMA))
                {
                foundComma = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (IsSame ("All", tkn) == DIFFERENT)
                    {
                    YvyraPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
                    return (ERROR);
                    }
                for (i=0; i<numCurrentDivisions; i++)
                    activeParts[i] = YES;
                expecting  = Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt > numCurrentDivisions)
                    {
                    YvyraPrint ("%s   Partition delimiter is too large\n", spacer);
                    return (ERROR);
                    }
                if (fromI == -1)
                    fromI = tempInt;
                else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
                    {
                    toJ = tempInt;
                    for (i=fromI-1; i<toJ; i++)
                        activeParts[i] = YES;
                    fromI = toJ = -1;
                    foundDash = NO;
                    }
                else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
                    {
                    activeParts[fromI-1] = YES;
                    fromI = tempInt;
                    foundComma = NO;
                    }
                    
                expecting  = Expecting(COMMA);
                expecting |= Expecting(DASH);
                expecting |= Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else
                return (ERROR);
            }
        /* set Ngammacat (numGammaCats) ************************************************************/
        else if (!strcmp(parmName, "Ngammacat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 2 && tempInt < MAX_RATE_CATS)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
                            {
                            modelParams[i].numGammaCats = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Ngammacat to %d\n", spacer, modelParams[i].numGammaCats);
                            else
                                YvyraPrint ("%s   Setting Ngammacat to %d for partition %d\n", spacer, modelParams[i].numGammaCats, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Ngammacat argument (should be between 2 and %d)\n", spacer, MAX_RATE_CATS);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nlnormcat (numLnormCats) ************************************************************/
        else if (!strcmp(parmName, "Nlnormcat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 2 && tempInt < MAX_RATE_CATS)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
                            {
                            modelParams[i].numLnormCats = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Nlnormcat to %d\n", spacer, modelParams[i].numLnormCats);
                            else
                                YvyraPrint ("%s   Setting Nlnormcat to %d for partition %d\n", spacer, modelParams[i].numLnormCats, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Nlnormcat argument (should be between 2 and %d)\n", spacer, MAX_RATE_CATS);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nmixtcat (numMixtCats) ************************************************************/
        else if (!strcmp(parmName, "Nmixtcat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 2 && tempInt < MAX_RATE_CATS)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
                            {
                            modelParams[i].numMixtCats = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Nmixtcat to %d\n", spacer, modelParams[i].numMixtCats);
                            else
                                YvyraPrint ("%s   Setting Nmixtcat to %d for partition %d\n", spacer, modelParams[i].numLnormCats, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Nmixtcat argument (should be between 2 and %d)\n", spacer, MAX_RATE_CATS);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Usegibbs (useGibbs) *************************************************************/
        else if (!strcmp(parmName, "Usegibbs"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            if (!strcmp(tempStr, "Yes"))
                                {
                                YvyraPrint ("%s   Downsampling of site rates ('usegibbs = yes') disabled temporarily because of conflict with likelihood calculators\n", spacer);
                                return (ERROR);
                                strcpy(modelParams[i].useGibbs, "Yes");
                                }
                            else
                                {
                                strcpy(modelParams[i].useGibbs, "No");
                                }

                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Usegibbs to %s (if applicable)\n", spacer, modelParams[i].useGibbs);
                            else
                                YvyraPrint ("%s   Setting Usegibbs to %s (if applicable) for partition %d\n", spacer, modelParams[i].useGibbs, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for Usegibbs (using Gibbs sampling of discrete gamma)\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }           
        /* set Gibbsfreq (gibbsFreq) ************************************************************/
        else if (!strcmp(parmName, "Gibbsfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 1 && tempInt <= 1000)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            modelParams[i].gibbsFreq = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Gibbsfreq to %d\n", spacer, modelParams[i].gibbsFreq);
                            else
                                YvyraPrint ("%s   Setting Gibbsfreq to %d for partition %d\n", spacer, modelParams[i].gibbsFreq, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Gibbsgammafreq argument (should be between 1 and 1000)\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set NumM10GammaCats (numM10GammaCats) ************************************************************/
        else if (!strcmp(parmName, "NumM10GammaCats"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 2 && tempInt < MAX_RATE_CATS)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
                            {
                            modelParams[i].numM10GammaCats = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting NumM10GammaCats to %d\n", spacer, modelParams[i].numM10GammaCats);
                            else
                                YvyraPrint ("%s   Setting NumM10GammaCats to %d for partition %d\n", spacer, modelParams[i].numM10GammaCats, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid NumM10GammaCats argument (should be between 2 and %d)\n", spacer, MAX_RATE_CATS);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set NumM10BetaCats (numM10BetaCats) ************************************************************/
        else if (!strcmp(parmName, "NumM10BetaCats"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 2 && tempInt < MAX_RATE_CATS)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
                            {
                            modelParams[i].numM10BetaCats = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting NumM10BetaCats to %d\n", spacer, modelParams[i].numM10BetaCats);
                            else
                                YvyraPrint ("%s   Setting NumM10BetaCats to %d for partition %d\n", spacer, modelParams[i].numM10BetaCats, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid NumM10GammaCats argument (should be between 2 and %d)\n", spacer, MAX_RATE_CATS);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nbetacat (numBetaCats) *****************************************************/
        else if (!strcmp(parmName, "Nbetacat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt >= 2 && tempInt < MAX_RATE_CATS)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
                            {
                            modelParams[i].numBetaCats = tempInt;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Nbetacat to %d\n", spacer, modelParams[i].numBetaCats);
                            else
                                YvyraPrint ("%s   Setting Nbetacat to %d for partition %d\n", spacer, modelParams[i].numBetaCats, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Nbetacat argument (should be between 2 and %d)\n", spacer, MAX_RATE_CATS);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Parsmodel (useParsModel) *******************************************************/
        else if (!strcmp(parmName, "Parsmodel"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            if (!strcmp(tempStr, "Yes"))
                                strcpy(modelParams[i].parsModel, "Yes");
                            else
                                strcpy(modelParams[i].parsModel, "No");

                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Parsmodel to %s\n", spacer, modelParams[i].parsModel);
                            else
                                YvyraPrint ("%s   Setting Parsmodel to %s for partition %d\n", spacer, modelParams[i].parsModel, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for using (so-called) parsimony model\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }           
        /* set Augment (augmentData) **********************************************************/
        else if (!strcmp(parmName, "Augment"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                            {
                            if (!strcmp(tempStr, "Yes"))
                                strcpy(modelParams[i].augmentData, "Yes");
                            else
                                strcpy(modelParams[i].augmentData, "No");

                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Augmentdata to %s\n", spacer, modelParams[i].augmentData);
                            else
                                YvyraPrint ("%s   Setting Augmentdata to %s for partition %d\n", spacer, modelParams[i].augmentData, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for data augmentation\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }           
        /* set Ploidy (ploidy) ***************************************************************/
        else if (!strcmp(parmName, "Ploidy"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) &&
                            (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                            {
                            strcpy(modelParams[i].ploidy, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting ploidy level to %s\n", spacer, modelParams[i].ploidy);
                            else
                                YvyraPrint ("%s   Setting ploidy level to %s for partition %d\n", spacer, modelParams[i].ploidy, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid ploidy level argument\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Rates (ratesModel) *************************************************************/
        else if (!strcmp(parmName, "Rates"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                            {
                            if (!strcmp(tempStr, "Adgamma") && (modelParams[i].dataType != DNA && modelParams[i].dataType != RNA && modelParams[i].dataType != PROTEIN))
                                {
                                /* we won't apply an adgamma model to anything but DNA, RNA, or PROTEIN data */
                                }
                            else if ((!strcmp(tempStr, "Propinv") ||  !strcmp(tempStr, "Invgamma")) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
                                {
                                /* we will not apply pinvar to standard or restriction site data */
                                }
                            else
                                {
                                strcpy(modelParams[i].ratesModel, tempStr);
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Rates to %s\n", spacer, modelParams[i].ratesModel);
                                else
                                    YvyraPrint ("%s   Setting Rates to %s for partition %d\n", spacer, modelParams[i].ratesModel, i+1);
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Rates argument\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Covarion (covarionModel) *******************************************************/
        else if (!strcmp(parmName, "Covarion"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) &&
                            (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                            {
                            strcpy(modelParams[i].covarionModel, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Covarion to %s\n", spacer, modelParams[i].covarionModel);
                            else
                                YvyraPrint ("%s   Setting Covarion to %s for partition %d\n", spacer, modelParams[i].covarionModel, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Rates argument\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Coding (missingType) ***********************************************************/
        else if (!strcmp(parmName, "Coding"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                for (i=0; i<numCurrentDivisions; i++)
                    modelParams[i].coding = ALL;
                
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(VERTICALBAR))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD))
                            {
                            if(!strcmp(tempStr, "All"))
                                {
                                modelParams[i].coding |= ALL; // Does nothing
                                }
                            else if(!strcmp(tempStr, "Nosingletons"))
                                {
                                modelParams[i].coding |= NOSINGLETONS;
                                }
                            else if(!strcmp(tempStr, "Variable"))
                                {
                                modelParams[i].coding |= VARIABLE;
                                }
                            else if(!strcmp(tempStr, "Informative"))
                                {
                                modelParams[i].coding |= INFORMATIVE;
                                }
                            else
                                {
                                if(modelParams[i].dataType != RESTRICTION)
                                    {
                                    YvyraPrint ("%s   Invalid coding for standard characters: %s\n", spacer, tempStr);
                                    return (ERROR);
                                    }
                                
                                if(!strcmp(tempStr, "Noabsencesites"))
                                    {
                                    modelParams[i].coding |= NOABSENCESITES;
                                    }
                                else if(!strcmp(tempStr, "Nopresencesites"))
                                    {
                                    modelParams[i].coding |= NOPRESENCESITES;
                                    }
                                else if(!strcmp(tempStr, "Nosingletonpresence"))
                                    {
                                    modelParams[i].coding |= NOSINGLETONPRESENCE;
                                    }
                                else if(!strcmp(tempStr, "Nosingletonabsence"))
                                    {
                                    modelParams[i].coding |= NOSINGLETONABSENCE;
                                    }
                                }
                            
                            CodingToString(modelParams[i].coding, modelParams[i].codingString);
                            
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Enabling Coding %s\n", spacer, tempStr);
                            else
                                YvyraPrint ("%s   Enabling Coding %s for partition %d\n", spacer, tempStr, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for missing patterns\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON) | Expecting(VERTICALBAR);
                }
            else
                return (ERROR);
            }
        /* set Statefrmod (statefreqModel) ************************************************************/  //SK
        else if (!strcmp(parmName, "Statefreqmodel") || !strcmp(parmName, "Statefrmod"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();

                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0))
                            {     
                            strcpy(modelParams[i].statefreqModel, tempStr);
                            modelParams[i].nStates = NumStates (i);
      
                            if (nApplied == 0 && numCurrentDivisions == 1) 
                                YvyraPrint ("%s   Setting Statefrmod to %s\n", spacer, modelParams[i].statefreqModel);
                            else  
                                YvyraPrint ("%s   Setting Statefrmod to %s for partition %d\n", spacer, modelParams[i].statefreqModel, i+1);
      
                            if (modelParams[i].dataType != RESTRICTION && strcmp(modelParams[i].statefreqModel, "Stationary"))
                                {     
                                YvyraPrint ("%s   Invalid setting for Statefrmod: non-stationary models only\n", spacer);
                                YvyraPrint ("%s   implemented for data type \"RESTRICTION\"\n", spacer);
                                return (ERROR);
                                }     
                            }
                        }
                   }
                else
                    {
                    YvyraPrint ("%s   Invalid setting for state frequency model\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* any additional setting goes here */

        else
            return (ERROR);
        }

    return (NO_ERROR);
}


int DoPlot (void)
{
    int             i, n, nHeaders, burnin, len, longestHeader, whichIsX, whichIsY, numPlotted;
    char            temp[100], **headerNames = NULL;
    SumpFileInfo    fileInfo;
    ParameterSample *parameterSamples;
    
#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return NO_ERROR;
#   endif
    
    /* initialize values */
    headerNames = NULL;
    nHeaders = 0;
    parameterSamples = NULL;
    
    /* tell user we are ready to go */
    YvyraPrint ("%s   Plotting parameters in file %s ...\n", spacer, plotParams.plotFileName);
    
    /* examine plot file */
    if (ExamineSumpFile (plotParams.plotFileName, &fileInfo, &headerNames, &nHeaders) == ERROR)
        return ERROR;
    
    /* Calculate burn in */
    burnin = fileInfo.firstParamLine - fileInfo.headerLine - 1;
    
    /* tell the user that everything is fine */
    YvyraPrint ("%s   Found %d parameter lines in file \"%s\"\n", spacer, fileInfo.numRows + burnin, plotParams.plotFileName);
    if (burnin > 0)
        YvyraPrint ("%s   Of the %d lines, %d of them will be summarized (starting at line %d)\n", spacer, fileInfo.numRows+burnin, fileInfo.numRows, fileInfo.firstParamLine);
    else
        YvyraPrint ("%s   All %d lines will be summarized (starting at line %d)\n", spacer, fileInfo.numRows, fileInfo.firstParamLine);
    YvyraPrint ("%s   (Only the last set of lines will be read, in case multiple\n", spacer);
    YvyraPrint ("%s   parameter blocks are present in the same file.)\n", spacer);
    
    /* allocate space to hold parameter information */
    if (AllocateParameterSamples (&parameterSamples, 1, fileInfo.numRows, fileInfo.numColumns) == ERROR)
        goto errorExit;
    
    /* Now we read the file for real. First, rewind file pointer to beginning of file... */
    if (ReadParamSamples (plotParams.plotFileName, &fileInfo, parameterSamples, 0) == ERROR)
        goto errorExit;
    
    /* get length of longest header */
    longestHeader = 9; /* 9 is the length of the word "parameter" (for printing table) */
    for (i=0; i<nHeaders; i++)
        {
        len = (int) strlen(headerNames[i]);
        if (len > longestHeader)
            longestHeader = len;
        }
    
    /* print x-y plot of parameter vs. generation */
    whichIsX = -1;
    for (i=0; i<nHeaders; i++)
        {
        if (IsSame (headerNames[i], "Gen") == SAME)
            whichIsX = i;
        }
    
    if (whichIsX < 0)
        {
        YvyraPrint ("%s   Could not find a column labelled \"Gen\" \n", spacer);
        goto errorExit;
        }
    
    numPlotted = 0;
    for (n=0; n<nHeaders; n++)
        {
        strcpy (temp, headerNames[n]);
        whichIsY = -1;
        if (!strcmp(plotParams.match, "Perfect"))
            {
            if (IsSame (temp, plotParams.parameter) == SAME)
                whichIsY = n;
            }
        else if (!strcmp(plotParams.match, "All"))
            {
            whichIsY = n;
            }
        else
            {
            if (IsSame (temp, plotParams.parameter) == CONSISTENT_WITH)
                whichIsY = n;
            }
        
        if (whichIsY >= 0 && whichIsX != whichIsY)
            {
            YvyraPrint ("\n%s   Rough trace plot of parameter %s:\n", spacer, headerNames[whichIsY]);
            if (PrintPlot (parameterSamples[whichIsX].values[0], parameterSamples[whichIsY].values[0], fileInfo.numRows) == ERROR)
                goto errorExit;
            numPlotted++;
            }
        }
    
    if (numPlotted == 0)
        {
        YvyraPrint ("%s   Did not find any parameters matching \"%s\" to plot\n", spacer, plotParams.parameter);
        }
    
    /* free memory */
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free(headerNames);
    FreeParameterSamples(parameterSamples);
    
    expecting = Expecting(COMMAND);
    
    return (NO_ERROR);
    
errorExit:
    
    /* free memory */
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free(headerNames);
    FreeParameterSamples(parameterSamples);
    
    expecting = Expecting(COMMAND);
    
    return (ERROR);
}


int DoPlotParm (char *parmName, char *tkn)
{
    int         tempI;
    YFlt      tempD;
    char        tempStr[100];
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before sumt can be used\n", spacer);
        return (ERROR);
        }
    
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            expecting  = Expecting(PARAMETER);
            expecting |= Expecting(SEMICOLON);
            }
        /* set Filename (plotParams.plotFileName) ***************************************************/
        else if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (plotParams.plotFileName, tkn);
                YvyraPrint ("%s   Setting plot filename to %s\n", spacer, plotParams.plotFileName);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Relburnin (chainParams.relativeBurnin) ********************************************************/
        else if (!strcmp(parmName, "Relburnin"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.relativeBurnin = YES;
                    else
                        chainParams.relativeBurnin = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for Relburnin\n", spacer);
                    return (ERROR);
                    }
                if (chainParams.relativeBurnin == YES)
                    YvyraPrint ("%s   Using relative burnin (a fraction of samples discarded).\n", spacer);
                else
                    YvyraPrint ("%s   Using absolute burnin (a fixed number of samples discarded).\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }
        /* set Burnin (chainParams.chainBurnIn) *******************************************************/
        else if (!strcmp(parmName, "Burnin"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                chainParams.chainBurnIn = tempI;
                YvyraPrint ("%s   Setting burnin to %d\n", spacer, chainParams.chainBurnIn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Burninfrac (chainParams.burninFraction) ************************************************************/
        else if (!strcmp(parmName, "Burninfrac"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                if (tempD < 0.01)
                    {
                    YvyraPrint ("%s   Burnin fraction too low (< 0.01)\n", spacer);
                    return (ERROR);
                    }
                if (tempD > 0.50)
                    {
                    YvyraPrint ("%s   Burnin fraction too high (> 0.50)\n", spacer);
                    return (ERROR);
                    }
                chainParams.burninFraction = tempD;
                YvyraPrint ("%s   Setting burnin fraction to %.2f\n", spacer, chainParams.burninFraction);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }
        /* set Parameter (plotParams.parameter) *******************************************************/
        else if (!strcmp(parmName, "Parameter"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (plotParams.parameter, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Parameter (plotParams.match) *******************************************************/
        else if (!strcmp(parmName, "Match"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    strcpy (plotParams.match, tempStr);
                else
                    return (ERROR);
                
                YvyraPrint ("%s   Setting plot matching to %s\n", spacer, plotParams.match);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        
        else
            return (ERROR);
        }
    
    return (NO_ERROR);
}


int DoPropset (void)
{
    YvyraPrint ("%s   Successfully set proposal parameters\n", spacer);
    
    return (NO_ERROR);
}


int DoPropsetParm (char *parmName, char *tkn)
{
    int                 i, j, k, nMatches, tempInt;
    YFlt              tempFloat;
    static MCMCMove     *mv = NULL;
    static YFlt       *theValue, theValueMin, theValueMax;
    static int          jump, runIndex, chainIndex;
    static char         *temp=NULL, *localTkn=NULL; /*freed at the end of the call*/
    static char         *tempName=NULL;         /*not freed at the end of the call*/

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before proposal parameters can be changed\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* we expect a move name with possible run and chain specification as follows:
               <move_name>$<tuning_param_name>(<run>,<chain>)=<number>   -- apply to run <run> and chain <chain>
               <move_name>$<tuning_param_name>(,<chain>)=<number>        -- apply to chain <chain> for all runs
               <move_name>$<tuning_param_name>(<run>,)=<number>          -- apply to all chains of run <run>
               <move_name>$prob(<run>,<chain>)=<number>                  -- change relative proposal probability
               <move_name>$targetrate(<run>,<chain>)=<number>            -- change target acc rate for autotuning

               the parsing is complicated by the fact that the move name can look something like:
               eTBR(Tau{all})
               eTBR(Tau{1,4,5})
               so we need to assemble the move name from several tokens that are parsed out separately;
               here we receive only the first part (before the left parenthesis)
            */
            
            /* copy to local move name */
            SafeStrcpy(&tempName, tkn);
            mv = NULL;
            foundComma = foundEqual = NO;
            expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
                Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
                Expecting(DOLLAR);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (mv == NULL)
            {
            /* we are still assembling the move name */
            SafeStrcat(&tempName, tkn);
            expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
                Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
                Expecting(DOLLAR);
            }
        else
            {
            /* we have a parameter name; now find the parameter name, case insensitive */
            SafeStrcpy(&localTkn, tkn);
            for (i=0; i<(int)strlen(localTkn); i++)
                localTkn[i] = tolower(localTkn[i]);
            nMatches = j = 0;
            for (i=0; i<mv->moveType->numTuningParams; i++)
                {
                SafeStrcpy(&temp, mv->moveType->shortTuningName[i]);
                for (k=0; k<(int)strlen(temp); k++)
                    temp[k] = tolower(temp[k]);
                if (strncmp(localTkn,temp,strlen(localTkn)) == 0)
                    {
                    j = i;
                    nMatches++;
                    }
                }
            if (strncmp(localTkn,"prob",strlen(localTkn)) == 0)
                {
                j = -1;
                nMatches++;
                }
            else if (strncmp(localTkn,"targetrate",strlen(localTkn)) == 0)
                {
                j = -2;
                nMatches++;
                }
            if (nMatches == 0)
                {
                YvyraPrint ("%s   Could not find move parameter to change '%s'\n", spacer, localTkn);  
                return (ERROR);
                }
            else if (nMatches > 1)
                {
                YvyraPrint ("%s   Several move parameters matched the abbreviated name '%s'\n", spacer, localTkn);
                return (ERROR);
                }
            
            if (j == -1)
                {
                theValue = mv->relProposalProb;
                theValueMin = 0.0;
                theValueMax = 1000.0;
                jump = 1;
                }
            else if (j == -2)
                {
                theValue = mv->targetRate;
                theValueMin = 0.10;
                theValueMax = 0.70;
                jump = 1;
                }
            else
                {
                theValue = &mv->tuningParam[0][j];
                theValueMin = mv->moveType->minimum[j];
                theValueMax = mv->moveType->maximum[j];
                jump = mv->moveType->numTuningParams;
                }
            chainIndex = -1;
            runIndex = -1;
            expecting = Expecting(LEFTPAR) | Expecting(EQUALSIGN);
            }
        }
    else if (expecting == Expecting(LEFTCURL) || expecting == Expecting(RIGHTCURL))
        {
        /* we are still assembling the move name */
        SafeStrcat (&tempName, tkn);
        expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
            Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
            Expecting(DOLLAR);
        }
    else if (expecting == Expecting(DOLLAR))
        {
        /* we know that the name is complete now; find the move by its name, 
           case insensitive */
        SafeStrcpy(&localTkn, tempName);
        j=(int)strlen(localTkn);
        for (i=0; i<j; i++)
            localTkn[i] = tolower(localTkn[i]);
            
        /* find the move */
        nMatches = j = 0;
        for (i=0; i<numApplicableMoves; i++)
            {
            mv = moves[i];
            SafeStrcpy(&temp,mv->name);
            for (k=0; k<(int)strlen(temp); k++)
                temp[k] = tolower(temp[k]);
            if (strncmp(temp,localTkn,strlen(localTkn)) == 0)
                {
                j = i;
                nMatches++;
                }
            }
        if (nMatches == 0)
            {
            YvyraPrint ("%s   Could not find move '%s'\n", spacer, localTkn);   
            return (ERROR);
            }
        else if (nMatches > 1)
            {
            YvyraPrint ("%s   Several moves matched the abbreviated name '%s'\n", spacer, localTkn);   
            return (ERROR);
            }
        else
            mv = moves[j];

        foundComma = foundEqual = NO;
        expecting = Expecting(ALPHA);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        if (mv == NULL)
            {
            /* we are still assembling the move name */
            SafeStrcat (&tempName, tkn);
            expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
                Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
                Expecting(DOLLAR);
            }
        else /* if (mv != NULL) */
            {
            /* we will be reading in run and chain indices */
            expecting = Expecting(NUMBER) | Expecting(COMMA);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (mv == NULL)
            {
            /* we are still assembling the move name */
            SafeStrcat (&tempName, tkn);
            expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
                Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
                Expecting(DOLLAR);
            }
        else if (foundEqual == YES)
            {
            sscanf (tkn, "%lf", &tempFloat);
            if (tempFloat < theValueMin || tempFloat > theValueMax)
                {
                YvyraPrint ("%s   The value is out of range (min = %lf; max = %lf)\n", spacer, theValueMin, theValueMax);
                return (ERROR);
                }
            if (runIndex == -1 && chainIndex == -1)
                {
                for (i=0; i<chainParams.numRuns; i++)
                    {
                    for (j=0; j<chainParams.numChains; j++)
                        {
                        *theValue = tempFloat;
                        theValue += jump;
                        }
                    }
                }
            else if (runIndex == -1 && chainIndex >= 0)
                {
                theValue += chainIndex*jump;
                for (i=0; i<chainParams.numRuns; i++)
                    {
                    *theValue = tempFloat;
                    theValue += chainParams.numChains*jump;
                    }
                }
            else if (runIndex >= 0 && chainIndex == -1)
                {
                theValue += runIndex*chainParams.numChains*jump;
                for (i=0; i<chainParams.numChains; i++)
                    {
                    *theValue = tempFloat;
                    theValue += jump;
                    }
                }
            else /* if (runIndex >= 0 && chainIndex >= 0) */
                {
                theValue[runIndex*chainParams.numChains*jump+chainIndex*jump] = tempFloat;
                }
            expecting = Expecting (PARAMETER) | Expecting(SEMICOLON);
            }
        else /* if (foundEqual == NO) */
            {
            sscanf (tkn, "%d", &tempInt);
            if (foundComma == NO)
                {
                if (tempInt <= 0 || tempInt > chainParams.numRuns)
                    {
                    YvyraPrint ("%s   Run index is out of range (min=1; max=%d)\n", spacer, chainParams.numRuns);
                    return (ERROR);
                    }
                runIndex = tempInt - 1;
                expecting = Expecting(COMMA);
                }
            else
                {
                if (tempInt <= 0 || tempInt > chainParams.numChains)
                    {
                    YvyraPrint ("%s   Chain index is out of range (min=1; max=%d)\n", spacer, chainParams.numChains);
                    return (ERROR);
                    }
                chainIndex = tempInt - 1;
                expecting = Expecting(RIGHTPAR);
                }
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        if (mv == NULL)
            {
            /* we are still assembling the move name */
            SafeStrcat (&tempName, tkn);
            expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
                Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
                Expecting(DOLLAR);
            }
        else
            {
            /* we will be reading in chain index, if present */
            foundComma = YES;
            expecting = Expecting(RIGHTPAR) | Expecting(NUMBER);
            }
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        if (mv == NULL)
            {
            /* we are still assembling the move name */
            SafeStrcat (&tempName, tkn);
            expecting = Expecting(LEFTCURL) | Expecting(RIGHTCURL) | Expecting(COMMA) |
                Expecting(LEFTPAR) | Expecting(RIGHTPAR) | Expecting(NUMBER) | Expecting(ALPHA) |
                Expecting(DOLLAR);
            }
        else
            expecting = Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        foundEqual = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    SAFEFREE (temp);
    SAFEFREE (localTkn);
    return (NO_ERROR);
}


int DoPrset (void)
{
    int         i, nApplied, lastActive=0;

    nApplied = NumActiveParts ();
    for (i=numCurrentDivisions-1; i>=0; i--)
        {
        if (activeParts[i] == YES)
            {
            lastActive = i;
            break;
            }
        }
            
    if (numCurrentDivisions == 1)
        YvyraPrint ("%s   Successfully set prior model parameters\n", spacer);
    else 
        {
        if (nApplied == numCurrentDivisions || nApplied == 0)
            {
            YvyraPrint ("%s   Successfully set prior model parameters to all\n", spacer);
            YvyraPrint ("%s   applicable data partitions \n", spacer);
            }
        else
            {
            YvyraPrint ("%s   Successfully set prior model parameters to\n", spacer);
            if (nApplied == 1)
                YvyraPrint ("%s   partition", spacer);
            else
                YvyraPrint ("%s   partitions", spacer);
            for (i=0; i<numCurrentDivisions; i++)
                {
                if (activeParts[i] == YES)
                    {
                    if (i == lastActive && nApplied > 1)
                        YvyraPrint (" and %d", i+1);
                    else
                        YvyraPrint (" %d", i+1);
                    if (nApplied > 2 && i != lastActive)
                        YvyraPrint (",");
                    }
                }
            YvyraPrint (" (if applicable)\n");
            }
        }
    
    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);

    return (NO_ERROR);
}


int DoPrsetParm (char *parmName, char *tkn)
{
    int         i, j, k, tempInt, nApplied, index, ns, flag=0;
    YFlt      tempD, sum;
    char        tempStr[100];

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
        return (ERROR);
        }
    if (inValidCommand == YES)
        {
        for (i=0; i<numCurrentDivisions; i++)
            activeParts[i] = NO;
        inValidCommand = NO;
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Applyto (Applyto) ***********************************************************/
        if (!strcmp(parmName, "Applyto"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(LEFTPAR);
            else if (expecting == Expecting(LEFTPAR))
                {
                for (i=0; i<numCurrentDivisions; i++)
                    activeParts[i] = NO;
                fromI = toJ = -1;
                foundDash = NO;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                if (fromI != -1)
                    activeParts[fromI-1] = YES;
#               if 0
                for (i=0; i<numCurrentDivisions; i++)
                    YvyraPrint ("%d ", activeParts[i]);
                YvyraPrint ("\n");
#               endif
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(COMMA))
                {
                foundComma = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (IsSame ("All", tkn) == DIFFERENT)
                    {
                    YvyraPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
                    return (ERROR);
                    }
                for (i=0; i<numCurrentDivisions; i++)
                    activeParts[i] = YES;
                expecting  = Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt > numCurrentDivisions)
                    {
                    YvyraPrint ("%s   Partition delimiter is too large\n", spacer);
                    return (ERROR);
                    }
                if (fromI == -1)
                    fromI = tempInt;
                else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
                    {
                    toJ = tempInt;
                    for (i=fromI-1; i<toJ; i++)
                        activeParts[i] = YES;
                    fromI = toJ = -1;
                    foundDash = NO;
                    }
                else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
                    {
                    activeParts[fromI-1] = YES;
                    fromI = tempInt;
                    foundComma = NO;
                    }
                    
                expecting  = Expecting(COMMA);
                expecting |= Expecting(DASH);
                expecting |= Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else
                return (ERROR);
            }
        /* set Tratiopr (tRatioPr) *********************************************************/
        /* set Revmatpr (revMatPr) *********************************************************/
        /* set Aarevmatpr (aaRevMatPr) *****************************************************/
        /* set Aamodelpr (aaModelPr) *******************************************************/
        /* set Revratepr (revSymDirPr) *****************************************************/
        else if (!strcmp(parmName, "Revratepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == ERROR)  /* we only allow symmetric dirichlet prior, so no need to store the value */
                    {
                    YvyraPrint ("%s   Invalid Revratepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        sscanf (tkn, "%lf", &tempD);
                        if (tempD <= 0.0)
                            {
                            YvyraPrint ("%s   Symmetric Dirichlet parameter must be positive\n", spacer);
                            return (ERROR);
                            }
                        modelParams[i].revMatSymDir = tempD;
                        if (nApplied == 0 && numCurrentDivisions == 1)
                            YvyraPrint ("%s   Setting Revratepr to Symmetric Dirichlet(%1.2lf)\n", spacer, modelParams[i].revMatSymDir);
                        else
                            YvyraPrint ("%s   Setting Revratepr to Symmetric Dirichlet(%1.2lf) for partition %d\n", spacer, modelParams[i].revMatSymDir, i+1);
                        expecting  = Expecting(RIGHTPAR);
                        }
                    }
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Omegapr (omegaPr) ***********************************************************/
        /* set Ny98omega1pr (ny98omega1pr) *************************************************/
        /* set Ny98omega3pr (ny98omega3pr) *************************************************/
        /* set M3omegapr (m3omegapr) *******************************************************/
        /* set M10betapr (m10betapr) *******************************************************/
        /* set M10gammapr (m10gammapr) *****************************************************/
        /* set Codoncatfreqs (codonCatFreqPr) **********************************************/
        /* set Shapepr (shapePr) ***********************************************************/
        else if (!strcmp(parmName, "Shapepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    flag = 0;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) &&
                            (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN ||
                             modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD))
                            {
                            strcpy(modelParams[i].shapePr, tempStr);
                            flag = 1;
                            }
                        }
                    if (flag == 0)
                        {
                        YvyraPrint ("%s   Warning: %s can be set only for partition containing data of at least one of following type:\n", spacer, parmName);
                        YvyraPrint ("%s       DNA, RNA, PROTEIN, RESTRICTION, STANDARD. Currently there is no active partition with such data.\n", spacer);
                        return (ERROR);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Shapepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) &&
                        (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN ||
                         modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD))
                        {
                        if (!strcmp(modelParams[i].shapePr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].shapeUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].shapeUni[0] >= modelParams[i].shapeUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].shapeUni[1] > MAX_SHAPE_PARAM)
                                    {
                                    YvyraPrint ("%s   Upper value for uniform cannot be greater than %1.2lf\n", spacer, MAX_SHAPE_PARAM);
                                    return (ERROR);
                                    }
                                if (modelParams[i].shapeUni[0] < MIN_SHAPE_PARAM)
                                    {
                                    YvyraPrint ("%s   Lower value for uniform cannot be less than %1.2lf\n", spacer, MIN_SHAPE_PARAM);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Shapepr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Shapepr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].shapePr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].shapeExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Shapepr to Exponential(%1.2lf)\n", spacer, modelParams[i].shapeExp);
                            else
                                YvyraPrint ("%s   Setting Shapepr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].shapeExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].shapePr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].shapeFix = tempD;
                            if (modelParams[i].shapeFix > MAX_SHAPE_PARAM)
                                {
                                YvyraPrint ("%s   Shape parameter cannot be greater than %1.2lf\n", spacer, MAX_SHAPE_PARAM);
                                return (ERROR);
                                }
                            if (modelParams[i].shapeFix < MIN_SHAPE_PARAM)
                                {
                                YvyraPrint ("%s   Shape parameter cannot be less than %1.2lf\n", spacer, MIN_SHAPE_PARAM);
                                return (ERROR);
                                }
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Shapepr to Fixed(%1.2lf)\n", spacer, modelParams[i].shapeFix);
                            else
                                YvyraPrint ("%s   Setting Shapepr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].shapeFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Pinvarpr (pInvarPr) *********************************************************/
        else if (!strcmp(parmName, "Pinvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    flag = 0;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) &&
                            (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                            {
                            strcpy(modelParams[i].pInvarPr, tempStr);
                            flag = 1;
                            }
                        }
                    if (flag == 0)
                        {
                        YvyraPrint ("%s   Warning: %s can be set only for partition containing data of at least one of following type:\n", spacer, parmName);
                        YvyraPrint ("%s            DNA, RNA, PROTEIN. Currently there is no active partition with such data.\n", spacer);
                        return (ERROR);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Pinvarpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) &&
                        (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                        {
                        if (!strcmp(modelParams[i].pInvarPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].pInvarUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].pInvarUni[0] >= modelParams[i].pInvarUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].pInvarUni[1] > 1.0)
                                    {
                                    YvyraPrint ("%s   Upper value for uniform should be less than or equal to 1.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Pinvarpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Pinvarpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].pInvarPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD > 1.0)
                                {
                                YvyraPrint ("%s   Value for Pinvar should be in the interval (0, 1)\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].pInvarFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Pinvarpr to Fixed(%1.2lf)\n", spacer, modelParams[i].pInvarFix);
                            else
                                YvyraPrint ("%s   Setting Pinvarpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].pInvarFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Ratecorrpr (adGammaCorPr) ***************************************************/
        else if (!strcmp(parmName, "Ratecorrpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    flag = 0;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
                            {
                            strcpy(modelParams[i].adGammaCorPr, tempStr);
                            flag = 1;
                            }
                        }
                    if (flag == 0)
                        {
                        YvyraPrint ("%s   Warning: %s can be set only for partition containing either DNA or RNA data.\n", spacer, parmName);
                        YvyraPrint ("%s       Currently there is no active partition with such data.\n", spacer);
                        return (ERROR);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Ratecorrpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                foundDash = NO;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting(NUMBER) | Expecting(DASH);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
                        {
                        if (!strcmp(modelParams[i].adGammaCorPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (foundDash == YES)
                            tempD *= -1.0;
                            modelParams[i].adgCorrUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].adgCorrUni[0] >= modelParams[i].adgCorrUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].adgCorrUni[1] > 1.0)
                                    {
                                    YvyraPrint ("%s   Upper value for uniform should be less than or equal to 1.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].adgCorrUni[0] < -1.0)
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than or equal to -1.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Ratecorrpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].adgCorrUni[0], modelParams[i].adgCorrUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Ratecorrpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].adgCorrUni[0], modelParams[i].adgCorrUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].adGammaCorPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (foundDash == YES)
                                tempD *= -1.0;
                            if (tempD > 1.0 || tempD < -1.0)
                                {
                                YvyraPrint ("%s   Value for Ratecorrpr should be in the interval (-1, +1)\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].adgCorrFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Ratecorrpr to Fixed(%1.2lf)\n", spacer, modelParams[i].adgCorrFix);
                            else
                                YvyraPrint ("%s   Setting Ratecorrpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].adgCorrFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                foundDash = NO;
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER) | Expecting(DASH);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Ratepr (ratePr) *************************************************************/
        else if (!strcmp(parmName, "Ratepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                            {
                            if (!strcmp(tempStr,"Variable"))
                                strcpy(modelParams[i].ratePr, "Dirichlet");
                            else
                                strcpy(modelParams[i].ratePr, tempStr);
                            modelParams[i].ratePrDir = 1.0;
                            if (!strcmp(tempStr,"Variable") || !strcmp(tempStr,"Fixed"))
                                {
                                if (tempStr[0]=='V')
                                    strcat (tempStr," [Dirichlet(..,1,..)]");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Ratepr to %s\n", spacer, tempStr);
                                else
                                    YvyraPrint ("%s   Setting Ratepr to %s for partition %d\n", spacer, tempStr, i+1);
                                if (tempStr[0]=='V')
                                    strcpy (tempStr,"Variable");
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Ratepr argument\n", spacer);
                    return (ERROR);
                    }
                if (!strcmp(tempStr,"Fixed") || !strcmp(tempStr,"Variable"))
                    expecting  = Expecting(PARAMETER) | Expecting(SEMICOLON);
                else
                    expecting = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting (NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                /* find next partition to fill in */
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    if ((activeParts[i] == YES || nApplied == 0) && numVars[i] == 0)
                        break;
                if (i == numCurrentDivisions)
                    {
                    YvyraPrint ("%s   Could not find first ratemultiplier partition\n", spacer);
                    return (ERROR);
                    }
                numVars[i] = 1;
                /* read in the parameter */
                sscanf (tkn, "%lf", &tempD);
                if (tempD < ALPHA_MIN || tempD > ALPHA_MAX)
                    {
                    YvyraPrint ("%s   Ratemultiplier Dirichlet parameter %lf out of range\n", spacer, tempD);
                    return (ERROR);
                    }
                /* set the parameter */
                modelParams[i].ratePrDir = tempD;               
                /* check if all partitions have been filled in */
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && numVars[i] == 0)
                        break;
                    }
                /* set expecting accordingly so that we know what should be coming next */
                if (i == numCurrentDivisions)
                    expecting = Expecting (RIGHTPAR);
                else
                    expecting = Expecting (COMMA);
                }
            else if (expecting == Expecting (COMMA))
                expecting = Expecting (NUMBER);
            else if (expecting == Expecting (RIGHTPAR))
                {
                /* print message */
                for (i=j=0; i<numCurrentDivisions; i++)
                    {
                    if (numVars[i] == 1)
                        {
                        j++;
                        if (j == 1)
                            {
                            YvyraPrint ("%s   Setting Ratepr to Dirichlet(%1.2f",
                                spacer, modelParams[i].ratePrDir);
                            }
                        else
                            YvyraPrint (",%1.2f", modelParams[i].ratePrDir);
                        }
                    }
                if (numCurrentDivisions == 1)
                    YvyraPrint (")\n");
                else
                    {
                    YvyraPrint (") for partition");
                    if (j > 1)
                        YvyraPrint ("s");
                    for (i=k=0; i<numCurrentDivisions; i++)
                        {
                        if (numVars[i] == 1)
                            {
                            k++;
                            if (k == j && j > 1)
                                YvyraPrint (", and %d", i+1);
                            else if (k == 1)
                                YvyraPrint (" %d", i+1);
                            else
                                YvyraPrint (", %d", i+1);
                            }
                        }
                    YvyraPrint ("\n");
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Generatepr (generatePr) *****************************************************/
        else if (!strcmp(parmName, "Generatepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                            {
                            if (!strcmp(tempStr,"Variable"))
                                strcpy(modelParams[i].generatePr, "Dirichlet");
                            else
                                strcpy(modelParams[i].generatePr, tempStr);
                            modelParams[i].generatePrDir = 1.0;
                            if (!strcmp(tempStr,"Variable") || !strcmp(tempStr,"Fixed"))
                                {
                                if (tempStr[0]=='V')
                                    strcat (tempStr," [Dirichlet(..,1,..)]");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Generatepr to %s\n", spacer, tempStr);
                                else
                                    YvyraPrint ("%s   Setting Generatepr to %s for partition %d\n", spacer, tempStr, i+1);
                                if (tempStr[0]=='V')
                                    strcpy (tempStr,"Variable");
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Generatepr argument\n", spacer);
                    return (ERROR);
                    }
                if (!strcmp(tempStr,"Fixed") || !strcmp(tempStr,"Variable"))
                    expecting  = Expecting(PARAMETER) | Expecting(SEMICOLON);
                else
                    expecting = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting (NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                /* find next partition to fill in */
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    if ((activeParts[i] == YES || nApplied == 0) && numVars[i] == 0)
                        break;
                if (i == numCurrentDivisions)
                    {
                    YvyraPrint ("%s   Could not find first generate multiplier partition\n", spacer);
                    return (ERROR);
                    }
                numVars[i] = 1;
                /* read in the parameter */
                sscanf (tkn, "%lf", &tempD);
                if (tempD < ALPHA_MIN || tempD > ALPHA_MAX)
                    {
                    YvyraPrint ("%s   Generate multiplier Dirichlet parameter %lf out of range\n", spacer, tempD);
                    return (ERROR);
                    }
                /* set the parameter */
                modelParams[i].generatePrDir = tempD;
                /* check if all partitions have been filled in */
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && numVars[i] == 0)
                        break;
                    }
                /* set expecting accordingly so that we know what should be coming next */
                if (i == numCurrentDivisions)
                    expecting = Expecting (RIGHTPAR);
                else
                    expecting = Expecting (COMMA);
                }
            else if (expecting == Expecting (COMMA))
                expecting = Expecting (NUMBER);
            else if (expecting == Expecting (RIGHTPAR))
                {
                /* print message */
                for (i=j=0; i<numCurrentDivisions; i++)
                    {
                    if (numVars[i] == 1)
                        {
                        j++;
                        if (j == 1)
                            {
                            YvyraPrint ("%s   Setting Generatepr to Dirichlet(%1.2f",
                                          spacer, modelParams[i].generatePrDir);
                            }
                        else
                            YvyraPrint (",%1.2f", modelParams[i].generatePrDir);
                        }
                    }
                if (numCurrentDivisions == 1)
                    YvyraPrint (")\n");
                else
                    {
                    YvyraPrint (") for partition");
                    if (j > 1)
                        YvyraPrint ("s");
                    for (i=k=0; i<numCurrentDivisions; i++)
                        {
                        if (numVars[i] == 1)
                            {
                            k++;
                            if (k == j && j > 1)
                                YvyraPrint (", and %d", i+1);
                            else if (k == 1)
                                YvyraPrint (" %d", i+1);
                            else
                                YvyraPrint (", %d", i+1);
                            }
                        }
                    YvyraPrint ("\n");
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Covswitchpr (covSwitchPr) ***************************************************/
        else if (!strcmp(parmName, "Covswitchpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    flag = 0;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) &&
                            (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                            {
                            strcpy(modelParams[i].covSwitchPr, tempStr);
                            flag = 1;
                            }
                        }
                    if (flag == 0)
                        {
                        YvyraPrint ("%s   Warning: %s can be set only for partition containing data of at least one of following type:\n", spacer, parmName);
                        YvyraPrint ("%s       DNA, RNA, PROTEIN. Currently there is no active partition with such data.\n", spacer);
                        return (ERROR);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Covswitchpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) &&
                        (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                        {
                        if (!strcmp(modelParams[i].covSwitchPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].covswitchUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].covswitchUni[0] >= modelParams[i].covswitchUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Covswitchpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Covswitchpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].covSwitchPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].covswitchExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Covswitchpr to Exponential(%1.2lf)\n", spacer, modelParams[i].covswitchExp);
                            else
                                YvyraPrint ("%s   Setting Covswitchpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].covswitchExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].covSwitchPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].covswitchFix[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Covswitchpr to Fixed(%1.4lf,%1.4lf)\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[1]);
                                else
                                    YvyraPrint ("%s   Setting Covswitchpr to Fixed(%1.4lf,%1.4lf) for partition %d\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Symdirihyperpr (piSymDirPr for standard characters) *************************/
        else if (!strcmp(parmName, "Symdirihyperpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                foundBeta = NO;
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (foundBeta == NO)
                    {
                    /* expecting to see Uniform, Exponential, or Fixed */
                    if (IsArgValid(tkn, tempStr) == NO_ERROR)
                        {
                        nApplied = NumActiveParts ();
                        flag = 0;
                        for (i=0; i<numCurrentDivisions; i++)
                            {
                            if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
                                {
                                strcpy(modelParams[i].symPiPr, tempStr);
                                flag = 1;
                                }
                            }
                        if (flag == 0)
                            {
                            YvyraPrint ("%s   Warning: %s can be set only for partition containing data", spacer, parmName);
                            YvyraPrint ("  of at least one of the following type: STANDARD, RESTRICTION.");
                            YvyraPrint ("Currently there is no active partition with such data. ");
                            return (ERROR);
                            }
                        }
                    else
                        {
                        YvyraPrint ("%s   Invalid Symdirihyperpr argument\n", spacer);
                        return (ERROR);
                        }
                    expecting  = Expecting(LEFTPAR);
                    for (i=0; i<numCurrentDivisions; i++)
                        numVars[i] = 0;
                    foundBeta = YES;    
                    }   
                else
                    {
                    /* expecting infinity */
                    if (IsSame("Infinity", tkn) == SAME || IsSame("Infinity", tkn) == CONSISTENT_WITH)
                        {
                        nApplied = NumActiveParts ();
                        for (i=0; i<numCurrentDivisions; i++)
                            {
                            if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
                                {
                                if (!strcmp(modelParams[i].symPiPr, "Fixed"))
                                    {
                                    modelParams[i].symBetaFix = -1;
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("%s   Setting Symdirihyperpr to Beta(Infinity)\n", spacer);
                                    else
                                        YvyraPrint ("%s   Setting Symdirihyperpr to Beta(Infinity) for partition %d\n", spacer, i+1);
                                    expecting  = Expecting(RIGHTPAR);
                                    }
                                else
                                    {
                                    YvyraPrint ("%s   Problem setting Symdirihyperpr\n", spacer);
                                    return (ERROR);
                                    }
                                }
                            }
                        expecting  = Expecting(RIGHTPAR);
                        }
                    }       
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
                        {
                        if (!strcmp(modelParams[i].symPiPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].symBetaFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Symdirihyperpr to Fixed(%1.2lf)\n", spacer, modelParams[i].symBetaFix);
                            else
                                YvyraPrint ("%s   Setting Symdirihyperpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].symBetaFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].symPiPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].symBetaExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Symdirihyperpr to Exponential(%1.2lf)\n", spacer, modelParams[i].symBetaExp);
                            else
                                YvyraPrint ("%s   Setting Symdirihyperpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].symBetaExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].symPiPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].symBetaUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)    
                                expecting = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].symBetaUni[0] >= modelParams[i].symBetaUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Symdirihyperpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Symdirihyperpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else
                            {
                            YvyraPrint ("%s   Problem setting Symdirihyperpr\n", spacer);
                            return (ERROR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Statefreqpr (stateFreqPr) ***************************************************/
        else if (!strcmp(parmName, "Statefreqpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsSame ("Equal", tkn) == DIFFERENT && IsSame ("Empirical", tkn) == DIFFERENT)
                    {
                    /* the user wants to specify a dirichlet or fixed prior */
                    if (IsArgValid(tkn, tempStr) == NO_ERROR)
                        {
                        nApplied = NumActiveParts ();
                        for (i=0; i<numCurrentDivisions; i++)
                            {
                            if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                                strcpy(modelParams[i].stateFreqPr, tempStr);
                            }
                        /* if (flag == 0)
                            {
                            YvyraPrint ("%s   Warning: %s can be set only for partition containing CONTINUOUS data.\
                            Currently there is no active partition with such data. ", spacer, parmName);
                            return (ERROR);
                            } */
                        }
                    else
                        {
                        YvyraPrint ("%s   Invalid Statefreqpr argument\n", spacer);
                        return (ERROR);
                        }
                    /* TODO: Here we set flat dirichlet parameters */
                    expecting  = Expecting(LEFTPAR);
                    }
                else
                    {
                    /* the user wants equal or empirical state frequencies */
                    nApplied = NumActiveParts ();
                    if (IsSame ("Equal", tkn) == SAME || IsSame ("Equal", tkn) == CONSISTENT_WITH)
                        {
                        for (i=0; i<numCurrentDivisions; i++)
                            {
                            if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                                strcpy(modelParams[i].stateFreqsFixType, "Equal");
                            }
                        }
                    else if (IsSame ("Empirical", tkn) == SAME || IsSame ("Empirical", tkn) == CONSISTENT_WITH)
                        {
                        for (i=0; i<numCurrentDivisions; i++)
                            {
                            if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                                strcpy(modelParams[i].stateFreqsFixType, "Empirical");
                            }
                        }
                    else
                        {
                        YvyraPrint ("%s   Invalid Statefreqpr delimiter\n", spacer);
                        return (ERROR);
                        }
                    expecting  = Expecting(RIGHTPAR);
                    }
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                tempNumStates = 0;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                sscanf (tkn, "%lf", &tempD);
                tempStateFreqs[tempNumStates++] = tempD;
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                        if (!strcmp(modelParams[i].stateFreqPr,"Fixed"))
                            strcpy(modelParams[i].stateFreqsFixType, "User");
                    }
                expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
                        {
                        ns = NumStates(i);
                        if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
                            {
                            if (tempNumStates == 1)
                                {
                                for (j=0; j<ns; j++)
                                    modelParams[i].stateFreqsDir[j] = tempStateFreqs[0] / ns;
                                YvyraPrint ("%s   Setting Statefreqpr to Dirichlet(", spacer);
                                for (j=0; j<ns; j++)
                                    {
                                    YvyraPrint ("%1.2lf", modelParams[i].stateFreqsDir[j]);
                                    if (j == ns - 1)
                                        YvyraPrint (")");
                                    else
                                        YvyraPrint (","); 
                                    }   
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("\n");
                                else
                                    YvyraPrint (" for partition %d\n", i+1); 
                                modelParams[i].numDirParams = ns;
                                }
                            else
                                {
                                if (tempNumStates != ns)
                                    {
                                    YvyraPrint ("%s   Found %d dirichlet parameters but expecting %d\n", spacer, tempNumStates, modelParams[i].nStates);
                                    return (ERROR);
                                    }
                                else
                                    {
                                    modelParams[i].numDirParams = ns;
                                    for (j=0; j<ns; j++)
                                        modelParams[i].stateFreqsDir[j] = tempStateFreqs[j];
                                    YvyraPrint ("%s   Setting Statefreqpr to Dirichlet(", spacer);
                                    for (j=0; j<ns; j++)
                                        {
                                        YvyraPrint ("%1.2lf", modelParams[i].stateFreqsDir[j]);
                                        if (j == ns - 1)
                                            YvyraPrint (")");
                                        else
                                            YvyraPrint (","); 
                                        }   
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("\n");
                                    else
                                        YvyraPrint (" for partition %d\n", i+1); 
                                    }
                                }
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed"))
                            {
                            if (tempNumStates == 0)
                                {
                                if (!strcmp(modelParams[i].stateFreqsFixType, "Equal"))
                                    YvyraPrint ("%s   Setting Statefreqpr to Fixed(Equal)", spacer);
                                else if (!strcmp(modelParams[i].stateFreqsFixType, "Empirical"))
                                    YvyraPrint ("%s   Setting Statefreqpr to Fixed(Empirical)", spacer);
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("\n");
                                else
                                    YvyraPrint (" for partition %d\n", i+1); 
                                }
                            else 
                                {
                                if (tempNumStates == ns)
                                    {
                                    sum = 0.0;
                                    for (j=0; j<ns; j++)
                                        sum += tempStateFreqs[j];
                                    if (AreDoublesEqual (sum, (YFlt) 1.0, (YFlt) 0.001) == NO)
                                        {
                                        YvyraPrint ("%s   State frequencies do not sum to 1.0\n", spacer);
                                        return (ERROR);
                                        }
                                    strcpy(modelParams[i].stateFreqsFixType, "User");
                                    for (j=0; j<ns; j++)
                                        modelParams[i].stateFreqsFix[j] = tempStateFreqs[j];
                                    YvyraPrint ("%s   Setting Statefreqpr to Fixed(", spacer);
                                    for (j=0; j<ns; j++)
                                        {
                                        YvyraPrint ("%1.2lf", modelParams[i].stateFreqsFix[j]);
                                        if (j == ns - 1)
                                            YvyraPrint (")");
                                        else
                                            YvyraPrint (","); 
                                        }   
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("\n");
                                    else
                                        YvyraPrint (" for partition %d\n", i+1); 
                                    }
                                else
                                    {
                                    YvyraPrint ("%s   Found %d state frequencies but expecting %d\n", spacer, tempNumStates, modelParams[i].nStates);
                                    return (ERROR);
                                    }
                                }
                                
                            }
                        }
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Rootfreqpr (rootFreqPr) *****************************************************/
        else if (!strcmp(parmName, "Rootfreqpr"))
            {     
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {  /*for root frequencies, only allow dirichlet and fixed (nb,nb,...) prior, no equal or empirical*/ 
      
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {     
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {     
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == RESTRICTION)
                            {     
                            strcpy(modelParams[i].rootFreqPr, tempStr);
                            flag=1;
                            }     
                        }     
                    if( flag == 0) 
                        {     
                        YvyraPrint ("%s   Warning: %s can only be set for partition containing RESTRICTION data.", spacer, parmName);
                        return (ERROR);
                        }     
                    }     
                else  
                    {     
                        YvyraPrint ("%s   Invalid Rootfreqpr argument\n", spacer);
                        return (ERROR);
                    }         
                // TODO: Here we set flat dirichlet parameters (SK: ??)
                expecting  = Expecting(LEFTPAR);
                }     
            else if (expecting == Expecting(LEFTPAR))
                { 
                tempNumStates = 0;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                sscanf (tkn, "%lf", &tempD);
                tempStateFreqs[tempNumStates++] = tempD;

                expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        ns = NumStates(i);
                        if (!strcmp(modelParams[i].rootFreqPr,"Dirichlet"))
                            {
                            if (tempNumStates == 1)
                                {
                                for (j=0; j<ns; j++)
                                    modelParams[i].rootFreqsDir[j] = tempStateFreqs[0] / ns;

                                YvyraPrint ("%s   Setting Rootfreqpr to Dirichlet(", spacer);

                                for (j=0; j<ns; j++)
                                    {
                                    YvyraPrint ("%1.2lf", modelParams[i].rootFreqsDir[j]);
                                    if (j == ns - 1)
                                        YvyraPrint (")");
                                    else
                                        YvyraPrint (",");
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("\n");
                                else
                                    YvyraPrint (" for partition %d\n", i+1);
                                modelParams[i].numDirParamsRoot = ns;
                                }
                            else
                                {
                                if (tempNumStates != ns)
                                    {
                                    YvyraPrint ("%s   Found %d dirichlet parameters but expecting %d\n", spacer, tempNumStates, modelParams[i].nStates);
                                    return (ERROR);
                                    }
                                else
                                    {
                                    modelParams[i].numDirParamsRoot = ns;
                                    for (j=0; j<ns; j++)
                                        modelParams[i].rootFreqsDir[j] = tempStateFreqs[j];
                                    YvyraPrint ("%s   Setting Rootfreqpr to Dirichlet(", spacer);
                                    for (j=0; j<ns; j++)
                                        {
                                        YvyraPrint ("%1.2lf", modelParams[i].rootFreqsDir[j]);
                                        if (j == ns - 1)
                                            YvyraPrint (")");
                                        else
                                            YvyraPrint (",");
                                        }
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("\n");
                                    else
                                        YvyraPrint (" for partition %d\n", i+1);
                                    }
                                }
                            }
                        else if (!strcmp(modelParams[i].rootFreqPr,"Fixed"))
                            {
                            if (tempNumStates == 0)
                                {
                                YvyraPrint ("%s   Rootfreqpr=Fixed requires specification of state frequencies at the root,", spacer);
                                YvyraPrint ("%s   e.g., rootfreqpr=fixed(0.5,0.5). ", spacer);
                                return (ERROR);
                                }
                            else
                                {
                                if (tempNumStates == ns)
                                    {
                                    sum = 0.0;
                                    for (j=0; j<ns; j++)
                                        sum += tempStateFreqs[j];
                                    if (AreDoublesEqual (sum, (YFlt) 1.0, (YFlt) 0.001) == NO)
                                        {
                                        YvyraPrint ("%s   Root state frequencies do not sum to 1.0\n", spacer);
                                        return (ERROR);
                                        }
                                    for (j=0; j<ns; j++)
                                        modelParams[i].rootFreqsFix[j] = tempStateFreqs[j];
                                    YvyraPrint ("%s   Setting Rootfreqpr to Fixed(", spacer);
                                    for (j=0; j<ns; j++)
                                        {
                                        YvyraPrint ("%1.2lf", modelParams[i].rootFreqsFix[j]);
                                        if (j == ns - 1)
                                            YvyraPrint (")");
                                        else
                                            YvyraPrint (",");
                                        }
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("\n");
                                    else
                                        YvyraPrint (" for partition %d\n", i+1);
                                    }
                                else
                                    {
                                    YvyraPrint ("%s   Found %d root state frequencies but expecting %d\n", spacer, tempNumStates, modelParams[i].nStates);
                                    return (ERROR);
                                    }
                                }
                            }
                        }
                    }
                    expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Topologypr (topologyPr) *****************************************************/
        else if (!strcmp(parmName, "Topologypr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                foundEqual = YES;
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (foundEqual == YES)
                    {
                    if (IsArgValid(tkn, tempStr) == NO_ERROR)
                        {
                        /* set topology prior */
                        nApplied = NumActiveParts ();
                        for (i=0; i<numCurrentDivisions; i++)
                            {
                            if (activeParts[i] == YES || nApplied == 0)
                                {
                                strcpy(modelParams[i].topologyPr, tempStr);
                                /* erase previous constraints, if any */
                                for (j=0; j<numDefinedConstraints; j++)
                                    modelParams[i].activeConstraints[j] = NO;
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Topologypr to %s\n", spacer, modelParams[i].topologyPr);
                                else
                                    YvyraPrint ("%s   Setting Topologypr to %s for partition %d\n", spacer, modelParams[i].topologyPr, i+1);
                                /* adjust branch length prior if necessary */
                                if (strcmp(modelParams[i].topologyPr,"Fixed") != 0 && strcmp(modelParams[i].brlensPr,"Fixed") == 0)
                                    {
                                    YvyraPrint ("%s   Resetting Brlenspr to default\n", spacer);
                                    if (strcmp(modelParams[i].clockPr,"Clock") == 0)
                                        strcpy(modelParams[i].brlensPr, "Uniform");
                                    else
                                        strcpy(modelParams[i].brlensPr, defaultModel.brlensPr);
                                    }
                                }
                            }
                        }
                    else
                        {
                        YvyraPrint ("%s   Invalid Topologypr argument\n", spacer);
                        return (ERROR);
                        }
                    /* make sure we know what to do next */
                    if (!strcmp(tempStr, "Uniform") || !strcmp(tempStr, "Speciestree"))
                        expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                    else
                        expecting = Expecting(LEFTPAR);
                    foundEqual = NO;
                    }
                else
                    {
                    /* find out whether we need a tree name or constraint name */
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        if (activeParts[i] == YES || nApplied == 0)
                            break;

                    if (foundDash == YES)   /* we must be collecting constraint numbers */
                        {
                        YvyraPrint ("%s   Expecting a number\n", spacer);
                        return (ERROR);
                        }
                    if (!strcmp(modelParams[i].topologyPr,"Constraints"))
                        {
                        /* find constraint number */
                        if (CheckString (constraintNames, numDefinedConstraints, tkn, &index) == ERROR)
                            {
                            YvyraPrint ("%s   Could not find constraint named %s\n", spacer, tkn);
                            return (ERROR);
                            }
                        tempActiveConstraints[index] = YES;
                        expecting = Expecting(RIGHTPAR);
                        expecting |= Expecting(COMMA);
                        }
                    else
                        {
                        /* find tree number */
                        if (GetUserTreeFromName (&index, tkn) == ERROR)
                            {
                            YvyraPrint ("%s   Could not set fixed topology from user tree '%s'\n", spacer, tkn);
                            return (ERROR);
                            }
                        fromI = index + 1;        /* fromI is used to hold the index of the user tree, 1-based */
                        expecting = Expecting(RIGHTPAR);
                        }
                    }
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                for (i=0; i<numDefinedConstraints; i++)
                    tempActiveConstraints[i] = NO;
                fromI = toJ = -1;
                foundDash = foundComma = NO;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(NUMBER))
                {
                /* find out whether we need a tree number or constraint number */
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    if (activeParts[i] == YES || nApplied == 0)
                        break;

                if (!strcmp(modelParams[i].topologyPr,"Constraints"))
                    {
                    if (numDefinedConstraints == 0)
                        {
                        YvyraPrint ("%s   No constraints have been defined\n", spacer);
                        return (ERROR);
                        }
                    sscanf (tkn, "%d", &tempInt);
                    if (tempInt > numDefinedConstraints)
                        {
                        YvyraPrint ("%s   Constraint number is too large\n", spacer);
                        return (ERROR);
                        }
                    if (fromI == -1)
                        {
                        if (foundDash == YES)
                            {
                            YvyraPrint ("%s   Unexpected dash\n", spacer);
                            return (ERROR);
                            }
                        fromI = tempInt;
                        tempActiveConstraints[fromI-1] = YES;
                        }
                    else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
                        {
                        toJ = tempInt;
                        //for (i=fromI-1; i<toJ; i++)
                        for (i=fromI; i<toJ; i++)
                            tempActiveConstraints[i] = YES;
                        fromI = toJ = -1;
                        foundDash = NO;
                        }
                    else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
                        {
                        fromI = tempInt;
                        tempActiveConstraints[fromI-1] = YES;
                        foundComma = NO;
                        }
                    expecting  = Expecting(COMMA);
                    expecting |= Expecting(DASH);
                    expecting |= Expecting(RIGHTPAR);
                    }
                else /* if (!strcmp(modelParams[i].topologyPr,"Fixed")) */
                    {
                    if (numUserTrees == 0)
                        {
                        YvyraPrint ("%s   No user trees have been defined\n", spacer);
                        return (ERROR);
                        }
                    sscanf (tkn, "%d", &tempInt);
                    if (tempInt > numUserTrees)
                        {
                        YvyraPrint ("%s   Tree number is too large\n", spacer);
                        return (ERROR);
                        }
                    if (tempInt < 1)
                        {
                        YvyraPrint ("%s   Tree number is too small\n", spacer);
                        return (ERROR);
                        }
                    fromI = tempInt;
                    expecting = Expecting(RIGHTPAR);    /* only one tree number acceptable */
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                foundComma = YES;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                /* find out whether we need a tree number or constraint number(s) */
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    if (activeParts[i] == YES || nApplied == 0)
                        break;

                if (!strcmp(modelParams[i].topologyPr,"Constraints"))
                    {
                    /* set constraints */
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            modelParams[i].numActiveConstraints = 0;
                            for (j=0; j<numDefinedConstraints; j++)
                                {
                                if (tempActiveConstraints[j] == YES)
                                    {
                                    modelParams[i].activeConstraints[j] = YES;
                                    modelParams[i].numActiveConstraints++;
                                    }
                                else
                                    modelParams[i].activeConstraints[j] = NO;
                                }
                            if (modelParams[i].numActiveConstraints == 0)
                                {
                                YvyraPrint ("%s   No constraints have been defined\n", spacer);
                                return (ERROR);
                                }
                            }
                        }
                    }
                else /* if (!strcmp(modelParams[i].topologyPr,"Constraints")) */
                    {
                    /* set start tree index */
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            modelParams[i].topologyFix = fromI-1;
                        }
                    }
#               if 0
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    YvyraPrint ("%4d -- ", i+1);
                    for (j=0; j<numDefinedConstraints; j++)
                        YvyraPrint (" %d", modelParams[i].activeConstraints[j]);
                    YvyraPrint ("\n");
                    }
#               endif               
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nodeagepr (nodeAgePr) *******************************************************/
        else if (!strcmp(parmName, "Nodeagepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].nodeAgePr, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Nodeagepr to %s\n", spacer, modelParams[i].nodeAgePr);
                            else
                                YvyraPrint ("%s   Setting Nodeagepr to %s for partition %d\n", spacer, modelParams[i].nodeAgePr, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Nodeagepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Brlenspr (brlensPr) *********************************************************/
        else if (!strcmp(parmName, "Brlenspr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = NO;
                foundEqual = YES;
                foundLeftPar = NO;
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (foundEqual == YES)
                    {
                    if (IsArgValid(tkn, tempStr) == NO_ERROR)
                        {
                        strcpy (colonPr, tempStr);
                        nApplied = NumActiveParts ();
                        for (i=0; i<numCurrentDivisions; i++)
                            if (activeParts[i] == YES || nApplied == 0)
                                strcpy(modelParams[i].brlensPr, tempStr);
                        }
                    else
                        {
                        YvyraPrint ("%s   Invalid Brlenspr argument\n", spacer);
                        return (ERROR);
                        }
                    foundEqual = NO;
                    if (!strcmp(colonPr,"Fixed"))
                        expecting = Expecting(LEFTPAR);
                    else
                        expecting = Expecting(COLON);
                    }
                else if (foundLeftPar == YES)
                    {
                    /*process argument of fixed() prior*/
                    /* find tree number */
                    if (GetUserTreeFromName (&tempInt, tkn) == ERROR)
                        {
                        YvyraPrint ("%s   Could not set fixed branch lengths from the user tree '%s'\n", spacer, tkn);
                        return (ERROR);
                        }
                    fromI = tempInt + 1;        /* fromI is used to hold the index of the user tree, 1-based */
                    expecting = Expecting(RIGHTPAR);
                    foundLeftPar = NO;
                    }
                else
                    {
                    if (!strcmp(colonPr, "Unconstrained"))
                        {
                        /* have unconstrained branch lengths, which we expect to have a uniform or exponential distribution */
                        nApplied = NumActiveParts ();
                        if (IsSame ("Uniform", tkn) == SAME || IsSame ("Uniform", tkn) == CONSISTENT_WITH)
                            {
                            for (i=0; i<numCurrentDivisions; i++)
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].unconstrainedPr, "Uniform");
                            }
                        else if (IsSame ("Exponential", tkn) == SAME || IsSame ("Exponential", tkn) == CONSISTENT_WITH)
                            {
                            for (i=0; i<numCurrentDivisions; i++)
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].unconstrainedPr, "Exponential");
                            }
                        else if (IsSame ("GammaDirichlet", tkn) == SAME || IsSame ("GammaDirichlet", tkn) == CONSISTENT_WITH)
                            {
                            for (i=0; i<numCurrentDivisions; i++)
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].unconstrainedPr, "GammaDir");
                            }
                        else if (IsSame ("invGamDirichlet", tkn) == SAME || IsSame ("invGamDirichlet", tkn) == CONSISTENT_WITH)
                            {
                            for (i=0; i<numCurrentDivisions; i++)
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].unconstrainedPr, "invGamDir");
                            }
                        else if (IsSame ("twoExponential", tkn) == SAME || IsSame ("twoExponential", tkn) == CONSISTENT_WITH)
                            {
                            for (i=0; i<numCurrentDivisions; i++)
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].unconstrainedPr, "twoExp");
                            }
                        else
                            {
                            YvyraPrint ("%s   Do not understand %s\n", spacer, tkn);
                            return (ERROR);
                            }
                        expecting  = Expecting(LEFTPAR);
                        }
                    else if (!strcmp(colonPr, "Clock"))
                        {
                        /* otherwise we have a clock constraint and expect uniform, birthdeath, coalescence or fixed prior */
                        nApplied = NumActiveParts ();
                        if (IsSame ("Uniform", tkn) == SAME || IsSame ("Uniform", tkn) == CONSISTENT_WITH)
                            {
                            strcpy (clockPr, "Uniform");
                            for (i=0; i<numCurrentDivisions; i++)
                                {
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].clockPr, "Uniform");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Uniform\n", spacer);
                                else
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Uniform for partition %d\n", spacer, i+1);
                                }
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        else if (IsSame ("Birthdeath", tkn) == SAME || IsSame ("Birthdeath", tkn) == CONSISTENT_WITH)
                            {
                            strcpy (clockPr, "Birthdeath");
                            for (i=0; i<numCurrentDivisions; i++)
                                {
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].clockPr, "Birthdeath");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Birthdeath\n", spacer);
                                else
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Birthdeath for partition %d\n", spacer, i+1);
                                }
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        else if (IsSame ("Coalescence", tkn) == SAME || IsSame ("Coalescence", tkn) == CONSISTENT_WITH)
                            {
                            strcpy (clockPr, "Coalescence");
                            for (i=0; i<numCurrentDivisions; i++)
                                {
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].clockPr, "Coalescence");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Coalescence\n", spacer);
                                else
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Coalescence for partition %d\n", spacer, i+1);
                                }
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        else if (IsSame ("Speciestreecoalescence", tkn) == SAME || IsSame ("Speciestreecoalescence", tkn) == CONSISTENT_WITH)
                            {
                            strcpy (clockPr, "Speciestreecoalescence");
                            for (i=0; i<numCurrentDivisions; i++)
                                {
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].clockPr, "Speciestreecoalescence");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Speciestreecoalescence\n", spacer);
                                else
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Speciestreecoalescence for partition %d\n", spacer, i+1);
                                }
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        else if (IsSame ("Fossilization", tkn) == SAME || IsSame ("Fossilization", tkn) == CONSISTENT_WITH)
                            {
                            strcpy (clockPr, "Fossilization");
                            for (i=0; i<numCurrentDivisions; i++)
                                {
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].clockPr, "Fossilization");
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Fossilization\n", spacer);
                                else
                                    YvyraPrint ("%s   Setting Brlenspr to Clock:Fossilization for partition %d\n", spacer, i+1);
                                }
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        else if (IsSame ("Fixed", tkn) == SAME || IsSame ("Fixed", tkn) == CONSISTENT_WITH)
                            {
                            strcpy (clockPr, "Fixed");
                            for (i=0; i<numCurrentDivisions; i++)
                                {
                                if (activeParts[i] == YES || nApplied == 0)
                                    strcpy(modelParams[i].clockPr, "Fixed");
                                }
                            expecting = Expecting(LEFTPAR);     /* Proceed with tree name */
                            }
                        else
                            {
                            YvyraPrint ("%s   Do not understand %s\n", spacer, tkn);
                            return (ERROR);
                            }
                        }
                    else
                        {
                        YvyraPrint ("%s   Do not understand %s\n", spacer, tkn);
                        return (ERROR);
                        }
                    }
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                foundLeftPar = YES;
                expecting  = Expecting(NUMBER);
                if (!strcmp(colonPr,"Fixed") || (!strcmp(colonPr,"Clock") && !strcmp(clockPr,"Fixed")))
                    {
                    expecting |= Expecting(ALPHA);
                    }
                else
                    {
                    for (i=0; i<numCurrentDivisions; i++)
                        numVars[i] = 0;
                    }
                }
            else if (expecting == Expecting(NUMBER))
                {
                if (!strcmp(colonPr, "Unconstrained"))
                    {
                    /* have unconstrained branch lengths */
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            if (!strcmp(modelParams[i].unconstrainedPr,"Uniform"))
                                {
                                sscanf (tkn, "%lf", &tempD);
                                modelParams[i].brlensUni[numVars[i]++] = tempD;
                                if (numVars[i] == 1)
                                    expecting  = Expecting(COMMA);
                                else
                                    {
                                    if (modelParams[i].brlensUni[0] > 0.000001)
                                        {
                                        YvyraPrint ("%s   Lower value for uniform must equal 0.0\n", spacer);
                                        return (ERROR);
                                        }
                                    if (modelParams[i].brlensUni[0] >= modelParams[i].brlensUni[1])
                                        {
                                        YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                        return (ERROR);
                                        }
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("%s   Setting Brlenspr to Unconstrained:Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brlensUni[0], modelParams[i].brlensUni[1]);
                                    else
                                        YvyraPrint ("%s   Setting Brlenspr to Unconstrained:Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brlensUni[0], modelParams[i].brlensUni[1], i+1);
                                    expecting  = Expecting(RIGHTPAR);
                                    }
                                }
                            else if (!strcmp(modelParams[i].unconstrainedPr,"Exponential"))
                                {
                                sscanf (tkn, "%lf", &tempD);
                                modelParams[i].brlensExp = tempD;
                                if (modelParams[i].brlensExp <= 0.0)
                                    {
                                    YvyraPrint ("%s   Value for exponential must > 0.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brlenspr to Unconstrained:Exponential(%1.2lf)\n", spacer, modelParams[i].brlensExp);
                                else
                                    YvyraPrint ("%s   Setting Brlenspr to Unconstrained:Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].brlensExp, i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            else if (!strcmp(modelParams[i].unconstrainedPr,"GammaDir"))
                                {
                                sscanf (tkn, "%lf", &tempD);
                                modelParams[i].brlensDir[numVars[i]++] = tempD;
                                if (numVars[i] == 1 || numVars[i] == 2 || numVars[i] == 3)
                                    expecting  = Expecting(COMMA);
                                else
                                    {
                                    if (modelParams[i].brlensDir[0] <= 0.0 || modelParams[i].brlensDir[1] <= 0.0 || modelParams[i].brlensDir[2] <= 0.0 || modelParams[i].brlensDir[3] <= 0.0)
                                        {
                                        YvyraPrint ("%s   alphaT & betaT & a & c for GammaDir prior must > 0.0\n", spacer);
                                        return (ERROR);
                                        }
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("%s   Setting Brlenspr to Unconstrained:GammaDir(%1.2lf,%1.4lf,%1.2lf,%1.2lf)\n", spacer, modelParams[i].brlensDir[0], modelParams[i].brlensDir[1], modelParams[i].brlensDir[2], modelParams[i].brlensDir[3]);
                                    else
                                        YvyraPrint ("%s   Setting Brlenspr to Unconstrained:GammaDir(%1.2lf,%1.4lf,%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brlensDir[0], modelParams[i].brlensDir[1], modelParams[i].brlensDir[2], modelParams[i].brlensDir[3], i+1);
                                    expecting  = Expecting(RIGHTPAR);
                                    }
                                }
                            else if (!strcmp(modelParams[i].unconstrainedPr,"invGamDir"))
                                {
                                sscanf (tkn, "%lf", &tempD);
                                modelParams[i].brlensDir[numVars[i]++] = tempD;
                                if (numVars[i] == 1 || numVars[i] == 2 || numVars[i] == 3)
                                    expecting  = Expecting(COMMA);
                                else
                                    {
                                    if (modelParams[i].brlensDir[0] <= 2.0 || modelParams[i].brlensDir[1] <= 0.0 || modelParams[i].brlensDir[2] <= 0.0 || modelParams[i].brlensDir[3] <= 0.0)
                                        {
                                        YvyraPrint ("%s   alphaT must > 2.0, betaT & a & c for invGamDir prior must > 0.0\n", spacer);
                                        return (ERROR);
                                        }
                                    if (nApplied == 0 && numCurrentDivisions == 1)
                                        YvyraPrint ("%s   Setting Brlenspr to Unconstrained:invGamDir(%1.2lf,%1.4lf,%1.2lf,%1.2lf)\n", spacer, modelParams[i].brlensDir[0], modelParams[i].brlensDir[1], modelParams[i].brlensDir[2], modelParams[i].brlensDir[3]);
                                    else
                                        YvyraPrint ("%s   Setting Brlenspr to Unconstrained:invGamDir(%1.2lf,%1.4lf,%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brlensDir[0], modelParams[i].brlensDir[1], modelParams[i].brlensDir[2], modelParams[i].brlensDir[3], i+1);
                                    expecting  = Expecting(RIGHTPAR);
                                    }
                                }
                            else if (!strcmp(modelParams[i].unconstrainedPr,"twoExp"))
                                {
                                    sscanf (tkn, "%lf", &tempD);
                                    modelParams[i].brlens2Exp[numVars[i]++] = tempD;
                                    if (numVars[i] == 1)
                                        expecting  = Expecting(COMMA);
                                    else
                                    {
                                        if (modelParams[i].brlens2Exp[0] <= 0.0 || modelParams[i].brlens2Exp[1] <= 0.0)
                                            {
                                            YvyraPrint ("%s   Values for the two exponentials must > 0.0\n", spacer);
                                            return (ERROR);
                                            }
                                        if (nApplied == 0 && numCurrentDivisions == 1)
                                            YvyraPrint ("%s   Setting Brlenspr to Unconstrained:twoExp(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brlens2Exp[0], modelParams[i].brlens2Exp[1]);
                                        else
                                            YvyraPrint ("%s   Setting Brlenspr to Unconstrained:twoExp(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brlens2Exp[0], modelParams[i].brlens2Exp[1], i+1);
                                        expecting  = Expecting(RIGHTPAR);
                                    }
                                }
                            }
                        }
                    }
                else if (!strcmp(colonPr,"Fixed") || !strcmp(colonPr,"Clock"))
                    {
                    sscanf (tkn, "%d", &tempInt);
                    if (tempInt < 1 || tempInt > numUserTrees)
                        {
                        YvyraPrint ("%s   Tree needs to be in the range %d to %d\n", spacer, 1, numUserTrees);
                        return (ERROR);
                        }
                    fromI = tempInt;
                    expecting = Expecting(RIGHTPAR);
                    }
                foundLeftPar = NO;
                }
            else if (expecting == Expecting(COLON))
                {
                expecting  = Expecting(ALPHA);
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                if (!strcmp(colonPr,"Fixed") || (!strcmp(colonPr,"Clock") && !strcmp(clockPr,"Fixed")))
                    {
                    /* index of a tree which set up branch lengths*/
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            modelParams[i].brlensFix = fromI-1;
                        if (!strcmp(colonPr,"Fixed"))
                            {
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Brlenspr to Fixed(%s)\n", spacer, userTree[fromI-1]->name);
                            else
                                YvyraPrint ("%s   Setting Brlenspr to Fixed(%s) for partition %d\n", spacer, userTree[fromI-1]->name, i+1);
                            }
                        else
                            {
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Brlenspr to Fixed(%s)\n", spacer, userTree[fromI-1]->name);
                            else
                                YvyraPrint ("%s   Setting Brlenspr to Fixed(%s) for partition %d\n", spacer, userTree[fromI-1]->name, i+1);
                            }
                        }
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Speciationpr (speciationPr) *************************************************/
        else if (!strcmp(parmName, "Speciationpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].speciationPr, tempStr);
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Speciationpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].speciationPr,"Uniform"))
                            {
                            modelParams[i].speciationUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].speciationUni[0] >= modelParams[i].speciationUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Speciationpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].speciationUni[0], modelParams[i].speciationUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Speciationpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].speciationUni[0], modelParams[i].speciationUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].speciationPr,"Exponential"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Exponential parameter must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].speciationExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Speciationpr to Exponential(%1.2lf)\n", spacer, modelParams[i].speciationExp);
                            else
                                YvyraPrint ("%s   Setting Speciationpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].speciationExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].speciationPr,"Fixed"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Net speciation rate must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].speciationFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Speciationpr to Fixed(%1.2lf)\n", spacer, modelParams[i].speciationFix);
                            else
                                YvyraPrint ("%s   Setting Speciationpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].speciationFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Extinctionpr (extinctionPr) *************************************************/
        else if (!strcmp(parmName, "Extinctionpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].extinctionPr, tempStr);
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Extinctionpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].extinctionPr,"Beta"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Beta parameter must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].extinctionBeta[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Extinctionpr to Beta(%1.2lf,%1.2lf)\n", spacer, modelParams[i].extinctionBeta[0], modelParams[i].extinctionBeta[1]);
                                else
                                    YvyraPrint ("%s   Setting Extinctionpr to Beta(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].extinctionBeta[0], modelParams[i].extinctionBeta[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].extinctionPr,"Exponential"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Exponential parameter must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].extinctionExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Extinctionpr to Exponential(%1.2lf)\n", spacer, modelParams[i].extinctionExp);
                            else
                                YvyraPrint ("%s   Setting Extinctionpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].extinctionExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].extinctionPr,"Fixed"))
                            {
                            if (tempD <= 0.0 || tempD >= 1.0)
                                {
                                YvyraPrint ("%s   Relative extinction rate must be in range (0,1)\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].extinctionFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Extinctionpr to Fixed(%1.2lf)\n", spacer, modelParams[i].extinctionFix);
                            else
                                YvyraPrint ("%s   Setting Extinctionpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].extinctionFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Fossilizationpr (fossilizationPr) */
        else if (!strcmp(parmName, "Fossilizationpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].fossilizationPr, tempStr);
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Fossilization argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].fossilizationPr,"Beta"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Beta parameter must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].fossilizationBeta[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Fossilizationpr to Beta(%1.2lf,%1.2lf)\n", spacer, modelParams[i].fossilizationBeta[0], modelParams[i].fossilizationBeta[1]);
                                else
                                    YvyraPrint ("%s   Setting Fossilizationpr to Beta(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].fossilizationBeta[0], modelParams[i].fossilizationBeta[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].fossilizationPr,"Exponential"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Exponential parameter must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].fossilizationExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Fossilizationpr to Exponential(%1.2lf)\n", spacer, modelParams[i].fossilizationExp);
                            else
                                YvyraPrint ("%s   Setting Fossilizationpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].fossilizationExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].fossilizationPr,"Fixed"))
                            {
                            if (tempD < 0.0 || tempD >= 1.0)
                                {
                                YvyraPrint ("%s   Relative fossilization rate must be in the range (0,1)\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].fossilizationFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Fossilizationpr to Fixed(%1.2lf)\n", spacer, modelParams[i].fossilizationFix);
                            else
                                YvyraPrint ("%s   Setting Fossilizationpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].fossilizationFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set SampleStrat (sampleStrat) ***************************************************/
        else if (!strcmp(parmName, "Samplestrat"))  // prset samplestrat = random 3: 250 100 60 [fossil sampling], 1: 150 [birth], 2: 200 100 [death];
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].sampleStrat, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting SampleStrat to %s\n", spacer, modelParams[i].sampleStrat);
                            else
                                YvyraPrint ("%s   Setting SampleStrat to %s for partition %d\n", spacer, modelParams[i].sampleStrat, i+1);
                            if (!strcmp(modelParams[i].sampleStrat,"Random") || !strcmp(modelParams[i].sampleStrat,"Diversity"))
                                {
                                foundFSNum[i]  = foundBSNum[i]  = foundDSNum[i]  = NO;
                                foundFSTime[i] = foundBSTime[i] = foundDSTime[i] = NO;
                                modelParams[i].fossilSamplingNum = 0;
                                modelParams[i].birthRateShiftNum = 0;
                                modelParams[i].deathRateShiftNum = 0;
                                numVars[i] = 0;
                                expecting  = Expecting(NUMBER);
                                expecting |= Expecting(PARAMETER);
                                expecting |= Expecting(SEMICOLON);
                                }
                            else
                                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Samplestrat argument\n", spacer);
                    return (ERROR);
                    }
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (foundFSNum[i] == NO)  // fossil sampling rate shift times
                            {
                            sscanf (tkn, "%d", &tempInt);
                            if (tempInt < 0 || tempInt > 99)
                                {
                                YvyraPrint ("%s   Number of fossil sampling rate shift must be between 0 and 100\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].fossilSamplingNum = tempInt;
                            foundFSNum[i] = YES;
                            if (tempInt == 0)
                                {
                                foundFSTime[i] = YES;
                                expecting  = Expecting(COMMA);
                                expecting |= Expecting(PARAMETER);
                                expecting |= Expecting(SEMICOLON);
                                }
                            else
                                expecting = Expecting(COLON);
                            }
                        else if (foundFSTime[i] == NO)
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Time of fossil sampling rate shift must be > 0.\n", spacer);
                                return (ERROR);
                                }
                            if (numVars[i] > 0 && modelParams[i].fossilSamplingTime[numVars[i]-1] < tempD)
                                {
                                YvyraPrint ("%s   Time of fossil sampling rate shift must be in decreasing order\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].fossilSamplingTime[numVars[i]] = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting %d fossil sampling rate shift time to %1.2lf\n", spacer, numVars[i]+1,
                                              modelParams[i].fossilSamplingTime[numVars[i]]);
                            else
                                YvyraPrint ("%s   Setting %d fossil sampling rate shift time to %1.2lf for partition %d\n", spacer, numVars[i]+1,
                                              modelParams[i].fossilSamplingTime[numVars[i]], i+1);
                            numVars[i]++;
                            expecting = Expecting(NUMBER);
                            if (numVars[i] == modelParams[i].fossilSamplingNum) {
                                foundFSTime[i] = YES;
                                numVars[i] = 0;
                                expecting  = Expecting(COMMA);
                                expecting |= Expecting(PARAMETER);
                                expecting |= Expecting(SEMICOLON);
                                }
                            }
                        else if (foundBSNum[i] == NO)  // birth rate shift times
                            {
                            sscanf (tkn, "%d", &tempInt);
                            if (tempInt < 0 || tempInt > 99)
                                {
                                YvyraPrint ("%s   Number of birth rate shift must be between 0 and 100\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].birthRateShiftNum = tempInt;
                            foundBSNum[i] = YES;
                            if (tempInt == 0)
                                {
                                foundBSTime[i] = YES;
                                expecting  = Expecting(COMMA);
                                expecting |= Expecting(PARAMETER);
                                expecting |= Expecting(SEMICOLON);
                                }
                            else
                                expecting = Expecting(COLON);
                            }
                        else if (foundBSTime[i] == NO)
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Time of birth rate shift must be > 0.\n", spacer);
                                return (ERROR);
                                }
                            if (numVars[i] > 0 && modelParams[i].birthRateShiftTime[numVars[i]-1] < tempD)
                                {
                                YvyraPrint ("%s   Time of birth rate shift must be in decreasing order\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].birthRateShiftTime[numVars[i]] = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting %d birth rate shift time to %1.2lf\n", spacer, numVars[i]+1,
                                              modelParams[i].birthRateShiftTime[numVars[i]]);
                            else
                                YvyraPrint ("%s   Setting %d birth rate shift time to %1.2lf for partition %d\n", spacer, numVars[i]+1,
                                              modelParams[i].birthRateShiftTime[numVars[i]], i+1);
                            numVars[i]++;
                            expecting = Expecting(NUMBER);
                            if (numVars[i] == modelParams[i].birthRateShiftNum) {
                                foundBSTime[i] = YES;
                                numVars[i] = 0;
                                expecting  = Expecting(COMMA);
                                expecting |= Expecting(PARAMETER);
                                expecting |= Expecting(SEMICOLON);
                                }
                            }
                        else if (foundDSNum[i] == NO)  // death rate shift times
                            {
                            sscanf (tkn, "%d", &tempInt);
                            if (tempInt < 0 || tempInt > 99)
                                {
                                YvyraPrint ("%s   Number of death rate shift must be between 0 and 100\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].deathRateShiftNum = tempInt;
                            foundDSNum[i] = YES;
                            if (tempInt == 0)
                                {
                                foundDSTime[i] = YES;
                                expecting  = Expecting(COMMA);
                                expecting |= Expecting(PARAMETER);
                                expecting |= Expecting(SEMICOLON);
                                }
                            else
                                expecting = Expecting(COLON);
                            }
                        else if (foundDSTime[i] == NO)
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Time of death rate shift must be > 0.\n", spacer);
                                return (ERROR);
                                }
                            if (numVars[i] > 0 && modelParams[i].deathRateShiftTime[numVars[i]-1] < tempD)
                                {
                                YvyraPrint ("%s   Time of death rate shift must be in decreasing order\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].deathRateShiftTime[numVars[i]] = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting %d death rate shift time to %1.2lf\n", spacer, numVars[i]+1,
                                              modelParams[i].deathRateShiftTime[numVars[i]]);
                            else
                                YvyraPrint ("%s   Setting %d death rate shift time to %1.2lf for partition %d\n", spacer, numVars[i]+1,
                                              modelParams[i].deathRateShiftTime[numVars[i]], i+1);
                            numVars[i]++;
                            expecting = Expecting(NUMBER);
                            if (numVars[i] == modelParams[i].deathRateShiftNum)
                                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COLON) || expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else
                {
                return (ERROR);
                }
            }
        /* set Sampleprob (sampleProb) *****************************************************/
        else if (!strcmp(parmName, "Sampleprob"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0))
                        {
                        sscanf (tkn, "%lf", &tempD);
                        if (tempD < 0.0 || tempD > 1.0)
                            {
                            YvyraPrint ("%s   Sampleprob should be in range [0,1]\n", spacer);
                            return (ERROR);
                            }
                        modelParams[i].sampleProb = tempD;
                        if (nApplied == 0 && numCurrentDivisions == 1)
                            YvyraPrint ("%s   Setting Sampleprob to %1.8lf\n", spacer, modelParams[i].sampleProb);
                        else
                            YvyraPrint ("%s   Setting Sampleprob to %1.8lf for partition %d\n", spacer, modelParams[i].sampleProb, i+1);
                        }
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Treeagepr (treeAgePr) *******************************************************/
        else if (!strcmp(parmName, "Treeagepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].treeAgePr.name, tempStr);
                            if (!strcmp(tempStr,"Fixed"))
                                modelParams[i].treeAgePr.prior = fixed;
                            else if (!strcmp(tempStr,"Uniform"))
                                modelParams[i].treeAgePr.prior = uniform;
                            else if (!strcmp(tempStr,"Offsetexponential"))
                                modelParams[i].treeAgePr.prior = offsetExponential;
                            else if (!strcmp(tempStr,"Truncatednormal"))
                                modelParams[i].treeAgePr.prior = truncatedNormal;
                            else if (!strcmp(tempStr,"Lognormal"))
                                modelParams[i].treeAgePr.prior = logNormal;
                            else if (!strcmp(tempStr,"Offsetlognormal"))
                                modelParams[i].treeAgePr.prior = offsetLogNormal;
                            else if (!strcmp(tempStr,"Gamma"))
                                modelParams[i].treeAgePr.prior = standardGamma;
                            else if (!strcmp(tempStr,"Offsetgamma"))
                                modelParams[i].treeAgePr.prior = offsetGamma;
                            }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Treeagepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    if (activeParts[i] == YES || nApplied == 0)
                        strcat(modelParams[i].treeAgePr.name, "(");
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                sprintf (tempStr, "%1.2lf", tempD);
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (numVars[i] == 0 && tempD < 0.0)
                            {
                            if (modelParams[i].treeAgePr.prior == uniform ||
                                modelParams[i].treeAgePr.prior == offsetExponential ||
                                modelParams[i].treeAgePr.prior == truncatedNormal ||
                                modelParams[i].treeAgePr.prior == offsetLogNormal ||
                                modelParams[i].treeAgePr.prior == offsetGamma)
                                YvyraPrint ("%s   Minimum, offset or truncation point must be nonnegative\n", spacer);
                            else if (modelParams[i].treeAgePr.prior == fixed)
                                YvyraPrint ("%s   Fixed age must be nonnegative\n", spacer);
                            else
                                YvyraPrint ("%s   Mean must be nonnegative\n", spacer);
                            break;
                            }
                        else if (numVars[i] == 1)
                            {
                            if (modelParams[i].treeAgePr.prior == uniform && tempD <= modelParams[i].treeAgePr.priorParams[0])
                                {
                                YvyraPrint ("%s   Max of uniform distribution must be larger than min\n", spacer);
                                break;
                                }
                            else if ((modelParams[i].treeAgePr.prior == standardGamma || modelParams[i].treeAgePr.prior == logNormal) && tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Standard deviation must be positive\n", spacer);
                                break;
                                }
                            else if ((modelParams[i].treeAgePr.prior == offsetExponential ||
                                      modelParams[i].treeAgePr.prior == offsetGamma ||
                                      modelParams[i].treeAgePr.prior == offsetLogNormal) && tempD <= modelParams[i].treeAgePr.priorParams[0])
                                {
                                YvyraPrint ("%s   Mean must be larger than offset\n", spacer);
                                break;
                                }
                            }
                        else if (numVars[i] == 2 && tempD <= 0.0)
                            {
                            YvyraPrint ("%s   Standard deviation must be positive\n", spacer);
                            break;
                            }
                        modelParams[i].treeAgePr.priorParams[numVars[i]++] = tempD;
                        sprintf (tempStr, "%1.2lf", tempD);
                        strcat(modelParams[i].treeAgePr.name, tempStr);
                        if (modelParams[i].treeAgePr.prior == fixed || numVars[i] == 3)
                            expecting = Expecting(RIGHTPAR);
                        else if (numVars[i] == 1)
                            expecting = Expecting(COMMA);
                        else if (modelParams[i].treeAgePr.prior == standardGamma ||
                                 modelParams[i].treeAgePr.prior == uniform ||
                                 modelParams[i].treeAgePr.prior == offsetExponential ||
                                 modelParams[i].treeAgePr.prior == logNormal)
                            expecting = Expecting(RIGHTPAR);
                        else
                            expecting = Expecting(COMMA);
                        }
                    }
                if (i < numCurrentDivisions)
                    {
                    /* An error occurred. Reset calibrations and bail out */
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].treeAgePr.name, "Gamma(1.00,1.00)");
                            modelParams[i].treeAgePr.prior = standardGamma;
                            modelParams[i].treeAgePr.priorParams[0] = 1.0;
                            modelParams[i].treeAgePr.priorParams[1] = 1.0;
                            modelParams[i].treeAgePr.priorParams[2] = -1.0;
                            modelParams[i].treeAgePr.min = 0.0;
                            modelParams[i].treeAgePr.max = POS_INFINITY;
                            }
                        }
                    return (ERROR);
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    if (activeParts[i] == YES || nApplied == 0)
                        strcat(modelParams[i].treeAgePr.name, ",");
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        strcat(modelParams[i].treeAgePr.name, ")");
                        YvyraPrint ("%s   Setting Treeagepr to %s\n", spacer, modelParams[i].treeAgePr.name);

                        if (modelParams[i].treeAgePr.prior == fixed)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbFix;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioFix;
                            modelParams[i].treeAgePr.min           = modelParams[i].treeAgePr.priorParams[0];
                            modelParams[i].treeAgePr.max           = modelParams[i].treeAgePr.priorParams[0];
                            }
                        else if (modelParams[i].treeAgePr.prior == uniform)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbUniform;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioUniform;
                            modelParams[i].treeAgePr.min           = modelParams[i].treeAgePr.priorParams[0];
                            modelParams[i].treeAgePr.max           = modelParams[i].treeAgePr.priorParams[1];
                            }
                        else if (modelParams[i].treeAgePr.prior == offsetExponential)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbOffsetExponential_Param_Offset_Mean;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioOffsetExponential_Param_Offset_Mean;
                            modelParams[i].treeAgePr.min           = modelParams[i].treeAgePr.priorParams[0];
                            modelParams[i].treeAgePr.max           = POS_INFINITY;
                            }
                        else if (modelParams[i].treeAgePr.prior == truncatedNormal)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbTruncatedNormal_Param_Trunc_Mean_Sd;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioTruncatedNormal_Param_Trunc_Mean_Sd;
                            modelParams[i].treeAgePr.min           = modelParams[i].treeAgePr.priorParams[0];
                            modelParams[i].treeAgePr.max           = POS_INFINITY;
                            }
                        else if (modelParams[i].treeAgePr.prior == logNormal)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbLognormal_Param_Mean_Sd;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioLognormal_Param_Mean_Sd;
                            modelParams[i].treeAgePr.min           = 0.0;
                            modelParams[i].treeAgePr.max           = POS_INFINITY;
                            }
                        else if (modelParams[i].treeAgePr.prior == offsetLogNormal)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbOffsetLognormal_Param_Offset_Mean_Sd;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioOffsetLognormal_Param_Offset_Mean_Sd;
                            modelParams[i].treeAgePr.min           = modelParams[i].treeAgePr.priorParams[0];
                            modelParams[i].treeAgePr.max           = POS_INFINITY;
                            }
                        else if (modelParams[i].treeAgePr.prior == standardGamma)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbGamma_Param_Mean_Sd;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioGamma_Param_Mean_Sd;
                            modelParams[i].treeAgePr.min           = 0.0;
                            modelParams[i].treeAgePr.max           = POS_INFINITY;
                            }
                        else if (modelParams[i].treeAgePr.prior == offsetGamma)
                            {
                            modelParams[i].treeAgePr.LnPriorProb   = &LnPriorProbOffsetGamma_Param_Offset_Mean_Sd;
                            modelParams[i].treeAgePr.LnPriorRatio  = &LnProbRatioOffsetGamma_Param_Offset_Mean_Sd;
                            modelParams[i].treeAgePr.min           = modelParams[i].treeAgePr.priorParams[0];
                            modelParams[i].treeAgePr.max           = POS_INFINITY;
                            }
                        }
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Clockratepr (clockRatePr) ***************************************************/
        else if (!strcmp(parmName, "Clockratepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                foundDash = NO;
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].clockRatePr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Clockratepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                expecting |= Expecting(DASH);   /* negative numbers possible */
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                if (foundDash == YES)
                    {
                    foundDash = NO;
                    tempD *= -1.0;
                    }
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].clockRatePr,"Normal"))
                            {
                            if (tempD <= 0.0)
                                {
                                if (numVars[i] == 0)
                                    YvyraPrint ("%s   Mean of the normal must be positive\n", spacer);
                                else if (numVars[i] == 1)
                                    YvyraPrint ("%s   Standard deviation of the normal must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].clockRateNormal[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 || numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Clockratepr to Normal(%1.6lf,%1.6lf)\n", spacer, modelParams[i].clockRateNormal[0], modelParams[i].clockRateNormal[1]);
                                else
                                    YvyraPrint ("%s   Setting Clockratepr to Normal(%1.6lf,%1.6lf) for partition %d\n", spacer, modelParams[i].clockRateNormal[0], modelParams[i].clockRateNormal[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].clockRatePr,"Lognormal"))
                            {
                            modelParams[i].clockRateLognormal[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 || numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Clockratepr to Lognormal(%1.2lf,%1.2lf)\n", spacer, modelParams[i].clockRateLognormal[0], modelParams[i].clockRateLognormal[1]);
                                else
                                    YvyraPrint ("%s   Setting Clockratepr to Lognormal(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].clockRateLognormal[0], modelParams[i].clockRateLognormal[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].clockRatePr,"Exponential"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Rate of the exponential must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].clockRateExp = tempD;
                            if (nApplied == 0 || numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Clockratepr to Exponential(%1.2lf)\n", spacer, modelParams[i].clockRateExp);
                            else
                                YvyraPrint ("%s   Setting Clockratepr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].clockRateExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].clockRatePr,"Gamma"))
                            {
                            if (tempD <= 0.0)
                                {
                                if (numVars[i] == 0)
                                    YvyraPrint ("%s   Shape of the gamma must be positive\n", spacer);
                                else if (numVars[i] == 1)
                                    YvyraPrint ("%s   Rate (inverse scale) of the gamma must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].clockRateGamma[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 || numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Clockratepr to Gamma(%1.2lf,%1.2lf)\n", spacer, modelParams[i].clockRateGamma[0], modelParams[i].clockRateGamma[1]);
                                else
                                    YvyraPrint ("%s   Setting Clockratepr to Gamma(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].clockRateGamma[0], modelParams[i].clockRateGamma[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].clockRatePr,"Fixed"))
                            {
                            if (tempD <= 0.0)
                                {
                                YvyraPrint ("%s   Fixed clock rate must be positive\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].clockRateFix = tempD;
                            if (nApplied == 0 || numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Clockratepr to Fixed(%1.6lf)\n", spacer, modelParams[i].clockRateFix);
                            else
                                YvyraPrint ("%s   Setting Clockratepr to Fixed(%1.6lf) for partition %d\n", spacer, modelParams[i].clockRateFix, i+1);
                            for (k=0; k<numGlobalChains; k++)
                                {
                                if (UpdateClockRate(tempD, k) == ERROR) 
                                    return (ERROR);
                                }
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Clockvarpr (clockVarPr) *****************************************************/
        else if (!strcmp(parmName, "Clockvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (strcmp(tempStr, "Bm") == 0)
                        strcpy(tempStr, "TK02");
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].clockVarPr, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Clockvarpr to %s\n", spacer, modelParams[i].clockVarPr);
                            else
                                YvyraPrint ("%s   Setting Clockvarpr to %s for partition %d\n", spacer, modelParams[i].clockVarPr, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Clockvarpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Growthpr (growthPr) *********************************************************/
        else if (!strcmp(parmName, "Growthpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                isNegative = NO;
                }
            else if (expecting == Expecting(ALPHA))
                {
                isNegative = NO;
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].growthPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Growthpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER) | Expecting(DASH);
                }
            else if (expecting == Expecting(DASH))
                {
                expecting  = Expecting(NUMBER);
                isNegative = YES;
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) &&
                        (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                        {
                        if (!strcmp(modelParams[i].growthPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (isNegative == YES)
                                tempD *= -1.0;
                            modelParams[i].growthUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].growthUni[0] >= modelParams[i].growthUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Growthpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].growthUni[0], modelParams[i].growthUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Growthpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].growthUni[0], modelParams[i].growthUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].growthPr,"Normal"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (isNegative == YES)
                                tempD *= -1.0;
                            modelParams[i].growthNorm[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].growthNorm[1] < 0.0)
                                    {
                                    YvyraPrint ("%s   Variance for normal distribution should be greater than zero\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Growthpr to Normal(%1.2lf,%1.2lf)\n", spacer, modelParams[i].growthNorm[0], modelParams[i].growthNorm[1]);
                                else
                                    YvyraPrint ("%s   Setting Growthpr to Normal(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].growthNorm[0], modelParams[i].growthNorm[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].growthPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].growthExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Growthpr to Exponential(%1.2lf)\n", spacer, modelParams[i].growthExp);
                            else
                                YvyraPrint ("%s   Setting Growthpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].growthExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].growthPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (isNegative == YES)
                                tempD *= -1.0;
                            modelParams[i].growthFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Growthpr to Fixed(%1.2lf)\n", spacer, modelParams[i].growthFix);
                            else
                                YvyraPrint ("%s   Setting Growthpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].growthFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                isNegative = NO;
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Popsizepr (popSizePr) *******************************************************/
        else if (!strcmp(parmName, "Popsizepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].popSizePr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Popsizepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) &&
                        (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
                        {
                        if (!strcmp(modelParams[i].popSizePr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].popSizeUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                {
                                if (tempD <= 0.0)
                                    {
                                    YvyraPrint ("%s   Lower value for Uniform must be a positive number\n", spacer);
                                    return (ERROR);
                                    }
                                expecting  = Expecting(COMMA);
                                }
                            else
                                {
                                if (modelParams[i].popSizeUni[0] >= modelParams[i].popSizeUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for Uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Popsizepr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].popSizeUni[0], modelParams[i].popSizeUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Popsizepr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].popSizeUni[0], modelParams[i].popSizeUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].popSizePr,"Lognormal"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].popSizeLognormal[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                {
                                expecting  = Expecting(COMMA);
                                }
                            else
                                {
                                if (modelParams[i].popSizeLognormal[1] <= 0.0)
                                    {
                                    YvyraPrint ("%s   Standard deviation of Lognormal must be a positive number\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Popsizepr to Lognormal(%1.2lf,%1.2lf)\n", spacer, modelParams[i].popSizeLognormal[0], modelParams[i].popSizeLognormal[1]);
                                else
                                    YvyraPrint ("%s   Setting Popsizepr to Lognormal(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].popSizeLognormal[0], modelParams[i].popSizeLognormal[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].popSizePr,"Normal"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].popSizeNormal[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                {
                                if (modelParams[i].popSizeNormal[0] <= 0.0)
                                    {
                                    YvyraPrint ("%s   Mean of Truncated Normal must be a positive number\n", spacer);
                                    return (ERROR);
                                    }
                                expecting  = Expecting(COMMA);
                                }
                            else
                                {
                                if (modelParams[i].popSizeNormal[1] <= 0.0)
                                    {
                                    YvyraPrint ("%s   Standard deviation of Truncated Normal must be a positive number\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Popsizepr to Truncated Normal(%1.2lf,%1.2lf)\n", spacer, modelParams[i].popSizeNormal[0], modelParams[i].popSizeNormal[1]);
                                else
                                    YvyraPrint ("%s   Setting Popsizepr to Truncated Normal(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].popSizeNormal[0], modelParams[i].popSizeNormal[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].popSizePr,"Gamma"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].popSizeGamma[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                {
                                if (modelParams[i].popSizeGamma[0] <= 0.0)
                                    {
                                    YvyraPrint ("%s   Shape (alpha) of Gamma must be a positive number\n", spacer);
                                    return (ERROR);
                                    }
                                expecting  = Expecting(COMMA);
                                }
                            else
                                {
                                if (modelParams[i].popSizeGamma[1] <= 0.0)
                                    {
                                    YvyraPrint ("%s   Rate (beta) of Gamma must be a positive number\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Popsizepr to Gamma(%1.2lf,%1.2lf)\n", spacer, modelParams[i].popSizeGamma[0], modelParams[i].popSizeGamma[1]);
                                else
                                    YvyraPrint ("%s   Setting Popsizepr to Gamma(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].popSizeGamma[0], modelParams[i].popSizeGamma[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].popSizePr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (AreDoublesEqual(tempD, 0.0, ETA)==YES)
                                {
                                YvyraPrint ("%s   Popsizepr cannot be fixed to 0.0\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].popSizeFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Popsizepr to Fixed(%1.5lf)\n", spacer, modelParams[i].popSizeFix);
                            else
                                YvyraPrint ("%s   Setting Popsizepr to Fixed(%1.5lf) for partition %d\n", spacer, modelParams[i].popSizeFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Popvarpr (popVarPr) *********************************************************/
        else if (!strcmp(parmName, "Popvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].popVarPr, tempStr);
                        
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Popvarpr to %s\n", spacer, modelParams[i].popVarPr);
                            else
                                YvyraPrint ("%s   Setting Popvarpr to %s for partition %d\n", spacer, modelParams[i].popVarPr, i+1);
                            }
                        }
                    expecting  = Expecting(PARAMETER) | Expecting(SEMICOLON);
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Popvarpr argument\n", spacer);
                    return (ERROR);
                    }
                }
            else
                return (ERROR);
            }
        /* set Compound Poisson Process rate prior (cppRatePr) *****************************/
        else if (!strcmp(parmName, "Cppratepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].cppRatePr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Cppratepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].cppRatePr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].cppRateExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Cppratepr to Exponential(%1.2lf)\n", spacer, modelParams[i].cppRateExp);
                            else
                                YvyraPrint ("%s   Setting Cppratepr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].cppRateExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].cppRatePr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < CPPLAMBDA_MIN || tempD > CPPLAMBDA_MAX)
                                {
                                YvyraPrint ("%s   CPP rate must be in the range %f - %f\n", spacer, CPPLAMBDA_MIN, CPPLAMBDA_MAX);
                                return (ERROR);
                                }
                            modelParams[i].cppRateFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Cppratepr to Fixed(%1.2lf)\n", spacer, modelParams[i].cppRateFix);
                            else
                                YvyraPrint ("%s   Setting Cppratepr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].cppRateFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Compound Poisson Process rate multiplier standard deviation (log scale) *****/
        else if (!strcmp(parmName, "Cppmultdevpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].cppMultDevPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Cppmultdevpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].cppMultDevPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < POS_MIN || tempD > POS_INFINITY)
                                {
                                YvyraPrint ("%s   The standard deviation of rate multipliers must be in the range %f - %f\n", spacer, POS_MIN, POS_INFINITY);
                                return (ERROR);
                                }
                            modelParams[i].cppMultDevFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Cppmultdevpr to Fixed(%1.2lf)\n", spacer, modelParams[i].cppMultDevFix);
                            else
                                YvyraPrint ("%s   Setting Cppmultdevpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].cppMultDevFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set prior for variance of lognormal of autocorrelated rates (tk02varPr) *********/
        else if (!strcmp(parmName, "TK02varpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].tk02varPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid TK02varpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].tk02varPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].tk02varUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].tk02varUni[0] >= modelParams[i].tk02varUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting TK02varpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].tk02varUni[0], modelParams[i].tk02varUni[1]);
                                else
                                    YvyraPrint ("%s   Setting TK02varpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].tk02varUni[0], modelParams[i].tk02varUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].tk02varPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].tk02varExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting TK02varpr to Exponential(%1.2lf)\n", spacer, modelParams[i].tk02varExp);
                            else
                                YvyraPrint ("%s   Setting TK02varpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].tk02varExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].tk02varPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < TK02VAR_MIN || tempD > TK02VAR_MAX)
                                {
                                YvyraPrint ("%s   Ratevar (nu) must be in the range %f - %f\n", spacer, TK02VAR_MIN, TK02VAR_MAX);
                                return (ERROR);
                                }
                            modelParams[i].tk02varFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting TK02varpr to Fixed(%1.2lf)\n", spacer, modelParams[i].tk02varFix);
                            else
                                YvyraPrint ("%s   Setting TK02varpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].tk02varFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set prior for variance of white noise rates (wnvarPr) ***************************/
        else if (!strcmp(parmName, "WNvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].wnvarPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid WNvarpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].wnvarPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].wnvarUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].wnvarUni[0] >= modelParams[i].wnvarUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting WNvarpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].wnvarUni[0], modelParams[i].wnvarUni[1]);
                                else
                                    YvyraPrint ("%s   Setting WNvarpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].wnvarUni[0], modelParams[i].wnvarUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].wnvarPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].wnvarExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting WNvarpr to Exponential(%1.2lf)\n", spacer, modelParams[i].wnvarExp);
                            else
                                YvyraPrint ("%s   Setting WNvarpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].wnvarExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].wnvarPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < WNVAR_MIN || tempD > WNVAR_MAX)
                                {
                                YvyraPrint ("%s   WNvar must be in the range %f - %f\n", spacer, WNVAR_MIN, WNVAR_MAX);
                                return (ERROR);
                                }
                            modelParams[i].wnvarFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting WNvarpr to Fixed(%1.2lf)\n", spacer, modelParams[i].wnvarFix);
                            else
                                YvyraPrint ("%s   Setting WNvarpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].wnvarFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set prior for variance of lognormal of independent rates (ilnvarPr) *************/
        else if (!strcmp(parmName, "ILNvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].ilnvarPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid ILNvarpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].ilnvarPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].ilnvarUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].ilnvarUni[0] >= modelParams[i].ilnvarUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting ILNvarpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].ilnvarUni[0], modelParams[i].ilnvarUni[1]);
                                else
                                    YvyraPrint ("%s   Setting ILNvarpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].ilnvarUni[0], modelParams[i].ilnvarUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].ilnvarPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].ilnvarExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting ILNvarpr to Exponential(%1.2lf)\n", spacer, modelParams[i].ilnvarExp);
                            else
                                YvyraPrint ("%s   Setting ILNvarpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].ilnvarExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].ilnvarPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < INDRVAR_MIN || tempD > INDRVAR_MAX)
                                {
                                YvyraPrint ("%s   ILNvar must be in the range %f - %f\n", spacer, INDRVAR_MIN, INDRVAR_MAX);
                                return (ERROR);
                                }
                            modelParams[i].ilnvarFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting ILNvarpr to Fixed(%1.2lf)\n", spacer, modelParams[i].ilnvarFix);
                            else
                                YvyraPrint ("%s   Setting ILNvarpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].ilnvarFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set prior for variance of independent branch rate gamma distribution (igrvarPr) */
        else if (!strcmp(parmName, "IGRvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].igrvarPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid IGRvarpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].igrvarPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].igrvarUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].igrvarUni[0] >= modelParams[i].igrvarUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting IGRvarpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].igrvarUni[0], modelParams[i].igrvarUni[1]);
                                else
                                    YvyraPrint ("%s   Setting IGRvarpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].igrvarUni[0], modelParams[i].igrvarUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].igrvarPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].igrvarExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting IGRvarpr to Exponential(%1.2lf)\n", spacer, modelParams[i].igrvarExp);
                            else
                                YvyraPrint ("%s   Setting IGRvarpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].igrvarExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].igrvarPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < INDRVAR_MIN || tempD > INDRVAR_MAX)
                                {
                                YvyraPrint ("%s   IGRvar must be in the range %f - %f\n", spacer, INDRVAR_MIN, INDRVAR_MAX);
                                return (ERROR);
                                }
                            modelParams[i].igrvarFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting IGRvarpr to Fixed(%1.2lf)\n", spacer, modelParams[i].igrvarFix);
                            else
                                YvyraPrint ("%s   Setting IGRvarpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].igrvarFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set prior for variance of mixed relaxed clock model (mixedvarPr) ****************/
        else if (!strcmp(parmName, "Mixedvarpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            strcpy(modelParams[i].mixedvarPr, tempStr);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Mixedvarpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if (activeParts[i] == YES || nApplied == 0)
                        {
                        if (!strcmp(modelParams[i].mixedvarPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].mixedvarUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].mixedvarUni[0] >= modelParams[i].mixedvarUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Mixedvarpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].mixedvarUni[0], modelParams[i].mixedvarUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Mixedvarpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].mixedvarUni[0], modelParams[i].mixedvarUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].mixedvarPr,"Exponential"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].mixedvarExp = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Mixedvarpr to Exponential(%1.2lf)\n", spacer, modelParams[i].mixedvarExp);
                            else
                                YvyraPrint ("%s   Setting Mixedvarpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].mixedvarExp, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        else if (!strcmp(modelParams[i].mixedvarPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (tempD < INDRVAR_MIN || tempD > INDRVAR_MAX)
                                {
                                YvyraPrint ("%s   Mixedvar must be in the range %f - %f\n", spacer, INDRVAR_MIN, INDRVAR_MAX);
                                return (ERROR);
                                }
                            modelParams[i].mixedvarFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Mixedvarpr to Fixed(%1.2lf)\n", spacer, modelParams[i].mixedvarFix);
                            else
                                YvyraPrint ("%s   Setting Mixedvarpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].mixedvarFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Browncorrpr (brownCorrPr) ***************************************************/
        else if (!strcmp(parmName, "Browncorrpr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    flag = 0;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
                            {
                            strcpy(modelParams[i].brownCorrPr, tempStr);
                            flag = 1;
                            }
                        }
                    if (flag == 0)
                        {
                        YvyraPrint ("%s   Warning: %s can be set only for partition containing CONTINUOUS data.\
                            Currently there is no active partition with such data.\n", spacer, parmName);
                        return (ERROR);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Browncorrpr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                foundDash = NO;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting(NUMBER) | Expecting(DASH);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
                        {
                        if (!strcmp(modelParams[i].brownCorrPr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (foundDash == YES)
                                tempD *= -1.0;
                            modelParams[i].brownCorrUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].brownCorrUni[0] >= modelParams[i].brownCorrUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].brownCorrUni[1] > 1.0)
                                    {
                                    YvyraPrint ("%s   Upper value for uniform should be less than or equal to 1.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].brownCorrUni[0] < -1.0)
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than or equal to -1.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Browncorrpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownCorrUni[0], modelParams[i].brownCorrUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Browncorrpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brownCorrUni[0], modelParams[i].brownCorrUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].brownCorrPr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            if (foundDash == YES)
                                tempD *= -1.0;
                            if (tempD > 1.0 || tempD < -1.0)
                                {
                                YvyraPrint ("%s   Value for Browncorrpr should be in the interval (-1, +1)\n", spacer);
                                return (ERROR);
                                }
                            modelParams[i].brownCorrFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Browncorrpr to Fixed(%1.2lf)\n", spacer, modelParams[i].brownCorrFix);
                            else
                                YvyraPrint ("%s   Setting Browncorrpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].brownCorrFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                foundDash = NO;
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER) | Expecting(DASH);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Brownscalepr (brownScalePr) *************************************************/
        else if (!strcmp(parmName, "Brownscalepr"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    flag = 0;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
                            {
                            strcpy(modelParams[i].brownScalePr, tempStr);
                            flag = 1;
                            }
                        }
                    if (flag == 0)
                            {
                            YvyraPrint ("%s   Warning: %s can be set only for partition containing CONTINUOUS data.\
                            Currently there is no active partition with such data. ", spacer, parmName);
                            return (ERROR);
                            }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid Brownscalepr argument\n", spacer);
                    return (ERROR);
                    }
                expecting  = Expecting(LEFTPAR);
                for (i=0; i<numCurrentDivisions; i++)
                    numVars[i] = 0;
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                nApplied = NumActiveParts ();
                for (i=0; i<numCurrentDivisions; i++)
                    {
                    if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
                        {
                        if (!strcmp(modelParams[i].brownScalePr,"Uniform"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].brownScaleUni[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (modelParams[i].brownScaleUni[0] <= 0.0)
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than 0.0\n", spacer);
                                    return (ERROR);
                                    }
                                if (modelParams[i].brownScaleUni[0] >= modelParams[i].brownScaleUni[1])
                                    {
                                    YvyraPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
                                    return (ERROR);
                                    }
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brownscalepr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownScaleUni[0], modelParams[i].brownScaleUni[1]);
                                else
                                    YvyraPrint ("%s   Setting Brownscalepr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brownScaleUni[0], modelParams[i].brownScaleUni[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].brownScalePr,"Gamma"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].brownScaleGamma[numVars[i]++] = tempD;
                            if (numVars[i] == 1)
                                expecting  = Expecting(COMMA);
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Setting Brownscalepr to Gamma(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownScaleGamma[0], modelParams[i].brownScaleGamma[1]);
                                else
                                    YvyraPrint ("%s   Setting Brownscalepr to Gamma(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brownScaleGamma[0], modelParams[i].brownScaleGamma[1], i+1);
                                expecting  = Expecting(RIGHTPAR);
                                }
                            }
                        else if (!strcmp(modelParams[i].brownScalePr,"Fixed"))
                            {
                            sscanf (tkn, "%lf", &tempD);
                            modelParams[i].brownScaleFix = tempD;
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting Brownscalepr to Fixed(%1.2lf)\n", spacer, modelParams[i].brownScaleFix);
                            else
                                YvyraPrint ("%s   Setting Brownscalepr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].brownScaleFix, i+1);
                            expecting  = Expecting(RIGHTPAR);
                            }
                        }
                    }
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }

        else
            return (ERROR);
        }

    return (NO_ERROR);
}


int DoQuit (void)
{
    int         i;
    char        tempName[100];
    
    /* free information for model and matrix */
    FreeModel ();
    FreeMatrix ();
    
    SafeFclose (&logFileFp);
    logToFile = NO;
    
    /* check to see if any memory has not been freed */
    for (i=0; i<NUM_ALLOCS; i++)
        {
        if (memAllocs[i] == YES)
            {
            YvyraPrint ("   WARNING: Memory (%d) has not been freed\n", i);
            if (mode == INTERACTIVE && quitOnError == NO)
                {
                YvyraPrint ("%s   Hit return key to continue  ", spacer);
                fflush (stdin);
                if (fgets (tempName, 100, stdin) == NULL)
                    {
                    printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
                    }
                }
            }
        }
    
    /* free modelIndicatorParams and modelElementNames */
    for (i=0; i<203; i++)
        free (modelElementNames[1][i]);
    for (i=0; i<3; i++)
        free (modelElementNames[i]);
    free (modelElementNames);
    free (modelIndicatorParams);
    
    YvyraPrint ("   Quitting program\n\n");
    
    /* If we quit while reading a yvyra block, then we need to make certain
     that we return a NO_ERROR_QUIT so we can break out of DoExecute cleanly,
     and dealloc "s" there. */
    if (inYvyraBlock == YES)
        {
        inYvyraBlock = NO;
        return (NO_ERROR_QUIT);
        }
    
    return (NO_ERROR);
}


int DoReport (void)
{
    /* TODO: smart update */
    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);
    return (NO_ERROR);
}


int DoReportParm (char *parmName, char *tkn)
{
    int         i, tempInt, nApplied;
    char        tempStr[100];

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before the report settings can be altered\n", spacer);
        return (ERROR);
        }
    if (inValidCommand == YES)
        {
        for (i=0; i<numCurrentDivisions; i++)
            activeParts[i] = NO;
        inValidCommand = NO;
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Applyto (Applyto) *************************************************************/
        if (!strcmp(parmName, "Applyto"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(LEFTPAR);
            else if (expecting == Expecting(LEFTPAR))
                {
                for (i=0; i<numCurrentDivisions; i++)
                    activeParts[i] = NO;
                fromI = toJ = -1;
                foundDash = NO;
                expecting = Expecting(NUMBER) | Expecting(ALPHA);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                if (fromI != -1)
                    activeParts[fromI-1] = YES;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(COMMA))
                {
                foundComma = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (IsSame ("All", tkn) == DIFFERENT)
                    {
                    YvyraPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
                    return (ERROR);
                    }
                for (i=0; i<numCurrentDivisions; i++)
                    activeParts[i] = YES;
                expecting  = Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt > numCurrentDivisions)
                    {
                    YvyraPrint ("%s   Partition delimiter is too large\n", spacer);
                    return (ERROR);
                    }
                if (fromI == -1)
                    fromI = tempInt;
                else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
                    {
                    toJ = tempInt;
                    for (i=fromI-1; i<toJ; i++)
                        activeParts[i] = YES;
                    fromI = toJ = -1;
                    foundDash = NO;
                    }
                else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
                    {
                    activeParts[fromI-1] = YES;
                    fromI = tempInt;
                    foundComma = NO;
                    }
                    
                expecting  = Expecting(COMMA);
                expecting |= Expecting(DASH);
                expecting |= Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else
                return (ERROR);
            }
        /* set report format of tratio ***************************************************/
        else if (!strcmp(parmName, "Tratio"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    tempInt = NO;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        /* check that data type is correct; we do not know yet if the user will specify
                        a nst=2 model so we cannot check that tratio is an active parameter */
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
                            {
                            strcpy(modelParams[i].tratioFormat, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting transition/transversion rate ratio (tratio) format to %s\n", spacer, modelParams[i].tratioFormat);
                            else
                                YvyraPrint ("%s   Setting transition/transversion rate ratio (tratio) format to %s for partition %d\n", spacer, modelParams[i].tratioFormat, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid transition/transversion rate ratio (tratio) format \n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set report format of revmat ***************************************************/
        else if (!strcmp(parmName, "Revmat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    tempInt = NO;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        /* check that data type is correct; we do not know yet if the user will specify
                           a nst=6 model so we cannot check that revmat is an active parameter */
                        if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
                            {
                            strcpy(modelParams[i].revmatFormat, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting reversible rate matrix (revmat) format to %s\n", spacer, modelParams[i].revmatFormat);
                            else
                                YvyraPrint ("%s   Setting reversible rate matrix (revmat) format to %s for partition %d\n", spacer, modelParams[i].revmatFormat, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid reversible rate matrix (revmat) format \n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set report format of ratemult *************************************************/
        else if (!strcmp(parmName, "Ratemult"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    tempInt = NO;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        /* we do not know yet if the user will specify variable rates across partitions 
                           so only check that we have more than one partition in the model */
                        if ((activeParts[i] == YES || nApplied == 0) && numCurrentDivisions > 1)
                            {
                            strcpy(modelParams[i].ratemultFormat, tempStr);
                            if (nApplied == 0 && numCurrentDivisions == 1)
                                YvyraPrint ("%s   Setting rate multiplier (ratemult) format to %s\n", spacer, modelParams[i].ratemultFormat);
                            else
                                YvyraPrint ("%s   Setting rate multiplier (ratemult) format to %s for partition %d\n", spacer, modelParams[i].ratemultFormat, i+1);
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid rate multiplier (ratemult) format \n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set report format of tree ***************************************************/
        else if (!strcmp(parmName, "Tree"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    tempInt = NO;
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        strcpy(modelParams[i].treeFormat, tempStr);
                        if (nApplied == 0 && numCurrentDivisions == 1)
                            YvyraPrint ("%s   Setting tree report format to %s\n", spacer, modelParams[i].treeFormat);
                        else
                            YvyraPrint ("%s   Setting tree report format to %s for partition %d\n", spacer, modelParams[i].treeFormat, i+1);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid tree report format \n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set inferancstates ***********************************************************/
        else if (!strcmp(parmName, "Ancstates"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy(modelParams[i].inferAncStates,tempStr);
                            if (!strcmp(tempStr,"Yes"))
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Reporting ancestral states (if applicable)\n", spacer);
                                else
                                    YvyraPrint ("%s   Reporting ancestral states for partition %d (if applicable)\n", spacer, i+1);
                                }
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Not reporting ancestral states\n", spacer);
                                else
                                    YvyraPrint ("%s   Not reporting ancestral states for partition %d\n", spacer, i+1);
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid ancstates option\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set inferSiteRates ***************************************************************/
        else if (!strcmp(parmName, "Siterates"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy (modelParams[i].inferSiteRates, tempStr);
                            if (!strcmp(tempStr,"Yes"))
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Reporting site rates (if applicable)\n", spacer);
                                else
                                    YvyraPrint ("%s   Reporting site rates for partition %d (if applicable)\n", spacer, i+1);
                                }
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Not reporting site rates\n", spacer);
                                else
                                    YvyraPrint ("%s   Not reporting site rates for partition %d\n", spacer, i+1);
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid siterates option\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set inferSiteLikes ***************************************************************/
        else if (!strcmp(parmName, "Sitelikes"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy (modelParams[i].inferSiteLikes, tempStr);
                            if (!strcmp(tempStr,"Yes"))
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Reporting per-site log-likelihoods\n", spacer);
                                else
                                    YvyraPrint ("%s   Reporting per-site log-likelihoods for partition %d\n", spacer, i+1);
                                }
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Not reporting per-site log-likelihoods\n", spacer);
                                else
                                    YvyraPrint ("%s   Not reporting per-site log-likelihoods for partition %d\n", spacer, i+1);
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid sitelikes option\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set inferpossel *************************************************************/
        else if (!strcmp(parmName, "Possel"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy (modelParams[i].inferPosSel, tempStr);
                            if (!strcmp(tempStr, "Yes"))
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Reporting positive selection (if applicable)\n", spacer);
                                else
                                    YvyraPrint ("%s   Reporting positive selection for partition %d (if applicable)\n", spacer, i+1);
                                }
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Not reporting positive selection\n", spacer);
                                else
                                    YvyraPrint ("%s   Not reporting positive selection for partition %d\n", spacer, i+1);
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid possel option\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set inferSiteOmegas *************************************************************/
        else if (!strcmp(parmName, "Siteomega"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting (ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    nApplied = NumActiveParts ();
                    for (i=0; i<numCurrentDivisions; i++)
                        {
                        if (activeParts[i] == YES || nApplied == 0)
                            {
                            strcpy (modelParams[i].inferSiteOmegas, tempStr);
                            if (!strcmp(tempStr, "Yes"))
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Reporting site omega values (if applicable)\n", spacer);
                                else
                                    YvyraPrint ("%s   Reporting site omega values for partition %d (if applicable)\n", spacer, i+1);
                                }
                            else
                                {
                                if (nApplied == 0 && numCurrentDivisions == 1)
                                    YvyraPrint ("%s   Not reporting site omega values\n", spacer);
                                else
                                    YvyraPrint ("%s   Not reporting site omega values for partition %d\n", spacer, i+1);
                                }
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid siteomega option\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else
            {
            return (ERROR);
            }
        }

    return (NO_ERROR);      
}


int DoStartvals (void)
{
    YvyraPrint ("%s   Successfully set starting values\n", spacer);
    /*
    for (i=0; i<numParams; i++)
        assert (IsTreeConsistent(&params[i], 0, 0) == YES);
    */
    return (NO_ERROR);
}


int DoStartvalsParm (char *parmName, char *tkn)
{
    int                 i, j, k, nMatches, tempInt=0, treeIndex, chainId, ret;
    YFlt              tempFloat=0.0, *value, *subValue;
    Tree                *theTree, *usrTree;
    PolyTree            *thePolyTree;
    YFlt              minRate=0.0, maxRate=0.0, clockRate;
    static Param        *param = NULL;
    static YFlt       *theValue, theValueMin, theValueMax;
    static int          useSubvalues, useStdStateFreqs, useIntValues, numExpectedValues, nValuesRead, runIndex, chainIndex, foundName, foundDash;
    static char         *temp=NULL, *tempName=NULL;

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before starting values can be changed\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* we expect a parameter name with possible run and chain specification as follows:
               <param_name>(<run>,<chain>)=(<number>,...)   -- apply to run <run> and chain <chain>
               <param_name>(,<chain>)=(<number>,...)        -- apply to chain <chain> for all runs
               <param_name>(<run>,)=(<number>,...)          -- apply to run <run> for all chains

               topology and branch length parameters are specified like
               <param_name>(<run>,<chain>)=<tree_name>|<Newick_tree_spec>

               parameter names will often be followed by partition specifiers like:
               pi{all}
               pinvar{1,4,5}
               so we need to assemble the parameter name from several tokens that are parsed out separately;
               here we receive only the first part (before the left curly, if present)
            */
            
            /* copy to local parameter name */
            SafeStrcpy (&tempName, tkn);
            param = NULL;
            runIndex = chainIndex = -1;
            useSubvalues = NO;
            useIntValues = NO;
            useStdStateFreqs = NO;
            foundComma = foundEqual = foundName = foundDash = NO;
            expecting = Expecting(LEFTCURL) | Expecting(LEFTPAR) |  Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (param == NULL)
            {
            /* we are still assembling the parameter name */
            SafeStrcat (&tempName, tkn);
            expecting = Expecting(RIGHTCURL) | Expecting(COMMA) | Expecting(LEFTPAR) | Expecting(EQUALSIGN);
            }
        else
            {
            /* we have a tree name; now find the tree based on name, case insensitive */
            if (GetUserTreeFromName (&treeIndex, tkn) == ERROR || treeIndex == -1)
                {
                YvyraPrint ("%s   Error in finding user tree\n", spacer);
                return (ERROR);
                }

            /* set the tree parameter */
            for (i=0; i<chainParams.numRuns; i++)
                {
                if (runIndex != -1 && i != runIndex)
                    continue;
                for (j=0; j<chainParams.numChains; j++)
                    {
                    if (chainIndex != -1 && j != chainIndex)
                        continue;
                    chainId = i*chainParams.numChains + j;
                    if (param->paramType != P_POPSIZE)
                        {
                        if (param->paramType == P_TOPOLOGY || param->paramType == P_BRLENS || param->paramType == P_SPECIESTREE)
                            {
                            /* topology or brlens or speciestree params */
                            theTree = GetTree (param, chainId, 0);
                            usrTree = GetTree (param, chainId, 1); /* use as scratch space */
                            }
                        else
                            {
                            /* relaxed clock params */
                            theTree = GetTree (modelSettings[param->relParts[0]].brlens, chainId, 0);
                            usrTree = GetTree (modelSettings[param->relParts[0]].brlens, chainId, 1);
                            }
                        CopyToTreeFromTree (usrTree, theTree);
                        if (param->paramType == P_SPECIESTREE)
                            thePolyTree = AllocatePolyTree(numSpecies);
                        else
                            thePolyTree = AllocatePolyTree (numTaxa);
                        CopyToPolyTreeFromPolyTree (thePolyTree, userTree[treeIndex]);
                        if (param->paramType == P_SPECIESTREE)
                            {
                            ResetIntNodeIndices(thePolyTree);
                            }
                        else
                            {
                            PrunePolyTree (thePolyTree);
                            ResetTipIndices(thePolyTree);
                            }
                        RandResolve (NULL, thePolyTree, &globalSeed, theTree->isRooted);
                        GetPolyDownPass(thePolyTree);
                        ResetIntNodeIndices(thePolyTree);
                        if (param->paramType == P_SPECIESTREE)
                            {
                            ret=CopyToSpeciesTreeFromPolyTree (usrTree, thePolyTree);
                            }
                        else
                            ret=CopyToTreeFromPolyTree (usrTree, thePolyTree);
                        FreePolyTree (thePolyTree);
                        if (ret==ERROR)
                            return ERROR;
                        }
                    else
                        {
                        /* param->paramType == P_POPSIZE */
                        theTree = GetTree (modelSettings[param->relParts[0]].speciesTree, chainId, 0);
                        usrTree = GetTree (modelSettings[param->relParts[0]].speciesTree, chainId, 1);
                        CopyToTreeFromTree (usrTree, theTree);
                        thePolyTree = AllocatePolyTree(numSpecies);
                        CopyToPolyTreeFromPolyTree (thePolyTree, userTree[treeIndex]);
                        ResetIntNodeIndices(thePolyTree);
                        RandResolve (NULL, thePolyTree, &globalSeed, theTree->isRooted);
                        CopyToSpeciesTreeFromPolyTree (usrTree, thePolyTree);
                        FreePolyTree (thePolyTree);
                        }
                    if (param->paramType == P_TOPOLOGY)
                        {
                        if (theTree->checkConstraints == YES && CheckSetConstraints (usrTree) == ERROR)
                            {
                            YvyraPrint ("%s   Could not set the constraints for topology parameter '%s'\n", spacer, param->name);
                            return (ERROR);
                            }
                        if (ResetTopologyFromTree (theTree, usrTree) == ERROR)
                            {
                            YvyraPrint ("%s   Could not set the topology parameter '%s'\n", spacer, param->name);
                            return (ERROR);
                            }
                        if (theTree->checkConstraints == YES && CheckSetConstraints (theTree)==ERROR)
                            {
                            YvyraPrint ("%s   Could not set the constraints for topology parameter '%s'\n", spacer, param->name);
                            return (ERROR);
                            }
                        FillTopologySubParams (param, chainId, 0, &globalSeed);
                        //YvyraPrint ("%s   Branch lengths and relaxed clock subparameters of a parameter '%s' are reset.\n", spacer, param->name);
                        /* FillSpeciesTreeParams removed: species tree not supported */
                        //assert (IsTreeConsistent(param, chainId, 0) == YES);
                        }
                    else if (param->paramType == P_BRLENS)
                        {
                        if (usrTree->allDownPass[0]->length == 0.0 && param->paramId != BRLENS_CLOCK_FOSSIL)
                            {
                            YvyraPrint ("%s   User tree '%s' does not have branch lengths so it cannot be used in setting parameter '%s'\n", spacer, userTree[treeIndex]->name, param->name);
                            return (ERROR);
                            }
                        if (AreTopologiesSame (theTree, usrTree) == NO)
                            {
                            YvyraPrint ("%s   Topology of user tree '%s' wrong in setting parameter '%s'\n", spacer, userTree[treeIndex]->name, param->name);
                            return (ERROR);
                            }
                        //assert (IsTreeConsistent(param, chainId, 0) == YES);
                        /* reset node depths to ensure that non-dated tips have node depth 0.0 */
                        /* if (usrTree->isClock == YES)
                            SetNodeDepths(usrTree);  */
                        if (ResetBrlensFromTree (theTree, usrTree) == ERROR)
                            {
                            YvyraPrint ("%s   Could not set parameter '%s' from user tree '%s'\n", spacer, param->name, userTree[treeIndex]->name);
                            return (ERROR);
                            }
                        if (theTree->isClock == YES && modelParams[theTree->relParts[0]].treeAgePr.prior == fixed)
                            {
                            // if (!strcmp(modelParams[theTree->relParts[0]].clockPr,"Uniform") ||
                            //    !strcmp(modelParams[theTree->relParts[0]].clockPr,"Birthdeath") ||
                            //    !strcmp(modelParams[theTree->relParts[0]].clockPr,"Fossilization"));
                            // We cannot check the consistency of root age and clock rate here because they can be set in any order.
                            // Defer this check to the CheckModel fxn, called just before starting the chain.
                            }
                        /* the test will find suitable clock rate and ages of nodes in theTree */
                        if (theTree->isClock == YES && IsClockSatisfied (theTree,0.001) == NO)
                            {
                            YvyraPrint ("%s   Non-calibrated tips are not at the same level after setting up starting tree branch lengths(%s) from user tree '%s'.\n",
                                          spacer, param->name, userTree[treeIndex]->name);
                            ShowNodes(theTree->root,0,YES);
                            return (ERROR);
                            }
                        if (theTree->isCalibrated == YES && IsCalibratedClockSatisfied (theTree, &minRate,&maxRate, 0.001) == NO)
                            {
                            YvyraPrint ("%s   Problem setting calibrated tree parameters\n", spacer);
                            return (ERROR);
                            }
                        if (theTree->isCalibrated == YES && !strcmp(modelParams[theTree->relParts[0]].clockRatePr, "Fixed"))
                            {
                            clockRate = modelParams[theTree->relParts[0]].clockRateFix;
                            if ((clockRate < minRate && AreDoublesEqual (clockRate, minRate , 0.001) == NO) ||
                                (clockRate > maxRate && AreDoublesEqual (clockRate, maxRate , 0.001) == NO))
                                {
                                YvyraPrint ("%s   Fixed branch lengths do not satisfy fixed clockrate\n", spacer);
                                return (ERROR);
                                }
                            }
                        theTree->fromUserTree=YES;
                        
                        FillBrlensSubParams (param, chainId, 0);
                        //YvyraPrint ("%s   Relaxed clock subparameters of a parameter '%s' are reset.\n", spacer, param->name);
                        //assert (IsTreeConsistent(param, chainId, 0) == YES);
                        /* FillSpeciesTreeParams removed: species tree not supported */
                        //assert (IsTreeConsistent(param, chainId, 0) == YES);
                        }
                    else if (param->paramType == P_CPPEVENTS || param->paramType == P_TK02BRANCHRATES || param->paramType == P_ILNBRANCHRATES ||
                             param->paramType == P_IGRBRANCHRATES || param->paramType == P_MIXEDBRCHRATES || param->paramType == P_WNBRANCHRATES)
                        {
                        if (theTree->isCalibrated == YES && theTree->fromUserTree == NO)
                            { /* if theTree is not set from user tree then we can not guarantee that branch lengths will stay the same
                                 by the time we start mcmc run because of clockrate adjustment. */
                            YvyraPrint ("%s    Set starting values for branch lengths first before setting starting values of relaxed parameters!\n", spacer);
                            return (ERROR);
                            }
                        if (theTree->isCalibrated == NO && IsClockSatisfied (usrTree, 0.001) == NO) // user tree is not calibrated so do not check it if calibration is in place
                            {
                            YvyraPrint ("%s   Branch lengths of the user tree '%s' do not satisfy clock in setting parameter '%s'\n", spacer, userTree[treeIndex], param->name);
                            ShowNodes(usrTree->root,0,YES);
                            return (ERROR);
                            }
                        if (AreTopologiesSame (theTree, usrTree) == NO)
                            {
                            YvyraPrint ("%s   Topology of user tree '%s' is wrong in setting parameter '%s'\n", spacer, userTree[treeIndex]->name, param->name);
                            return (ERROR);
                            }
                        if (SetRelaxedClockParam (param, chainId, 0, userTree[treeIndex]) == ERROR)
                            {
                            YvyraPrint ("%s   Could not set parameter '%s' from user tree '%s'\n", spacer, param->name, userTree[treeIndex]->name);
                            return (ERROR);
                            }
                        //assert (IsTreeConsistent(param, chainId, 0) == YES);
                        }
                    else if (param->paramType == P_POPSIZE)
                        {
                        if (AreTopologiesSame (theTree, usrTree) == NO)
                            {
                            YvyraPrint ("%s   Topology of user tree '%s' is wrong in setting parameter '%s'\n", spacer, userTree[treeIndex]->name, param->name);
                            return (ERROR);
                            }
                        if (SetPopSizeParam (param, chainId, 0, userTree[treeIndex]) == ERROR)
                            {
                            YvyraPrint ("%s   Could not set parameter '%s' from user tree '%s'\n", spacer, param->name, userTree[treeIndex]->name);
                            return (ERROR);
                            }
                        }
                    /* P_SPECIESTREE handling removed: not supported in yvyra */
                    }
                }
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(LEFTCURL))
        {
        /* we are still assembling the parameter name */
        SafeStrcat (&tempName, tkn);
        expecting = Expecting(NUMBER) | Expecting(ALPHA) | Expecting(LEFTPAR) | Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(RIGHTCURL))
        {
        /* we are still assembling the parameter name */
        SafeStrcat (&tempName, tkn);
        foundComma = NO; /*! if there was a comma in the partition part, we should reset this variable. Otherwise we can't parse something like A{1,2}(3,4) */
        expecting = Expecting(LEFTPAR) | Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        if (foundEqual == NO)
            {
            foundName = YES;    /* we found the name */
            /* we will be reading in run and chain indices */
            expecting = Expecting(NUMBER) | Expecting(COMMA);
            }
        else
            {
            expecting = Expecting(NUMBER);
            expecting |= Expecting(DASH);
            }
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundName == NO)
            {
            /* we are still assembling the parameter name */
            SafeStrcat (&tempName, tkn);
            expecting = Expecting(COMMA) | Expecting(LEFTPAR) | Expecting(RIGHTCURL) | Expecting(EQUALSIGN);
            }
        else if (foundEqual == YES)
            {
            theValueMin = param->min;
            theValueMax = param->max;

            /* we are reading in a parameter value */
            if (0) /* P_OMEGA removed */
                {
                /* continue with subvalues */
                nValuesRead = 0;
                numExpectedValues = param->nSubValues/2;
                useSubvalues = YES;
                theValueMin = ETA;
                theValueMax = 1.0;
                }
            if (0) /* P_OMEGA removed */
                {
                /* continue with subvalues */
                nValuesRead = 0;
                numExpectedValues = param->nSubValues/2;
                useSubvalues = YES;
                theValueMin = ETA;
                theValueMax = 1.0;
                }
            if (param->nIntValues > 0 && nValuesRead==numExpectedValues && useIntValues == NO)
                {
                /* continue with intValues */
                nValuesRead = 0;
                numExpectedValues = param->nIntValues;
                useIntValues = YES;
                }
            if (param->paramType==P_PI && modelSettings[param->relParts[0]].dataType == STANDARD && param->paramId != SYMPI_EQUAL
                && nValuesRead==numExpectedValues && useStdStateFreqs == NO)
                {
                /* we have read alpha_symdir, continue with multistate char state freqs */
                nValuesRead = 0;
                numExpectedValues = param->nStdStateFreqs;
                if (param->hasBinaryStd == YES)
                    numExpectedValues -= 2 * modelSettings[param->relParts[0]].numBetaCats;
                useStdStateFreqs = YES;
                theValueMin = ETA;
                theValueMax = 1.0;
                }
            nValuesRead++;
            if (nValuesRead > numExpectedValues)
                {
                if (param->nIntValues > 0)   
                    YvyraPrint ("%s   Only %d values were expected for parameter '%s'\n", spacer, param->nValues+param->nIntValues, param->name);
                else
                    YvyraPrint ("%s   Only %d values were expected for parameter '%s'\n", spacer, numExpectedValues, param->name);
                return (ERROR);
                }
            if (useIntValues == YES)
                sscanf (tkn, "%d", &tempInt);
            else
                sscanf (tkn, "%lf", &tempFloat);
            if (foundDash == YES)
                {
                if (useIntValues == YES)
                    tempInt = -tempInt;
                else
                    tempFloat = -tempFloat;
                foundDash = NO;
                }
            if (useIntValues == NO && (tempFloat < theValueMin || tempFloat > theValueMax))
                {
                YvyraPrint ("%s   The value is out of range (min = %lf; max = %lf)\n", spacer, theValueMin, theValueMax);
                return (ERROR);
                }
            for (i=0; i<chainParams.numRuns; i++)
                {
                if (runIndex != -1 && runIndex != i)
                    continue;
                for (j=0; j<chainParams.numChains; j++)
                    {
                    if (chainIndex != -1 && chainIndex != j)
                        continue;
                    if (useIntValues == YES)
                        {
                        GetParamIntVals (param, i*chainParams.numChains+j, 0)[nValuesRead-1] = tempInt;
                        }
                    else
                        {
                        if (useSubvalues == NO && useStdStateFreqs == NO)
                            theValue = GetParamVals (param, i*chainParams.numChains+j, 0);
                        else if (useSubvalues == YES)
                            theValue = GetParamSubVals (param, i*chainParams.numChains+j, 0);
                        else if (useStdStateFreqs == YES)
                            {
                            theValue = GetParamStdStateFreqs (param, i*chainParams.numChains+j, 0);
                            if (param->hasBinaryStd == YES)
                                theValue += 2 * modelSettings[param->relParts[0]].numBetaCats;
                            }
                        else
                            return (ERROR);
                        if (param->paramType == P_CLOCKRATE)
                            {
                            if (UpdateClockRate(tempFloat, i*chainParams.numChains+j) == ERROR) 
                                {
                                return (ERROR);
                                }
                            }
                        theValue[nValuesRead-1] = tempFloat;
                        }
                    }
                }
            expecting = Expecting (COMMA) | Expecting(RIGHTPAR);
            }
        else /* if (foundEqual == NO) */
            {
            sscanf (tkn, "%d", &tempInt);
            if (foundComma == NO)
                {
                if (tempInt <= 0 || tempInt > chainParams.numRuns)
                    {
                    YvyraPrint ("%s   Run index is out of range (min=1; max=%d)\n", spacer, chainParams.numRuns);
                    return (ERROR);
                    }
                runIndex = tempInt - 1;
                expecting = Expecting(COMMA);
                }
            else
                {
                if (tempInt <= 0 || tempInt > chainParams.numChains)
                    {
                    YvyraPrint ("%s   Chain index is out of range (min=1; max=%d)\n", spacer, chainParams.numChains);
                    return (ERROR);
                    }
                chainIndex = tempInt - 1;
                foundComma = NO;
                expecting = Expecting(RIGHTPAR);
                }
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        if (foundEqual == YES)
            {
            /* we expect another parameter value */
            expecting = Expecting(NUMBER);
            }
        else /* if (foundEqual == NO) */
            {
            /* we will be reading in chain index, if present */
            foundComma = YES;
            expecting = Expecting(RIGHTPAR) | Expecting(NUMBER); 
            /* if the comma is in a list of partitions (so between { and }) we have to add the comma to the parameter name */
            if (param == NULL && strchr(tempName, '}')==NULL && strchr(tempName, '{')!=NULL) 
              SafeStrcat (&tempName, ",");
            }
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        if (foundEqual == NO)
            {
            /* this is the end of the run and chain specification */
            expecting = Expecting(EQUALSIGN);
            }
        else /* if (foundEqual == YES) */
            {
            /* this is the end of the parameter values */
            if (nValuesRead != numExpectedValues)
                {
                YvyraPrint ("%s   Expected %d values but only found %d values for parameter '%s'\n", spacer, numExpectedValues, nValuesRead, param->name);
                return (ERROR);
                }
            /* Post processing needed for some parameters */
            /* FIXME: param is NULL here (from clang static analyzer) */
            if (param->paramType == P_SHAPE || param->paramType == P_CORREL)
                {
                for (i=0; i<chainParams.numRuns; i++)
                    {
                    if (runIndex != -1 && runIndex != i)
                        continue;
                    for (j=0; j<chainParams.numChains; j++)
                        {
                        if (chainIndex != -1 && chainIndex != j)
                            continue;
                        value = GetParamVals(param,i*chainParams.numChains+j,0);
                        subValue = GetParamSubVals(param,i*chainParams.numChains+j,0);
                        if (param->paramType == P_SHAPE && !strncmp(param->name, "Alpha", 5))
                            {
                            if (DiscreteGamma (subValue, value[0], value[0], param->nSubValues, 0) == ERROR)
                                return (ERROR);
                            }
                        else if (param->paramType == P_SHAPE && !strncmp(param->name, "Sigma", 5))
                            {
                            if( DiscreteLogNormal(subValue, value[0], param->nSubValues, 1) == ERROR)
                                return (ERROR);
                            }
                        else if (param->paramType == P_CORREL)
                            {
                            if (AutodGamma (subValue, value[0], (int)(sqrt(param->nSubValues) + 0.5)) == ERROR)
                                return (ERROR);
                            }
                        }
                    }
                }
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        foundEqual = YES;
        foundName = YES;

        /* we now know that the name is complete; try to find the parameter with this name (case insensitive) */
        /* FIXME: tempName is NULL? (from clang static analyzer) */
        for (i=0; i<(int)strlen(tempName); i++)
            tempName[i] = (char)(tolower(tempName[i]));
        
        /* first check exact matches */
        nMatches = j = 0;
        for (i=0; i<numParams; i++)
            {
            param = &params[i];
            SafeStrcpy (&temp, param->name);
            for (k=0; k<(int)(strlen(temp)); k++)
                temp[k] = (char)(tolower(temp[k]));
            if (strcmp(tempName,temp) == 0)
                {
                j = i;
                nMatches++;
                }
            }
        /* now check unambiguous abbreviation matches */
        if (nMatches == 0)
            {
            nMatches = j = 0;
            for (i=0; i<numParams; i++)
                {
                param = &params[i];
                SafeStrcpy (&temp, param->name);
                for (k=0; k<(int)strlen(temp); k++)
                    temp[k] = (char)(tolower(temp[k]));
                if (strncmp(tempName,temp,strlen(tempName)) == 0)
                    {
                    j = i;
                    nMatches++;
                    }
                }
            }

        if (nMatches == 0)
            {
            extern char *tokenP;
            YvyraPrint ("%s   Could not find parameter '%s': ignoring values\n", spacer, tempName);
            while (*tokenP && *tokenP++!=')') {}; 
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            return (*tokenP ? NO_ERROR:ERROR);
            }
        else if (nMatches > 1)
            {
            YvyraPrint ("%s   Several parameters matched the abbreviated name '%s'\n", spacer, tempName);
            return (ERROR);
            }
            
        param = &params[j];
        if (param->printParam == NO && !(param->paramType == P_TOPOLOGY && strcmp(modelParams[param->relParts[0]].topologyPr,"Fixed")!=0)
                                    && !(param->paramType == P_CPPEVENTS)
                                    && !(param->paramType == P_TK02BRANCHRATES)
                                    && !(param->paramType == P_WNBRANCHRATES)
                                    && !(param->paramType == P_IGRBRANCHRATES)
                                    && !(param->paramType == P_ILNBRANCHRATES)
                                    && !(param->paramType == P_MIXEDBRCHRATES)
                                    && !(param->paramType == P_POPSIZE && param->nValues > 1))
            {
            YvyraPrint ("%s   The parameter '%s' is fixed so the starting value cannot be set\n", spacer, param->name);
            return (ERROR);
            }
        if (param->paramType == P_BRLENS || param->paramType == P_TOPOLOGY || param->paramType == P_CPPEVENTS ||
            param->paramType == P_TK02BRANCHRATES || param->paramType == P_WNBRANCHRATES ||
            param->paramType == P_IGRBRANCHRATES || param->paramType == P_ILNBRANCHRATES || param->paramType == P_MIXEDBRCHRATES ||
            param->paramType == P_SPECIESTREE || (param->paramType == P_POPSIZE && param->nValues > 1))
            {
            /* all these parameters are set from a tree */
            expecting = Expecting (ALPHA);
            }
        else
            /* run of the mill character */
            {
            theValueMin = param->min;
            theValueMax = param->max;
            if ((param->paramType == P_PI && modelParams[param->relParts[0]].dataType != STANDARD) ||
                param->paramType == P_MIXTURE_RATES)
                {
                useSubvalues = YES;
                useIntValues = NO;
                numExpectedValues = param->nSubValues;
                }
            else if (param->nValues == 0 && param->nIntValues > 0)
                {
                useSubvalues = NO;
                useIntValues = YES;
                numExpectedValues = param->nIntValues;
                }
            else if (param->nValues > 0)
                {
                useSubvalues = NO;
                useIntValues = NO;
                numExpectedValues = param->nValues;
                }
            else
                {
                YvyraPrint ("%s   Not expecting any values for parameter '%s'\n", spacer, param->name);
                return (ERROR);
                }
            nValuesRead = 0;
            expecting = Expecting(LEFTPAR);
            }
        }
    else
        return (ERROR);

    SAFEFREE (temp);
    return (NO_ERROR);
}


int DoUnlink (void)
{
    int         i, j;
    
    YvyraPrint ("%s   Unlinking\n", spacer);
    
    /* update status of linkTable */
    for (j=0; j<NUM_LINKED; j++)
        {
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (tempLinkUnlink[j][i] == YES)
                {
                linkTable[j][i] = ++linkNum;
                }
            }
        }
    
#   if 0
    for (j=0; j<NUM_LINKED; j++)
        {
        YvyraPrint ("%s   ", spacer);
        for (i=0; i<numCurrentDivisions; i++)
            YvyraPrint ("%d", linkTable[j][i]);
        YvyraPrint ("\n");
        }
#   endif

    /* reinitialize the temporary table */
    for (j=0; j<NUM_LINKED; j++)
        for (i=0; i<numCurrentDivisions; i++)
            tempLinkUnlink[j][i] = NO;

    /* set up parameters and moves */
    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);

    return (NO_ERROR);
}


int DoShowMcmcTrees (void)
{
    int         run, chain, chainIndex, i;
    Tree        *t;

    for (run=0; run<chainParams.numRuns; run++)
        {
        for (chain=0; chain<chainParams.numChains; chain++)
            {
            chainIndex = run*chainParams.numChains + chain;
            for (i=0; i<numTrees; i++)
                {
                t = GetTreeFromIndex (i, chainIndex, 0);
                if (t->isRooted == YES)
                    YvyraPrint ("\n   Tree '%s' [rooted]:\n\n", t->name);
                else
                    YvyraPrint ("\n   Tree '%s' [unrooted]:\n\n", t->name);
                if (ShowTree (t) == ERROR)
                    return (ERROR);
                else
                    YvyraPrint ("\n");
                }
            }
        }

    return (NO_ERROR);
}


int DoShowModel (void)
{
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
        return (ERROR);
        }

    if (ShowModel() == ERROR)
        return (ERROR);

    return (NO_ERROR);
}


int DoShowMoves (void)
{
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before moves can be assigned\n", spacer);
        return (ERROR);
        }

    YvyraPrint ("%s   Moves that will be used by MCMC sampler (rel. proposal prob. > 0.0):\n\n", spacer);
    if (ShowMoves(YES) == ERROR)
        return (ERROR);

    if (showmovesParams.allavailable == YES)
        {
        YvyraPrint ("%s   Other available moves (rel. proposal prob. = 0.0):\n\n", spacer);
        if (ShowMoves(NO) == ERROR)
            return (ERROR);
        }
    else
        YvyraPrint ("%s   Use 'Showmoves allavailable=yes' to see a list of all available moves\n", spacer);

    return (NO_ERROR);
}


int DoShowmovesParm (char *parmName, char *tkn)
{
    char    tempStr[100];
    
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Allavailable **********************************************************/
        if (!strcmp(parmName, "Allavailable"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        showmovesParams.allavailable = YES;
                    else
                        showmovesParams.allavailable = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for allavailable\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else
            return (ERROR);

        }
        
    return (NO_ERROR);
}


int DoShowParams (void)
{
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before model parameters can be shown\n", spacer);
        return (ERROR);
        }

    if (ShowParameters(YES, YES, YES) == ERROR)
        return (ERROR);

    return (NO_ERROR);
}


/*------------------------------------------------------------------------
|
|   FillNormalParams: Allocate and fill in non-tree parameters
|
-------------------------------------------------------------------------*/
int FillNormalParams (RandLong *seed, int fromChain, int toChain)
{
    int         i, j, k, chn, tempInt, *intValue;
    YFlt      *bs, *value, *subValue, scaler;
    Tree        *tree;
    Param       *p;
    ModelInfo   *m;
    ModelParams *mp;

    /* fill in values for nontree params for state 0 of chains */
    for (chn=fromChain; chn<toChain; chn++)
        {
        for (k=0; k<numParams; k++)
            {
            p  = &params[k];
            mp = &modelParams[p->relParts[0]];
            m  = &modelSettings[p->relParts[0]];

            /* find model settings and nStates, pInvar, invar cond likes */

            value = GetParamVals (p, chn, 0);
            subValue = GetParamSubVals (p, chn, 0);
            intValue = GetParamIntVals (p, chn, 0);

            if (p->paramType == P_PI)
                {
                /* Fill in state frequencies ****************************************************************************/
                /* note that standard chars are mainly dealt with in ProcessStdChars in mcmc.c */
                if (p->paramId == SYMPI_UNI || p->paramId == SYMPI_UNI_MS)
                    value[0] = 1.0;
                else if (p->paramId == SYMPI_EXP || p->paramId == SYMPI_EXP_MS)
                    value[0] = 1.0;
                else if (p->paramId == SYMPI_FIX || p->paramId == SYMPI_FIX_MS)
                    value[0] = mp->symBetaFix;

                else if (p->paramId == PI_DIR)
                    {
                    if (mp->numDirParams != mp->nStates && mp->numDirParams != 0)
                        {
                        YvyraPrint ("%s   Mismatch between number of dirichlet parameters (%d) and the number of states (%d)\n", spacer, mp->numDirParams, m->numStates);
                        return ERROR;
                        }

                    /* if user has not set dirichlet parameters, go with default */
                    /* overall variance equals number of states */
                    if (mp->numDirParams == 0)
                        for (i=0; i<mp->nStates; i++)
                            value[i] = mp->stateFreqsDir[i] = 1.0;
                    else
                        for (i=0; i<m->numStates; i++)
                            value[i] = mp->stateFreqsDir[i];

                    /* now fill in subvalues */
                    for (i=0; i<m->numStates; i++)
                        subValue[i] = 1.0 / mp->nStates;
                    }

                else if (p->paramId == PI_USER)
                    {
                    for (i=0; i<m->numStates; i++)
                        subValue[i] =  mp->stateFreqsFix[i];
                    }
                    
                else if (p->paramId == PI_FIXED)
                    {
                    /* AA model Pi arrays removed: molecular models not supported in yvyra */
                    }

                else if (p->paramId == PI_EMPIRICAL)
                    {
                    if (GetEmpiricalFreqs (p->relParts, p->nRelParts) == ERROR)
                        return (ERROR);
                    for (i=0; i<mp->nStates; i++)
                        subValue[i] = empiricalFreqs[i];
                    }

                else if (p->paramId == PI_EQUAL)
                    {
                    for (i=0; i<mp->nStates; i++)
                        subValue[i] = 1.0 / mp->nStates;
                    }
                /* fill for the directional model of evolution for state frequencies */
                else if (p->paramId == DIRPI_DIRxDIR || p->paramId == DIRPI_DIRxFIXED || p->paramId == DIRPI_FIXEDxDIR || p->paramId == DIRPI_MIX)
                    {     
                    /* check if both Dirichlet priors are set correctly, if applicable */
                    if (p->paramId == DIRPI_DIRxDIR || p->paramId == DIRPI_DIRxFIXED || p->paramId == DIRPI_MIX)   
                        {
                        if (mp->numDirParams != mp->nStates && mp->numDirParams != 0)
                            { 
                            YvyraPrint ("%s   Mismatch between number of dirichlet parameters for the equillibrium frequencies (%d) and the number of states (%d)\n", spacer, mp->numDirParams, m->numStates);    
                            return ERROR;
                            }
                        }         
                    if (p->paramId == DIRPI_DIRxDIR || p->paramId == DIRPI_FIXEDxDIR || p->paramId == DIRPI_MIX)       
                        { 
                        if (mp->numDirParamsRoot != mp->nStates && mp->numDirParamsRoot != 0)
                            { 
                            YvyraPrint ("%s   Mismatch between number of dirichlet parameters for the root frequencies (%d) and the number of states (%d)\n", spacer, mp->numDirParamsRoot, m->numStates);        
                            return ERROR;
                            } 
                        }         
                    /* now fill in values and subvalues */
                    if (p->paramId == DIRPI_DIRxDIR || p->paramId == DIRPI_DIRxFIXED || p->paramId == DIRPI_MIX)   
                        { 
                        /* if user has not set dirichlet parameters, go with default */
                        /* overall variance equals number of states */
                        if (mp->numDirParams == 0)
                            for (i=0; i<mp->nStates; i++)
                                value[i] = mp->stateFreqsDir[i] = 1.0;
                        else
                            for (i=0; i<m->numStates; i++)
                                value[i] = mp->stateFreqsDir[i];

                        /* now fill in subvalues for equilibrium frequencies */
                        for (i=0; i<m->numStates; i++)
                                subValue[i] =  (1.0 / mp->nStates);
                        /* and the fixed root freqs, if any */
                        if (p->paramId == DIRPI_DIRxFIXED)
                            {
                            for (i=0; i<m->numStates; i++)
                                subValue[i+m->numStates] =  mp->rootFreqsFix[i];
                            }
                        }
                    if (p->paramId == DIRPI_DIRxDIR || p->paramId == DIRPI_FIXEDxDIR || p->paramId == DIRPI_MIX)
                        {
                        /* now fill in the root freqs, 2nd half of vector */
                        if (mp->numDirParamsRoot == 0)
                            for (i=0; i<m->numStates; i++)
                                value[i+m->numStates] = mp->rootFreqsDir[i] = 1.0;
                        else
                            for (i=0; i<m->numStates; i++)
                                value[i+m->numStates] = mp->rootFreqsDir[i];

                        for (i=0; i<m->numStates; i++)
                                subValue[i+m->numStates] =  (1.0 / mp->nStates);

                        if (p->paramId == DIRPI_FIXEDxDIR)
                            {
                            for (i=0; i<m->numStates; i++)
                                subValue[i] =  mp->stateFreqsFix[i];
                            }
                        }
                    }
                else if (p->paramId == DIRPI_FIXEDxFIXED)
                    {
                    for (i=0; i<m->numStates; i++)
                        {
                        subValue[i] =  mp->stateFreqsFix[i];
                        subValue[i+m->numStates] =  mp->rootFreqsFix[i];
                        }
                    }  
                }
            else if (p->paramType == P_MIXTURE_RATES)
                {
                /* Fill in rates of site rate mixture *******************************************************************/
                
                /* We use value array for dirichlet prior parameters. We use a flat prior, so this will be a series of 1.0 values */
                for (i=0; i<m->numRateCats; ++i)
                    value[i] = 1.0;
                
                /* Now fill in subvalues by setting them to be equal, Note that we use rates and not rate proportions. */
                for (i=0; i<m->numRateCats; ++i)
                    subValue[i] = 1.0;
                }
            else if (p->paramType == P_SHAPE)
                {
                /* Fill in shape values ********************************************************************************/
                /* first get hyperprior */
                if (p->paramId == SHAPE_UNI)
                    {
                    value[0] = 1.0;
                    if (value[0] < mp->shapeUni[0] || value[0] > mp->shapeUni[1])
                        value[0] = mp->shapeUni[0] + (mp->shapeUni[1] - mp->shapeUni[0]) *  0.5;
                    }
                else if (p->paramId == SHAPE_EXP)
                    value[0] = 1.0;  // was 100.0
                else if (p->paramId == SHAPE_FIX)
                    value[0] = mp->shapeFix;
                /* now fill in rates */
                if (!strcmp(mp->ratesModel, "LNorm"))
                    {
                    if (DiscreteLogNormal (subValue, value[0], mp->numGammaCats, 1) == ERROR)
                        return (ERROR);
                    }
                else  /* gamma rate */
                    {
                    if (DiscreteGamma (subValue, value[0], value[0], mp->numGammaCats, 0) == ERROR)
                        return (ERROR);
                    }
                }
            else if (p->paramType == P_PINVAR)
                {
                /* Fill in pInvar ***************************************************************************************/
                if (p->paramId == PINVAR_UNI)
                    value[0] = 0.0;

                else if (p->paramId == PINVAR_FIX)
                    value[0] =  mp->pInvarFix;
                }
            else if (p->paramType == P_CORREL)
                {
                /* Fill in correlation parameter of adgamma model *******************************************************/
                if (p->paramId == CORREL_UNI)
                    value[0] = 0.0;

                else if (p->paramId == CORREL_FIX)
                    value[0] =  mp->adgCorrFix;
                
                /* Fill in correlation matrices */
                AutodGamma (subValue, value[0], mp->numGammaCats);
                }
            else if (p->paramType == P_SWITCH)
                {
                /* switchRates for covarion model removed: molecular models not supported *****************************/
                for (j=0; j<2; j++)
                    {
                    if (p->paramId == SWITCH_UNI)
                        value[j] = RandomNumber(seed) * (mp->covswitchUni[1] - mp->covswitchUni[0]) + mp->covswitchUni[0];

                    else if (p->paramId == SWITCH_EXP)
                        value[j] = -(1.0/mp->covswitchExp) * log(RandomNumber(seed));

                    else if (p->paramId == SWITCH_FIX)
                        value[j] = mp->covswitchFix[j];
                    }
                }
            else if (p->paramType == P_RATEMULT)
                {
                /* Fill in division rates *****************************************************************************/
                for (j=0; j<p->nValues; j++)
                    {
                    value[j] = 1.0;
                    /* fill in more info about the divisions if this is a true rate multiplier
                       and not a base rate */
                    if (p->nSubValues > 0)
                        {
                        /* num uncompressed chars */
                        subValue[j] = (modelSettings[p->relParts[j]].numUncompressedChars);
                        /* Dirichlet parameters */
                        subValue[p->nValues + j] = modelParams[p->relParts[j]].ratePrDir;
                        }
                    }
                }
            else if (p->paramType == P_SPECRATE)
                {
                /* Fill in speciation rates *****************************************************************************/
                for (j=0; j<p->nValues; j++)
                    {
                    if (p->paramId == SPECRATE_FIX)
                        value[j] = mp->speciationFix;
                    else
                        value[j] = 0.3;
                    }
                }
            else if (p->paramType == P_EXTRATE)
                {
                /* Fill in extinction rates *****************************************************************************/
                for (j=0; j<p->nValues; j++)
                    {
                    if (p->paramId == EXTRATE_FIX)
                        value[j] = mp->extinctionFix;
                    else
                        value[j] = 0.2;
                    }
                }
            else if (p->paramType == P_FOSLRATE)
                {
                /* Fill in fossilization rates */
                for (j=0; j<p->nValues; j++)
                    {
                    if (p->paramId == FOSLRATE_FIX)
                        value[j] = mp->fossilizationFix;
                    else
                        value[j] = 0.1;
                    }
                }
            else if (p->paramType == P_GROWTH)
                {
                /* Fill in growth rate **********************************************************************************/
                if (p->paramId == GROWTH_FIX)
                    value[0] = mp->growthFix;
                else
                    value[0] = 1.0;
                }
            else if (p->paramType == P_POPSIZE)
                {
                /* Fill in population size ******************************************************************************/
                for (j=0; j<p->nValues; j++)
                    {
                    if (p->paramId == POPSIZE_UNI)
                        value[j] = RandomNumber(seed) * (mp->popSizeUni[1] - mp->popSizeUni[0]) + mp->popSizeUni[0];
                    else if (p->paramId == POPSIZE_LOGNORMAL)
                        value[j] = exp(mp->popSizeLognormal[0]);
                    else if (p->paramId == POPSIZE_NORMAL)
                        value[j] = mp->popSizeNormal[0];
                    else if (p->paramId == POPSIZE_GAMMA)
                        value[j] = mp->popSizeGamma[0] / mp->popSizeGamma[1];
                    else if (p->paramId == POPSIZE_FIX)
                        value[j] = mp->popSizeFix;
                    }
                }
            else if (p->paramType == P_BMCORR)
                {
                /* Fill in correlation parameter for brownian motion ****************************************************/
                if (p->paramId == BMCORR_FIX)
                    value[0] = mp->brownCorrFix;
                else
                    value[0] = 0.0;
                }
            else if (p->paramType == P_BMSIGMA)
                {
                /* Fill in sigma parameter for brownian motion **********************************************************/
                if (p->paramId == BMSIGMA_FIX)
                    value[0] = mp->brownScaleFix;
                else if (p->paramId == BMSIGMA_UNI)
                    value[0] = RandomNumber(seed) * (mp->brownScaleUni[1] - mp->brownScaleUni[0]) + mp->brownScaleUni[0];
                else if (p->paramId == BMSIGMA_GAMMA)
                    value[0] = mp->brownScaleGamma[0] / mp->brownScaleGamma[1];
                }
            else if (p->paramType == P_CPPRATE)
                {
                /* Fill in lambda (cpp rate) ********************************************************************************************/
                if (p->paramId == CPPRATE_EXP)
                    value[0] = -(1.0/mp->cppRateExp) * log(RandomNumber(seed));
                else if (p->paramId == CPPRATE_FIX)
                    value[0] = mp->cppRateFix;
                }
            else if (p->paramType == P_CPPMULTDEV)
                {
                /* Fill in log standard deviation (for relaxed clock rate multiplier) ***********************************************************/
                if (p->paramId == CPPMULTDEV_FIX)
                    value[0] = mp->cppMultDevFix;
                }
            else if (p->paramType == P_CPPEVENTS)
                {
                /* We fill in these when we fill in tree params **************************************************************************/
                }
            else if (p->paramType == P_TK02VAR)
                {
                /* Fill in variance of relaxed clock lognormal **************************************************************************/
                if (p->paramId == TK02VAR_UNI)
                    value[0] = RandomNumber(seed) * (mp->tk02varUni[1] - mp->tk02varUni[0]) + mp->tk02varUni[0];
                else if (p->paramId == TK02VAR_EXP)
                    value[0] = 2.0/(mp->tk02varExp);
                else if (p->paramId == TK02VAR_FIX)
                    value[0] = mp->tk02varFix;
                }
            else if (p->paramType == P_WNVAR)
                {
                /* Fill in variance of relaxed clock white noise **************************************************************************/
                if (p->paramId == WNVAR_UNI)
                    value[0] = RandomNumber(seed) * (mp->wnvarUni[1] - mp->wnvarUni[0]) + mp->wnvarUni[0];
                else if (p->paramId == WNVAR_EXP)
                    value[0] = 1.0/(mp->wnvarExp);
                else if (p->paramId == WNVAR_FIX)
                    value[0] = mp->wnvarFix;
                }
            else if (p->paramType == P_ILNVAR)
                {
                /* Fill in variance of relaxed clock lognormal **************************************************************************/
                if (p->paramId == ILNVAR_UNI)
                    value[0] = RandomNumber(seed) * (mp->ilnvarUni[1] - mp->ilnvarUni[0]) + mp->ilnvarUni[0];
                else if (p->paramId == ILNVAR_EXP)
                    value[0] = 1.0/(mp->ilnvarExp);
                else if (p->paramId == ILNVAR_FIX)
                    value[0] = mp->ilnvarFix;
                }
             else if (p->paramType == P_IGRVAR)
                {
                /* Fill in variance of relaxed clock gamma      **************************************************************************/
                if (p->paramId == IGRVAR_UNI)
                    value[0] = RandomNumber(seed) * (mp->igrvarUni[1] - mp->igrvarUni[0]) + mp->igrvarUni[0];
                else if (p->paramId == IGRVAR_EXP)
                    value[0] = 1.0/(mp->igrvarExp);
                else if (p->paramId == IGRVAR_FIX)
                    value[0] = mp->igrvarFix;
                }
            else if (p->paramType == P_MIXEDVAR)
                {
                /* Fill in variance of mixed relaxed clock      **************************************************************************/
                if (p->paramId == MIXEDVAR_UNI)
                    value[0] = RandomNumber(seed) * (mp->mixedvarUni[1] - mp->mixedvarUni[0]) + mp->mixedvarUni[0];
                else if (p->paramId == MIXEDVAR_EXP)
                    value[0] = 1.0/(mp->mixedvarExp);
                else if (p->paramId == MIXEDVAR_FIX)
                    value[0] = mp->mixedvarFix;
                }
            else if (p->paramType == P_TK02BRANCHRATES || p->paramType == P_WNBRANCHRATES ||
                     p->paramType == P_ILNBRANCHRATES || p->paramType == P_IGRBRANCHRATES)
                {
                /* We fill in these when we fill in tree params **************************************************************************/
                }
            else if (p->paramType == P_MIXEDBRCHRATES)
                {
                /* initialize the mixed relaxed clock model to TK02 or ILN */
                intValue[0] = (RandomNumber(seed) <0.5) ? RCL_IGR : RCL_ILN;
                /* We fill in the rest when we fill in tree params **************************************************************************/
                }
            else if (p->paramType == P_CLOCKRATE)
                {
                /* Fill in base rate of molecular clock **************************************************************************/
                if (p->paramId == CLOCKRATE_FIX)
                    value[0] = mp->clockRateFix;
                else if (p->paramId == CLOCKRATE_NORMAL)
                    value[0] = mp->clockRateNormal[0];
                else if (p->paramId == CLOCKRATE_LOGNORMAL)
                    value[0] = exp(mp->clockRateLognormal[0]);
                else if (p->paramId == CLOCKRATE_EXP)
                    value[0] = 1.0/(mp->clockRateExp);
                else if (p->paramId == CLOCKRATE_GAMMA)
                    value[0] = mp->clockRateGamma[0]/mp->clockRateGamma[1];
                }
            }   /* next param */
        }   /* next chain */

    return NO_ERROR;
}

    
int FillRelPartsString (Param *p, char **relPartString)
{
    int         i, n, filledString;
    char        *tempStr;
    int             tempStrSize=50;

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        YvyraPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    if (numCurrentDivisions == 1)
        {
        filledString = NO;
        SafeStrcpy (relPartString, "");
        }
    else
        {
        filledString = YES;
        if (p->nRelParts == numCurrentDivisions)
            {
            SafeStrcpy (relPartString, "{all}");
            }
        else
            {
            SafeStrcpy (relPartString, "{");
            for (i=n=0; i<p->nRelParts; i++)
                {
                n++;
                SafeSprintf (&tempStr, &tempStrSize, "%d", p->relParts[i] + 1);
                SafeStrcat (relPartString, tempStr);
                if (n < p->nRelParts)
                    SafeStrcat (relPartString, ",");
                }
            SafeStrcat (relPartString, "}");
            }
        }
    free (tempStr);
    return (filledString);
}


/*--------------------------------------------------------------
 |
 |  FillStdStateFreqs: fills stationary frequencies for standard data divisions of chains in range [chfrom, chto)
 |
 ---------------------------------------------------------------*/
void FillStdStateFreqs (int chfrom, int chto, RandLong *seed)
{
    int     chn, n, i, j, k, b, c, nb, index;
    YFlt  *subValue, sum, symDir[MAX_STD_STATES];
    Param   *p;
    
    for (chn=chfrom; chn<chto; chn++)
        {
        for (k=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->paramType != P_PI || modelParams[p->relParts[0]].dataType != STANDARD)
                continue;
            subValue = GetParamStdStateFreqs (p, chn, 0);
            if (p->paramId == SYMPI_EQUAL)
                {
                for (n = index = 0; n < MAX_STD_STATES-1; n++)
                    {
                    for (i=0; i<p->nRelParts; i++)
                        if (modelSettings[p->relParts[i]].isTiNeeded[n] == YES)
                            break;
                    if (i < p->nRelParts)
                        {
                        for (j=0; j<(n+2); j++)
                            {
                            subValue[index++] = 1.0 / (n+2);
                            }
                        }
                    }
                for (n = MAX_STD_STATES-1; n < 2*MAX_STD_STATES-3; n++)
                    {
                    for (i=0; i<p->nRelParts; i++)
                        if (modelSettings[p->relParts[i]].isTiNeeded[n] == YES)
                            break;
                    if (i < p->nRelParts)
                        {
                        for (j = 0; j < n-MAX_STD_STATES+4; j++)
                            {
                            subValue[index++] = 1.0 / (n-MAX_STD_STATES+4);
                            }
                        }
                    }
                /* fill equal state frequencies for USERTYPE characters */
                for (i=0; i<p->nRelParts; i++)
                    {
                    ModelInfo *mi = &modelSettings[p->relParts[i]];
                    for (c=0; c<mi->numChars; c++)
                        {
                        if (mi->cType[c] == USERTYPE)
                            {
                            n = mi->nStates[c];
                            for (j=0; j<n; j++)
                                subValue[mi->bsIndex[c] + j] = 1.0 / n;
                            }
                        }
                    }
                }
            
            /* Deal with transition asymmetry for standard characters */
            /* First, fill in stationary frequencies for beta categories if needed; */
            /* discard category frequencies (assume equal) */
            if (p->paramId == SYMPI_FIX || p->paramId == SYMPI_UNI || p->paramId == SYMPI_EXP ||
                p->paramId == SYMPI_FIX_MS || p->paramId == SYMPI_UNI_MS || p->paramId == SYMPI_EXP_MS)
                {
                if (p->hasBinaryStd == YES)
                    {
                    nb=modelParams[p->relParts[0]].numBetaCats;
                    BetaBreaks (p->values[0], p->values[0], subValue, nb);
                    b = 2*nb;
                    for (i=b-2; i>0; i-=2)
                        {
                        subValue[i] = subValue[i/2];
                        }
                    for (i=1; i<b; i+=2)
                        {
                        subValue[i] = (1.0 - subValue[i-1]);
                        }
                    subValue += (2 * nb);
                    }
                
                /* Then fill in state frequencies for multistate chars, one set for each */
                for (i=0; i<MAX_STD_STATES; i++)
                    symDir[i] = p->values[0];
                
                for (c=0; c<p->nSympi; c++)
                    {
                    /* now fill in subvalues */
                    DirichletRandomVariable (symDir, subValue, p->sympinStates[c], seed);
                    sum = 0.0;
                    for (i=0; i<p->sympinStates[c]; i++)
                        {
                        if (subValue[i] < 0.0001)
                            subValue[i] = 0.0001;
                        sum += subValue[i];
                        }
                    for (i=0; i<modelParams[p->relParts[0]].nStates; i++)
                        subValue[i] /= sum;
                    subValue += p->sympinStates[c];
                    }
                }
            }   /* next parameter */
        }   /* next chain */
}


/* FillTopologySubParams: Fill subparams (brlens) for a topology */
int FillTopologySubParams (Param *param, int chn, int state, RandLong *seed)
{
    int         i,returnVal;
    Tree        *tree, *tree1;
    Param       *q;
    YFlt      clockRate;
    PolyTree    *sourceTree;
    YFlt      minRate=0.0, maxRate=0.0;

    tree = GetTree (param, chn, state);
    
    for (i=1; i<param->nSubParams; i++)
        {
        q = param->subParams[i];
        tree1 = GetTree (q, chn, state);
        if (CopyToTreeFromTree(tree1, tree) == ERROR)
            return (ERROR);
        }
    for (i=0; i<param->nSubParams; i++)
        {
        q = param->subParams[i];
        tree = GetTree (q, chn, state);
        if (q->paramId == BRLENS_FIXED || q->paramId == BRLENS_CLOCK_FIXED)
            {
            if (param->paramId == TOPOLOGY_NCL_FIXED ||
                param->paramId == TOPOLOGY_NCL_FIXED_HOMO ||
                param->paramId == TOPOLOGY_NCL_FIXED_HETERO ||
                param->paramId == TOPOLOGY_CL_FIXED  ||
                param->paramId == TOPOLOGY_RCL_FIXED ||
                param->paramId == TOPOLOGY_CCL_FIXED ||
                param->paramId == TOPOLOGY_RCCL_FIXED||
                param->paramId == TOPOLOGY_FIXED)
                {
                sourceTree = AllocatePolyTree(numTaxa);
                CopyToPolyTreeFromPolyTree (sourceTree, userTree[modelParams[q->relParts[0]].brlensFix]);
                PrunePolyTree (sourceTree);
                ResetTipIndices (sourceTree);
                ResetIntNodeIndices (sourceTree);
                if (tree->isRooted != sourceTree->isRooted)
                    {
                    YvyraPrint ("%s   Cannot set fixed branch lengths because of mismatch in rootedness", spacer);
                    FreePolyTree (sourceTree);
                    return (ERROR);
                    }
                if (CopyToTreeFromPolyTree(tree,sourceTree) == ERROR)
                    {
                    YvyraPrint ("%s   Problem setting fixed branch lengths", spacer);
                    FreePolyTree (sourceTree);
                    return (ERROR);
                    }
                FreePolyTree (sourceTree);
                if (tree->isClock == YES && IsClockSatisfied(tree, 0.001) == NO)
                    {
                    YvyraPrint ("%s   Fixed branch lengths do not satisfy clock", spacer);
                    return (ERROR);
                    }
                if (tree->isCalibrated == YES && IsCalibratedClockSatisfied (tree,&minRate,&maxRate, 0.001) == NO)
                    {
                    YvyraPrint ("%s   Fixed branch lengths do not satisfy calibrations", spacer);
                    return (ERROR);
                    }
                if (tree->isCalibrated == YES && !strcmp(modelParams[tree->relParts[0]].clockRatePr, "Fixed"))
                    {
                    clockRate = modelParams[tree->relParts[0]].clockRateFix;
                    if ((clockRate < minRate && AreDoublesEqual (clockRate, minRate , 0.001) == NO) ||
                        (clockRate > maxRate && AreDoublesEqual (clockRate, maxRate , 0.001) == NO))
                        {
                        YvyraPrint ("%s   Fixed branch lengths do not satisfy fixed clockrate", spacer);
                        return (ERROR);
                        }
                    }

                tree->fromUserTree=YES;
                returnVal = NO_ERROR;
                }
            else
                {
                YvyraPrint ("%s   Fixed branch lengths can only be used for a fixed topology\n", spacer);
                return (ERROR);
                }
            }
        else if (tree->isCalibrated == YES ||
                 (tree->isClock == YES && (!strcmp(modelParams[tree->relParts[0]].clockPr,"Uniform") ||
                                           !strcmp(modelParams[tree->relParts[0]].clockPr,"Birthdeath") ||
                                           !strcmp(modelParams[tree->relParts[0]].clockPr,"Fossilization"))))
            {
            assert (tree->isClock == YES);
            clockRate = *GetParamVals(modelSettings[tree->relParts[0]].clockRate, chn, state);
            returnVal = InitCalibratedBrlens (tree, clockRate, seed);
            if (IsClockSatisfied (tree, 0.001) == NO)
                {
                YvyraPrint ("%s   Branch lengths of the tree does not satisfy clock\n",  spacer);
                return (ERROR);
                }
            tree->fromUserTree=NO;
            }
        else if (tree->isClock == YES)
            returnVal = InitClockBrlens (tree);
        else
            returnVal = InitBrlens (tree, 0.02);

        if (returnVal == ERROR)
            return (ERROR);

        if (FillBrlensSubParams (q, chn, state) == ERROR)
            return (ERROR);
        }

    return (NO_ERROR);
}


/* FillBrlensSubParams: Fill any relaxed clock subparams of a brlens param */
int FillBrlensSubParams (Param *param, int chn, int state)
{
    int         i, j, *nEvents;
    YFlt      *brlen, *branchRate, **position, **rateMult;
    Tree        *tree;
    TreeNode    *p;
    Param       *q;

    tree = GetTree (param, chn, state);
    
    for (i=0; i<param->nSubParams; i++)
        {
        q = param->subParams[i];
        if (q->paramType == P_CPPEVENTS)
            {
            nEvents = q->nEvents[2*chn+state];
            position = q->position[2*chn+state];
            rateMult = q->rateMult[2*chn+state];
            brlen = GetParamSubVals (q, chn, state);
            for (j=0; j<tree->nNodes-1; j++)
                {
                p = tree->allDownPass[j];
                if (nEvents[p->index] != 0)
                    {
                    free (position[p->index]);
                    position[p->index] = NULL;
                    free (rateMult[p->index]);
                    rateMult[p->index] = NULL;
                    }
                nEvents[p->index] = 0;
                // assert (j==tree->nNodes-2 || fabs(p->length - (p->anc->nodeDepth - p->nodeDepth)) < 0.000001);
                brlen[p->index] = p->length;
                }
            }
        else if (q->paramType == P_TK02BRANCHRATES || q->paramType == P_WNBRANCHRATES ||
                 q->paramType == P_IGRBRANCHRATES || q->paramType == P_ILNBRANCHRATES || q->paramType == P_MIXEDBRCHRATES)
            {
            branchRate = GetParamVals (q, chn, state);
            brlen = GetParamSubVals (q, chn, state);
            for (j=0; j<tree->nNodes-1; j++)
                {
                p = tree->allDownPass[j];
                branchRate[p->index] = 1.0;
                brlen[p->index] = p->length;
                }
            }
        }

    return (NO_ERROR);
}


/* Note: In PruneConstraintPartitions() we can not rely on specific rooting of a tree since different partitions
   may theoretically have different clock models, while constraints apply to all partitions/trees */
int PruneConstraintPartitions (void)
{
    int             i, j, constraintId, nLongsNeeded;
    
    nLongsNeeded = (numLocalTaxa - 1) / nBitsInALong + 1;

    for (constraintId=0; constraintId<numDefinedConstraints; constraintId++)
        {
        definedConstraintPruned[constraintId] = (BitsLong *) SafeRealloc ((void *)definedConstraintPruned[constraintId], nLongsNeeded*sizeof(BitsLong));
        if (!definedConstraintPruned[constraintId])
            {
            YvyraPrint ("%s   Problems allocating constraintPartition in PruneConstraintPartitions", spacer);
            return (ERROR);
            }

        /* initialize bits in partition to add; get rid of deleted taxa in the process */
        ClearBits(definedConstraintPruned[constraintId], nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            if (IsBitSet(i, definedConstraint[constraintId]) == YES)
                SetBit(j, definedConstraintPruned[constraintId]);
            j++;
            }
        assert (j == numLocalTaxa);

        if (definedConstraintsType[constraintId] == PARTIAL)
            {
            definedConstraintTwoPruned[constraintId] = (BitsLong *) SafeRealloc ((void *)definedConstraintTwoPruned[constraintId], nLongsNeeded*sizeof(BitsLong));
            if (!definedConstraintTwoPruned[constraintId])
                {
                YvyraPrint ("%s   Problems allocating constraintPartition in PruneConstraintPartitions", spacer);
                return (ERROR);
                }

            /* initialize bits in partition to add; get rid of deleted taxa in the process */
            ClearBits(definedConstraintTwoPruned[constraintId], nLongsNeeded);
            for (i=j=0; i<numTaxa; i++)
                {
                if (taxaInfo[i].isDeleted == YES)
                    continue;
                if (IsBitSet(i, definedConstraintTwo[constraintId]) == YES)
                    SetBit(j, definedConstraintTwoPruned[constraintId]);
                j++;
                }
            assert (j == numLocalTaxa);
            }
        else if (definedConstraintsType[constraintId] == NEGATIVE || (definedConstraintsType[constraintId] == HARD))
            {
            /* Here we create definedConstraintTwoPruned[constraintId] which is complemente of definedConstraintPruned[constraintId] */
            definedConstraintTwoPruned[constraintId] = (BitsLong *) SafeRealloc ((void *)definedConstraintTwoPruned[constraintId], nLongsNeeded*sizeof(BitsLong));
            if (!definedConstraintTwoPruned[constraintId])
                {
                YvyraPrint ("%s   Problems allocating constraintPartition in PruneConstraintPartitions", spacer);
                return (ERROR);
                }

            /* initialize bits in partition to add; get rid of deleted taxa in the process */
            ClearBits(definedConstraintTwoPruned[constraintId], nLongsNeeded);
            for (i=j=0; i<numTaxa; i++)
                {
                if (taxaInfo[i].isDeleted == YES)
                    continue;
                if (IsBitSet(i, definedConstraint[constraintId]) == NO)
                    SetBit(j, definedConstraintTwoPruned[constraintId]);
                j++;
                }
            assert (j == numLocalTaxa);         
            }
    }
    
    return NO_ERROR;
}


int DoesTreeSatisfyConstraints(Tree *t)
{
    int         i, k, numTaxa, nLongsNeeded;
    TreeNode    *p;
    int         CheckFirst, CheckSecond; /*Flag indicating whether corresponding set(first/second) of partial constraint has to be checked*/
#   if defined (DEBUG_CONSTRAINTS)
    int         locks_count=0;
#   endif

    if (t->checkConstraints == NO)
        return YES;
    /* get some handy numbers */
    numTaxa = t->nNodes - t->nIntNodes - (t->isRooted == YES ? 1 : 0);
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;

    if (t->bitsets == NULL)
        {
        AllocateTreePartitions(t);
        }
    else
        {
        ResetTreePartitions(t);  /*Inefficient function, rewrite faster version*/
        }
#   if defined (DEBUG_CONSTRAINTS)
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (p->isLocked == YES)
            {
            if (IsUnionEqThird (definedConstraintPruned[p->lockID], definedConstraintPruned[p->lockID], p->partition, nLongsNeeded) == NO && IsUnionEqThird (definedConstraintTwoPruned[p->lockID], definedConstraintTwoPruned[p->lockID], p->partition, nLongsNeeded) == NO)
                {
                printf ("DEBUG ERROR: Locked node does not represent right partition. \n");
                return ABORT;
                }
            else
                {
                locks_count++;
                }
            }
        }

    if (locks_count != t->nLocks)
        {
        printf ("DEBUG ERROR: locks_count:%d should be locks_count:%d\n", locks_count, t->nLocks);
        return ABORT;
        }
#   endif

    for (k=0; k<numDefinedConstraints; k++)
        {
#   if defined (DEBUG_CONSTRAINTS)
        if (t->constraints[k] == YES && definedConstraintsType[k] == HARD)
            {
            if (t->isRooted == YES)
                {
                CheckFirst = YES;
                CheckSecond = NO; 
                }
            else
                {
                /*exactly one of next two will be YES*/
                CheckFirst = IsBitSet(localOutGroup, definedConstraintPruned[k])==YES ? NO : YES;
                CheckSecond = IsBitSet(localOutGroup, definedConstraintTwoPruned[k])==YES ? NO : YES;
                assert ((CheckFirst^CheckSecond)==1);
                }

            for (i=0; i<t->nIntNodes; i++)
                {
                p = t->intDownPass[i];
                if (p->anc != NULL)
                    {
                    if (CheckFirst==YES &&  IsPartNested(definedConstraintPruned[k], p->partition, nLongsNeeded) && IsPartNested(p->partition,definedConstraintPruned[k], nLongsNeeded))
                        break;
                    if (CheckSecond==YES &&  IsPartNested(definedConstraintTwoPruned[k], p->partition, nLongsNeeded) && IsPartNested(p->partition, definedConstraintTwoPruned[k], nLongsNeeded))
                        break;
                    }
                }

            if (i==t->nIntNodes)
                {
                printf ("DEBUG ERROR: Hard constraint is not satisfied. \n");
                return ABORT;
                //assert (0);
                }
            }
#   endif

        if (t->constraints[k] == NO || definedConstraintsType[k] == HARD)
            continue;

        if (definedConstraintsType[k] == PARTIAL)
            {
            /* alternative way
            if (t->isRooted == NO && !IsBitSet(localOutGroup, definedConstraintPruned[k]))
                {
                m = FirstTaxonInPartition (constraintPartition, nLongsNeeded);
                for (i=0; t->nodes[i].index != m; i++)
                    ;
                p = &t->nodes[i];

                p=p->anc;
                while (!IsPartNested(definedConstraintPruned[k], p->partition, nLongsNeeded))
                    p=p->anc;

                if (IsSectionEmpty(definedConstraintTwoPruned[k], p->partition, nLongsNeeded))
                    continue;
                }
            */
            if (t->isRooted == YES)
                {
                CheckFirst = YES;
                CheckSecond = NO; /* In rooted case even if we have a node with partition fully containing second set and not containing the first set it would not satisfy the constraint */
                }
            else
                {
                if (NumBits(definedConstraintPruned[k], nLongsNeeded) == 1 || NumBits(definedConstraintTwoPruned[k], nLongsNeeded) == 1)
                    continue;
                /*one or two of the next two statements will be YES*/
                CheckFirst = IsBitSet(localOutGroup, definedConstraintPruned[k])==YES ? NO : YES;
                CheckSecond = IsBitSet(localOutGroup, definedConstraintTwoPruned[k])==YES ? NO : YES;
                assert ((CheckFirst|CheckSecond)==1);
                }
            for (i=0; i<t->nIntNodes; i++)
                {
                p = t->intDownPass[i];
                if (p->anc != NULL)
                    { 
                    if (CheckFirst== YES && IsPartNested(definedConstraintPruned[k], p->partition, nLongsNeeded) && IsSectionEmpty(definedConstraintTwoPruned[k], p->partition, nLongsNeeded))
                        break;
                    if (CheckSecond==YES && IsPartNested(definedConstraintTwoPruned[k], p->partition, nLongsNeeded) && IsSectionEmpty(definedConstraintPruned[k], p->partition, nLongsNeeded))
                        break;
                    }
                }
            if (i==t->nIntNodes)
                return NO;
            }
        else
            {
            assert (definedConstraintsType[k] == NEGATIVE);
            if (t->isRooted == YES)
                {
                CheckFirst = YES;
                CheckSecond = NO; 
                }
            else
                {
                /*exactly one of next two will be YES*/
                CheckFirst = IsBitSet(localOutGroup, definedConstraintPruned[k])==YES ? NO : YES;
                CheckSecond = IsBitSet(localOutGroup, definedConstraintTwoPruned[k])==YES ? NO : YES;
                assert ((CheckFirst^CheckSecond)==1);
                }

            for (i=0; i<t->nIntNodes; i++)
                {
                p = t->intDownPass[i];
                if (p->anc != NULL)
                    {
                    if (CheckFirst==YES && AreBitfieldsEqual(definedConstraintPruned[k], p->partition, nLongsNeeded))
                        break;
                    if (CheckSecond==YES && AreBitfieldsEqual(definedConstraintTwoPruned[k], p->partition, nLongsNeeded))
                        break;
                    }
                }
            if (i!=t->nIntNodes)
                return NO;
            }
        }
    
    return YES;
}


/*------------------------------------------------------------------
|
|   FillTreeParams: Fill in trees and branch lengths
|                   
|   Note: Should be run after FillNormalParams because
|   clockrate needs to be set if calibrated tree needs
|   to be filled.
|
------------------------------------------------------------------*/
int FillTreeParams (RandLong *seed, int fromChain, int toChain)
{
    int         i, k, chn, nTaxa, tmp;
    Param       *p, *q;
    Tree        *tree;
    PolyTree    *constraintTree;
    PolyTree    *constraintTreeRef;

    if (PruneConstraintPartitions() == ERROR)
        return ERROR;

    /* Build starting trees for state 0 */
    for (chn=fromChain; chn<toChain; chn++)
        {
        for (k=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->paramType == P_TOPOLOGY)
                {
                q = p->subParams[0];
                tree = GetTree (q, chn, 0);
                if (tree->isRooted == YES)
                    nTaxa = tree->nNodes - tree->nIntNodes - 1;
                else
                    nTaxa = tree->nNodes - tree->nIntNodes;
                
                /* fixed topology */
                if (p->paramId == TOPOLOGY_NCL_FIXED ||
                    p->paramId == TOPOLOGY_RNCL_FIXED ||
                    p->paramId == TOPOLOGY_NCL_FIXED_HOMO ||
                    p->paramId == TOPOLOGY_NCL_FIXED_HETERO ||
                    p->paramId == TOPOLOGY_CL_FIXED  ||
                    p->paramId == TOPOLOGY_RCL_FIXED ||
                    p->paramId == TOPOLOGY_CCL_FIXED ||
                    p->paramId == TOPOLOGY_RCCL_FIXED||
                    p->paramId == TOPOLOGY_FIXED     ||
                    p->paramId == TOPOLOGY_PARSIMONY_FIXED)
                    {
                    constraintTree = AllocatePolyTree (numTaxa);
                    CopyToPolyTreeFromPolyTree (constraintTree, userTree[modelParams[p->relParts[0]].topologyFix]);
                    PrunePolyTree (constraintTree);
                    ResetTipIndices(constraintTree);
                    ResetIntNodeIndices(constraintTree);
                    RandResolve (NULL, constraintTree, seed, constraintTree->isRooted);
                    if (tree->nIntNodes != constraintTree->nIntNodes)
                        {
                        if (tree->isRooted != constraintTree->isRooted)
                            {
                            YvyraPrint ("%s   Could not fix topology because user tree '%s' differs in rootedness with the model tree.\n", spacer,
                                          userTree[modelParams[p->relParts[0]].topologyFix]->name);
                            YvyraPrint ("%s   The user tree %s is%srooted, while expected model tree is%srooted.\n", spacer,
                                          userTree[modelParams[p->relParts[0]].topologyFix]->name, (constraintTree->isRooted?" ":" not "), (tree->isRooted?" ":" not "));
                            YvyraPrint ("%s   Check brlenspr is set correctly before fixing topology.\n", spacer);
                            }
                        else
                            YvyraPrint ("%s   Could not fix topology because user tree '%s' is not fully resolved.\n",
                                                spacer, userTree[modelParams[p->relParts[0]].topologyFix]->name);
                        FreePolyTree (constraintTree);
                        return (ERROR);
                        }
                    if (CopyToTreeFromPolyTree(tree, constraintTree) == ERROR)
                        {
                        YvyraPrint ("%s   Could not fix topology according to user tree '%s'\n", spacer, userTree[modelParams[p->relParts[0]].topologyFix]->name);
                        FreePolyTree (constraintTree);
                        return (ERROR);
                        }
                    FreePolyTree (constraintTree);
                    }
                /* constrained topology */
                else if (tree->nConstraints > 0)
                    {
                    constraintTreeRef = AllocatePolyTree (nTaxa);
                    if (!constraintTreeRef)
                        return (ERROR);
                    if (BuildConstraintTree (tree, constraintTreeRef, localTaxonNames) == ERROR)
                        {
                        FreePolyTree (constraintTreeRef);
                        return (ERROR);
                        }
                    if (AllocatePolyTreePartitions (constraintTreeRef) == ERROR)
                        return (ERROR);

                    constraintTree = AllocatePolyTree (nTaxa);
                    if (!constraintTree)
                        return (ERROR);
                    if (AllocatePolyTreePartitions (constraintTree) == ERROR)
                        return (ERROR);

                    for (i=0;i<100;i++)
                        {
                        CopyToPolyTreeFromPolyTree(constraintTree,constraintTreeRef);
                        tmp = RandResolve (tree, constraintTree, &globalSeed, tree->isRooted);
                        if (tmp != NO_ERROR)
                            {
                            if (tmp  == ERROR)
                                {
                                FreePolyTree (constraintTree);
                                return (ERROR);
                                }
                            else
                                {   
                                assert (tmp  == ABORT);
                                continue;
                                }
                            }
                   
                        CopyToTreeFromPolyTree(tree, constraintTree);
                        if (DoesTreeSatisfyConstraints(tree)==YES)
                            break;
                        }
#   if defined (DEBUG_CONSTRAINTS)
                    if (theTree->checkConstraints == YES && CheckConstraints (tree) == ERROR)
                        {
                        printf ("Error in constraints of starting tree\n");
                        getchar();
                        }
#   endif
                    FreePolyTree (constraintTree);
                    FreePolyTree (constraintTreeRef);
                    if (i==100)
                        {
                        YvyraPrint ("%s   Could not build a starting tree satisfying all constraints\n", spacer);                     
                        return (ERROR);
                        }
                    }
                /* random topology */
                else
                    {
                    if (tree->isRooted == YES)
                        {
                        if (BuildRandomRTopology (tree, &globalSeed) == ERROR)
                            return (ERROR);
                        }
                    else
                        {
                        if (BuildRandomUTopology (tree, &globalSeed) == ERROR)
                            return (ERROR);
                        if (MoveCalculationRoot (tree, localOutGroup) == ERROR)
                            return (ERROR);
                        }
                    }
                
                if (LabelTree (tree, localTaxonNames) == ERROR)
                    return (ERROR);
                if (q == p)
                    continue;   /* this is a parsimony tree without branch lengths */
                if (InitializeTreeCalibrations (tree) == ERROR)
                    return (ERROR);
                if (FillTopologySubParams(p, chn, 0, seed)== ERROR)
                    return (ERROR);
                }
            }
        }

    /* FillSpeciesTreeParams removed: species tree not supported in yvyra */

    return (NO_ERROR);
}


void FreeCppEvents (Param *p)
{
    int i, j;
    
    if (p->paramType != P_CPPEVENTS)
        return;

    if (p->nEvents != NULL)
        {
        free (p->nEvents[0]);
        free (p->nEvents);
        for (i=0; i<numGlobalChains; i++)
            {
            for (j=0; j<2*numLocalTaxa; j++)
                {
                free (p->position[2*i][j]);
                free (p->rateMult[2*i][j]);
                free (p->position[2*i+1][j]);
                free (p->rateMult[2*i+1][j]);
                }
            }
        free (p->position[0]);
        free (p->position);
        free (p->rateMult[0]);
        free (p->rateMult);
        p->nEvents = NULL;
        p->position = NULL;
        p->rateMult = NULL;
        }
}


int FreeModel (void)
{
    int             i;
    Param          *p;

    if (memAllocs[ALLOC_MODEL] == YES)
        {
        for (i = 0; i < numCurrentDivisions; i++)
            free (modelParams[i].activeConstraints);

        free (modelParams);
        free (modelSettings);
        memAllocs[ALLOC_MODEL] = NO;
        }

    if (memAllocs[ALLOC_MOVES] == YES)
        {
        for (i = 0; i < numApplicableMoves; i++)
            FreeMove (moves[i]);

        SAFEFREE (moves);
        numApplicableMoves = 0;
        memAllocs[ALLOC_MOVES] = NO;
        }

    if (memAllocs[ALLOC_COMPMATRIX] == YES)
        {
        free (compMatrix);
        memAllocs[ALLOC_COMPMATRIX] = NO;
        }

    if (memAllocs[ALLOC_NUMSITESOFPAT] == YES)
        {
        free (numSitesOfPat);
        memAllocs[ALLOC_NUMSITESOFPAT] = NO;
        }

    if (memAllocs[ALLOC_COMPCOLPOS] == YES)
        {
        free (compColPos);
        memAllocs[ALLOC_COMPCOLPOS] = NO;
        }

    if (memAllocs[ALLOC_COMPCHARPOS] == YES)
        {
        free (compCharPos);
        memAllocs[ALLOC_COMPCHARPOS] = NO;
        }

    if (memAllocs[ALLOC_ORIGCHAR] == YES)
        {
        free (origChar);
        memAllocs[ALLOC_ORIGCHAR] = NO;
        }

    if (memAllocs[ALLOC_STDTYPE] == YES)
        {
        free (stdType);
        memAllocs[ALLOC_STDTYPE] = NO;
        }

    if (memAllocs[ALLOC_STDSTATEFREQS] == YES)
        {
        SAFEFREE (stdStateFreqs);
        memAllocs[ALLOC_STDSTATEFREQS] = NO;
        }

    if (memAllocs[ALLOC_PARAMVALUES] == YES)
        {
        for (i = 0; i < numParams; i++)
            {
            p = &params[i];

            if (p->paramType == P_CPPEVENTS)
                FreeCppEvents (p);
            }

        SAFEFREE (paramValues);
        SAFEFREE (intValues);
        paramValsRowSize = intValsRowSize = 0;
        memAllocs[ALLOC_PARAMVALUES] = NO;
        }

    if (memAllocs[ALLOC_PARAMS] == YES)
        {
        for (i = 0; i < numParams; i++)
            {
            SAFEFREE (params[i].name);

            if (params[i].paramHeader)
                SAFEFREE (params[i].paramHeader);
            }

        SAFEFREE (params);
        SAFEFREE (relevantParts);
        numParams = 0;
        memAllocs[ALLOC_PARAMS] = NO;
        }

    if (memAllocs[ALLOC_MCMCTREES] == YES)
        {
            /* FIXME: Trees needs to be deallocated, but I can't figure
                out how many there are.  The loop below tries to free
                unallocated memory...
             */
            /*
        for (i = 0; i < numParams; ++i)
            {
            p = &params[i];

            for (j = 0; j < numGlobalChains; ++j)
                FreeTree (GetTree (p, j, 0));
            }
            */

        SAFEFREE (mcmcTree);
        SAFEFREE (subParamPtrs);
        memAllocs[ALLOC_MCMCTREES] = NO;
        }

    if (memAllocs[ALLOC_SYMPIINDEX] == YES)
        {
        free (sympiIndex);
        memAllocs[ALLOC_SYMPIINDEX] = NO;
        }

    if (memAllocs[ALLOC_LOCTAXANAMES] == YES)
        {
        free (localTaxonNames);
        memAllocs[ALLOC_LOCTAXANAMES] = NO;
        }

    if (memAllocs[ALLOC_LOCALTAXONCALIBRATION] == YES)
        {
        free (localTaxonCalibration);
        memAllocs[ALLOC_LOCALTAXONCALIBRATION] = NO;
        }

    return (NO_ERROR);
}


void FreeMove (MCMCMove *mv)
{
    free (mv->tuningParam[0]);
    free (mv->tuningParam);
    free (mv->relProposalProb);
    free (mv->nAccepted);
    free (mv->name);
    free (mv);
}


/* Compute empirical state freq are return it in global array empiricalFreqs[] */
int GetEmpiricalFreqs (int *relParts, int nRelParts)
{
    (void)relParts; (void)nRelParts;
    YvyraPrint ("%s   Empirical state frequencies not supported in yvyra (use equal or Dirichlet)\n", spacer);
    return (ERROR);
}


int GetNumDivisionChars (void)
{
    int         c, d, n;
    ModelInfo   *m;

    /* count number of characters in each division */
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        
        n = 0;
        for (c=0; c<numChar; c++)
            {
            if (charInfo[c].isExcluded == NO && partitionId[c][partitionNum] == d+1)
                n++;
            }
        m->numUncompressedChars = n;
        }
    
    return (NO_ERROR);
}


int *GetParamIntVals (Param *parm, int chain, int state)
{
    return parm->intValues + (2 * chain + state) * intValsRowSize;
}


YFlt  *GetParamStdStateFreqs (Param *parm, int chain, int state)
{
    return parm->stdStateFreqs + (2 * chain + state) * stdStateFreqsRowSize;
}


YFlt  *GetParamSubVals (Param *parm, int chain, int state)
{
    return parm->subValues + (2 * chain + state) * paramValsRowSize;
}


YFlt  *GetParamVals (Param *parm, int chain, int state)
{
    return parm->values + (2 * chain + state) * paramValsRowSize;
}
void GetPossibleRestrictionSites (int resSiteCode, int *sites)
{
    int     m;
    
    for (m=0; m<2; m++)
        sites[m] = 0;
        
    if (resSiteCode == 1)
        sites[0] = 1;
    else if (resSiteCode == 2)
        sites[1] = 1;
    else
        sites[0] = sites[1] = 1;

#   if 0
    printf ("%2d -- ", aaCode);
    for (m=0; m<20; m++)
        printf ("%d", aa[m]);
    printf ("\n");
#   endif
}


Tree *GetTree (Param *parm, int chain, int state)
{
    return mcmcTree[parm->treeIndex + ((2 * chain + state) * numTrees)];
}


Tree *GetTreeFromIndex (int index, int chain, int state)
{
    return mcmcTree[index + ((2 * chain + state) * numTrees)];
}


/*-----------------------------------------------------------
|
|   GetUserTreeFromName: Do case-insensitive search for user
|      tree, return index if match, -1 and ERROR if error
|
------------------------------------------------------------*/
int GetUserTreeFromName (int *index, char *treeName)
{
    int     i, j, k, nMatches;
    char    localName[100], temp[100];

    (*index) = -1;  /* appropriate return if no match */

    if ((int)strlen(treeName) > 99)
        {
        YvyraPrint ("%s   Too many characters in tree name\n", spacer);
        return (ERROR);
        }

    strcpy (localName, treeName);
    for (i=0; i<(int)strlen(localName); i++)
        localName[i] = tolower(localName[i]);

    nMatches = j = 0;
    for (i=0; i<numUserTrees; i++)
        {
        strcpy (temp, userTree[i]->name);
        for (k=0; k<(int)strlen(temp); k++)
            temp[k] = tolower(temp[k]);
        if (strcmp(localName,temp) == 0)
            {
            j = i;
            nMatches++;
            }
        }
    if (nMatches==0)
        {
        for (i=0; i<numUserTrees; i++)
            {
            strcpy (temp, userTree[i]->name);
            for (k=0; k<(int)strlen(temp); k++)
                temp[k] = tolower(temp[k]);
            if (strncmp(localName,temp,strlen(localName)) == 0)
                {
                j = i;
                nMatches++;
                }
            }
        }
    if (nMatches == 0)
        {
        YvyraPrint ("%s   Could not find tree '%s'\n", spacer, localName);  
        return (ERROR);
        }
    else if (nMatches > 1)
        {
        YvyraPrint ("%s   Several trees matched the abbreviated name '%s'\n", spacer, localName);
        return (ERROR);
        }
    else
        {
        (*index) = j;
        return (NO_ERROR);
        }
}


/*----------------------------------------------------------------------
 |
 |   InitializeChainTrees: 'Constructor' for chain trees
 |
 -----------------------------------------------------------------------*/
int InitializeChainTrees (Param *p, int from, int to, int isRooted)
{
    int             i, st, isCalibrated, isClock, nTaxa;
    Tree           *tree, **treeHandle;
    Model          *mp;

    mp = &modelParams[p->relParts[0]];

    if (p->paramType == P_SPECIESTREE)
        nTaxa = numSpecies;
    else
        nTaxa = numLocalTaxa;

    /* figure out whether the trees are clock */
    if (!strcmp (mp->brlensPr, "Clock"))
        isClock = YES;
    else
        isClock = NO;

    /* figure out whether the trees are calibrated */
    if (!strcmp (mp->brlensPr, "Clock")
            && (strcmp (mp->nodeAgePr, "Calibrated") == 0
                || strcmp (mp->clockRatePr, "Fixed") != 0
                || (strcmp (mp->clockRatePr, "Fixed") == 0
                    && AreDoublesEqual (mp->clockRateFix, 1.0, 1E-6) == NO)))
        isCalibrated = YES;
    else
        isCalibrated = NO;

    /* allocate space for and construct the trees */
    /* NOTE: The memory allocation scheme used here must match GetTree
             and GetTreeFromIndex */
    for (i = from; i < to; i++)
        {
        treeHandle = mcmcTree + p->treeIndex + 2 * i * numTrees;

        if (*treeHandle)
            free (*treeHandle);

        if ((*treeHandle = AllocateTree (nTaxa)) == NULL)
            {
            YvyraPrint ("%s   Problem allocating MCMC trees\n",
                          spacer);
            return (ERROR);
            }

        treeHandle = mcmcTree + p->treeIndex + (2 * i + 1) * numTrees;

        if (*treeHandle)
            free (*treeHandle);

        if ((*treeHandle = AllocateTree (nTaxa)) == NULL)
            {
            YvyraPrint ("%s   Problem allocating MCMC trees\n",
                          spacer);
            return (ERROR);
            }
        }

    /* initialize the trees */
    for (i = from; i < to; i++)
        {
        for (st = 0; st < 2; st++)
            {
            tree = GetTree (p, i, st);

            if (numTrees > 1)
                sprintf (tree->name, "mcmc.tree%d_%d",
                         p->treeIndex + 1, i + 1);
            else        /* if (numTrees == 1) */
                sprintf (tree->name, "mcmc.tree_%d", i + 1);

            tree->nRelParts = p->nRelParts;
            tree->relParts = p->relParts;
            tree->isRooted = isRooted;
            tree->isClock = isClock;
            tree->isCalibrated = isCalibrated;

            if (p->paramType == P_SPECIESTREE)
                {
                tree->nNodes = 2 * numSpecies;
                tree->nIntNodes = numSpecies - 1;
                }
            else if (tree->isRooted == YES)
                {
                tree->nNodes = 2 * numLocalTaxa;
                tree->nIntNodes = numLocalTaxa - 1;
                }
            else        /* if (tree->isRooted == NO) */
                {
                tree->nNodes = 2 * numLocalTaxa - 2;
                tree->nIntNodes = numLocalTaxa - 2;
                }

            if (p->checkConstraints == YES)
                {
                tree->checkConstraints = YES;
                tree->nLocks = NumInformativeHardConstraints (mp);
                tree->nConstraints = mp->numActiveConstraints;      /* nConstraints is number of constraints to check */
                tree->constraints = mp->activeConstraints;
                }
            else
                {
                tree->checkConstraints = NO;
                tree->nConstraints = tree->nLocks = 0;
                tree->constraints = NULL;
                }
            }
        }

    return (NO_ERROR);
}


int InitializeLinks (void)
{
    int         i, j;
    
    linkNum = 0;
    for (i=0; i<NUM_LINKED; i++)
        {
        for (j=0; j<numCurrentDivisions; j++)
            linkTable[i][j] = linkNum;
        }

    return (NO_ERROR);
}


/* InitializeTreeCalibrations: Set calibrations for tree nodes */
int InitializeTreeCalibrations (Tree *t)
{
    int         i;
    TreeNode    *p;
    
    if (t->isCalibrated == NO)
        return (NO_ERROR);
    
    /* Set tip calibrations */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->right == NULL && localTaxonCalibration[p->index]->prior != unconstrained)
            {
            p->isDated = YES;
            p->calibration = localTaxonCalibration[p->index];
            p->age = p->calibration->min;
            }
        else if (p->left == NULL && p->right == NULL)
            {
            p->isDated = NO;
            p->calibration = NULL;
            p->age = -1.0;
            } 
        }

    /* Initialize interior calibrations */
    if (CheckSetConstraints(t) == ERROR)
        return (ERROR);

    return (NO_ERROR);
}


int IsApplicable (Param *param)
{
    if (param == NULL)
        return NO;

    return YES;
}


int IsApplicable_ThreeTaxaOrMore (Param *param)
{
    if (param->paramType != P_TOPOLOGY )
        return NO;
        
    if (LargestMovableSubtree (param) >= 3)
        return YES;
    else
        return NO;
}


int IsApplicable_FourTaxaOrMore (Param *param)
{
    if (param->paramType != P_TOPOLOGY )
        return NO;
    
    if (LargestMovableSubtree (param) >= 4)
        return YES;
    else
        return NO;
}


int IsApplicable_FiveTaxaOrMore (Param *param)
{
    if (param->paramType != P_TOPOLOGY )
        return NO;
    
    if (LargestMovableSubtree (param) >= 5)
        return YES;
    else
        return NO;
}


int IsApplicable_TreeAgeMove (Param *param)
{
    Tree        *t;
    TreeNode    *p;

    if (param == NULL)
        return NO;

    if (param->paramType != P_BRLENS)
        return NO;
    
    t = GetTree (param, 0, 0);

    p = t->root->left;
    if (p->isDated == NO)
        return NO;
    if (p->calibration->prior == fixed)
        return NO;
    else
        return YES;
}


int IsApplicable_AncestralFossil (Param *param)
{
    ModelParams *mp = &modelParams[param->relParts[0]];

    if (!strcmp(mp->sampleStrat, "FossilTip"))
        return NO;
    else if (mp->fossilizationFix == 0.0)
        return NO;
    else  /* fossils may be ancestors of other fossils or extant species */
        return YES;
}


int IsModelSame (int whichParam, int part1, int part2, int *isApplic1, int *isApplic2)
{
    int         i, isSame, isFirstNucleotide, isSecondNucleotide, isFirstProtein, isSecondProtein, nDiff, temp1, temp2;

    isSame = YES;
    *isApplic1 = YES;
    *isApplic2 = YES;

    /* We cannot rely on SetModelInfo being called first so we need to be smart in figuring out model data type by looking at several model params */
    isFirstNucleotide = isSecondNucleotide = NO;
    if ((modelParams[part1].dataType == DNA || modelParams[part1].dataType == RNA) && strcmp(modelParams[part1].nucModel,"Protein") != 0)
        isFirstNucleotide = YES;
    if ((modelParams[part2].dataType == DNA || modelParams[part2].dataType == RNA) && strcmp(modelParams[part2].nucModel,"Protein") != 0)
        isSecondNucleotide = YES;       
    isFirstProtein = isSecondProtein = NO;
    if (modelParams[part1].dataType == PROTEIN || ((modelParams[part1].dataType == DNA || modelParams[part1].dataType == RNA) && !strcmp(modelParams[part1].nucModel,"Protein")))
        isFirstProtein = YES;
    if (modelParams[part2].dataType == PROTEIN || ((modelParams[part2].dataType == DNA || modelParams[part2].dataType == RNA) && !strcmp(modelParams[part2].nucModel,"Protein")))
        isSecondProtein = YES;      
    
    else if (whichParam == P_PI)
        {
        /* Check the state frequencies for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; /* part1 has a parsimony model and state frequencies do not apply */
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; /* part2 has a parsimony model and state frequencies do not apply */

        /* Check that the data are not CONTINUOUS for partitions 1 and 2 */
        if (modelParams[part1].dataType == CONTINUOUS)
            *isApplic1 = NO; /* state frequencies do not make sense for part1 */
        if (modelParams[part2].dataType == CONTINUOUS)
            *isApplic2 = NO; /* state frequencies do not make sense for part2 */
            
        /* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
        if (isFirstNucleotide != isSecondNucleotide)
            isSame = NO; /* data are not both nucleotide or both not nucleotide */
        else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
            isSame = NO; /* data are not the same */

        /* Check that the model structure is the same for both partitions */
        if (strcmp(modelParams[part1].nucModel, modelParams[part2].nucModel))
            isSame = NO; /* the nucleotide models are different */
        if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel) && !(!strcmp(modelParams[part1].nucModel, "Codon") && !strcmp(modelParams[part2].nucModel, "Codon")))
            isSame = NO; /* the models have different covarion struture */
 
        /* If both partitions have nucmodel=codon, then we also have to make certain that the same genetic code is used. */
        if (!strcmp(modelParams[part1].nucModel, "Codon") && !strcmp(modelParams[part2].nucModel, "Codon"))
            {
            if (strcmp(modelParams[part1].geneticCode, modelParams[part2].geneticCode))
                isSame = NO; /* the models have different genetic codes */
            }
        
        /* Let's see if the prior is the same. */
        if (modelParams[part1].dataType == STANDARD && modelParams[part2].dataType == STANDARD)
            {
            /* The data are morphological (STANDARD). The state frequencies are specified by a
               symmetric beta distribution, the parameter of which needs to be the same to apply to both
               partitions. Note that symPiPr = -1 is equivalent to setting the variance to 0.0. */
            /* In principle, the two partitions can be linked but if we need to take the possibility of
             * asymmetry in transition rates into account, then we do not know at this stage whether we
             * will have to sample from the state frequencies or integrate them out. Since these two
             * model levels are both represented by the same data structure, we do not allow partitions
             * to be linked unless symPiPr == -1, that is, if state frequencies are equal so that they
             * do not have to be sampled under any conditions. */
            if (strcmp(modelParams[part1].symPiPr,"Fixed") != 0 || strcmp(modelParams[part2].symPiPr,"Fixed") != 0)
                {
                isSame = NO;
                }
            else if (modelParams[part1].symBetaFix != -1 ||  modelParams[part2].symBetaFix != -1)
                {
                isSame = NO;
                }
            }
        if (modelSettings[part1].dataType == PROTEIN && modelSettings[part2].dataType == PROTEIN)
            {
            /* We are dealing with protein data. */
            if (strcmp(modelParams[part1].aaModelPr, modelParams[part2].aaModelPr) == 0)
                {
                if (strcmp(modelParams[part1].aaModelPr, "Fixed") == 0)
                    {
                    /* only have a single, fixed, amino acid rate matrix */
                    if (strcmp(modelParams[part1].aaModel, modelParams[part2].aaModel) != 0)
                        isSame = NO; /* we have different amino acid models, and the state frequencies must be different */
                    /* if we have an equalin model or Gtr model, then we need to check the prior on the state frequencies */
                    if (!strcmp(modelParams[part1].aaModel, "Equalin") || !strcmp(modelParams[part1].aaModel, "Gtr"))
                        {
                        if (!strcmp(modelParams[part1].stateFreqPr, modelParams[part2].stateFreqPr))
                            {
                            /* the prior form is the same */
                            if (!strcmp(modelParams[part1].stateFreqPr, "Dirichlet")) /* both prior models must be dirichlet */
                                {
                                for (i=0; i<modelParams[part1].nStates; i++)
                                    if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (YFlt) 0.00001) == NO)
                                        isSame = NO; /* the dirichlet parameters are different */
                                }
                            else /* both prior models must be fixed */
                                {
                                if (!strcmp(modelParams[part1].stateFreqsFixType, modelParams[part2].stateFreqsFixType))
                                    {
                                    /* if (!strcmp(modelParams[part1].stateFreqsFixType, "Empirical"))
                                        isSame = NO;     Even though it is unlikely that the empirical values for both partitions are exactly the same, we will
                                                         allow isSame to equal YES. This means pooled base frequencies are used to determine the empirical
                                                         base frequencies. The user can still unlink this parameter. */
                                    if (!strcmp(modelParams[part1].stateFreqsFixType, "User"))
                                        {
                                        for (i=0; i<modelParams[part1].nStates; i++)
                                            if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (YFlt) 0.00001) == NO)
                                                isSame = NO; /* the user-specified base frequencies are different */
                                        }
                                    /* if the frequencies were both fixed to "equal", we are golden, and they are the same */
                                    }
                                else
                                    isSame = NO; /* the fixed parameters must be the same. The only other possibility is that the
                                                    user specified equal or empirical for one partition and then specified specific
                                                    numbers (user) for the other _and_ happened to set the user values to the equal
                                                    or empirical values. We ignore this possibility. */
                                }
                            }
                        else
                            isSame = NO;
                        }
                    }
                else
                    {
                    /* averaging over models */
                    /* P_AAMODEL link check removed: molecular models not supported */
                    }
                }
            }
        else
            {
            /* Otherwise, we are dealing with RESTRICTION or NUCLEOTIDE data. The dirichlet should be the same
               for both partitions. */
            if (!strcmp(modelParams[part1].stateFreqPr, modelParams[part2].stateFreqPr))
                {
                /* the prior form is the same */
                if (!strcmp(modelParams[part1].stateFreqPr, "Dirichlet")) /* both prior models must be dirichlet */
                    {
                    for (i=0; i<modelParams[part1].nStates; i++)
                        if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (YFlt) 0.00001) == NO)
                            isSame = NO; /* the dirichlet parameters are different */
                    }
                else /* both prior models must be fixed */
                    {
                    if (!strcmp(modelParams[part1].stateFreqsFixType, modelParams[part2].stateFreqsFixType))
                        {
                        /* if (!strcmp(modelParams[part1].stateFreqsFixType, "Empirical"))
                            isSame = NO;     Even though it is unlikely that the empirical values for both partitions are exactly the same, we will
                                             allow isSame to equal YES. This means pooled base frequencies are used to determine the empirical
                                             base frequencies. The user can still unlink this parameter. */
                        if (!strcmp(modelParams[part1].stateFreqsFixType, "User"))
                            {
                            for (i=0; i<modelParams[part1].nStates; i++)
                                if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (YFlt) 0.00001) == NO)
                                    isSame = NO; /* the user-specified base frequencies are different */
                            }
                        /* if the frequencies were both fixed to "equal", we are golden, and they are the same */
                        }
                    else
                        isSame = NO; /* the fixed parameters must be the same. The only other possibility is that the
                                        user specified equal or empirical for one partition and then specified specific
                                        numbers (user) for the other _and_ happened to set the user values to the equal
                                        or empirical values. We ignore this possibility. */
                    }
                }
            else
                isSame = NO;

            /* also check if both partitions use Stationary or Directional model - in the latter case, check if priors on root frequencies match */
            if (!strcmp(modelParams[part1].statefreqModel, modelParams[part2].statefreqModel))
                {     
                if (!strcmp(modelParams[part1].statefreqModel, "Directional"))
                    {     
                    if (!strcmp(modelParams[part1].rootFreqPr, modelParams[part2].rootFreqPr)) 
                        {     
                        /* the prior form is the same */
                        if (!strcmp(modelParams[part1].rootFreqPr, "Dirichlet")) /* both prior models must be dirichlet */
                            {     
                            for (i=0; i<modelParams[part1].nStates; i++)
                                if (AreDoublesEqual (modelParams[part1].rootFreqsDir[i], modelParams[part2].rootFreqsDir[i], (YFlt) 0.00001) == NO)
                                    isSame = NO; /* the dirichlet parameters are different */
                            }     
                        else /* in this case both prior models must be fixed to user values */
                            {     
                            for (i=0; i<modelParams[part1].nStates; i++)
                                if (AreDoublesEqual (modelParams[part1].rootFreqsDir[i], modelParams[part2].rootFreqsDir[i], (YFlt) 0.00001) == NO)
                                    isSame = NO; /* the user-specified state frequencies are different */
                            }     
                        }     
                    else  
                        isSame = NO;
                    }
                }
            else
                isSame = NO;
            }

        /* Check to see if the state frequencies are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the state frequencies are inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_MIXTURE_RATES)
        {
        /* Check the mixture rate parameter for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; /* part1 has a parsimony model and mixture rate parameter does not apply */
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; /* part2 has a parsimony model and mixture rate parameter does not apply */
        
        /* Check that the data are not CONTINUOUS for partitions 1 and 2 */
        if (modelParams[part1].dataType == CONTINUOUS)
            *isApplic1 = NO; /* the mixture rate parameter does not make sense for part1 */
        if (modelParams[part2].dataType == CONTINUOUS)
            *isApplic2 = NO; /* the mixture rate parameter does not make sense for part2 */
        
        /* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
        if (isFirstNucleotide != isSecondNucleotide)
            isSame = NO; /* data are not both nucleotide */
        else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
            isSame = NO; /* data are not the same */
        
        /* Let's check that the mixture rate parameter is relevant for the two partitions */
        if (strcmp(modelParams[part1].ratesModel, "Kmixture") != 0)
            *isApplic1 = NO; /* the rate mixture parameter does not apply to part1 */
        if (strcmp(modelParams[part2].ratesModel, "Kmixture") != 0)
            *isApplic2 = NO; /* the rate mixture parameter does not apply to part2 */
        
        /* We may have a nucleotide model. Make certain the models are not of type codon. */
        if (!strcmp(modelParams[part1].nucModel, "Codon"))
            *isApplic1 = NO; /* we have a codon model for part1, and a shape parameter does not apply */
        if (!strcmp(modelParams[part2].nucModel, "Codon"))
            *isApplic2 = NO;/* we have a codon model for part2, and a shape parameter does not apply */
        
        /* Check that the model structure is the same for both partitions */
        if ((!strcmp(modelParams[part1].nucModel, "4by4") || !strcmp(modelParams[part1].nucModel, "Doublet")) && !strcmp(modelParams[part2].nucModel, "Codon"))
            isSame = NO; /* the nucleotide models are incompatible with the same rate mixture parameter */
        if ((!strcmp(modelParams[part2].nucModel, "4by4") || !strcmp(modelParams[part2].nucModel, "Doublet")) && !strcmp(modelParams[part1].nucModel, "Codon"))
            isSame = NO; /* the nucleotide models are incompatible with the same rate mixture parameter */
        
        /* if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
         isSame = NO; */ /* the models have different covarion struture */
        /* NOTE: Perhaps we should allow the possiblity that the shape parameter is the same for the case
         where one partition has a covarion model but the other does not and both datatypes are the same. */
        
        /* Check that the number of rate components is the same */
        if (modelParams[part1].numMixtCats != modelParams[part2].numMixtCats)
            isSame = NO; /* the number of rate components is not the same, so we cannot set the parameter to be equal for both partitions */
        
        /* Check that the priors are the same. */
        /* For now we only allow a flat Dirichlet prior, so this is not needed */
        
        /* Check to see if the rate mixture parameter is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the rate mixture parameter is inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_SHAPE)
        {
        /* Check the shape parameter for partitions 1 and 2 (this applies to the lnorm as well as various gamma models of rate variation across sites) */
        
        /* Check if the model is parsimony for either partition */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; /* part1 has a parsimony model and shape parameter does not apply */
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; /* part2 has a parsimony model and shape parameter does not apply */

        /* Check that the data are not CONTINUOUS for partitions 1 and 2 */
        if (modelParams[part1].dataType == CONTINUOUS)
            *isApplic1 = NO; /* the shape parameter does not make sense for part1 */
        if (modelParams[part2].dataType == CONTINUOUS)
            *isApplic2 = NO; /* the shape parameter does not make sense for part2 */

        /* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
        if (isFirstNucleotide != isSecondNucleotide)
            isSame = NO; /* data are not both nucleotide */
        else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
            isSame = NO; /* data are not the same */

        /* Let's check that the shape parameter is even relevant for the two partitions */
        if (!strcmp(modelParams[part1].ratesModel, "Equal") || !strcmp(modelParams[part1].ratesModel, "Propinv")
            || !strcmp(modelParams[part1].ratesModel, "Kmixture"))
            *isApplic1 = NO; /* the shape parameter does not make sense for part1 */
        if (!strcmp(modelParams[part2].ratesModel, "Equal") || !strcmp(modelParams[part2].ratesModel, "Propinv")
            || !strcmp(modelParams[part2].ratesModel, "Kmixture"))
            *isApplic2 = NO; /* the shape parameter does not make sense for part2 */

        /* Check that the model is either lnorm or gamma for both partitions */
        if (!strcmp(modelParams[part1].ratesModel, "Lnorm") && strcmp(modelParams[part1].ratesModel, "Lnorm") != 0)
            isSame = NO;    /* if the first is lnorm, the second must be lnorm */
        if (strcmp(modelParams[part1].ratesModel, "Lnorm") != 0 && !strcmp(modelParams[part1].ratesModel, "Lnorm"))
            isSame = NO;    /* if the first is not lnorm, the second cannot be lnorm */

        /* We may have a nucleotide model. Make certain the models are not of type codon. */
        if (!strcmp(modelParams[part1].nucModel, "Codon"))
            *isApplic1 = NO; /* we have a codon model for part1, and a shape parameter does not apply */
        if (!strcmp(modelParams[part2].nucModel, "Codon"))
            *isApplic2 = NO;/* we have a codon model for part2, and a shape parameter does not apply */

        /* Check that the model structure is the same for both partitions */
        if ((!strcmp(modelParams[part1].nucModel, "4by4") || !strcmp(modelParams[part1].nucModel, "Doublet")) && !strcmp(modelParams[part2].nucModel, "Codon"))
            isSame = NO; /* the nucleotide models are incompatible with the same shape parameter */
        if ((!strcmp(modelParams[part2].nucModel, "4by4") || !strcmp(modelParams[part2].nucModel, "Doublet")) && !strcmp(modelParams[part1].nucModel, "Codon"))
            isSame = NO; /* the nucleotide models are incompatible with the same shape parameter */
        /* if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
            isSame = NO; */ /* the models have different covarion struture */
        /* NOTE: Perhaps we should allow the possiblity that the shape parameter is the same for the case
                 where one partition has a covarion model but the other does not and both datatypes are the same. */

        /* Check that the number of rate categories is the same */
        if (!strcmp(modelParams[part1].ratesModel, "Lnorm") && (modelParams[part1].numLnormCats != modelParams[part2].numLnormCats))
            isSame = NO; /* the number of rate categories is not the same, so we cannot set the parameter to be equal for both partitions */
        else if (modelParams[part1].numGammaCats != modelParams[part2].numGammaCats)
            isSame = NO; /* the number of rate categories is not the same, so we cannot set the parameter to be equal for both partitions */

        /* Check that the priors are the same. */
        if (!strcmp(modelParams[part1].shapePr,"Uniform") && !strcmp(modelParams[part2].shapePr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].shapeUni[0], modelParams[part2].shapeUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].shapeUni[1], modelParams[part2].shapeUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].shapePr,"Exponential") && !strcmp(modelParams[part2].shapePr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].shapeExp, modelParams[part2].shapeExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].shapePr,"Fixed") && !strcmp(modelParams[part2].shapePr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].shapeFix, modelParams[part2].shapeFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
        
        /* Check to see if the shape parameter is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the shape parameter is inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_PINVAR)
        {
        /* Check the proportion of invariable sites parameter for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; /* part1 has a parsimony model and proportion of invariable sites parameter does not apply */
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; /* part2 has a parsimony model and proportion of invariable sites parameter does not apply */

        /* Check that the data are not CONTINUOUS for partitions 1 and 2 */
        if (modelParams[part1].dataType == CONTINUOUS)
            *isApplic1 = NO; /* the proportion of invariable sites parameter does not make sense for part1 */
        if (modelParams[part2].dataType == CONTINUOUS)
            *isApplic2 = NO; /* the proportion of invariable sites parameter does not make sense for part2 */

        /* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
        if (isFirstNucleotide != isSecondNucleotide)
            isSame = NO; /* data are not both nucleotide */
        else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
            isSame = NO; /* data are not the same */

        /* Let's check that proportion of invariable sites parameter is even relevant for the two partitions */
        if (!strcmp(modelParams[part1].ratesModel, "Equal") || !strcmp(modelParams[part1].ratesModel, "Gamma") ||
            !strcmp(modelParams[part1].ratesModel, "LNorm") || !strcmp(modelParams[part1].ratesModel, "Adgamma") ||
            !strcmp(modelParams[part1].ratesModel, "Kmixture"))
            *isApplic1 = NO; /* the proportion of invariable sites parameter does not make sense for part1 */
        if (!strcmp(modelParams[part2].ratesModel, "Equal") || !strcmp(modelParams[part2].ratesModel, "Gamma") ||
            !strcmp(modelParams[part2].ratesModel, "LNorm") || !strcmp(modelParams[part2].ratesModel, "Adgamma") ||
            !strcmp(modelParams[part1].ratesModel, "Kmixture"))
            *isApplic2 = NO; /* the proportion of invariable sites parameter does not make sense for part2 */
            
        /* It is not sensible to have a covarion model and a proportion of invariable sites */
        if (!strcmp(modelParams[part1].covarionModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].covarionModel, "Yes"))
            *isApplic2 = NO;
        
        /* We have a nucleotide model. Make certain the models are not of type codon. */
        if (!strcmp(modelParams[part1].nucModel, "Codon"))
            *isApplic1 = NO; /* we have a codon model for part1, and a proportion of invariable sites parameter does not apply */
        if (!strcmp(modelParams[part2].nucModel, "Codon"))
            *isApplic2 = NO;/* we have a codon model for part2, and a proportion of invariable sites parameter does not apply */

        /* Check that the model structure is the same for both partitions */
        if (strcmp(modelParams[part1].nucModel, modelParams[part2].nucModel))
            isSame = NO; /* the nucleotide models are different */
        if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
            isSame = NO; /* the models have different covarion struture */
        
        /* check the priors */
        if (!strcmp(modelParams[part1].pInvarPr,"Uniform") && !strcmp(modelParams[part2].pInvarPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].pInvarUni[0], modelParams[part2].pInvarUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].pInvarUni[1], modelParams[part2].pInvarUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].pInvarPr,"Fixed") && !strcmp(modelParams[part2].pInvarPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].pInvarFix, modelParams[part2].pInvarFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
        
        /* Check to see if the switching rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the switching rates are inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_CORREL)
        {
        /* Check the autocorrelation parameter for gamma rates on partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; /* part1 has a parsimony model and autocorrelation parameter does not apply */
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; /* part2 has a parsimony model and autocorrelation parameter does not apply */

        /* Check that the data are either DNA, RNA, or PROTEIN for partitions 1 and 2 */
        if (modelSettings[part1].dataType != DNA && modelSettings[part1].dataType != RNA && modelSettings[part1].dataType != PROTEIN)
            *isApplic1 = NO; /* the switching rates do not make sense for part1 */
        if (modelSettings[part2].dataType != DNA && modelSettings[part2].dataType != RNA && modelSettings[part2].dataType != PROTEIN)
            *isApplic2 = NO; /* the switching rates do not make sense for part2 */
            
        /* Now, check that the data are the same (i.e., both nucleotide or both amino acid). */
        if (isFirstNucleotide != isSecondNucleotide)
            isSame = NO; /* one or the other is nucleotide, so they cannot be the same */
        else if (modelSettings[part1].dataType != modelSettings[part2].dataType && isFirstNucleotide == NO)
            isSame = NO; /* data are not both nucleotide or both amino acid */

        /* Let's check that autocorrelation parameter is even relevant for the two partitions */
        if (strcmp(modelParams[part1].ratesModel, "Adgamma"))
            *isApplic1 = NO; /* the autocorrelation parameter does not make sense for part1 */
        if (strcmp(modelParams[part2].ratesModel, "Adgamma"))
            *isApplic2 = NO; /* the autocorrelation parameter does not make sense for part2 */

        /* Assuming that we have a nucleotide model, make certain the models are not of type codon. */
        if (!strcmp(modelParams[part1].nucModel, "Codon"))
            *isApplic1 = NO; /* we have a codon model for part1, and a autocorrelation parameter does not apply */
        if (!strcmp(modelParams[part2].nucModel, "Codon"))
            *isApplic2 = NO; /* we have a codon model for part2, and a autocorrelation parameter does not apply */
        
        /* Check that the model structure is the same for both partitions */
        if (strcmp(modelParams[part1].nucModel, modelParams[part2].nucModel))
            isSame = NO; /* the nucleotide models are different */
        if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
            isSame = NO; /* the models have different covarion struture */

        /* Check the priors for both partitions. */
        if (!strcmp(modelParams[part1].adGammaCorPr,"Uniform") && !strcmp(modelParams[part2].adGammaCorPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].adgCorrUni[0], modelParams[part2].adgCorrUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].adgCorrUni[1], modelParams[part2].adgCorrUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].adGammaCorPr,"Fixed") && !strcmp(modelParams[part2].adGammaCorPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].adgCorrFix, modelParams[part2].adgCorrFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

        /* Check to see if the switching rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the switching rates are inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_SWITCH)
        {
        /* Check the covarion switching rates on partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; /* part1 has a parsimony model and switching rates do not apply */
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; /* part2 has a parsimony model and switching rates do not apply */
        
        /* Check that the data are either DNA, RNA, or PROTEIN for partitions 1 and 2 */
        if (modelSettings[part1].dataType != DNA && modelSettings[part1].dataType != RNA && modelSettings[part1].dataType != PROTEIN)
            *isApplic1 = NO; /* the switching rates do not make sense for part1 */
        if (modelSettings[part2].dataType != DNA && modelSettings[part2].dataType != RNA && modelSettings[part2].dataType != PROTEIN)
            *isApplic2 = NO; /* the switching rates do not make sense for part2 */
            
        /* Now, check that the data are the same (i.e., both nucleotide or both amino acid). */
        if (isFirstNucleotide != isSecondNucleotide)
            isSame = NO; /* one or the other is nucleotide, so they cannot be the same */
        else if (modelSettings[part1].dataType != modelSettings[part2].dataType && isFirstNucleotide == NO)
            isSame = NO; /* data are not both nucleotide or both amino acid */

        /* Lets check that covarion model has been selected for partitions 1 and 2 */
        if (!strcmp(modelParams[part1].covarionModel, "No"))
            *isApplic1 = NO; /* the switching rates do not make sense for part1 */
        if (!strcmp(modelParams[part2].covarionModel, "No"))
            *isApplic2 = NO; /* the switching rates do not make sense for part2 */

        /* If we have a nucleotide model make certain the models are not of type codon or doublet. */
        if (!strcmp(modelParams[part1].nucModel, "Codon") || !strcmp(modelParams[part1].nucModel, "Doublet"))
            *isApplic1 = NO; /* we have a codon model for part1, and a covarion switch parameter does not apply */
        if (!strcmp(modelParams[part2].nucModel, "Codon") || !strcmp(modelParams[part2].nucModel, "Doublet"))
            *isApplic2 = NO; /* we have a codon model for part2, and a covarion switch parameter does not apply */

        /* Check that the priors are the same. */
        if (!strcmp(modelParams[part1].covSwitchPr,"Uniform") && !strcmp(modelParams[part2].covSwitchPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].covswitchUni[0], modelParams[part2].covswitchUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].covswitchUni[1], modelParams[part2].covswitchUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].covSwitchPr,"Exponential") && !strcmp(modelParams[part2].covSwitchPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].covswitchExp, modelParams[part2].covswitchExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].covSwitchPr,"Fixed") && !strcmp(modelParams[part2].covSwitchPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].covswitchFix[0], modelParams[part2].covswitchFix[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].covswitchFix[1], modelParams[part2].covswitchFix[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

        /* Check to see if the switching rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the switching rates are inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_RATEMULT)
        {
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
           silly and doesn't take this information into account) and a rate multiplier is nonsensical. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
        
        /* Check that the branch lengths are at least proportional. */
        if (IsModelSame (P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;

        /* See if the rate prior is the same for the partitions */
        if (strcmp(modelParams[part1].ratePr, modelParams[part2].ratePr) != 0)
            isSame = NO;

        /* Check to see if rate multipliers are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; 
            
        }
    else if (whichParam == P_TOPOLOGY)
        {
        /* Check the topology for partitions 1 and 2. */
        
        /* If the prior is different, then the topologies cannot be the same. */
        if (strcmp(modelParams[part1].topologyPr, modelParams[part2].topologyPr))
            isSame = NO;

        /* If both partitions have topologies constrained, then we need to make certain that the constraints are the same. */
        /* This also guarantees that any calibrations will be the same. */
        if (!strcmp(modelParams[part1].topologyPr, "Constraints") && !strcmp(modelParams[part2].topologyPr, "Constraints"))
            {
            if (modelParams[part1].numActiveConstraints != modelParams[part2].numActiveConstraints)
                isSame = NO;
            else
                {
                nDiff = 0;
                for (i=0; i<numDefinedConstraints; i++)
                    if (modelParams[part1].activeConstraints[i] != modelParams[part2].activeConstraints[i])
                        nDiff++;
                if (nDiff != 0)
                    isSame = NO;
                }
            }
        }
    else if (whichParam == P_BRLENS)
        {
        /* Check the branch lengths for partitions 1 and 2. */

        /* First, if the topologies are different, the same branch lengths cannot apply. */
        if (IsModelSame (P_TOPOLOGY, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_TOPOLOGY][part1] != linkTable[P_TOPOLOGY][part2])
            isSame = NO;

        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
           silly and doesn't take this information into account). */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check to see if the branch lengths are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;

        /* Make sure the branch lengths have the same priors and are not unlinked */
        if (*isApplic1 == YES && *isApplic2 == YES)
            {
            /* We are dealing with real branch lengths (not parsimony) for both partitions */
            
            /* The branch length prior should be the same */
            if (strcmp(modelParams[part1].brlensPr, modelParams[part2].brlensPr))
                isSame = NO;
                
            /* if both partitions have unconstrained brlens, then we need to check that the priors on the branch lengths are the same */
            if (!strcmp(modelParams[part1].brlensPr, "Unconstrained") && !strcmp(modelParams[part2].brlensPr, "Unconstrained"))
                {
                if (strcmp(modelParams[part1].unconstrainedPr, modelParams[part2].unconstrainedPr))
                    isSame = NO;
                else
                    {
                    if (!strcmp(modelParams[part1].unconstrainedPr, "Uniform"))
                        {
                        if (AreDoublesEqual (modelParams[part1].brlensUni[0], modelParams[part2].brlensUni[0], (YFlt) 0.00001) == NO)
                            isSame = NO;
                        if (AreDoublesEqual (modelParams[part1].brlensUni[1], modelParams[part2].brlensUni[1], (YFlt) 0.00001) == NO)
                            isSame = NO;
                        }
                    else if (!strcmp(modelParams[part1].unconstrainedPr, "Exponential"))
                        {
                        if (AreDoublesEqual (modelParams[part1].brlensExp, modelParams[part2].brlensExp, (YFlt) 0.00001) == NO)
                            isSame = NO;
                        }
                    else if (!strcmp(modelParams[part1].unconstrainedPr, "twoExp"))
                        {
                            if (AreDoublesEqual (modelParams[part1].brlens2Exp[0], modelParams[part2].brlens2Exp[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].brlens2Exp[1], modelParams[part2].brlens2Exp[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                        }
                    else
                        {
                        if (AreDoublesEqual (modelParams[part1].brlensDir[0], modelParams[part2].brlensDir[0], (YFlt) 0.00001) == NO)  
                            isSame = NO;
                        if (AreDoublesEqual (modelParams[part1].brlensDir[1], modelParams[part2].brlensDir[1], (YFlt) 0.00001) == NO) 
                            isSame = NO;
                        if (AreDoublesEqual (modelParams[part1].brlensDir[2], modelParams[part2].brlensDir[2], (YFlt) 0.00001) == NO)  
                            isSame = NO;
                        if (AreDoublesEqual (modelParams[part1].brlensDir[3], modelParams[part2].brlensDir[3], (YFlt) 0.00001) == NO) 
                            isSame = NO;
                        }   
                    }
                }
            
            /* if both partitions have clock brlens, then we need to check that the priors on the clock are the same */
            if (!strcmp(modelParams[part1].brlensPr, "Clock") && !strcmp(modelParams[part2].brlensPr, "Clock"))
                {
                if (strcmp(modelParams[part1].clockPr, modelParams[part2].clockPr))
                    isSame = NO;
                else
                    {
                    if (!strcmp(modelParams[part1].clockPr, "Birthdeath"))
                        {
                        if (!strcmp(modelParams[part1].speciationPr,"Uniform") && !strcmp(modelParams[part2].speciationPr,"Uniform"))
                            {
                            if (AreDoublesEqual (modelParams[part1].speciationUni[0], modelParams[part2].speciationUni[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].speciationUni[1], modelParams[part2].speciationUni[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].speciationPr,"Exponential") && !strcmp(modelParams[part2].speciationPr,"Exponential"))
                            {
                            if (AreDoublesEqual (modelParams[part1].speciationExp, modelParams[part2].speciationExp, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].speciationPr,"Fixed") && !strcmp(modelParams[part2].speciationPr,"Fixed"))
                            {
                            if (AreDoublesEqual (modelParams[part1].speciationFix, modelParams[part2].speciationFix, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else
                            isSame = NO;

                        if (!strcmp(modelParams[part1].extinctionPr,"Beta") && !strcmp(modelParams[part2].extinctionPr,"Beta"))
                            {
                            if (AreDoublesEqual (modelParams[part1].extinctionBeta[0], modelParams[part2].extinctionBeta[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].extinctionBeta[1], modelParams[part2].extinctionBeta[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].extinctionPr,"Fixed") && !strcmp(modelParams[part2].extinctionPr,"Fixed"))
                            {
                            if (AreDoublesEqual (modelParams[part1].extinctionFix, modelParams[part2].extinctionFix, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else
                            isSame = NO;

                        if (AreDoublesEqual (modelParams[part1].sampleProb, modelParams[part2].sampleProb, 0.00001) == NO)
                            isSame = NO;
                        if (strcmp(modelParams[part1].sampleStrat,modelParams[part2].sampleStrat))
                            isSame = NO;
                        }
                    else if (!strcmp(modelParams[part1].clockPr, "Coalescence") || !strcmp(modelParams[part1].clockPr, "Speciestreecoalescence"))
                        {
                        if (!strcmp(modelParams[part1].popSizePr,"Uniform") && !strcmp(modelParams[part2].popSizePr,"Uniform"))
                            {
                            if (AreDoublesEqual (modelParams[part1].popSizeUni[0], modelParams[part2].popSizeUni[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].popSizeUni[1], modelParams[part2].popSizeUni[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].popSizePr,"Lognormal") && !strcmp(modelParams[part2].popSizePr,"Lognormal"))
                            {
                            if (AreDoublesEqual (modelParams[part1].popSizeLognormal[0], modelParams[part2].popSizeLognormal[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].popSizeLognormal[1], modelParams[part2].popSizeLognormal[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].popSizePr,"Normal") && !strcmp(modelParams[part2].popSizePr,"Normal"))
                            {
                            if (AreDoublesEqual (modelParams[part1].popSizeNormal[0], modelParams[part2].popSizeNormal[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].popSizeNormal[1], modelParams[part2].popSizeNormal[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].popSizePr,"Gamma") && !strcmp(modelParams[part2].popSizePr,"Gamma"))
                            {
                            if (AreDoublesEqual (modelParams[part1].popSizeGamma[0], modelParams[part2].popSizeGamma[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].popSizeGamma[1], modelParams[part2].popSizeGamma[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].popSizePr,"Fixed") && !strcmp(modelParams[part2].popSizePr,"Fixed"))
                            {
                            if (AreDoublesEqual (modelParams[part1].popSizeFix, modelParams[part2].popSizeFix, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else
                            isSame = NO;

                        if (strcmp(modelParams[part1].ploidy, modelParams[part2].ploidy) != 0)
                            isSame = NO;
                        }
                    if (strcmp(modelParams[part1].clockPr, "Uniform") == 0 && strcmp(modelParams[part1].nodeAgePr, "Calibrated") != 0)
                        {
                        if (modelParams[part1].treeAgePr.prior != modelParams[part2].treeAgePr.prior)
                            isSame = NO;
                        if (modelParams[part1].treeAgePr.prior == fixed)
                            {
                            if (AreDoublesEqual (modelParams[part1].treeAgePr.priorParams[0], modelParams[part2].treeAgePr.priorParams[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (modelParams[part1].treeAgePr.prior == offsetLogNormal ||
                            modelParams[part1].treeAgePr.prior == truncatedNormal ||
                            modelParams[part1].treeAgePr.prior == offsetGamma)
                            {
                            for (i=0; i<3; i++)
                                {
                                if (AreDoublesEqual (modelParams[part1].treeAgePr.priorParams[i], modelParams[part2].treeAgePr.priorParams[i], (YFlt) 0.00001) == NO)
                                    isSame = NO;
                                }
                            }
                        else
                            {
                            for (i=0; i<2; i++)
                                {
                                if (AreDoublesEqual (modelParams[part1].treeAgePr.priorParams[i], modelParams[part2].treeAgePr.priorParams[i], (YFlt) 0.00001) == NO)
                                    isSame = NO;
                                }
                            }
                        }
                    else if (strcmp(modelParams[part1].clockPr, "Fossilization") == 0)
                        {
                        if (modelParams[part1].treeAgePr.prior != modelParams[part2].treeAgePr.prior)
                            isSame = NO;
                        if (modelParams[part1].treeAgePr.prior == fixed)
                            {
                            if (AreDoublesEqual (modelParams[part1].treeAgePr.priorParams[0], modelParams[part2].treeAgePr.priorParams[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (modelParams[part1].treeAgePr.prior == offsetLogNormal ||
                            modelParams[part1].treeAgePr.prior == truncatedNormal ||
                            modelParams[part1].treeAgePr.prior == offsetGamma)
                            {
                            for (i=0; i<3; i++)
                                {
                                if (AreDoublesEqual (modelParams[part1].treeAgePr.priorParams[i], modelParams[part2].treeAgePr.priorParams[i], (YFlt) 0.00001) == NO)
                                    isSame = NO;
                                }
                            }
                        else
                            {
                            for (i=0; i<2; i++)
                                {
                                if (AreDoublesEqual (modelParams[part1].treeAgePr.priorParams[i], modelParams[part2].treeAgePr.priorParams[i], (YFlt) 0.00001) == NO)
                                    isSame = NO;
                                }
                            }
                        
                        if (!strcmp(modelParams[part1].speciationPr,"Uniform") && !strcmp(modelParams[part2].speciationPr,"Uniform"))
                            {
                            if (AreDoublesEqual (modelParams[part1].speciationUni[0], modelParams[part2].speciationUni[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].speciationUni[1], modelParams[part2].speciationUni[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].speciationPr,"Exponential") && !strcmp(modelParams[part2].speciationPr,"Exponential"))
                            {
                            if (AreDoublesEqual (modelParams[part1].speciationExp, modelParams[part2].speciationExp, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].speciationPr,"Fixed") && !strcmp(modelParams[part2].speciationPr,"Fixed"))
                            {
                            if (AreDoublesEqual (modelParams[part1].speciationFix, modelParams[part2].speciationFix, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else
                            isSame = NO;
                        
                        if (!strcmp(modelParams[part1].extinctionPr,"Beta") && !strcmp(modelParams[part2].extinctionPr,"Beta"))
                            {
                            if (AreDoublesEqual (modelParams[part1].extinctionBeta[0], modelParams[part2].extinctionBeta[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].extinctionBeta[1], modelParams[part2].extinctionBeta[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].extinctionPr,"Exponential") && !strcmp(modelParams[part2].extinctionPr,"Exponential"))
                            {
                            if (AreDoublesEqual (modelParams[part1].extinctionExp, modelParams[part2].extinctionExp, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].extinctionPr,"Fixed") && !strcmp(modelParams[part2].extinctionPr,"Fixed"))
                            {
                            if (AreDoublesEqual (modelParams[part1].extinctionFix, modelParams[part2].extinctionFix, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else
                            isSame = NO;
                        
                        if (!strcmp(modelParams[part1].fossilizationPr,"Beta") && !strcmp(modelParams[part2].fossilizationPr,"Beta"))
                            {
                            if (AreDoublesEqual (modelParams[part1].fossilizationBeta[0], modelParams[part2].fossilizationBeta[0], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            if (AreDoublesEqual (modelParams[part1].fossilizationBeta[1], modelParams[part2].fossilizationBeta[1], (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].fossilizationPr,"Exponential") && !strcmp(modelParams[part2].fossilizationPr,"Exponential"))
                            {
                            if (AreDoublesEqual (modelParams[part1].fossilizationExp, modelParams[part2].fossilizationExp, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else if (!strcmp(modelParams[part1].fossilizationPr,"Fixed") && !strcmp(modelParams[part2].fossilizationPr,"Fixed"))
                            {
                            if (AreDoublesEqual (modelParams[part1].fossilizationFix, modelParams[part2].fossilizationFix, (YFlt) 0.00001) == NO)
                                isSame = NO;
                            }
                        else
                            isSame = NO;
                        
                        if (AreDoublesEqual (modelParams[part1].sampleProb, modelParams[part2].sampleProb, 0.00001) == NO)
                            isSame = NO;
                        if (strcmp(modelParams[part1].sampleStrat,modelParams[part2].sampleStrat))
                            isSame = NO;
                        }
                    }

                /* if the same clock prior, we need to check calibrations */
                if (strcmp(modelParams[part1].nodeAgePr,modelParams[part2].nodeAgePr) != 0)
                    isSame = NO;
                
                /* If fixed clock brlens, check if the brlens come from the same tree */
                if (!strcmp(modelParams[part1].clockPr, "Fixed") && !strcmp(modelParams[part2].clockPr, "Fixed"))
                    {
                    if (modelParams[part1].brlensFix != modelParams[part2].brlensFix)
                        isSame = NO;
                    }
                }
            /* If fixed brlens, check if the brlens come from the same tree */
            if (!strcmp(modelParams[part1].brlensPr, "Fixed") && !strcmp(modelParams[part2].brlensPr, "Fixed"))
                {
                if (modelParams[part1].brlensFix != modelParams[part2].brlensFix)
                    isSame = NO;
                }
            }
        }
    else if (whichParam == P_SPECRATE)
        {
        /* Check the speciation rates for partitions 1 and 2. */

        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
           silly and doesn't take this information into account) and a speciation rate cannot be estimated. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
            
        /* Check that the branch length prior is a clock:birthdeath or clock:fossilization for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        if (strcmp(modelParams[part1].clockPr, "Birthdeath") != 0 && strcmp(modelParams[part1].clockPr, "Fossilization") != 0)
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockPr, "Birthdeath") != 0 && strcmp(modelParams[part2].clockPr, "Fossilization") != 0)
            *isApplic2 = NO;
        
        /* Now, check that the prior on the speciation rates are the same. */
        if (!strcmp(modelParams[part1].speciationPr,"Uniform") && !strcmp(modelParams[part2].speciationPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].speciationUni[0], modelParams[part2].speciationUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].speciationUni[1], modelParams[part2].speciationUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].speciationPr,"Exponential") && !strcmp(modelParams[part2].speciationPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].speciationExp, modelParams[part2].speciationExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].speciationPr,"Fixed") && !strcmp(modelParams[part2].speciationPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].speciationFix, modelParams[part2].speciationFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if the speciation rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_EXTRATE)
        {
        /* Check the extinction rates for partitions 1 and 2. */

        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
           silly and doesn't take this information into account) and a extinction rate cannot be estimated. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
            
        /* Check that the branch length prior is a clock:birthdeath or clock:fossilization for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        if (strcmp(modelParams[part1].clockPr, "Birthdeath")!= 0 && strcmp(modelParams[part1].clockPr, "Fossilization") != 0)
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockPr, "Birthdeath")!= 0 && strcmp(modelParams[part2].clockPr, "Fossilization") != 0)
            *isApplic2 = NO;
        
        /* Now, check that the prior on the extinction rates are the same. */
        if (!strcmp(modelParams[part1].extinctionPr,"Beta") && !strcmp(modelParams[part2].extinctionPr,"Beta"))
            {
            if (AreDoublesEqual (modelParams[part1].extinctionBeta[0], modelParams[part2].extinctionBeta[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].extinctionBeta[1], modelParams[part2].extinctionBeta[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].extinctionPr,"Exponential") && !strcmp(modelParams[part2].extinctionPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].extinctionExp, modelParams[part2].extinctionExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].extinctionPr,"Fixed") && !strcmp(modelParams[part2].extinctionPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].extinctionFix, modelParams[part2].extinctionFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if the extinction rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_FOSLRATE)
        {
        /* Check the fossilization rates for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;
        
        /* Check that the branch length prior is a clock:fossilization for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        if (strcmp(modelParams[part1].clockPr, "Fossilization") != 0)
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockPr, "Fossilization") != 0)
            *isApplic2 = NO;
        
        /* Now, check that the prior on the fossilization rates are the same. */
        if (!strcmp(modelParams[part1].fossilizationPr,"Beta") && !strcmp(modelParams[part2].fossilizationPr,"Beta"))
            {
            if (AreDoublesEqual (modelParams[part1].fossilizationBeta[0], modelParams[part2].fossilizationBeta[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].fossilizationBeta[1], modelParams[part2].fossilizationBeta[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].fossilizationPr,"Exponential") && !strcmp(modelParams[part2].fossilizationPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].fossilizationExp, modelParams[part2].fossilizationExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].fossilizationPr,"Fixed") && !strcmp(modelParams[part2].fossilizationPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].fossilizationFix, modelParams[part2].fossilizationFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if the fossilization rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_POPSIZE)
        {
        /* Check population size for partitions 1 and 2. */

        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
           silly and doesn't take this information into account) and population size cannot be estimated. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
            
        /* Check that the branch length prior is a clock:coalescence or clock:speciestreecoalescence for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        if (strcmp(modelParams[part1].clockPr, "Coalescence") != 0 && strcmp(modelParams[part1].clockPr, "Speciestreecoalescence") != 0)
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockPr, "Coalescence") != 0 && strcmp(modelParams[part2].clockPr, "Speciestreecoalescence") != 0)
            *isApplic2 = NO;
        
        /* Now, check that the prior on population size is the same. */
        if (!strcmp(modelParams[part1].popSizePr,"Uniform") && !strcmp(modelParams[part2].popSizePr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].popSizeUni[0], modelParams[part2].popSizeUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].popSizeUni[1], modelParams[part2].popSizeUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].popSizePr,"Lognormal") && !strcmp(modelParams[part2].popSizePr,"Lognormal"))
            {
            if (AreDoublesEqual (modelParams[part1].popSizeLognormal[0], modelParams[part2].popSizeLognormal[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].popSizeLognormal[1], modelParams[part2].popSizeLognormal[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].popSizePr,"Normal") && !strcmp(modelParams[part2].popSizePr,"Normal"))
            {
            if (AreDoublesEqual (modelParams[part1].popSizeNormal[0], modelParams[part2].popSizeNormal[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].popSizeNormal[1], modelParams[part2].popSizeNormal[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].popSizePr,"Gamma") && !strcmp(modelParams[part2].popSizePr,"Gamma"))
            {
            if (AreDoublesEqual (modelParams[part1].popSizeGamma[0], modelParams[part2].popSizeGamma[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].popSizeGamma[1], modelParams[part2].popSizeGamma[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].popSizePr,"Fixed") && !strcmp(modelParams[part2].popSizePr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].popSizeFix, modelParams[part2].popSizeFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if population size is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_GROWTH)
        {
        /* Check growth rate for partitions 1 and 2. */

        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
           silly and doesn't take this information into account) and growth rate cannot be estimated. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
            
        /* Check that the branch length prior is a clock:coalescence for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        if (strcmp(modelParams[part1].clockPr, "Coalescence"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockPr, "Coalescence"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on growth rate is the same. */
        if (!strcmp(modelParams[part1].growthPr,"Uniform") && !strcmp(modelParams[part2].growthPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].growthUni[0], modelParams[part2].growthUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].growthUni[1], modelParams[part2].growthUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].growthPr,"Exponential") && !strcmp(modelParams[part2].growthPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].growthExp, modelParams[part2].growthExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].growthPr,"Fixed") && !strcmp(modelParams[part2].growthPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].growthFix, modelParams[part2].growthFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if growth rate is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; 
        }
    else if (whichParam == P_BMCORR)
        {
        /* Check the correlation parameter for brownian motion 1 and 2. */
        
        /* Check that the data are either CONTINUOUS for partitions 1 and 2 */
        if (modelParams[part1].dataType != CONTINUOUS)
            *isApplic1 = NO; /* the correlation parameter does not make sense for part1 */
        if (modelParams[part2].dataType != CONTINUOUS)
            *isApplic2 = NO; /* the correlation parameter does not make sense for part2 */
            
        /* Now, check that the data are the same. */
        if (modelParams[part1].dataType != modelParams[part2].dataType)
            isSame = NO; /* data are not both continuous */

        /* Check the priors for both partitions. */
        if (!strcmp(modelParams[part1].brownCorrPr,"Uniform") && !strcmp(modelParams[part2].brownCorrPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].brownCorrUni[0], modelParams[part2].brownCorrUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].brownCorrUni[1], modelParams[part2].brownCorrUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].brownCorrPr,"Fixed") && !strcmp(modelParams[part2].brownCorrPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].brownCorrFix, modelParams[part2].brownCorrFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

        /* Check to see if the correlation parameters are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the correlation parameters are inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_BMSIGMA)
        {
        /* Check the sigma parameter for brownian motion 1 and 2. */
        
        /* Check that the data are either CONTINUOUS for partitions 1 and 2 */
        if (modelParams[part1].dataType != CONTINUOUS)
            *isApplic1 = NO; /* the sigma parameter does not make sense for part1 */
        if (modelParams[part2].dataType != CONTINUOUS)
            *isApplic2 = NO; /* the sigma parameter does not make sense for part2 */
            
        /* Now, check that the data are the same. */
        if (modelParams[part1].dataType != modelParams[part2].dataType)
            isSame = NO; /* data are not both continuous */

        /* Check the priors for both partitions. */
        if (!strcmp(modelParams[part1].brownScalePr,"Uniform") && !strcmp(modelParams[part2].brownScalePr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].brownScaleUni[0], modelParams[part2].brownScaleUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].brownScaleUni[1], modelParams[part2].brownScaleUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].brownScalePr,"Fixed") && !strcmp(modelParams[part2].brownScalePr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].brownScaleFix, modelParams[part2].brownScaleFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].brownScalePr,"Gamma") && !strcmp(modelParams[part2].brownScalePr,"Gamma"))
            {
            if (AreDoublesEqual (modelParams[part1].brownScaleGamma[0], modelParams[part2].brownScaleGamma[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].brownScaleGamma[1], modelParams[part2].brownScaleGamma[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

        /* Check to see if the sigma parameters are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO; /* if the sigma parameters are inapplicable for either partition, then the parameter cannot be the same */
        }
    else if (whichParam == P_CPPRATE)
        {
        /* Check cpp rate for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
        silly and doesn't take this information into account) and cpp rate cannot be estimated. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is cpp for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "Cpp"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "Cpp"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on cpp rate is the same. */
        if (!strcmp(modelParams[part1].cppRatePr,"Exponential") && !strcmp(modelParams[part2].cppRatePr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].cppRateExp, modelParams[part2].cppRateExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].cppRatePr,"Fixed") && !strcmp(modelParams[part2].cppRatePr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].cppRateFix, modelParams[part2].cppRateFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
    
        /* Check to see if cpp rate is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;    
        }
    else if (whichParam == P_CPPMULTDEV)
        {
        /* Check cpp multiplier deviation prior for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
        silly and doesn't take this information into account) and this parameter is inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is cpp for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "Cpp"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "Cpp"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on sigma is the same. */
        if (!strcmp(modelParams[part1].cppMultDevPr,"Fixed") && !strcmp(modelParams[part2].cppMultDevPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].cppMultDevFix, modelParams[part2].cppMultDevFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if cpp multiplier sigma is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;    
        }
    else if (whichParam == P_CPPEVENTS)
        {
        /* Check cpp events for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply (as parsimony is very
        silly and doesn't take this information into account) and cpp events are inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is cpp for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "Cpp"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "Cpp"))
            *isApplic2 = NO;
        
        /* Now, check that the cpp parameter is the same */
        if (IsModelSame (P_CPPRATE, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_CPPRATE][part1] != linkTable[P_CPPRATE][part2])
            isSame = NO;
    
        /* ... and that the psigamma parameter is the same */
        if (IsModelSame (P_CPPMULTDEV, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_CPPMULTDEV][part1] != linkTable[P_CPPRATE][part2])
            isSame = NO;
    
        /* Not same if branch lengths are not the same */
        if (IsModelSame(P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;

        /* Set isSame to NO if cpp events are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_TK02VAR)
        {
        /* Check prior for variance of rate autocorrelation for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
        silly and doesn't take this information into account) and ratevar is inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is tk02 for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "TK02"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "TK02"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on tk02 variance is the same. */
        if (!strcmp(modelParams[part1].tk02varPr,"Uniform") && !strcmp(modelParams[part2].tk02varPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].tk02varUni[0], modelParams[part2].tk02varUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].tk02varUni[1], modelParams[part2].tk02varUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].tk02varPr,"Exponential") && !strcmp(modelParams[part2].tk02varPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].tk02varExp, modelParams[part2].tk02varExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].tk02varPr,"Fixed") && !strcmp(modelParams[part2].tk02varPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].tk02varFix, modelParams[part2].tk02varFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if tk02 variance is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;    
        }
    else if (whichParam == P_TK02BRANCHRATES)
        {
        /* Check TK02 relaxed clock branch rates for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply (as parsimony is very
        silly and doesn't take this information into account) and tk02 relaxed clock branch rates are inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is tk02 for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "TK02"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "TK02"))
            *isApplic2 = NO;
        
        /* Now, check that the tk02 variance parameter is the same */
        if (IsModelSame (P_TK02VAR, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_TK02VAR][part1] != linkTable[P_TK02VAR][part2])
            isSame = NO;

        /* Not same if branch lengths are not the same */
        if (IsModelSame(P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;

        /* Set isSame to NO if tk02 branch rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_WNVAR)
        {
        /* Check prior for wn shape for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is wn for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "WN"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "WN"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on wn shape is the same. */
        if (!strcmp(modelParams[part1].wnvarPr,"Uniform") && !strcmp(modelParams[part2].wnvarPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].wnvarUni[0], modelParams[part2].wnvarUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].wnvarUni[1], modelParams[part2].wnvarUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].wnvarPr,"Exponential") && !strcmp(modelParams[part2].wnvarPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].wnvarExp, modelParams[part2].wnvarExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].wnvarPr,"Fixed") && !strcmp(modelParams[part2].wnvarPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].wnvarFix, modelParams[part2].wnvarFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if wn variance is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_WNBRANCHRATES)
        {
        /* Check WN relaxed clock branch rates for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is igr for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "WN"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "WN"))
            *isApplic2 = NO;
        
        /* Now, check that the wn shape parameter is the same */
        if (IsModelSame (P_WNVAR, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_WNVAR][part1] != linkTable[P_WNVAR][part2])
            isSame = NO;
    
        /* Not same if branch lengths are not the same */
        if (IsModelSame(P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;

        /* Set isSame to NO if igr branch rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_IGRVAR)
        {
        /* Check prior for igr shape for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
        silly and doesn't take this information into account) and igr shape is inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO; 
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is igr for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "IGR"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "IGR"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on igr shape is the same. */
        if (!strcmp(modelParams[part1].igrvarPr,"Uniform") && !strcmp(modelParams[part2].igrvarPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].igrvarUni[0], modelParams[part2].igrvarUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].igrvarUni[1], modelParams[part2].igrvarUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].igrvarPr,"Exponential") && !strcmp(modelParams[part2].igrvarPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].igrvarExp, modelParams[part2].igrvarExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].igrvarPr,"Fixed") && !strcmp(modelParams[part2].igrvarPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].igrvarFix, modelParams[part2].igrvarFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if igr variance is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;    
        }
    else if (whichParam == P_IGRBRANCHRATES)
        {
        /* Check IGR relaxed clock branch rates for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply (as parsimony is very
        silly and doesn't take this information into account) and igr relaxed clock branch rates are inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is igr for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "IGR"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "IGR"))
            *isApplic2 = NO;
        
        /* Now, check that the igr shape parameter is the same */
        if (IsModelSame (P_IGRVAR, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_IGRVAR][part1] != linkTable[P_IGRVAR][part2])
            isSame = NO;
    
        /* Not same if branch lengths are not the same */
        if (IsModelSame(P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;

        /* Set isSame to NO if igr branch rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_ILNVAR)
        {
        /* Check prior for iln variance for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
        silly and doesn't take this information into account) and iln variance is inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is iln for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "ILN"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "ILN"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on iln variance is the same. */
        if (!strcmp(modelParams[part1].ilnvarPr,"Uniform") && !strcmp(modelParams[part2].ilnvarPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].ilnvarUni[0], modelParams[part2].ilnvarUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].ilnvarUni[1], modelParams[part2].ilnvarUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].ilnvarPr,"Exponential") && !strcmp(modelParams[part2].ilnvarPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].ilnvarExp, modelParams[part2].ilnvarExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].ilnvarPr,"Fixed") && !strcmp(modelParams[part2].ilnvarPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].ilnvarFix, modelParams[part2].ilnvarFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if iln variance is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_ILNBRANCHRATES)
        {
        /* Check ILN relaxed clock branch rates for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply (as parsimony is very
        silly and doesn't take this information into account) and iln relaxed clock branch rates are inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Check that the clock rate prior is iln for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "ILN"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "ILN"))
            *isApplic2 = NO;
        
        /* Now, check that the iln shape parameter is the same */
        if (IsModelSame (P_ILNVAR, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_ILNVAR][part1] != linkTable[P_ILNVAR][part2])
            isSame = NO;
    
        /* Not same if branch lengths are not the same */
        if (IsModelSame(P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;

        /* Set isSame to NO if iln branch rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_MIXEDVAR)
        {
        /* Check prior for mixed var for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
         silly and doesn't take this information into account) and variance is inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        
        /* Check that the clock rate prior is mixed for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "Mixed"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "Mixed"))
            *isApplic2 = NO;
        
        /* Now, check that the prior on var mixed is the same. */
        if (!strcmp(modelParams[part1].mixedvarPr,"Uniform") && !strcmp(modelParams[part2].mixedvarPr,"Uniform"))
            {
            if (AreDoublesEqual (modelParams[part1].mixedvarUni[0], modelParams[part2].mixedvarUni[0], (YFlt) 0.00001) == NO)
                isSame = NO;
            if (AreDoublesEqual (modelParams[part1].mixedvarUni[1], modelParams[part2].mixedvarUni[1], (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].mixedvarPr,"Exponential") && !strcmp(modelParams[part2].mixedvarPr,"Exponential"))
            {
            if (AreDoublesEqual (modelParams[part1].mixedvarExp, modelParams[part2].mixedvarExp, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else if (!strcmp(modelParams[part1].mixedvarPr,"Fixed") && !strcmp(modelParams[part2].mixedvarPr,"Fixed"))
            {
            if (AreDoublesEqual (modelParams[part1].mixedvarFix, modelParams[part2].mixedvarFix, (YFlt) 0.00001) == NO)
                isSame = NO;
            }
        else
            isSame = NO;
        
        /* Check to see if variance is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_MIXEDBRCHRATES)
        {
        /* Check mixed relaxed clock branch rates for partitions 1 and 2. */
        
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply (as parsimony is very
         silly and doesn't take this information into account) and relaxed clock branch rates are inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO;
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;
        
        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;
        
        /* Check that the clock rate prior is mixed for both partitions */
        if (strcmp(modelParams[part1].clockVarPr, "Mixed"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].clockVarPr, "Mixed"))
            *isApplic2 = NO;
        
        /* Now, check that the var parameter is the same */
        if (IsModelSame (P_MIXEDVAR, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_MIXEDVAR][part1] != linkTable[P_MIXEDVAR][part2])
            isSame = NO;
        
        /* Not same if branch lengths are not the same */
        if (IsModelSame(P_BRLENS, part1, part2, &temp1, &temp2) == NO)
            isSame = NO;
        if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
            isSame = NO;
        
        /* Set isSame to NO if mixed branch rates are inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_CLOCKRATE)
        {
        /* Check base substitution rates of clock tree for partitions 1 and 2. */
    
        /* Check if the model is parsimony for either partition. If so, then branch lengths do not apply (as parsimony is very
        silly and doesn't take this information into account) and clock branch rates are inapplicable. */
        if (!strcmp(modelParams[part1].parsModel, "Yes"))
            *isApplic1 = NO; 
        if (!strcmp(modelParams[part2].parsModel, "Yes"))
            *isApplic2 = NO;

        /* Check that the branch length prior is clock for both partitions. */
        if (strcmp(modelParams[part1].brlensPr, "Clock"))
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].brlensPr, "Clock"))
            *isApplic2 = NO;

        /* Set isSame to NO if base substitution rate parameter is inapplicable for either partition. */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_SPECIESTREE)
        {
        /* Species tree; check that it is used in both partitions */
        if (strcmp(modelParams[part1].topologyPr, "Speciestree") != 0)
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].topologyPr, "Speciestree") != 0)
            *isApplic2 = NO;

        /* Not same if inapplicable to either partition */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else if (whichParam == P_GENETREERATE)
        {
        /* Gene tree rate; check that it is used in both partitions */
        if (strcmp(modelParams[part1].topologyPr, "Speciestree") != 0)
            *isApplic1 = NO;
        if (strcmp(modelParams[part2].topologyPr, "Speciestree") != 0)
            *isApplic2 = NO;

        /* Not same if inapplicable to either partition */
        if ((*isApplic1) == NO || (*isApplic2) == NO)
            isSame = NO;
        }
    else
        {
        YvyraPrint ("%s   Could not find parameter in IsModelSame\n", spacer);
        return (NO);
        }
    
    return (isSame);
}


int LargestMovableSubtree(Param *treeParam)
{
    int         i, j, k, a, nLongsNeeded, numPartitions, largestSubtree;
    BitsLong    **constraintPartition, *subtreePartition, *testPartition, *mask;
    ModelParams *mp;
    int         foundAllSpeciesPartition;

    mp = &modelParams[treeParam->relParts[0]];

    if (treeParam->paramType == P_SPECIESTREE)
        return numLocalTaxa;    /* no constraints allowed in species tree; set constraints in gene trees instead */
    
    /* This is difficult because we cannot rely on the tree being initialized.
       We need to retrieve the bitfields ourselves and figure out what they mean. */
    nLongsNeeded = ((numLocalTaxa - 1) / nBitsInALong) + 1;
    subtreePartition = (BitsLong *) SafeCalloc(3*nLongsNeeded, sizeof(BitsLong));
    constraintPartition = (BitsLong **) SafeCalloc (numDefinedConstraints+1, sizeof(BitsLong *));
    constraintPartition[0] = (BitsLong *) SafeCalloc ((numDefinedConstraints+1)*nLongsNeeded, sizeof(BitsLong));
    for (i=1; i<numDefinedConstraints+1; i++)
        constraintPartition[i] = constraintPartition[i-1] + nLongsNeeded;
    testPartition = subtreePartition + nLongsNeeded;
    mask = testPartition + nLongsNeeded;
    
    /* set mask (needed to take care of unused bits when flipping partitions) */
    for (i=0; i<numLocalTaxa; i++)
        SetBit (i, mask);
    
    /* retrieve partitions */
    numPartitions = 0;
    foundAllSpeciesPartition = NO;
    for (a=0; a<numDefinedConstraints; a++)
        {
        if (mp->activeConstraints[a] == NO || definedConstraintsType[a] != HARD)
            continue;
        
        /* set bits in partition under consideration */
        ClearBits(constraintPartition[numPartitions], nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            if (IsBitSet(i, definedConstraint[a]) == YES)
                SetBit(j, constraintPartition[numPartitions]);
            j++;
            }
        
        /* make sure outgroup is outside constrained partition (marked 0) */
        if (strcmp(mp->brlensPr,"Clock") != 0 && IsBitSet(localOutGroup, constraintPartition[numPartitions]) == YES)
            FlipBits(constraintPartition[numPartitions], nLongsNeeded, mask);
        
        /* skip partition if uninformative */
        k = NumBits(constraintPartition[numPartitions], nLongsNeeded);
        if (k == 0 || k == 1)
            continue;
        
        /* record if we hit an all-species partition */
        if (k == numLocalTaxa)
            foundAllSpeciesPartition = YES;
        
        numPartitions++;
        }
    
    /* Add all-species partition if not already present */
    if (foundAllSpeciesPartition == NO)
    {
        ClearBits(constraintPartition[numPartitions], nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
        {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            SetBit(j, constraintPartition[numPartitions]);
            j++;
        }
        numPartitions++;
    }
    
    /* Now we have all constraints. Calculate the movable subtree for each */
    largestSubtree = 0;
    for (i=0; i<numPartitions; i++)
        {
        CopyBits (subtreePartition, constraintPartition[i], nLongsNeeded);
        k = 0;
        for (j=0; j<numPartitions; j++)
            {
            if (j==i)
                continue;
            if (IsPartNested(constraintPartition[j], constraintPartition[i], nLongsNeeded))
                {
                k++;    /* add one for clade we are removing from subtreePartition */
                CopyBits (testPartition, constraintPartition[j], nLongsNeeded);
                FlipBits (testPartition, nLongsNeeded, mask);
                for (k=0; k<nLongsNeeded; k++)
                    subtreePartition[k] = subtreePartition[k] & testPartition[k];
                }
            }
        k += NumBits (subtreePartition, nLongsNeeded);  /* add remaining free tips in subtreePartition */
        /* add calculation root if an unrooted tree and we are dealing with the root partition */
        if (strcmp(mp->brlensPr,"Clock") != 0 && NumBits (constraintPartition[i], nLongsNeeded) == numLocalTaxa - 1)
            k++;
        if (k > largestSubtree)
            largestSubtree = k;
        }

    free(subtreePartition);
    free(constraintPartition[0]);
    free(constraintPartition);
   
    return largestSubtree;
}


int Link (void)
{
    int         i, j;
    
    for (j=0; j<NUM_LINKED; j++)
        {
        YvyraPrint ("%4d -- ", j+1);
        for (i=0; i<numCurrentDivisions; i++)
            YvyraPrint (" %2d", tempLinkUnlink[j][i]);
        YvyraPrint ("\n");
        }
        
    return (NO_ERROR);
}


int NumActiveParts (void)
{
    int     i, nApplied;
    
    nApplied = 0;
    for (i=0; i<numCurrentDivisions; i++)
        if (activeParts[i] == YES)
            nApplied++;

    return (nApplied);
}


int NumInformativeHardConstraints (ModelParams *mp)
{
    int             i, j, k, a, numInformativeHardConstraints, nLongsNeeded;
    BitsLong        *constraintPartition, *mask;
       
    numInformativeHardConstraints = 0;
    
    nLongsNeeded = ((numLocalTaxa - 1) / nBitsInALong) + 1;
    constraintPartition = (BitsLong *) SafeCalloc (2*nLongsNeeded, sizeof(BitsLong));
    if (!constraintPartition)
        {
            YvyraPrint ("%s   Problems allocating constraintPartition", spacer);
            return ERROR;
        }
    mask = constraintPartition + nLongsNeeded;

    /* set mask (needed to take care of unused bits when flipping partitions) */
    for (i=0; i<numLocalTaxa; i++)
        SetBit (i, mask);
        
    for (a=0; a<numDefinedConstraints; a++)
        {
        if (mp->activeConstraints[a] == NO || definedConstraintsType[a] != HARD)
            continue;
            
        /* set bits in partition to add */
        ClearBits(constraintPartition, nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
            {
                if (taxaInfo[i].isDeleted == YES)
                    continue;
                if (IsBitSet(i, definedConstraint[a]) == YES)
                    SetBit(j, constraintPartition);
                j++;
            }
            
        /* make sure outgroup is outside constrained partition (marked 0) */
        if (strcmp(mp->brlensPr,"Clock") != 0 && IsBitSet(localOutGroup, constraintPartition) == YES)
            FlipBits(constraintPartition, nLongsNeeded, mask);
            
        /* skip partition if uninformative */
        k = NumBits(constraintPartition, nLongsNeeded);
        if (k == 0 || k == 1)
            continue;

        numInformativeHardConstraints++;
        }
        
    return numInformativeHardConstraints;
}


int NumNonExcludedChar (void)
{
    int     i, n;
    
    /* count number of non-excluded characters */
    n = 0;
    for (i=0; i<numChar; i++)
        {
        if (charInfo[i].isExcluded == NO)
            {
            n++;
            }
        }
    
    return n;
}


int NumStates (int part)
{
    if (modelParams[part].dataType == STANDARD)
        return (MAX_STD_STATES);
    return (-1);
}


/*-----------------------------------------------------------------------
|
|   PrintCompMatrix: Print compressed matrix
|
------------------------------------------------------------------------*/
int PrintCompMatrix (void)
{
    int             i, j, k, c, d;
    ModelInfo       *m;
    ModelParams     *mp;
    char            tempName[100];
    char            (*whichChar)(int);

    if (!compMatrix)
        return ERROR;

    whichChar = &WhichStand;

    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        mp = &modelParams[d];

        whichChar = &WhichStand;

        YvyraPrint ("\nCompressed matrix for division %d\n\n", d+1);
        
        k = 66;
        if (mp->dataType == CONTINUOUS)
            k /= 4;

        for (c=m->compMatrixStart; c<m->compMatrixStop; c+=k)
            {
            for (i=0; i<numLocalTaxa; i++)
                {
                strcpy (tempName, localTaxonNames[i]);
                YvyraPrint ("%-10.10s   ", tempName);
                for (j=c; j<c+k; j++)
                    {
                    if (j >= m->compMatrixStop)
                        break;
                    if (mp->dataType == CONTINUOUS)
                        YvyraPrint ("%3d ", compMatrix[pos(i,j,compMatrixRowSize)]);
                    else
                        YvyraPrint ("%c", whichChar((int)compMatrix[pos(i,j,compMatrixRowSize)]));
                    }
                YvyraPrint ("\n");
                }
            YvyraPrint ("\nNo. sites    ");
            for (j=c; j<c+k; j++)
                {
                if (j >= m->compMatrixStop)
                    break;
                i = (int) numSitesOfPat[m->compCharStart + (0*numCompressedChars) + (j - m->compMatrixStart)/m->nCharsPerSite]; /* NOTE: We are printing the unadulterated site pat nums */
                if (i>9)
                    i = 'A' + i - 10;
                else
                    i = '0' + i;
                if (mp->dataType == CONTINUOUS)
                    YvyraPrint ("   %c ", i);
                else
                    {
                    if ((j-m->compMatrixStart) % m->nCharsPerSite == 0)
                        YvyraPrint ("%c", i);
                    else
                        YvyraPrint (" ");
                    }
                }
            YvyraPrint ("\nOrig. char   ");
            for (j=c; j<c+k; j++)
                {
                if (j >= m->compMatrixStop)
                    break;
                i = origChar[j];
                if (i>9)
                    i = '0' + (i % 10);
                else
                    i = '0' +i;
                if (mp->dataType == CONTINUOUS)
                    YvyraPrint ("   %c ", i);
                else
                    YvyraPrint ("%c", i);
                }

            if (mp->dataType == STANDARD && m->nStates != NULL)
                {
                YvyraPrint ("\nNo. states   ");
                for (j=c; j<c+k; j++)
                    {
                    if (j >= m->compMatrixStop)
                        break;
                    i = m->nStates[j-m->compCharStart];
                    YvyraPrint ("%d", i);
                    }
                YvyraPrint ("\nCharType     ");
                for (j=c; j<c+k; j++)
                    {
                    if (j >= m->compMatrixStop)
                        break;
                    i = m->cType[j-m->compMatrixStart];
                    if (i == ORD)
                        YvyraPrint ("%c", 'O');
                    else if (i == UNORD)
                        YvyraPrint ("%c", 'U');
                    else
                        YvyraPrint ("%c", 'I');
                    }
                YvyraPrint ("\ntiIndex      ");
                for (j=c; j<c+k; j++)
                    {
                    if (j >= m->compMatrixStop)
                        break;
                    i = m->tiIndex[j-m->compCharStart];
                    YvyraPrint ("%d", i % 10);
                    }
                YvyraPrint ("\nbsIndex      ");
                for (j=c; j<c+k; j++)
                    {
                    if (j >= m->compMatrixStop)
                        break;
                    i = m->bsIndex[j-m->compCharStart];
                    YvyraPrint ("%d", i % 10);
                    }
                }
            YvyraPrint ("\n\n");
            }
        YvyraPrint ("Press return to continue\n");
        getchar();
        }   /* next division */

    return NO_ERROR;
}


/*----------------------------------------------------------------------
|
|   PrintMatrix: Print data matrix
|
|
------------------------------------------------------------------------*/
int PrintMatrix (void)
{
    int             i, j=0, c, printWidth, nextColumn;

    if (!matrix)
        return ERROR;
    
    YvyraPrint ("\nData matrix\n\n");
    
    printWidth = 79;

    for (c=0; c<numChar; c=j)
        {
        for (i=0; i<numTaxa; i++)
            {
            YvyraPrint ("%-10.10s   ", taxaNames[i]);
            j = c;
            for (nextColumn=13; nextColumn < printWidth; nextColumn++)
                {
                if (j >= numChar)
                    break;
                if (charInfo[j].charType == CONTINUOUS && nextColumn < printWidth - 3)
                    break;
                if (charInfo[j].charType == CONTINUOUS)
                    {   
                    YvyraPrint ("%3d ", matrix[pos(i,j,numChar)]);
                    nextColumn += 3;
                    }
                else if (charInfo[j].charType == STANDARD)
                    YvyraPrint ("%c", WhichStand(matrix[pos(i,j,numChar)]));
                j++;
                }
            YvyraPrint ("\n");
            }
        YvyraPrint ("\n");
        }

    return NO_ERROR;
}


/*--------------------------------------------------------------
|
|   ProcessStdChars: process standard characters
|
---------------------------------------------------------------*/
int ProcessStdChars (RandLong *seed)
{
    int             c, d, i, j, k, n, ts, index, numStandardChars, origCharPos, *bsIndex;
    char            piHeader[30];
    ModelInfo       *m;
    ModelParams     *mp=NULL;
    Param           *p;

    /* set character type, no. states, ti index and bs index for standard characters */
    /* first calculate how many standard characters we have */
    numStandardChars = 0;
    for (d=0; d<numCurrentDivisions; d++)
        {
        mp = &modelParams[d];
        m = &modelSettings[d];

        if (mp->dataType != STANDARD)
            continue;

        numStandardChars += m->numChars;
        }
    
    /* return if there are no standard characters */
    if (numStandardChars == 0)
        return (NO_ERROR);

    /* we are still here so we have standard characters and need to deal with them */
    
    /* first allocate space for stdType, stateSize, tiIndex, bsIndex */
    if (memAllocs[ALLOC_STDTYPE] == YES)
        {
        free (stdType);
        stdType = NULL;
        memAllocs[ALLOC_STDTYPE] = NO;
        }
    stdType = (int *)SafeCalloc(4 * (size_t)numStandardChars, sizeof(int));
    if (!stdType)
        {
        YvyraPrint ("%s   Problem allocating stdType (%d ints)\n", 4 * numStandardChars);
        return ERROR;
        }
    memAllocs[ALLOC_STDTYPE] = YES;
    stateSize = stdType + numStandardChars;
    tiIndex = stateSize + numStandardChars;
    bsIndex = tiIndex + numStandardChars;

    /* then fill in stdType and stateSize, set pointers */
    /* also fill in isTiNeeded for each division and tiIndex for each character */
    for (d=j=0; d<numCurrentDivisions; d++)
        {
        mp = &modelParams[d];
        m = &modelSettings[d];
        
        if (mp->dataType != STANDARD)
            continue;

        m->cType = stdType + j;
        m->nStates = stateSize + j;
        m->userTypeIdx = (int *) SafeCalloc (m->numChars, sizeof(int));
        if (m->userTypeIdx == NULL)
            return ERROR;
        for (c=0; c<m->numChars; c++)
            m->userTypeIdx[c] = -1;
        m->hasUserType = NO;
        m->tiIndex = tiIndex + j;
        m->bsIndex = bsIndex + j;

        m->cijkLength = 0;
        m->nCijkParts = 0;
        for (c=0; c<m->numChars; c++)
            {
            if (origChar[c+m->compMatrixStart] < 0)
                {
                /* this is a dummy character — find its correction group */
                int g, dummyIdx, found = NO;
                dummyIdx = c + m->compCharStart;
                for (g = 0; g < m->numCorrGroups; g++)
                    {
                    int gStart = m->corrGroupDummyStart[g] - m->compCharStart;
                    int gEnd = gStart + m->corrGroupNStates[g];
                    if (c >= gStart && c < gEnd)
                        {
                        m->cType[c] = m->corrGroupCType[g];
                        m->nStates[c] = m->corrGroupNStates[g];
                        if (m->corrGroupCType[g] == USERTYPE)
                            {
                            m->userTypeIdx[c] = m->corrGroupUserTypeIdx[g];
                            m->hasUserType = YES;
                            }
                        found = YES;
                        break;
                        }
                    }
                if (found == NO)
                    {
                    /* Singleton or other dummy — default to UNORD binary */
                    m->cType[c] = UNORD;
                    m->nStates[c] = 2;
                    }
                }
            else
                {
                /* this is an ordinary character */
                m->cType[c] = charInfo[origChar[c + m->compMatrixStart]].ctype;
                m->nStates[c] = charInfo[origChar[c + m->compMatrixStart]].numStates;
                if (m->cType[c] == USERTYPE)
                    {
                    m->userTypeIdx[c] = charInfo[origChar[c + m->compMatrixStart]].userTypeIndex;
                    m->hasUserType = YES;
                    }
                }
            
            /* check ctype settings */
            if (m->nStates[c] < 2)
                {
                YvyraPrint ("%s   WARNING: Compressed character %d (original character %d) of division %d has \n", spacer, c+m->compCharStart,origChar[c+m->compCharStart]+1, d+1);
                YvyraPrint ("%s            less than two observed states; it will be assumed to have two states.\n", spacer);
                m->nStates[c] = 2;
                }
            if (m->nStates[c] > 6 && m->cType[c] != UNORD && m->cType[c] != USERTYPE)
                {
                YvyraPrint ("%s   Only unordered model supported for characters with more than 6 states\n", spacer);
                return ERROR;
                }
            if (m->nStates[c] == 2 && m->cType[c] == ORD)
                m->cType[c] = UNORD;
            if (m->cType[c] == IRREV)
                {
                YvyraPrint ("%s   Irreversible model not yet supported\n", spacer);
                return ERROR;
                }
            
            /* find max number of states */
            if (m->nStates[c] > m->numModelStates)
                m->numModelStates = m->nStates[c];

            /* update Cijk info */
            if (m->cType[c] == USERTYPE)
                {
                /* Usertype characters always need eigendecomposition */
                ts = m->nStates[c];
                m->cijkLength += (ts * ts * ts) + (2 * ts);
                m->nCijkParts++;
                }
            else if (strcmp(mp->symPiPr,"Fixed") != 0 || AreDoublesEqual(mp->symBetaFix, -1.0, 0.00001) == NO)
                {
                /* Asymmetry between stationary state frequencies -- we need one cijk and eigenvalue
                    set for each multistate character */
                if (m->nStates[c] > 2 && (m->cType[c] == UNORD || m->cType[c] == ORD))
                    {
                    ts = m->nStates[c];
                    m->cijkLength += (ts * ts * ts) + (2 * ts);
                    m->nCijkParts++;
                    }
                }

            /* set the ti probs needed */
            if (m->stateFreq->nValues == 0 || m->nStates[c] == 2)
                {
                if (m->cType[c] == UNORD)
                    m->isTiNeeded[m->nStates[c]-2] = YES;
                if (m->cType[c] == ORD)
                    m->isTiNeeded[m->nStates[c]+MAX_STD_STATES-4] = YES;
                if (m->cType[c] == IRREV)
                    m->isTiNeeded[m->nStates[c]+2*MAX_STD_STATES-5] = YES;
                }
            }

        /* set ti index for each compressed character first         */
        /* set bs index later (below)                               */

        /* set base index, valid for binary chars */
        for (c=0; c<m->numChars; c++)
            m->tiIndex[c] = 0;

        /* first adjust for unordered characters */
        for (k = 0; k < MAX_STD_STATES-1; k++)
            {
            if (m->isTiNeeded[k] == NO)
                continue;

            for (c=0; c<m->numChars; c++)
                {
                if (m->cType[c] != UNORD || m->nStates[c] > k + 2)
                    {
                    m->tiIndex[c] += (k + 2) * (k + 2) * m->numRateCats;
                    }
                }
            }

        /* second for ordered characters */
        for (k = MAX_STD_STATES-1; k < 2*MAX_STD_STATES-3; k++)
            {
            if (m->isTiNeeded [k] == NO)
                continue;

            for (c=0; c<m->numChars; c++)
                {
                if (m->cType[c] == ORD && m->nStates[c] > k-MAX_STD_STATES+4)
                    {
                    m->tiIndex[c] += (k-MAX_STD_STATES+4) * (k-MAX_STD_STATES+4) * m->numRateCats;
                    }
                }
            }

        /* add USERTYPE characters to tiIndex after the ORD block */
        {
        int tiOffset = 0;
        /* calculate total space used by UNORD and ORD chars */
        for (k = 0; k < MAX_STD_STATES-1; k++)
            if (m->isTiNeeded[k] == YES)
                tiOffset += (k + 2) * (k + 2) * m->numRateCats;
        for (k = MAX_STD_STATES-1; k < 2*MAX_STD_STATES-3; k++)
            if (m->isTiNeeded[k] == YES)
                tiOffset += (k-MAX_STD_STATES+4) * (k-MAX_STD_STATES+4) * m->numRateCats;
        for (c=0; c<m->numChars; c++)
            {
            if (m->cType[c] == USERTYPE)
                {
                m->tiIndex[c] = tiOffset;
                tiOffset += m->nStates[c] * m->nStates[c] * m->numRateCats;
                }
            }
        }

        /* finally take beta cats into account in tiIndex        */
        /* the beta cats will only be used for binary characters */
        /* multistate characters get their ti indices reset here */
        if (m->numBetaCats > 1 && m->isTiNeeded[0] == YES)
            {
            k = 4 * m->numBetaCats * m->numRateCats;   /* offset for binary character ti probs */
            for (c=0; c<m->numChars; c++)
                {
                if (m->nStates[c] > 2)
                    {
                    m->tiIndex[c] = k;
                    k += m->nStates[c] * m->nStates[c] * m->numRateCats;
                    }
                }
            }
        j += m->numChars;
        }
    
    /* deal with bsIndex */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType != P_PI || modelParams[p->relParts[0]].dataType != STANDARD)
            continue;
        p->nSympi = 0;
        p->hasBinaryStd = NO;
        for (i=0; i<p->nRelParts; i++)
            if (modelSettings[p->relParts[i]].isTiNeeded[0] == YES)
                break;
        if (i < p->nRelParts)
            p->hasBinaryStd = YES;
        if (p->paramId == SYMPI_EQUAL)
            {
            /* calculate the number of state frequencies needed */
            /* also set bsIndex appropriately                   */
            for (n = index = 0; n < MAX_STD_STATES-1; n++)
                {
                for (i=0; i<p->nRelParts; i++)
                    if (modelSettings[p->relParts[i]].isTiNeeded[n] == YES)
                        break;
                if (i < p->nRelParts)
                    {
                    for (i=0; i<p->nRelParts; i++)
                        {
                        m = &modelSettings[p->relParts[i]];
                        for (c=0; c<m->numChars; c++)
                            {
                            if (m->cType[c] != UNORD || m->nStates[c] > n + 2)
                                {
                                m->bsIndex[c] += (n + 2);
                                }
                            }
                        }
                    index += (n + 2);
                    }
                }
            for (n = MAX_STD_STATES-1; n < 2*MAX_STD_STATES-3; n++)
                {
                for (i=0; i<p->nRelParts; i++)
                    if (modelSettings[p->relParts[i]].isTiNeeded[n] == YES)
                        break;
                if (i < p->nRelParts)
                    {
                    for (i=0; i<p->nRelParts; i++)
                        {
                        m = &modelSettings[p->relParts[i]];
                        for (c=0; c<m->numChars; c++)
                            {
                            if (m->cType[c] == ORD && m->nStates[c] > n-MAX_STD_STATES+4)
                                {
                                m->bsIndex[c] += n-MAX_STD_STATES+4;
                                }
                            }
                        }
                    index += n-MAX_STD_STATES+4;
                    }
                }
            /* For USERTYPE characters with equal base freqs, add them to nSympi
               so they go through eigendecomposition path */
            for (i=0; i<p->nRelParts; i++)
                {
                m = &modelSettings[p->relParts[i]];
                for (c=0; c<m->numChars; c++)
                    {
                    if (m->cType[c] == USERTYPE)
                        {
                        m->bsIndex[c] = index;
                        index += m->nStates[c];
                        p->nSympi++;
                        }
                    }
                }
            p->nStdStateFreqs = index;
            }
        else
            {
            /* if not equal we need space for beta category frequencies */
            index = 0;
            if (p->hasBinaryStd == YES)
                index += (2 * modelSettings[p->relParts[0]].numBetaCats);
            /* as well as one set of frequencies for each multistate character */
            for (i=0; i<p->nRelParts; i++)
                {
                m = &modelSettings[p->relParts[i]];
                for (c=0; c<m->numChars; c++)
                    {
                    if (((m->nStates[c] > 2 && (m->cType[c] == UNORD || m->cType[c] == ORD)) || m->cType[c] == USERTYPE))
                        {
                        m->bsIndex[c] = index;
                        index += m->nStates[c];
                        p->nSympi++;
                        }
                    }
                }
            p->nStdStateFreqs = index;
            }
        }
    
    /* allocate space for sympiIndex, stdStateFreqs; then fill */
    /* first count number of sympis needed */
    for (k=n=i=0; k<numParams; k++)
        {
        p = &params[k];
        n += p->nSympi;     /* nsympi calculated above */
        }
    
    /* then allocate and fill in */
    if (n > 0)
        {
        if (memAllocs[ALLOC_SYMPIINDEX] == YES)
            {
            sympiIndex = (int *) SafeRealloc ((void *) sympiIndex, 3*n * sizeof (int));
            for (i=0; i<3*n; i++)
                sympiIndex[i] = 0;
            }
        else
            sympiIndex = (int *) SafeCalloc (3*n, sizeof (int));
        if (!sympiIndex)
            {
            YvyraPrint ("%s   Problem allocating sympiIndex\n", spacer);
            return (ERROR);
            }
        else
            memAllocs[ALLOC_SYMPIINDEX] = YES;
        
        /* set up sympi pointers and fill sympiIndex */
        for (k=i=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->nSympi > 0)
                {
                if (p->nValues > 0)
                    p->printParam = YES;    /* print even if fixed alpha_symdir because we sample state freqs and do not integrate them out */
                index = 0;
                p->sympiBsIndex = sympiIndex + i;
                p->sympinStates = sympiIndex + i + n;
                p->sympiCType = sympiIndex + i + (2 * n);
                p->sympiUserTypeIndex = (int *) SafeCalloc (p->nSympi, sizeof(int));
                for (j=0; j<p->nSympi; j++)
                    p->sympiUserTypeIndex[j] = -1;
                for (j=0; j<p->nRelParts; j++)
                    {
                    m = &modelSettings[p->relParts[j]];
                    for (c=0; c<m->numChars; c++)
                        {
                        if (((m->nStates[c] > 2 && (m->cType[c] == UNORD || m->cType[c] == ORD)) || m->cType[c] == USERTYPE))
                            {
                            p->sympinStates[index] = m->nStates[c];
                            p->sympiBsIndex[index] = m->bsIndex[c];
                            p->sympiCType[index] = m->cType[c];
                            if (m->cType[c] == USERTYPE)
                                p->sympiUserTypeIndex[index] = m->userTypeIdx[c];
                            origCharPos = origChar[m->compCharStart + c];
                            for (ts=0; ts<m->nStates[c]; ts++)
                                {
                                sprintf (piHeader, "\tpi_%d(%d)", origCharPos+1, ts);
                                SafeStrcat(&p->paramHeader, piHeader);
                                }
                            index++;
                            }
                        }
                    }
                assert (index == p->nSympi);
                i += p->nSympi;
                }
            }
        assert (i == n);
        }
    
    /* count space needed for state frequencies */
    for (k=n=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType != P_PI || modelParams[p->relParts[0]].dataType != STANDARD)
            continue;
        n += p->nStdStateFreqs;
        }
    
    stdStateFreqsRowSize = n;
    
    /* allocate space */
    if (memAllocs[ALLOC_STDSTATEFREQS] == YES)
        {
        free (stdStateFreqs);
        stdStateFreqs = NULL;
        memAllocs[ALLOC_STDSTATEFREQS] = NO;
        }
    stdStateFreqs = (YFlt *) SafeCalloc (n * 2 * numGlobalChains, sizeof (YFlt));
    if (!stdStateFreqs)
        {
        YvyraPrint ("%s   Problem allocating stdStateFreqs in ProcessStdChars\n", spacer);
        return (ERROR);
        }
    else
        memAllocs[ALLOC_STDSTATEFREQS] = YES;
    
    /* set pointers */
    for (k=n=0; k<numParams; k++)
        {
        p = &params[k];
        if (p->paramType != P_PI || modelParams[p->relParts[0]].dataType != STANDARD)
            continue;
        p->stdStateFreqs = stdStateFreqs + n;
        n += p->nStdStateFreqs;
        }
    
    FillStdStateFreqs (0 , numGlobalChains, seed);

    return (NO_ERROR);
}
int SetLocalTaxa (void)
{
    int         i, j;
    
    /* free memory if allocated */
    if (memAllocs[ALLOC_LOCTAXANAMES] == YES)
        {
        free (localTaxonNames);
        localTaxonNames = NULL;
        memAllocs[ALLOC_LOCTAXANAMES] = NO;
        }
    if (memAllocs[ALLOC_LOCALTAXONCALIBRATION] == YES)
        {
        free (localTaxonCalibration);
        localTaxonCalibration = NULL;
        memAllocs[ALLOC_LOCALTAXONCALIBRATION] = NO;
        }
    
    /* count number of non-excluded taxa */
    numLocalTaxa = 0;
    for (i=0; i<numTaxa; i++)
        {
        if (taxaInfo[i].isDeleted == NO)
            numLocalTaxa++;
        }
        
    /* allocate memory */
    localTaxonNames = (char **)SafeCalloc((size_t)numLocalTaxa, sizeof(char *));
    if (!localTaxonNames)
        return (ERROR);
    memAllocs[ALLOC_LOCTAXANAMES] = YES;

    localTaxonCalibration = (Calibration **)SafeCalloc((size_t)numLocalTaxa, sizeof(Calibration *));
    if (!localTaxonCalibration)
        return (ERROR);
    memAllocs[ALLOC_LOCALTAXONCALIBRATION] = YES;
        
    /* point to names and calibrations of non-excluded taxa */
    localOutGroup = 0;
    for (i=j=0; i<numTaxa; i++)
        {
        if (taxaInfo[i].isDeleted == NO)
            {
            localTaxonNames[j] = taxaNames[i];
            localTaxonCalibration[j] = &tipCalibration[i];
            if (i == outGroupNum)
                localOutGroup = j;
            j++;
            }
        }

#   if 0
    /* show non-excluded taxa */
    for (i=0; i<numLocalTaxa; i++)
        YvyraPrint ("%s   %4d %s\n", spacer, i+1, localTaxonNames[i]);
#   endif
        
    return (NO_ERROR);
}


/*----------------------------------------------------------------------------
|
|   SetModelDefaults: This function will set up model defaults in modelParams.
|       It will also initialize parameters and moves by calling SetUpAnalysis.
|
-----------------------------------------------------------------------------*/
int SetModelDefaults (void)
{
    int         j;

    YvyraPrint ("%s   Setting model defaults\n", spacer);
    YvyraPrint ("%s   Seed (for generating default start values) = %d\n", spacer, globalSeed);

    if (InitializeLinks () == ERROR)
        {
        YvyraPrint ("%s   Problem initializing link table\n", spacer);
        return (ERROR);
        }

    /* Check that models are allocated */
    if (memAllocs[ALLOC_MODEL] == NO)
        {
        YvyraPrint ("%s   Model not allocated in SetModelDefaults\n", spacer);
        return (ERROR);
        }

    /* model parameters */
    for (j=0; j<numCurrentDivisions; j++)
        {
        modelParams[j] = defaultModel;                      /* start with default settings */
        
        modelParams[j].dataType = DataType (j);             /* data type for partition                      */

        if (modelParams[j].dataType == STANDARD)
            {   /* set default ascertainment bias for partition */
            modelParams[j].coding = VARIABLE;
            strcpy(modelParams[j].codingString, "Variable"); 
            }
        else if (0) /* RESTRICTION removed */
            {
            modelParams[j].coding = NOABSENCESITES;
            strcpy(modelParams[j].codingString, "Noabsencesites");   
            }
        else
            {
            modelParams[j].coding = ALL;
            strcpy(modelParams[j].codingString, "All");
            }

        modelParams[j].nStates = NumStates (j);             /* number of states for partition             */

        if (numDefinedConstraints > 0)
            modelParams[j].activeConstraints = (int *) SafeCalloc((size_t)(numDefinedConstraints), sizeof(int));  /* allocate space for active constraints (yes/no) */
        }

    return (NO_ERROR);
}


/*----------------------------------------------------------------------------
|
|   SetModelInfo: This function will set up model info using model
|       params
|
-----------------------------------------------------------------------------*/
int SetModelInfo (void)
{
    int             i, j, chn, ts;
    ModelParams     *mp;
    ModelInfo       *m;
    
    /* wipe all model settings */
    inferSiteRates = NO;
    inferSiteLikes = NO;
    inferAncStates = NO;
    inferSiteOmegas = NO;

    for (i=0; i<numCurrentDivisions; i++)
        {
        m = &modelSettings[i];

        /* make certain that we set this intentionally to "NO" so we 
           calculate cijk information and calculate cond likes when needed */
        m->upDateCijk = YES;
        m->upDateCl = YES;
        m->upDateAll = YES;

        /* make certain that we start with a parsimony branch length of zero */
        for (j=0; j<MAX_CHAINS; j++)
            m->parsTreeLength[j*2] = m->parsTreeLength[j*2+1] = 0.0;

        m->tRatio = NULL;
        m->revMat = NULL;
        m->omega = NULL;
        m->stateFreq = NULL;
        m->mixtureRates = NULL;
        m->shape = NULL;
        m->pInvar = NULL;
        m->correlation = NULL;
        m->switchRates = NULL;
        m->rateMult = NULL;
        m->topology = NULL;
        m->brlens = NULL;
        m->speciationRates = NULL;
        m->extinctionRates = NULL;
        m->fossilizationRates = NULL;
        m->popSize = NULL;
        m->aaModel = NULL;
        m->cppRate = NULL;
        m->cppEvents = NULL;
        m->cppMultDev = NULL;
        m->tk02var = NULL;
        m->tk02BranchRates = NULL;
        m->wnvar = NULL;
        m->wnBranchRates = NULL;
        m->igrvar = NULL;
        m->igrBranchRates = NULL;
        m->ilnvar = NULL;
        m->ilnBranchRates = NULL;
        m->mixedvar = NULL;
        m->mixedBrchRates = NULL;
        m->clockRate = NULL;

        m->CondLikeDown = NULL;
        m->CondLikeRoot = NULL;
        m->CondLikeScaler = NULL;
        m->Likelihood = NULL;
        m->TiProbs = NULL;

        m->CondLikeUp = NULL;
        m->StateCode = NULL;
        m->PrintAncStates = NULL;
        m->PrintSiteRates = NULL;
        
        m->printPosSel = NO;
        m->printAncStates = NO;
        m->printSiteRates = NO;

        m->nStates = NULL;
        m->bsIndex = NULL;
        m->cType = NULL;
        m->tiIndex = NULL;

        m->gibbsGamma = NO;
        m->gibbsFreq = 0;

        m->parsimonyBasedMove = NO;

        /* likelihood calculator flags */
        m->useVec = VEC_NONE;                 /* use SIMD code for this partition?            */

#if defined (SSE_ENABLED)
        m->numVecChars = 0;
        m->numFloatsPerVec = 0;
        m->lnL_Vec = NULL;
        m->lnLI_Vec = NULL;
        m->clP_SSE = NULL;
#if defined (AVX_ENABLED)
        m->clP_AVX = NULL;
#endif
#endif
            
        /* set all memory pointers to NULL */
        m->parsSets = NULL;
        m->numParsSets = 0;
        m->parsNodeLens = NULL;
        m->numParsNodeLens = 0;

        m->condLikes = NULL;
        m->tiProbs = NULL;
        m->scalers = NULL;
        m->numCondLikes = 0;
        m->numTiProbs = 0;
        m->numScalers = 0;

        m->condLikeIndex = NULL;
        m->condLikeScratchIndex = NULL;
        m->tiProbsIndex = NULL;
        m->tiProbsScratchIndex = NULL;
        m->nodeScalerIndex = NULL;
        m->nodeScalerScratchIndex = NULL;
        m->siteScalerIndex = NULL;
        m->siteScalerScratchIndex = -1;

        m->cijks = NULL;
        m->nCijkParts = 0;
        m->cijkIndex = NULL;
        m->cijkScratchIndex = -1;
        }

    /* set state of all chains to zero */
    for (chn=0; chn<numGlobalChains; chn++)
        state[chn] = 0;

    /* fill in modelSettings info with some basic model characteristics */
    for (i=0; i<numCurrentDivisions; i++)
        {
        mp = &modelParams[i];
        m = &modelSettings[i];
        
        if (!strcmp(mp->nucModel,"Protein") && (mp->dataType == DNA || mp->dataType == RNA))
            m->dataType = PROTEIN;
        else
            m->dataType = mp->dataType;

        /* parsimony model? */
        if (!strcmp(mp->parsModel, "Yes"))
            m->parsModelId = YES;
        else
            m->parsModelId = NO;

        /* number of rate categories */
        if (activeParams[P_SHAPE][i] > 0)
            {
            if (!strcmp(mp->ratesModel, "Lnorm"))
                m->numRateCats = mp->numLnormCats;
            else
                m->numRateCats = mp->numGammaCats;
            }
        else if (activeParams[P_MIXTURE_RATES][i] > 0)
            m->numRateCats = mp->numMixtCats;
        else
            m->numRateCats = 1;

        /* number of beta categories */
        if (mp->dataType == STANDARD && !(AreDoublesEqual(mp->symBetaFix, -1.0, 0.00001) == YES && !strcmp(mp->symPiPr,"Fixed")))
            m->numBetaCats = mp->numBetaCats;
        else
            m->numBetaCats = 1;

        /* number of omega categories */
        if ((mp->dataType == DNA || mp->dataType == RNA) && (!strcmp(mp->omegaVar, "Ny98") || !strcmp(mp->omegaVar, "M3")) && !strcmp(mp->nucModel, "Codon"))
            {
            m->numOmegaCats = 3;
            m->numRateCats = 1; /* if we are here, then we cannot have gamma or beta variation */
            m->numBetaCats = 1;
            }
        else if ((mp->dataType == DNA || mp->dataType == RNA) && !strcmp(mp->omegaVar, "M10") && !strcmp(mp->nucModel, "Codon"))
            {
            m->numOmegaCats = mp->numM10BetaCats + mp->numM10GammaCats;
            m->numRateCats = 1; /* if we are here, then we cannot have gamma or beta variation */
            m->numBetaCats = 1;
            }
        else
            m->numOmegaCats = 1;

        /* number of transition matrices depends on numGammaCats, numBetaCats, and numOmegaCats */
        m->numTiCats = m->numRateCats * m->numBetaCats * m->numOmegaCats;

        /* TODO: check that numStates and numModelStates are set
            appropriately for codon and doublet models */

        /* number of observable states */
        if (m->dataType == STANDARD)
            m->numStates = 0;   /* zero, meaning variable */
        else if (!strcmp(mp->nucModel,"Protein") && (mp->dataType == DNA || mp->dataType == RNA))
            m->numStates = 20;
        else
            m->numStates = mp->nStates;
        
        /* number of model states including hidden ones */
        if ((mp->dataType == DNA || mp->dataType == RNA) && !strcmp (mp->covarionModel, "Yes") && !strcmp(mp->nucModel, "4by4"))
            m->numModelStates = mp->nStates * 2;
        else if (m->dataType == PROTEIN && !strcmp (mp->covarionModel, "Yes"))
            m->numModelStates = mp->nStates * 2;
        else if ((mp->dataType == DNA || mp->dataType == RNA) && !strcmp(mp->nucModel,"Protein") && !strcmp (mp->covarionModel, "Yes"))
            m->numModelStates = 20 * 2;
        else if (mp->dataType == CONTINUOUS)
            m->numModelStates = 0;
        else if (mp->dataType == STANDARD)
            {
            /* use max possible for now; we don't know what chars will be included */
            m->numModelStates = MAX_STD_STATES;
            }
        else
            m->numModelStates = m->numStates;
            
        /* Fill in some information for calculating cijk. We will use m->cijkLength to 
           figure out if we need to diagonalize Q to calculate transition probabilities.
           If cijkLength = 0, then we won't bother. We use cijkLength later in this function. */
        m->cijkLength = 0;
        m->nCijkParts = 1;
        if (m->dataType == PROTEIN)
            {
            ts = m->numModelStates;
            m->cijkLength = (ts * ts * ts) + (2 * ts);
            if (!strcmp (mp->covarionModel, "Yes"))
                {
                m->cijkLength *= m->numRateCats;
                m->nCijkParts = m->numRateCats;
                }
            }
        else if (m->dataType == STANDARD)
            {
            /* set to 0 for now, update in ProcessStdChars */
            m->nCijkParts = 0;
            }
        else if (m->dataType == DNA || m->dataType == RNA)
            {
            if (m->nucModelId == NUCMODEL_4BY4)
                {
                if (!strcmp (mp->covarionModel, "No") && m->nst != 6 && m->nst != 203)
                    m->cijkLength = 0;
                else
                    {
                    ts = m->numModelStates;
                    m->cijkLength = (ts * ts * ts) + (2 * ts);
                    }
                if (!strcmp (mp->covarionModel, "Yes"))
                    {
                    m->cijkLength *= m->numRateCats;
                    m->nCijkParts = m->numRateCats;
                    }
                }
            else if (m->nucModelId == NUCMODEL_DOUBLET)
                {
                ts = m->numModelStates;
                m->cijkLength = (ts * ts * ts) + (2 * ts);
                }
            else if (m->nucModelId == NUCMODEL_CODON)
                {
                ts = m->numModelStates;
                m->cijkLength = (ts * ts * ts) + (2 * ts);
                m->cijkLength *= m->numOmegaCats;
                m->nCijkParts = m->numOmegaCats;
                }
            else
                {
                YvyraPrint ("%s   BUG: Inconsistent model. Please report this as a bug.\n");
                return ERROR;
                }
            }

        /* check if we should calculate ancestral states */
        if (!strcmp(mp->inferAncStates,"Yes"))
            {
            if (m->dataType == PROTEIN && !strcmp(mp->covarionModel, "No"))
                m->printAncStates = YES;
            else if (m->dataType == DNA || m->dataType == RNA)
                {
                if (!strcmp(mp->nucModel,"4by4") && !strcmp(mp->covarionModel, "No"))
                    m->printAncStates = YES;
                if (!strcmp(mp->nucModel,"Doublet"))
                    m->printAncStates = YES;
                if (!strcmp(mp->nucModel,"Codon") && !strcmp(mp->omegaVar,"Equal"))
                    m->printAncStates = YES;
                }
            else if (m->dataType == STANDARD || m->dataType == RESTRICTION)
                m->printAncStates = YES;
            if (m->printAncStates == YES)
                inferAncStates = YES;
            else
                YvyraPrint ("%s   Print out of ancestral states is not applicable for partition %d.\n",spacer,i);
            }

        /* check if we should calculate site rates */
        if (!strcmp(mp->inferSiteRates,"Yes"))
            {
            if (m->numRateCats > 1)
                {
                m->printSiteRates = YES;
                inferSiteRates = YES;
                }
            }

        /* check if we should output per-site log-likelihoods */
        if (!strcmp(mp->inferSiteLikes,"Yes"))
            {
            m->printSiteLikes = YES;
            inferSiteLikes = YES;
            }

        /* check if we should calculate positive selection */
        if (!strcmp(mp->inferPosSel, "Yes"))
            {
            if (m->numOmegaCats > 1)
                {
                m->printPosSel = YES;
                inferPosSel = YES;
                }
            }

        /* check if we should calculate site omegas */
        if (!strcmp(mp->inferSiteOmegas, "Yes"))
            {
            if (m->numOmegaCats > 1)
                {
                m->printSiteOmegas = YES;
                inferSiteOmegas = YES;
                }
            }
        /* check if we should use gibbs sampling of gamma (don't use it with pinvar or invgamma) */
        if (!strcmp(mp->useGibbs,"Yes") && (m->numRateCats > 1))
            {
            if (m->dataType == DNA || m->dataType == RNA || m->dataType == PROTEIN)
                {
                if (activeParams[P_CORREL][i] <= 0 && m->printSiteRates == NO && activeParams[P_PINVAR][i] <= 0)
                    {
                    m->gibbsGamma = YES;
                    m->gibbsFreq = mp->gibbsFreq;
                    }
                }
            }
        }

    return (NO_ERROR);  
}


/*-----------------------------------------------------------------
|
|   SetModelParams: Set up parameter structs for all model
|       parameters, including trees
|
|----------------------------------------------------------------*/
int SetModelParams (void)
{
    int             c, i, j, k, n, n1, n2, *isPartTouched, numRelParts, nRelParts, areAllPartsParsimony,
                    nClockBrlens, nRelaxedBrlens, nCalibratedBrlens;
    char            tempCodon[15], tempMult[20], *tempStr, temp[30];
    char static     *partString=NULL; /* mad static to avoid possible memory leak on return ERROR if it would be introduced later */
    Param           *p;
    ModelParams     *mp;
    ModelInfo       *m;
    int             tempStrSize = 300;

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    isPartTouched = (int *) SafeMalloc ((size_t)numCurrentDivisions * sizeof (int));
    if (!tempStr || !isPartTouched)
        {
        YvyraPrint ("%s   Problem allocating tempString (%d) or isPartTouched (%d)\n", spacer,
            tempStrSize * sizeof(char), numCurrentDivisions * sizeof(int));
        return (ERROR);
        }

#   if defined DEBUG_SETCHAINPARAMS
    /* only for debugging */
    YFlt      lnPriorRatio = 0.0, lnProposalRatio = 0.0;
#   endif

    /* allocate space for parameters */
    if (memAllocs[ALLOC_PARAMS] == YES)
        {
        for (i=0; i<numParams; i++)
            {
            SAFEFREE (params[i].name);
            if (params[i].paramHeader)
                {
                free (params[i].paramHeader);
                params[i].paramHeader = NULL;
                }
            }
        free (params);
        free (relevantParts);
        params = NULL;
        relevantParts = NULL;
        memAllocs[ALLOC_PARAMS] = NO;
        }

    /* wipe all chain parameter information */
    numParams = 0;
    numTrees = 0;
    chainHasAdgamma = NO;

    /* figure out number of parameters */
    /* this relies on activeParams[j][i] being set to 1, 2, ..., numParams */
    /* which is taken care of in SetUpLinkTable () */
    nRelParts = 0;
    for (j=0; j<NUM_LINKED; j++)
        {
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (activeParams[j][i] > numParams)
                numParams = activeParams[j][i];
            if (activeParams[j][i] > 0)
                nRelParts++;
            }
        }

    params = (Param *) SafeMalloc (numParams * sizeof(Param));
    relevantParts = (int *) SafeMalloc (nRelParts * sizeof(int));
    if (!params || !relevantParts)
        {
        YvyraPrint ("%s   Problem allocating params and relevantParts\n", spacer);
        if (params)
            free (params);
        if (relevantParts)
            free (relevantParts);
        free (tempStr);
        return ERROR;
        }
    else
        memAllocs[ALLOC_PARAMS] = YES;

    /* fill in info on each parameter */
    nRelParts = 0;  /* now cumulative number of relevant partitions */
    for (k=0; k<numParams; k++)
        {
        p = &params[k];

        /* find affected partitions */
        numRelParts = 0;
        for (j=0; j<NUM_LINKED; j++)
            {
            for (i=0; i<numCurrentDivisions; i++)
                {
                if (activeParams[j][i] == k + 1)
                    {
                    numRelParts++;
                    isPartTouched[i] = YES;
                    }
                else
                    isPartTouched[i] = NO;
                }
            if (numRelParts > 0)
                break;
            }

        /* find pointer to modelParams and modelSettings of first relevant partition */
        /* this will be handy later on */
        for (i=0; i<numCurrentDivisions; i++)
            if (isPartTouched[i] == YES)
                break;
        mp = &modelParams[i];
        m  = &modelSettings[i];
        
        /* Set default min and max */
        p->min = p->max = NEG_INFINITY;

        /* Parameter nValues and nSubValues, which are needed for memory allocation
           are calculated for each case in the code below. nSympi, however, is
           only used for one special type of parameter and it therefore makes
           sense to initialize it to 0 here. The same applies to hasBinaryStd
           and nIntValues. To be safe, we set all to 0 here. */
        p->nValues = 0;
        p->nSubValues = 0;
        p->nIntValues = 0;
        p->nSympi = 0;
        p->hasBinaryStd = NO;
        
        /* should this parameter be printed to a file? */
        p->printParam = NO;

        /* set print subparams to 0 */
        p->nPrintSubParams = 0;
        
        /* check constraints for tree parameter ? */
        p->checkConstraints = NO;

        /* set index number of parameter */
        p->index = k;
        
        /* set prior function to NULL */
        p->LnPriorRatio = NULL;

        /* set prior paramters to NULL */
        p->priorParams = NULL;

        /* set affectsLikelihood to NO */
        p->affectsLikelihood = NO;

        /* set cpp event pointers to NULL */
        p->nEvents = NULL;
        p->position = NULL;
        p->rateMult = NULL;

        /* set header and name to NULL */
        p->paramHeader = NULL;
        p->name = NULL;

        /* set up relevant partitions */
        p->nRelParts = numRelParts;
        p->relParts = relevantParts + nRelParts;
        nRelParts += numRelParts;
        for (i=n=0; i<numCurrentDivisions; i++)
            if (isPartTouched[i] == YES)
                p->relParts[n++] = i;

        /* get partition descriptor */
        SafeStrcat(&partString,"");
        FillRelPartsString (p, &partString);
            
        /* set up information for parameter */
        if (j == P_PI)
            {
            /* Set up state frequencies *****************************************************************************/
            p->paramType = P_PI;
            p->min = 0.0;
            p->max = 1.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].stateFreq = p;

            if (mp->dataType == STANDARD)
                {
                p->paramTypeName = "Symmetric diricihlet distribution alpha_i parameter";
                SafeStrcat(&p->name, "Alpha_symdir");
                /* boundaries for alpha_i */
                p->min = ETA;
                p->max = POS_INFINITY;
                }
            else
                {
                p->paramTypeName = "Stationary state frequencies";
                SafeStrcat(&p->name, "Pi");
                }
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            /* and the number of values and subvalues needed */
            if (mp->dataType == STANDARD)
                {
                /* find number of model states */
                m->numModelStates = 2;
                for (c=0; c<numChar; c++)
                    {
                    for (i=0; i<p->nRelParts; i++)
                        {
                        if (partitionId[c][partitionNum] == p->relParts[i] + 1 && charInfo[c].numStates > m->numModelStates)
                            m->numModelStates = charInfo[c].numStates;
                        }
                    }
                for (i=0; i<p->nRelParts; i++)
                    modelSettings[p->relParts[i]].numModelStates = m->numModelStates;

                /* symmetric hyperprior with only one variable (0 if equal) */
                p->nValues = 1;     /* set to 0 below for the SYMPI_EQUAL model */
                if (!strcmp(mp->symPiPr,"Uniform"))
                    {
                    if (m->numModelStates > 2) 
                        p->paramId = SYMPI_UNI_MS;
                    else
                        p->paramId = SYMPI_UNI;
                    }
                else if (!strcmp(mp->symPiPr,"Exponential"))
                    {
                    if (m->numModelStates > 2)
                        p->paramId = SYMPI_EXP_MS;
                    else
                        p->paramId = SYMPI_EXP;
                    }
                else if (!strcmp(mp->symPiPr,"Fixed"))
                    {
                    if (AreDoublesEqual(mp->symBetaFix, -1.0, 0.00001) == YES)
                        {
                        p->paramId = SYMPI_EQUAL;
                        p->nValues = 0;
                        }
                    else
                        {
                        if (m->numModelStates > 2) 
                            p->paramId = SYMPI_FIX_MS;
                        else
                            p->paramId = SYMPI_FIX;
                        }
                    }
                p->nSubValues = 0;  /* store state frequencies in p->stdStateFreqs */
                if (p->nValues == 1 && strcmp(mp->symPiPr,"Fixed") != 0)
                    {
                    p->printParam = YES;
                    SafeStrcat (&p->paramHeader, "alpha_symdir");
                    SafeStrcat (&p->paramHeader, partString);
                    }
                /* further processing done in ProcessStdChars */
                }
            else
                {
                /* deal with all models except standard */
                /* no hyperprior or fixed to one value, set default to 0  */
                p->nValues = 0;
                /* first check if model is stationary; if not, we process it differently */
                if (!strcmp(mp->statefreqModel,"Stationary"))          
                    {
                    /* one subvalue for each state */
                    p->nSubValues = mp->nStates;    /* mp->nStates is set to 20 if DNA || RNA && nucmodel==PROTEIN */
                    if (!strcmp(mp->stateFreqPr, "Dirichlet"))
                        {
                        p->paramId = PI_DIR;
                        p->nValues = mp->nStates;
                        }
                    else if (!strcmp(mp->stateFreqPr, "Fixed") && !strcmp(mp->stateFreqsFixType,"User"))
                        p->paramId = PI_USER;
                    else if (!strcmp(mp->stateFreqPr, "Fixed") && !strcmp(mp->stateFreqsFixType,"Empirical"))
                        p->paramId = PI_EMPIRICAL;
                    else if (!strcmp(mp->stateFreqPr, "Fixed") && !strcmp(mp->stateFreqsFixType,"Equal"))
                        {
                        p->paramId = PI_EQUAL;
                        }
                    
                    if (m->dataType == PROTEIN)
                        {
                        if (!strcmp(mp->aaModelPr, "Fixed"))
                            {
                            if (!strcmp(mp->aaModel, "Poisson"))
                                p->paramId = PI_EQUAL;
                            else if (!strcmp(mp->aaModel, "Equalin") || !strcmp(mp->aaModel, "Gtr"))
                                {
                                /* p->paramId stays to what it was set to above */
                                }
                            else
                                p->paramId = PI_FIXED;
                            }
                        else
                            p->paramId = PI_FIXED;
                        }
                    
                    if (p->paramId == PI_DIR)
                        p->printParam = YES;
                    if (0) /* molecular datatypes not supported */
                        {
                        if (!strcmp(mp->nucModel, "4by4"))
                            {
                            sprintf (temp, "pi(%c)", '?');
                            SafeStrcat (&p->paramHeader,temp);
                            SafeStrcat (&p->paramHeader,partString);
                            for (n1=1; n1<4; n1++)
                                {
                                sprintf (temp, "\tpi(%c)", (char)(n1));
                                SafeStrcat (&p->paramHeader,temp);
                                SafeStrcat (&p->paramHeader,partString);
                                }
                            }
                        else if (!strcmp(mp->nucModel, "Doublet"))
                            {
                            (void)(tempCodon,0);
                            sprintf (temp, "pi(%s)", tempCodon);
                            SafeStrcat (&p->paramHeader,temp);
                            SafeStrcat (&p->paramHeader,partString);
                            for (n1=1; n1<16; n1++)
                                {
                                (void)(tempCodon,n1);
                                sprintf (temp, "\tpi(%s)", tempCodon);
                                SafeStrcat (&p->paramHeader,temp);
                                SafeStrcat (&p->paramHeader,partString);
                                }
                            }
                        else if (!strcmp(mp->nucModel, "Codon"))
                            {
                            for (c=0; c<p->nSubValues; c++)
                                {
                                if (mp->codonNucs[c][0] == 0)
                                    strcpy (tempCodon, "pi(A");
                                else if (mp->codonNucs[c][0] == 1)
                                    strcpy (tempCodon, "pi(C");
                                else if (mp->codonNucs[c][0] == 2)
                                    strcpy (tempCodon, "pi(G");
                                else
                                    strcpy (tempCodon, "pi(T");
                                if (mp->codonNucs[c][1] == 0)
                                    strcat (tempCodon, "A");
                                else if (mp->codonNucs[c][1] == 1)
                                    strcat (tempCodon, "C");
                                else if (mp->codonNucs[c][1] == 2)
                                    strcat (tempCodon, "G");
                                else
                                    strcat (tempCodon, "T");
                                if (mp->codonNucs[c][2] == 0)
                                    strcat (tempCodon, "A)");
                                else if (mp->codonNucs[c][2] == 1)
                                    strcat (tempCodon, "C)");
                                else if (mp->codonNucs[c][2] == 2)
                                    strcat (tempCodon, "G)");
                                else
                                    strcat (tempCodon, "T)");
                                if (c == 0)
                                    {
                                    SafeStrcat (&p->paramHeader, tempCodon);
                                    SafeStrcat (&p->paramHeader, partString);
                                    }
                                else
                                    {
                                    SafeStrcat (&p->paramHeader, "\t");
                                    SafeStrcat (&p->paramHeader, tempCodon);
                                    SafeStrcat (&p->paramHeader, partString);
                                    }
                                }
                            }
                        }
                    else if (m->dataType == PROTEIN)
                        {
                        if (FillRelPartsString (p, &partString) == YES)
                            {
                            SafeSprintf (&tempStr, &tempStrSize, "pi(Ala)%s\tpi(Arg)%s\tpi(Asn)%s\tpi(Asp)%s\tpi(Cys)%s\tpi(Gln)%s\tpi(Glu)%s\tpi(Gly)%s\tpi(His)%s\tpi(Ile)%s\tpi(Leu)%s\tpi(Lys)%s\tpi(Met)%s\tpi(Phe)%s\tpi(Pro)%s\tpi(Ser)%s\tpi(Thr)%s\tpi(Trp)%s\tpi(Tyr)%s\tpi(Val)%s",
                            partString, partString, partString, partString, partString, partString, partString, partString, partString, partString,
                            partString, partString, partString, partString, partString, partString, partString, partString, partString, partString);
                            SafeStrcat (&p->paramHeader, tempStr);
                            }
                        else
                            SafeStrcat (&p->paramHeader, "pi(Ala)\tpi(Arg)\tpi(Asn)\tpi(Asp)\tpi(Cys)\tpi(Gln)\tpi(Glu)\tpi(Gly)\tpi(His)\tpi(Ile)\tpi(Leu)\tpi(Lys)\tpi(Met)\tpi(Phe)\tpi(Pro)\tpi(Ser)\tpi(Thr)\tpi(Trp)\tpi(Tyr)\tpi(Val)");
                        }
                    else if (mp->dataType == RESTRICTION)
                        {
                        if (FillRelPartsString (p, &partString) == YES)
                            {
                            SafeSprintf (&tempStr, &tempStrSize, "pi(0)%s\tpi(1)%s", partString, partString);
                            SafeStrcat (&p->paramHeader, tempStr);
                            }
                        else
                            SafeStrcat (&p->paramHeader, "pi(0)\tpi(1)");
                        }
                    else
                        {
                        YvyraPrint ("%s   Unknown data type in SetModelParams\n", spacer);
                        }
                    }
                else /* if statefreqmodel != stationary */
                    {
                    if (mp->dataType != RESTRICTION)
                        {   
                        YvyraPrint ("%s   Error in SetModelParams: non-stationary models currently only implemented for RESTRICTION data type\n", spacer);
                        }
                    else
                        {
                        /* two subvalues for each state: stationary and root frequencies */
                        p->nSubValues = mp->nStates * 2;

                        if (!strcmp(mp->statefreqModel, "Mixed"))
                            {
                            p->paramId = DIRPI_MIX;
                            p->nValues = mp->nStates * 2;
                            }
                        else if (!strcmp(mp->stateFreqPr, "Dirichlet"))
                            {
                            if (!strcmp(mp->rootFreqPr, "Dirichlet"))
                                {
                                p->paramId = DIRPI_DIRxDIR;
                                }
                            else if (!strcmp(mp->rootFreqPr, "Fixed"))
                                {
                                p->paramId = DIRPI_DIRxFIXED;
                                }
                            else
                                {
                                YvyraPrint ("%s   Error in SetModelParams: unknown setting for Rootfreqpr\n", spacer); 
                                }
                            p->nValues = m->numModelStates  * 2; //SK: use whole nValues vector even if only one of the priors is Dirichlet, start OR root freqs
                            }
                        else if (!strcmp(mp->stateFreqPr, "Fixed"))
                            {
                            if (!strcmp(mp->rootFreqPr, "Dirichlet"))
                                {
                                p->paramId = DIRPI_FIXEDxDIR;
                                p->nValues = mp->nStates * 2;
                                }
                            else if (!strcmp(mp->rootFreqPr, "Fixed"))
                                {
                                p->paramId = DIRPI_FIXEDxFIXED;
                                p->nValues = 0;
                                }
                            else
                                {
                                YvyraPrint ("%s   Error in SetModelParams: unknown setting for Rootfreqpr\n", spacer);
                                }
                            }
                        else
                            {
                            YvyraPrint ("%s   Error in SetModelParams: unknown setting for Statfreqpr, choose either 'Dirichlet' or 'Fixed'. \n", spacer);
                            }

                        if (p->paramId != DIRPI_FIXEDxFIXED)
                            p->printParam = YES;

                        if (FillRelPartsString (p, &partString) == YES)
                            {
                            SafeSprintf(&tempStr, &tempStrSize, "pi(0)%s\tpi(1)%s\trootpi(0)%s\trootpi(1)%s", partString, partString, partString, partString);
                            SafeStrcat (&p->paramHeader, tempStr);
                            if (p->paramId == DIRPI_MIX)
                                {
                                SafeSprintf(&tempStr, &tempStrSize, "\tstatefrmod%s", partString);
                                SafeStrcat (&p->paramHeader, tempStr);
                                }
                            }
                        else
                            {
                            SafeStrcat (&p->paramHeader, "pi(0)\tpi(1)\trootpi(0)\trootpi(1)");
                            if (p->paramId == DIRPI_MIX)
                                SafeStrcat(&p->paramHeader, "\tstatefrmod");
                            }
                        }
                    }   // end nonstationary model 
                }
            }
        else if (j == P_MIXTURE_RATES)
            {
            /* Set up mixture of site rates *****************************************************************/
            p->paramType = P_MIXTURE_RATES;
            p->nValues = mp->numMixtCats;   /* used for the Dirichlet prior parameters */
            p->nSubValues = mp->numMixtCats;
            p->min = 1E-5;
            p->max = 100;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].mixtureRates = p;
            p->paramTypeName = "Rates of site rate mixture";
            SafeStrcpy(&p->name, "Mixturerates");
            SafeStrcat(&p->name, partString);
            for (i=0; i<p->nValues; ++i)
                {
                if (i!=0)
                    SafeStrcat(&p->paramHeader, "\t");
                SafeSprintf(&tempStr, &tempStrSize, "mixturerates%s[%d]", partString, i+1);
                SafeStrcat(&p->paramHeader, tempStr);
                }
            
            /* find the parameter x prior type */
            p->paramId = MIXTURE_RATES;
            
            /* always print */
            p->printParam = YES;
            }
        else if (j == P_SHAPE)
            {
            /* Set up shape parameter of gamma/lnorm *****************************************************************/
            p->paramType = P_SHAPE;
            p->nValues = 1;
            if(!strcmp(mp->ratesModel, "Lnorm"))
                p->nSubValues = mp->numLnormCats;
            else
                p->nSubValues = mp->numGammaCats;
            p->min = 1E-5;
            p->max = 100;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].shape = p;
            if(!strcmp(mp->ratesModel, "LNorm"))
                {
                p->paramTypeName = "SD of scaled lognormal distribution of site rates";
                SafeStrcat(&p->name, "Sigma");
                SafeStrcat(&p->paramHeader, "sigma");
                }
            else
                {
                p->paramTypeName = "Shape of scaled gamma distribution of site rates";
                SafeStrcat(&p->name, "Alpha");
                SafeStrcat(&p->paramHeader, "alpha");
                }
            SafeStrcat(&p->name, partString);
            SafeStrcat(&p->paramHeader, partString);

            /* find the parameter x prior type */
            mp = &modelParams[p->relParts[0]];
            if (!strcmp(mp->shapePr,"Uniform"))
                p->paramId = SHAPE_UNI;
            else if (!strcmp(mp->shapePr,"Exponential"))
                p->paramId = SHAPE_EXP;
            else
                p->paramId = SHAPE_FIX;
                
            if (p->paramId != SHAPE_FIX)
                p->printParam = YES;
            }
        else if (j == P_PINVAR)
            {
            /* Set up proportion of invariable sites ****************************************************************/
            p->paramType = P_PINVAR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 0.0;
            p->max = 1.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].pInvar = p;

            p->paramTypeName = "Proportion of invariable sites";
            SafeStrcat(&p->name, "Pinvar");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->pInvarPr,"Uniform"))
                p->paramId = PINVAR_UNI;
            else
                p->paramId = PINVAR_FIX;
                
            if (p->paramId != PINVAR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "pinvar");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_CORREL)
            {
            /* Set up correlation parameter of adgamma model ********************************************************/
            p->paramType = P_CORREL;
            p->nValues = 1;
            p->nSubValues = mp->numGammaCats * mp->numGammaCats;
            p->min = -1.0;
            p->max = 1.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].correlation = p;

            p->paramTypeName = "Autocorrelation of gamma distribution of site rates";
            SafeStrcat(&p->name, "Rho");
            SafeStrcat(&p->name, partString);
            chainHasAdgamma = YES;

            /* find the parameter x prior type */
            if (!strcmp(mp->adGammaCorPr,"Uniform"))
                p->paramId = CORREL_UNI;
            else
                p->paramId = CORREL_FIX;

            if (p->paramId != CORREL_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "rho");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_SWITCH)
            {
            /* Set up switchRates for covarion model ****************************************************************/
            p->paramType = P_SWITCH;
            p->nValues = 2;
            p->nSubValues = mp->numGammaCats * mp->numGammaCats;
            p->min = 0.0;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].switchRates = p;

            p->paramTypeName = "Switch rates of covarion model";
            SafeStrcat(&p->name, "s_cov");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->covSwitchPr,"Uniform"))
                p->paramId = SWITCH_UNI;
            else if (!strcmp(mp->covSwitchPr,"Exponential"))
                p->paramId = SWITCH_EXP;
            else
                p->paramId = SWITCH_FIX;
                
            if (p->paramId != SWITCH_FIX)
                p->printParam = YES;

            SafeStrcat (&p->paramHeader, "s(off->on)");
            SafeStrcat (&p->paramHeader, partString);
            SafeStrcat (&p->paramHeader, "\ts(on->off)");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_RATEMULT)
            {
            /* Set up rateMult for partition specific rates ***********************************************************/
            p->paramType = P_RATEMULT;
            if (!strcmp(mp->ratePr,"Fixed") || numCurrentDivisions == 1)
                {
                p->nValues = 1;
                p->nSubValues = 0;
                }
            else
                {
                p->nValues = p->nRelParts = numRelParts; /* keep scaled division rates in value                        */
                p->nSubValues = p->nValues * 2;          /* keep number of uncompressed chars for scaling in subValue  */
                }                                        /* also keep Dirichlet prior parameters here                  */
            p->min = POS_MIN;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].rateMult = p;

            p->paramTypeName = "Partition-specific rate multiplier";
            SafeStrcat(&p->name, "Ratemultiplier");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (p->nSubValues == 0)
                p->paramId = RATEMULT_FIX;
            else
                p->paramId = RATEMULT_DIR;

            if (p->paramId != RATEMULT_FIX)
                p->printParam = YES;
            for (i=0; i<numCurrentDivisions; i++)
                {
                if (isPartTouched[i] == YES)
                    {
                    sprintf (tempMult, "m{%d}", i+1);
                    if (i == 0)
                        SafeStrcat (&p->paramHeader, tempMult);
                    else
                        {
                        SafeStrcat (&p->paramHeader, "\t");
                        SafeStrcat (&p->paramHeader, tempMult);
                        }
                    }
                }
            }
        else if (j == P_TOPOLOGY)
            {
            /* Set up topology **************************************************************************************/
            p->paramType = P_TOPOLOGY;
            p->nValues = 0;
            p->nSubValues = 0;
            p->min = NEG_INFINITY;  /* NA */
            p->max = NEG_INFINITY;  /* NA */
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].topology = p;

            p->paramTypeName = "Topology";
            SafeStrcat(&p->name, "Tau");
            SafeStrcat(&p->name, partString);
                    
            /* check that the model is not parsimony for all of the relevant partitions */
            areAllPartsParsimony = YES;
            for (i=0; i<p->nRelParts; i++)
                {
                if (modelSettings[p->relParts[i]].parsModelId == NO)
                    areAllPartsParsimony = NO;
                }

            /* check if topology is clock or nonclock */
            nClockBrlens = 0;
            nCalibratedBrlens = 0;
            nRelaxedBrlens = 0;
            if (areAllPartsParsimony == NO)
                {
                for (i=0; i<p->nRelParts; i++)
                    {
                    if (!strcmp(modelParams[p->relParts[i]].brlensPr, "Clock"))
                        {
                        nClockBrlens++;
                        if (strcmp(modelParams[p->relParts[i]].clockVarPr,"Strict") != 0)  /* not strict */
                            nRelaxedBrlens++;
                        if (!strcmp(modelParams[p->relParts[i]].nodeAgePr,"Calibrated"))
                            nCalibratedBrlens++;
                        }
                    }
                }
            
            /* now find the parameter x prior type */
            if (areAllPartsParsimony == YES)
                {
                if (!strcmp(mp->topologyPr, "Uniform"))
                    p->paramId = TOPOLOGY_PARSIMONY_UNIFORM;
                else if (!strcmp(mp->topologyPr,"Constraints"))
                    p->paramId = TOPOLOGY_PARSIMONY_CONSTRAINED;
                else
                    p->paramId = TOPOLOGY_PARSIMONY_FIXED;
                /* For this case, we also need to set the brlens ptr of the relevant partitions
                   so that it points to the topology parameter, since the rest of the
                   program will try to access the tree through this pointer. In FillTreeParams,
                   we will make sure that a pure parsimony topology parameter contains a pointer
                   to the relevant tree (like a brlens parameter would normally) */
                for (i=0; i<p->nRelParts; i++)
                    modelSettings[p->relParts[i]].brlens = p;
                }
            else
                {
                /* we assume for now that there is only one branch length; we will correct this
                   later in AllocateTreeParams if there are more than one set of branch lengths,
                   which is only possible for non-clock trees */
                if (!strcmp(mp->topologyPr, "Speciestree"))
                    p->paramId = TOPOLOGY_SPECIESTREE;
                else if (mp->dataType == RESTRICTION && strcmp(mp->statefreqModel, "Stationary") != 0)
                    {
                    /* rooted non-clock trees (directional substitution models) */
                    if (!strcmp(mp->topologyPr, "Uniform"))
                        p->paramId = TOPOLOGY_RNCL_UNIFORM;
                    else if (!strcmp(mp->topologyPr, "Constraints"))
                        p->paramId = TOPOLOGY_RNCL_CONSTRAINED;
                    else
                        p->paramId = TOPOLOGY_RNCL_FIXED;
                    }
                else if (!strcmp(mp->topologyPr, "Uniform") && nClockBrlens == 0)
                    p->paramId = TOPOLOGY_NCL_UNIFORM_HOMO;
                else if (!strcmp(mp->topologyPr,"Constraints") && nClockBrlens == 0)
                    p->paramId = TOPOLOGY_NCL_CONSTRAINED_HOMO;
                else if (!strcmp(mp->topologyPr,"Fixed") && nClockBrlens == 0)
                    p->paramId = TOPOLOGY_NCL_FIXED_HOMO;
                /* all below have clock branch lengths */
                else if (!strcmp(mp->topologyPr, "Uniform") && nRelaxedBrlens == 0)
                    {
                    if (nCalibratedBrlens >= 1)
                        p->paramId = TOPOLOGY_CCL_UNIFORM;
                    else
                        p->paramId = TOPOLOGY_CL_UNIFORM;
                    }
                else if (!strcmp(mp->topologyPr,"Constraints") && nRelaxedBrlens == 0)
                    {
                    if (nCalibratedBrlens >= 1)
                        p->paramId = TOPOLOGY_CCL_CONSTRAINED;
                    else
                        p->paramId = TOPOLOGY_CL_CONSTRAINED;
                    }
                else if (!strcmp(mp->topologyPr,"Fixed") && nRelaxedBrlens == 0)
                    {
                    if (nCalibratedBrlens >= 1)
                        p->paramId = TOPOLOGY_CCL_FIXED;
                    else
                        p->paramId = TOPOLOGY_CL_FIXED;
                    }
                /* all below have relaxed clock branch lengths */
                else if (!strcmp(mp->topologyPr, "Uniform"))
                    {
                    if (nCalibratedBrlens >= 1)
                        p->paramId = TOPOLOGY_RCCL_UNIFORM;
                    else
                        p->paramId = TOPOLOGY_RCL_UNIFORM;
                    }
                else if (!strcmp(mp->topologyPr,"Constraints"))
                    {
                    if (nCalibratedBrlens >= 1)
                        p->paramId = TOPOLOGY_RCCL_CONSTRAINED;
                    else
                        p->paramId = TOPOLOGY_RCL_CONSTRAINED;
                    }
                else if (!strcmp(mp->topologyPr,"Fixed"))
                    {
                    if (nCalibratedBrlens >= 1)
                        p->paramId = TOPOLOGY_RCCL_FIXED;
                    else
                        p->paramId = TOPOLOGY_RCL_FIXED;
                    }
                }

            /* should we print the topology? */
            for (i=0; i<p->nRelParts; i++)
                if (strcmp(modelParams[p->relParts[i]].treeFormat,"Topology") != 0)
                    break;
            if (i == p->nRelParts)
                p->printParam = YES;
            else
                p->printParam = NO;
            /* but always print parsimony topology */
            if (areAllPartsParsimony == YES)
                p->printParam = YES;
            /* and never print fixed topologies */
            if (p->paramId == TOPOLOGY_RCL_FIXED ||
                p->paramId == TOPOLOGY_RCCL_FIXED ||
                p->paramId == TOPOLOGY_CL_FIXED ||
                p->paramId == TOPOLOGY_CCL_FIXED ||
                p->paramId == TOPOLOGY_RNCL_FIXED ||
                p->paramId == TOPOLOGY_NCL_FIXED ||
                p->paramId == TOPOLOGY_PARSIMONY_FIXED)
                p->printParam = NO;
            }
        else if (j == P_BRLENS)
            {
            /* Set up branch lengths ********************************************************************************/
            p->paramType = P_BRLENS;
            p->nValues = 0;
            p->nSubValues = 0;
            p->min = 0.0;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].brlens = p;

            p->paramTypeName = "Branch lengths";
            SafeStrcat(&p->name, "V");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (modelSettings[p->relParts[0]].parsModelId == YES)
                p->paramId = BRLENS_PARSIMONY;
            else if (!strcmp(mp->brlensPr, "Clock"))
                {
                if (!strcmp(mp->clockPr,"Uniform"))
                    p->paramId = BRLENS_CLOCK_UNI;
                else if (!strcmp(mp->clockPr,"Coalescence"))
                    p->paramId = BRLENS_CLOCK_COAL;
                else if (!strcmp(mp->clockPr,"Birthdeath"))
                    p->paramId = BRLENS_CLOCK_BD;
                else if (!strcmp(mp->clockPr,"Speciestreecoalescence"))
                    p->paramId = BRLENS_CLOCK_SPCOAL;
                else if (!strcmp(mp->clockPr,"Fossilization"))
                    p->paramId = BRLENS_CLOCK_FOSSIL;
                else if (!strcmp(mp->clockPr,"Fixed"))
                    p->paramId = BRLENS_CLOCK_FIXED;
                }
            else if (!strcmp(mp->brlensPr, "Unconstrained"))
                {
                if (!strcmp(mp->unconstrainedPr,"Uniform"))
                    p->paramId = BRLENS_UNI;
                else if (!strcmp(mp->unconstrainedPr,"Exponential"))
                    p->paramId = BRLENS_EXP;
                else if (!strcmp(mp->unconstrainedPr,"GammaDir"))
                    p->paramId = BRLENS_GamDir;
                else if (!strcmp(mp->unconstrainedPr,"invGamDir"))
                    p->paramId = BRLENS_iGmDir;
                else if (!strcmp(mp->unconstrainedPr,"twoExp"))
                    p->paramId = BRLENS_twoExp;
                }
            else if (!strcmp(mp->brlensPr,"Fixed"))
                {
                p->paramId = BRLENS_FIXED;
                }

            /* should we print the branch lengths? */
            for (i=0; i<p->nRelParts; i++)
                if (strcmp(modelParams[p->relParts[i]].treeFormat,"Topology") != 0)
                    break;
            if (i < p->nRelParts)
                p->printParam = YES;
            else
                p->printParam = NO;
            }
        else if (j == P_SPECIESTREE)
            {
            /* Set up species tree **********************************************************************************/
            p->paramType = P_SPECIESTREE;
            p->nValues = 0;
            p->nSubValues = 0;
            p->min = NEG_INFINITY;  /* NA */
            p->max = NEG_INFINITY;  /* NA */
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].speciesTree = p;

            p->paramTypeName = "Species tree";
            SafeStrcat(&p->name, "Spt");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            p->paramId = SPECIESTREE_UNIFORM;

            /* should we print the tree? */
            p->printParam = YES;
            }
        else if (j == P_SPECRATE)
            {
            /* Set up speciation rate *******************************************************************************/
            p->paramType = P_SPECRATE;
            if (!strcmp(mp->clockPr,"Fossilization"))
                p->nValues = mp->birthRateShiftNum +1;  // rate in each time interval
            else
                p->nValues = 1;
            p->nSubValues = 0;
            p->min = POS_MIN;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].speciationRates = p;

            p->paramTypeName = "Speciation rate";
            SafeStrcat(&p->name, "Net_speciation");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->speciationPr,"Uniform"))
                p->paramId = SPECRATE_UNI;
            else if (!strcmp(mp->speciationPr,"Exponential"))
                p->paramId = SPECRATE_EXP;
            else
                p->paramId = SPECRATE_FIX;

            if (p->paramId != SPECRATE_FIX)
                p->printParam = YES;
                
            if (p->nValues == 1)
                {
                SafeStrcat (&p->paramHeader, "net_speciation");
                SafeStrcat (&p->paramHeader, partString);
                }
            else for (i = 0; i < p->nValues; i++)
                {
                SafeSprintf(&tempStr, &tempStrSize, "\tnet_speciation_%d", i+1);
                SafeStrcat (&p->paramHeader, tempStr);
                SafeStrcat (&p->paramHeader, partString);
                }
            }
        else if (j == P_EXTRATE)
            {
            /* Set up extinction rates ******************************************************************************/
            p->paramType = P_EXTRATE;
            if (!strcmp(mp->clockPr,"Fossilization"))
                p->nValues = mp->deathRateShiftNum +1;  // rate in each time interval
            else
                p->nValues = 1;
            p->nSubValues = 0;
            p->min = 0.0;
            p->max = 1.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].extinctionRates = p;

            p->paramTypeName = "Extinction rate";
            SafeStrcat(&p->name, "Relative_extinction");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->extinctionPr,"Beta"))
                p->paramId = EXTRATE_BETA;
            else if (!strcmp(mp->extinctionPr,"Exponential"))
                p->paramId = EXTRATE_EXP;
            else
                p->paramId = EXTRATE_FIX;

            if (p->paramId != EXTRATE_FIX)
                p->printParam = YES;
                
            if (p->nValues == 1)
                {
                SafeStrcat (&p->paramHeader, "relative_extinction");
                SafeStrcat (&p->paramHeader, partString);
                }
            else for (i = 0; i < p->nValues; i++)
                {
                SafeSprintf(&tempStr, &tempStrSize, "\trelative_extinction_%d", i+1);
                SafeStrcat (&p->paramHeader, tempStr);
                SafeStrcat (&p->paramHeader, partString);
                }
            }
        else if (j == P_FOSLRATE)
            {
            /* Set up fossilization rates */
            p->paramType = P_FOSLRATE;
            p->nValues = mp->fossilSamplingNum +1;  // rate in each time interval
            p->nSubValues = 0;
            p->min = 0.0;
            p->max = 1.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].fossilizationRates = p;
            
            p->paramTypeName = "Fossilization rate";
            SafeStrcat(&p->name, "Relative_fossilization");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            if (!strcmp(mp->fossilizationPr,"Beta"))
                p->paramId = FOSLRATE_BETA;
            else if (!strcmp(mp->fossilizationPr,"Exponential"))
                p->paramId = FOSLRATE_EXP;
            else
                p->paramId = FOSLRATE_FIX;
            
            if (p->paramId != FOSLRATE_FIX)
                p->printParam = YES;
            
            if (p->nValues == 1)
                {
                SafeStrcat (&p->paramHeader, "relative_fossilization");
                SafeStrcat (&p->paramHeader, partString);
                }
            else for (i = 0; i < p->nValues; i++)
                {
                SafeSprintf(&tempStr, &tempStrSize, "\trelative_fossilization_%d", i+1);
                SafeStrcat (&p->paramHeader, tempStr);
                SafeStrcat (&p->paramHeader, partString);
                }
            }
        else if (j == P_POPSIZE)
            {
            /* Set up population size *******************************************************************************/
            p->paramType = P_POPSIZE;
            if (!strcmp(mp->topologyPr,"Speciestree") && !strcmp(mp->popVarPr, "Variable"))
                p->nValues = 2 * numSpecies;
            else
                p->nValues = 1;
            p->nSubValues = 0;
            p->min = POS_MIN;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].popSize = p;

            p->paramTypeName = "Coalescent process population size parameter";
            SafeStrcat(&p->name, "Popsize");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->popSizePr,"Uniform"))
                {
                p->paramId      = POPSIZE_UNI;
                p->LnPriorRatio = &LnProbRatioUniform;
                p->priorParams  = mp->popSizeUni;
                p->LnPriorProb  = &LnPriorProbUniform;
                }
            else if (!strcmp(mp->popSizePr,"Normal"))
                {
                p->paramId      = POPSIZE_NORMAL;
                p->LnPriorRatio = &LnProbRatioTruncatedNormal;
                p->priorParams  = mp->popSizeNormal;
                p->LnPriorProb  = &LnPriorProbTruncatedNormal;
                }
            else if (!strcmp(mp->popSizePr,"Lognormal"))
                {
                p->paramId      = POPSIZE_LOGNORMAL;
                p->LnPriorRatio = &LnProbRatioLognormal;
                p->priorParams  = mp->popSizeLognormal;
                p->LnPriorProb  = &LnPriorProbLognormal;
                }
            else if (!strcmp(mp->popSizePr,"Gamma"))
                {
                p->paramId      = POPSIZE_GAMMA;
                p->LnPriorRatio = &LnProbRatioGamma;
                p->priorParams  = mp->popSizeGamma;
                p->LnPriorProb  = &LnPriorProbGamma;
                }
            else
                {
                p->paramId      = POPSIZE_FIX;
                p->LnPriorRatio = NULL;
                p->priorParams  = &mp->popSizeFix;
                p->LnPriorProb  = &LnPriorProbFix;
                }

            if (p->paramId != POPSIZE_FIX && p->nValues == 1)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "popsize");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_GROWTH)
            {
            /* Set up growth rate ***********************************************************************************/
            p->paramType = P_GROWTH;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;  // TODO: growth rate can be negative (also change Move_Growth_M)
            p->max = 100.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].growthRate = p;

            p->paramTypeName = "Population growth rate";
            SafeStrcat(&p->name, "R_pop");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->growthPr,"Uniform"))
                {
                p->paramId = GROWTH_UNI;
                p->priorParams = mp->growthUni;
                p->min = mp->growthUni[0];
                p->max = mp->growthUni[1];
                p->LnPriorRatio = &LnProbRatioUniform;
                p->LnPriorProb = &LnPriorProbUniform;
                }
            else if (!strcmp(mp->growthPr,"Exponential"))
                {
                p->paramId = GROWTH_EXP;
                p->priorParams = &mp->growthExp;
                p->LnPriorRatio = &LnProbRatioExponential;
                p->LnPriorProb = &LnPriorProbExponential;
                }
            else if (!strcmp(mp->growthPr,"Normal"))
                {
                p->paramId = GROWTH_NORMAL;
                p->priorParams = mp->growthNorm;
                p->LnPriorRatio = &LnProbRatioNormal;
                p->LnPriorProb = &LnPriorProbNormal;
                }
            else
                {
                p->paramId = GROWTH_FIX;
                p->priorParams = &mp->growthFix;
                p->min = p->max = mp->growthFix;
                p->LnPriorRatio = &LnProbRatioFix;
                p->LnPriorProb = &LnPriorProbFix;
                }

            if (p->paramId != GROWTH_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "growthRate");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_BMCORR)
            {
            /* Set up correlation parameter for brownian motion */
            p->paramType = P_BMCORR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = -1.0;
            p->max = 1.0;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].brownCorr = p;

            p->paramTypeName = "Browian correlation";
            SafeStrcat(&p->name, "BM_corr");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->brownCorrPr,"Uniform"))
                {
                p->paramId = BMCORR_UNI;
                p->priorParams = mp->brownCorrUni;
                p->min = mp->brownCorrUni[0];
                p->max = mp->brownCorrUni[1];
                p->LnPriorRatio = &LnProbRatioUniform;
                p->LnPriorProb = &LnPriorProbUniform;
                }
            else
                {
                p->paramId = BMCORR_FIX;
                p->priorParams = &mp->brownCorrFix;
                p->min = p->max = mp->brownCorrFix;
                p->LnPriorRatio = &LnProbRatioFix;
                p->LnPriorProb = &LnPriorProbFix;
                }

            if (p->paramId != BMCORR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "BM_corr");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_BMSIGMA)
            {
            /* Set up sigma parameter for brownian motion */
            p->paramType = P_BMSIGMA;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].brownSigma = p;

            p->paramTypeName = "Browian sigma";
            SafeStrcat(&p->name, "BM_sigma");
            SafeStrcat(&p->name, partString);

            /* find the parameter x prior type */
            if (!strcmp(mp->brownScalePr,"Uniform"))
                {
                p->paramId = BMSIGMA_UNI;
                p->priorParams = mp->brownScaleUni;
                p->min = mp->brownScaleUni[0];
                p->max = mp->brownScaleUni[1];
                p->LnPriorRatio = &LnProbRatioUniform;
                p->LnPriorProb = &LnPriorProbUniform;
                }
            else if (!strcmp(mp->brownScalePr,"Gamma"))
                {
                p->paramId = BMSIGMA_GAMMA;
                p->priorParams = mp->brownScaleGamma;
                p->LnPriorRatio = &LnProbRatioGamma;
                p->LnPriorProb = &LnPriorProbGamma;
                }
            else
                {
                p->paramId = BMSIGMA_FIX;
                p->priorParams = &mp->brownScaleFix;
                p->min = p->max = mp->brownScaleFix;
                p->LnPriorRatio = &LnProbRatioFix;
                p->LnPriorProb = &LnPriorProbFix;
                }

            if (p->paramId != BMSIGMA_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "BM_sigma");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_CPPRATE)
            {
            /* Set up cpprate ***************************************************************************************/
            p->paramType = P_CPPRATE;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 0.0;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].cppRate = p;

            p->paramTypeName = "Rate of rate-multiplying compound Poisson process";
            SafeStrcat(&p->name, "Lambda_cpp");
            SafeStrcat(&p->name, partString);
                    
            /* find the parameter x prior type */
            if (!strcmp(mp->cppRatePr,"Exponential"))
                p->paramId = CPPRATE_EXP;
            else
                p->paramId = CPPRATE_FIX;
            
            if (p->paramId != CPPRATE_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "cppRate");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_CPPMULTDEV)
            {
            /* Set up sigma of cpp rate multipliers *****************************************************************/
            p->paramType = P_CPPMULTDEV;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 0.0;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].cppMultDev = p;

            p->paramTypeName = "Standard deviation (log) of CPP rate multipliers";
            SafeStrcat(&p->name, "Sigma_cpp");
            SafeStrcat(&p->name, partString);
                    
            /* find the parameter x prior type */
            if (!strcmp(mp->cppMultDevPr,"Fixed"))
                p->paramId = CPPMULTDEV_FIX;
            
            if (p->paramId != CPPMULTDEV_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "sigma_cpp");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_CPPEVENTS)
            {
            /* Set up cpp events parameter **************************************************************************/
            p->paramType = P_CPPEVENTS;
            p->nValues = 0;
            p->nSubValues = 2*numLocalTaxa;     /* keep effective branch lengths here (for all nodes to be on the safe side) */
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].cppEvents = p;

            p->paramTypeName = "Events of rate-multiplying compound Poisson process";
            SafeStrcat(&p->name, "CppEvents");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            p->paramId = CPPEVENTS;
            
            /* should we print values to .p file? */
            p->printParam = NO;
            
            SafeStrcat (&p->paramHeader, "cppEvents");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_TK02VAR)
            {
            /* Set up tk02 relaxed clock variance parameter *********************************************************/
            p->paramType = P_TK02VAR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].tk02var = p;

            p->paramTypeName = "Variance of lognormal distribution of branch rates";
            SafeStrcat(&p->name, "TK02var");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            if (!strcmp(mp->tk02varPr,"Uniform"))
                p->paramId = TK02VAR_UNI;
            else if (!strcmp(mp->tk02varPr,"Exponential"))
                p->paramId = TK02VAR_EXP;
            else
                p->paramId = TK02VAR_FIX;
            
            if (p->paramId != TK02VAR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "tk02var");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_TK02BRANCHRATES)
            {
            /* Set up tk02 relaxed clock rates parameter ************************************************************/
            p->paramType = P_TK02BRANCHRATES;
            p->nValues = 2*numLocalTaxa;     /* use to hold the branch rates; we need one rate for the root */
            p->nSubValues = 2*numLocalTaxa;  /* use to hold the effective branch lengths */
            p->min = RATE_MIN;
            p->max = RATE_MAX;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].tk02BranchRates = p;
            
            p->paramTypeName = "Branch rates of TK02 relaxed clock";
            SafeStrcat(&p->name, "TK02Brlens");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            p->paramId = TK02BRANCHRATES;
            
            /* should we print values to .p file? */
            p->printParam = NO;

            SafeStrcat (&p->paramHeader, "tk02_brlens");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_WNVAR)
            {
            /* Set up white noise relaxed clock variance parameter **************************************************/
            p->paramType = P_WNVAR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].wnvar = p;

            p->paramTypeName = "Variance of WN model branch rates";
            SafeStrcat(&p->name, "WNvar");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            if (!strcmp(mp->wnvarPr,"Uniform"))
                p->paramId = WNVAR_UNI;
            else if (!strcmp(mp->wnvarPr,"Exponential"))
                p->paramId = WNVAR_EXP;
            else
                p->paramId = WNVAR_FIX;
            
            if (p->paramId != WNVAR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "wnvar");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_WNBRANCHRATES)
            {
            /* Set up white noise relaxed clock rates parameter *****************************************************/
            p->paramType = P_WNBRANCHRATES;
            p->nValues = 2*numLocalTaxa;     /* use to hold the branch rates; we need one rate for the root */
            p->nSubValues = 2*numLocalTaxa;  /* use to hold the effective branch lengths */
            p->min = RATE_MIN;
            p->max = RATE_MAX;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].igrBranchRates = p;
            
            p->paramTypeName = "Branch lengths of WN relaxed clock";
            SafeStrcat(&p->name, "WNBrlens");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            p->paramId = WNBRANCHRATES;
            
            /* should we print values to .p file? */
            p->printParam = NO;

            SafeStrcat (&p->paramHeader, "wn_brlens");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_IGRVAR)
            {
            /* Set up igr relaxed clock variance parameter **********************************************************/
            p->paramType = P_IGRVAR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].igrvar = p;

            p->paramTypeName = "Variance of IGR model branch rates";
            SafeStrcat(&p->name, "IGRvar");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            if (!strcmp(mp->igrvarPr,"Uniform"))
                p->paramId = IGRVAR_UNI;
            else if (!strcmp(mp->igrvarPr,"Exponential"))
                p->paramId = IGRVAR_EXP;
            else
                p->paramId = IGRVAR_FIX;
            
            if (p->paramId != IGRVAR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "igrvar");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_IGRBRANCHRATES)
            {
            /* Set up igr relaxed clock rates parameter *************************************************************/
            p->paramType = P_IGRBRANCHRATES;
            p->nValues = 2*numLocalTaxa;     /* use to hold the branch rates; we need one rate for the root */
            p->nSubValues = 2*numLocalTaxa;  /* use to hold the effective branch lengths */
            p->min = RATE_MIN;
            p->max = RATE_MAX;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].igrBranchRates = p;
            
            p->paramTypeName = "Branch lengths of IGR relaxed clock";
            SafeStrcat(&p->name, "IgrBrlens");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            p->paramId = IGRBRANCHRATES;
            
            /* should we print values to .p file? */
            p->printParam = NO;

            SafeStrcat (&p->paramHeader, "igr_brlens");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_ILNVAR)
            {
            /* Set up iln relaxed clock variance parameter **********************************************************/
            p->paramType = P_ILNVAR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].ilnvar = p;

            p->paramTypeName = "Variance of ILN model branch rates";
            SafeStrcat(&p->name, "ILNvar");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            if (!strcmp(mp->ilnvarPr,"Uniform"))
                p->paramId = ILNVAR_UNI;
            else if (!strcmp(mp->ilnvarPr,"Exponential"))
                p->paramId = ILNVAR_EXP;
            else
                p->paramId = ILNVAR_FIX;
            
            if (p->paramId != ILNVAR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "ilnvar");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_ILNBRANCHRATES)
            {
            /* Set up iln relaxed clock rates parameter *************************************************************/
            p->paramType = P_ILNBRANCHRATES;
            p->nValues = 2*numLocalTaxa;     /* use to hold the branch rates; we need one rate for the root */
            p->nSubValues = 2*numLocalTaxa;  /* use to hold the effective branch lengths */
            p->min = RATE_MIN;
            p->max = RATE_MAX;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].ilnBranchRates = p;
            
            p->paramTypeName = "Branch lengths of ILN relaxed clock";
            SafeStrcat(&p->name, "IlnBrlens");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            p->paramId = ILNBRANCHRATES;
            
            /* should we print values to .p file? */
            p->printParam = NO;

            SafeStrcat (&p->paramHeader, "iln_brlens");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_MIXEDVAR)
            {
            /* Set up mixed relaxed clock variance parameter ********************************************************/
            p->paramType = P_MIXEDVAR;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = 1E-6;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].mixedvar = p;
            
            p->paramTypeName = "Variance shared for mixed relaxed clock model";
            SafeStrcat(&p->name, "Mixedvar");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            if (!strcmp(mp->mixedvarPr,"Uniform"))
                p->paramId = MIXEDVAR_UNI;
            else if (!strcmp(mp->mixedvarPr,"Exponential"))
                p->paramId = MIXEDVAR_EXP;
            else
                p->paramId = MIXEDVAR_FIX;
            
            if (p->paramId != MIXEDVAR_FIX)
                p->printParam = YES;
            SafeStrcat (&p->paramHeader, "mixedvar");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_MIXEDBRCHRATES)
            {
            /* Set up mixed relaxed clock rates parameter ***********************************************************/
            p->paramType = P_MIXEDBRCHRATES;
            p->nValues = 2*numLocalTaxa;     /* use to hold the branch rates; we need one rate for the root */
            p->nSubValues = 2*numLocalTaxa;  /* use to hold the effective branch lengths */
            p->nIntValues = 1;               /* use to hold the model indicator: IGR or ILN */
            p->min = RATE_MIN;
            p->max = RATE_MAX;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].mixedBrchRates = p;
            
            p->paramTypeName = "Branch rates of mixed relaxed clock";
            SafeStrcat(&p->name, "MixedBrlens");
            SafeStrcat(&p->name, partString);
            
            /* find the parameter x prior type */
            p->paramId = MIXEDBRCHRATES;
            
            /* how to print model indicator (0 or 1) to .p file? */
            p->printParam = NO;
            
            SafeStrcat (&p->paramHeader, "mixed_brlens");
            SafeStrcat (&p->paramHeader, partString);
            }
        else if (j == P_CLOCKRATE)
            {
            /* Set up clockRate *************************************************************************************/
            p->paramType = P_CLOCKRATE;
            p->nValues = 1;
            p->nSubValues = 0;
            p->min = POS_MIN;
            p->max = POS_INFINITY;
            for (i=0; i<numCurrentDivisions; i++)
                if (isPartTouched[i] == YES)
                    modelSettings[i].clockRate = p;
    
            p->paramTypeName = "Base rate of clock";
            SafeStrcat(&p->name, "Clockrate");
            SafeStrcat(&p->name, partString);

            /* parameter does affect likelihoods */
            p->affectsLikelihood = YES;
            
            /* find the parameter x prior type */
            if (!strcmp(mp->clockRatePr,"Normal"))
                {
                p->paramId      = CLOCKRATE_NORMAL;
                p->LnPriorRatio = &LnProbRatioTruncatedNormal;
                p->priorParams  = mp->clockRateNormal;
                p->LnPriorProb  = &LnPriorProbTruncatedNormal;
                }
            else if (!strcmp(mp->clockRatePr,"Lognormal"))
                {
                p->paramId      = CLOCKRATE_LOGNORMAL;
                p->LnPriorRatio = &LnProbRatioLognormal;
                p->priorParams  = mp->clockRateLognormal;
                p->LnPriorProb  = &LnPriorProbLognormal;
                }
            else if (!strcmp(mp->clockRatePr,"Exponential"))
                {
                p->paramId      = CLOCKRATE_EXP;
                p->LnPriorRatio = &LnProbRatioExponential;
                p->priorParams  = &mp->clockRateExp;
                p->LnPriorProb  = &LnPriorProbExponential;
                }
            else if (!strcmp(mp->clockRatePr,"Gamma"))
                {
                p->paramId      = CLOCKRATE_GAMMA;
                p->LnPriorRatio = &LnProbRatioGamma;
                p->priorParams  = mp->clockRateGamma;
                p->LnPriorProb  = &LnPriorProbGamma;
                }
            else
                {
                p->paramId      = CLOCKRATE_FIX;
                p->LnPriorRatio = NULL;
                p->priorParams  = &mp->clockRateFix;
                p->LnPriorProb  = &LnPriorProbFix;
                }
                
            SafeStrcat (&p->paramHeader, "clockrate");
            SafeStrcat (&p->paramHeader, partString);
            if (p->paramId != CLOCKRATE_FIX)
                p->printParam = YES;
            }
        }
    free (tempStr);
    free (isPartTouched);
    SAFEFREE (partString);

    return (NO_ERROR);
}


/*----------------------------------------------------------------------------
|
|   SetMoves: This function will set up the applicable moves that could
|       potentially be used in updating the model parameters
|
-----------------------------------------------------------------------------*/
int SetMoves (void)
{
    int         i, j, k, moveIndex;
    Param       *param;
    
    /* free up previous moves if any */
    if (memAllocs[ALLOC_MOVES] == YES)
        {
        for (i=0; i<numApplicableMoves; i++)
            FreeMove (moves[i]);
        free (moves);
        moves = NULL;
        memAllocs[ALLOC_MOVES] = NO;
        }

    /* set up applicable moves                                   */
    /* each combination of moveType and param is a separate move */
    
    /* first count applicable moves */
    numApplicableMoves = 0;
    for (k=0; k<numParams; k++)
        {
        param = &params[k];
        for (i=0; i<NUM_MOVE_TYPES; i++)
            {
            if (moveTypes[i].level > userLevel)
                continue;
            if (moveTypes[i].isApplicable(param) == NO)
                continue;
            for (j=0; j<moveTypes[i].nApplicable; j++)
                {
                if (moveTypes[i].applicableTo[j] == param->paramId)
                    {
                    numApplicableMoves++;
                    break;
                    }
                }
            }
        }

    /* then allocate space for move pointers */
    moves = (MCMCMove **) SafeMalloc (numApplicableMoves * sizeof (MCMCMove *));
    if (!moves && numApplicableMoves > 0)   /* in the middle of defining a model, numApplicableMoves can be 0 */
        {
        YvyraPrint ("%s   Problem allocating moves\n", spacer);
        return (ERROR);
        }
    memAllocs[ALLOC_MOVES] = YES;

    /* finally allocate space for and set move defaults */
    moveIndex = 0;
    for (k=0; k<numParams; k++)
        {
        param = &params[k];
        for (i=0; i<NUM_MOVE_TYPES; i++)
            {   
            if (moveTypes[i].level > userLevel)
                continue;
            if (moveTypes[i].isApplicable(param) == NO)
                continue;
            for (j=0; j<moveTypes[i].nApplicable; j++)
                {
                if (moveTypes[i].applicableTo[j] == param->paramId)
                    {
                    if ((moves[moveIndex] = AllocateMove (&moveTypes[i], param)) == NULL)
                        break;
                    else
                        {
                        moves[moveIndex]->parm = param;
                        moveIndex++;
                        break;
                        }
                    }
                }
            }
        }

    if (moveIndex < numApplicableMoves)
        {
        for (i=0; i<moveIndex; i++)
            FreeMove (moves[i]);
        free (moves);
        memAllocs[ALLOC_MOVES] = NO;
        YvyraPrint ("%s   Problem setting moves\n", spacer);
        return (ERROR);
        }

    return (NO_ERROR);
}


/** SetPopSizeParam: Set population size values for a species tree from an input tree */
int SetPopSizeParam (Param *param, int chn, int state, PolyTree *pt)
{
    int         i, j, k, nLongsNeeded;
    YFlt      *values;
    Tree        *speciesTree;
    PolyNode    *pp;
    TreeNode    *p=NULL;

    nLongsNeeded = 1 + (pt->nNodes - pt->nIntNodes - 1) / nBitsInALong;

    /* Get pointer to values to be set */
    values = GetParamVals (param, chn, state);

    /* Get species tree */
    speciesTree = GetTree (modelSettings[param->relParts[0]].speciesTree, chn, state);

    /* Set them based on index of matching partitions */
    AllocatePolyTreePartitions(pt);
    AllocateTreePartitions(speciesTree);
    for (i=0; i<pt->nNodes; i++)
        {
        pp = pt->allDownPass[i];
        for (j=0; j<speciesTree->nNodes-1; j++)
            {
            p = speciesTree->allDownPass[j];
            for (k=0; k<nLongsNeeded; k++)
                {
                if (pp->partition[k] != p->partition[k])
                    break;
                }
            if (k == nLongsNeeded)
                break;
            }
        if (j == speciesTree->nNodes - 1)
            {
            YvyraPrint ("%s   Non-matching partitions when setting population size parameter", spacer);
            FreePolyTreePartitions(pt);
            FreeTreePartitions(speciesTree);
            return (ERROR);
            }
        values[p->index] = pt->popSize[pp->index];
    }

    FreePolyTreePartitions(pt);
    FreeTreePartitions(speciesTree);
    
    return (NO_ERROR);
}


/* SetRelaxedClockParam: set values for a relaxed clock param from an input tree */
int SetRelaxedClockParam (Param *param, int chn, int state, PolyTree *pt)
{
    int         i, j, k, *nEvents=NULL, *nEventsP=NULL, nLongsNeeded, isEventSet;
    YFlt       *effectiveBranchLengthP=NULL, *branchRate=NULL,
                **position=NULL, **rateMult=NULL, **positionP=NULL, **rateMultP=NULL,
                baseRate;
    Tree        *t;
    PolyNode    *pp;
    TreeNode    *p=NULL, *q;

    nLongsNeeded = 1 + (numLocalTaxa - 1) / nBitsInALong;

    /* set pointers to the right set of values */
    isEventSet = NO;
    if (param->paramType == P_CPPEVENTS)
        {
        /* find the right event set */
        for (i=0; i<pt->nESets; i++)
            {
            if (!strcmp(param->name,pt->eSetName[i]))
                break;
            }
        if (i == pt->nESets)
            {
            for (i=0; i<pt->nBSets; i++)
                if (!strcmp(param->name,pt->bSetName[i]))
                    break;

            if (i == pt->nBSets)
                return (NO_ERROR);
            else
                isEventSet = NO;
            }
        else
            isEventSet = YES;

        if (isEventSet == YES)
            {
            nEventsP  = pt->nEvents[i];
            positionP = pt->position[i];
            rateMultP = pt->rateMult[i];
            }
        else
            effectiveBranchLengthP = pt->effectiveBrLen[i];

        nEvents  = param->nEvents[2*chn+state];
        position = param->position[2*chn+state];
        rateMult = param->rateMult[2*chn+state];
        }
    else if (param->paramType == P_TK02BRANCHRATES || param->paramType == P_WNBRANCHRATES ||
             param->paramType == P_IGRBRANCHRATES || param->paramType == P_ILNBRANCHRATES || param->paramType == P_MIXEDBRCHRATES)
        {
        /* find the right effective branch length set */
        for (i=0; i<pt->nBSets; i++)
            if (!strcmp(param->name,pt->bSetName[i]))
                break;
        if (i == pt->nBSets)
            return (NO_ERROR);

        effectiveBranchLengthP = pt->effectiveBrLen[i];
        branchRate = GetParamVals (param, chn, state);
        }

    t = GetTree (param, chn, state);
    AllocatePolyTreePartitions (pt);
    AllocateTreePartitions (t);
    
    for (i=pt->nNodes-1; i>=0; i--)
        {
        pp = pt->allDownPass[i];
        for (j=0; j<t->nNodes; j++)
            {
            p = t->allDownPass[j];
            for (k=0; k<nLongsNeeded; k++)
                if (p->partition[k] != pp->partition[k])
                    break;
            if (k == nLongsNeeded)
                break;  /* match */
            }
        if (param->paramType == P_CPPEVENTS)
            {
            if (isEventSet == NO)
                {
                if (nEvents[p->index] != 1)
                    {
                    position[p->index] = (YFlt *) SafeRealloc ((void *) position[p->index], 1*sizeof(YFlt));
                    rateMult[p->index] = (YFlt *) SafeRealloc ((void *) rateMult[p->index], 1*sizeof(YFlt));
                    nEvents [p->index] = 1;
                    }
                position[p->index][0] = 0.5;
                if (p->anc->anc == NULL)
                    rateMult[p->index][0] = 1.0;
                else
                    {
                    baseRate = 1.0;
                    q = p->anc;
                    while (q->anc != NULL)
                        {
                        baseRate *= rateMult[q->index][0];
                        q = q->anc;
                        }
                    rateMult[p->index][0] = 2.0 * effectiveBranchLengthP[pp->index] / (p->length * baseRate) - 1.0;
                    }
                }
            else
                {
                if (nEvents[p->index] != nEventsP[pp->index])
                    {
                    if (nEventsP[pp->index] == 0)
                        {
                        free (position[p->index]);
                        free (rateMult[p->index]);
                        }
                    else
                        {
                        position[p->index] = (YFlt *) SafeRealloc ((void *) position[p->index], nEventsP[pp->index]*sizeof(YFlt));
                        rateMult[p->index] = (YFlt *) SafeRealloc ((void *) rateMult[p->index], nEventsP[pp->index]*sizeof(YFlt));
                        }
                    nEvents[p->index] = nEventsP[pp->index];
                    }
                for (j=0; j<nEventsP[pp->index]; j++)
                    {
                    position[p->index][j] = positionP[pp->index][j];
                    rateMult[p->index][j] = rateMultP[pp->index][j];
                    }
                }
            }
        else if (param->paramType == P_TK02BRANCHRATES)
            {
            if (p->anc->anc == NULL)
                branchRate[p->index] = 1.0;
            else if (p->length > 0.0)
                branchRate[p->index] = 2.0 * effectiveBranchLengthP[pp->index] / p->length - branchRate[p->anc->index];
            else
                branchRate[p->index] = branchRate[p->anc->index];
            }
        else if (param->paramType == P_ILNBRANCHRATES || param->paramType == P_IGRBRANCHRATES ||
                 param->paramType == P_MIXEDBRCHRATES || param->paramType == P_WNBRANCHRATES)
            {
            if (p->length > 0.0)
                branchRate[p->index] = effectiveBranchLengthP[pp->index] / p->length;
            else
                branchRate[p->index] = branchRate[p->anc->index];
            }
        }

    if (param->paramType == P_CPPEVENTS)
        {
        if (UpdateCppEvolLengths (param, t->root->left, chn) == ERROR)
            return (ERROR);
        }
    else if (param->paramType == P_TK02BRANCHRATES)
        {
        if (UpdateTK02EvolLengths (param, t, chn) == ERROR)
            return (ERROR);
        }
    else if (param->paramType == P_ILNBRANCHRATES || param->paramType == P_IGRBRANCHRATES ||
             param->paramType == P_MIXEDBRCHRATES || param->paramType == P_WNBRANCHRATES)
        {
        if (UpdateIndBrachLengths (param, t, chn) == ERROR)
            return (ERROR);
        }

    FreePolyTreePartitions (pt);
    FreeTreePartitions (t);

    return (NO_ERROR);
}


/*------------------------------------------------------------------------
|
|   SetUpAnalysis: Set parameters and moves
|
------------------------------------------------------------------------*/
int SetUpAnalysis (RandLong *seed)
{
    setUpAnalysisSuccess=NO;

    /* calculate number of characters and taxa */
    numLocalChar = NumNonExcludedChar ();

    /* we are checking later to make sure no partition is without characters */
    SetLocalTaxa ();
    if (numLocalTaxa <= 2)
        {
        YvyraPrint ("%s   There must be at least two included taxa, now there is %s\n", spacer,
            numLocalTaxa == 0 ? "none" : "only one");
        return (ERROR);
        }

    /* calculate number of global chains */
    numGlobalChains = chainParams.numRuns * chainParams.numChains;

    /* Set up link table */
    if (SetUpLinkTable () == ERROR)
        return (ERROR);
    
    /* Set up model info */
    if (SetModelInfo() == ERROR) 
        return (ERROR);
    
    /* Calculate number of (uncompressed) characters for each division */
    if (GetNumDivisionChars() == ERROR)
        return (ERROR);

    /* Compress data and calculate some things needed for setting up params. */
    if (CompressData() == ERROR)
        return (ERROR);

    /* Add dummy characters, if needed. */
    if (AddDummyChars() == ERROR)
        return (ERROR);

    /* Set up parameters for the chain. */
    if (SetModelParams () == ERROR)
        return (ERROR);
    
    /* Allocate normal params */
    if (AllocateNormalParams () == ERROR)
        return (ERROR);
    
    /* Allocate tree params */
    if (AllocateTreeParams () == ERROR)
        return (ERROR);
    
    /* Set default number of trees for sumt to appropriate number */
    sumtParams.numTrees = numTrees;

    /* Fill in normal parameters */
    if (FillNormalParams (seed, 0, numGlobalChains) == ERROR) 
        return (ERROR);

    /* Process standard characters (calculates bsIndex, tiIndex, and more). */
    if (ProcessStdChars(seed) == ERROR)
        return (ERROR);

    /* Fill in trees */
    if (FillTreeParams (seed, 0, numGlobalChains) == ERROR)
        return (ERROR);
    
    /* Set the applicable moves that could be used by the chain. */
    if (SetMoves () == ERROR)
        return (ERROR);
    
    setUpAnalysisSuccess=YES;
    
    return (NO_ERROR);
}


int SetUpLinkTable (void)
{
    int         i, j, k, m, paramCount, isApplicable1, isApplicable2,
                isFirst, isSame;
    int         *check, *modelId;

    check = (int *) SafeMalloc (2 * (size_t)numCurrentDivisions * sizeof (int));
    if (!check)
        {
        YvyraPrint ("%s   Problem allocating check (%d)\n", spacer, 2 * numCurrentDivisions * sizeof(int));
        return (ERROR);
        }
    modelId = check + numCurrentDivisions;

    for (j=0; j<NUM_LINKED; j++)
        for (i=0; i<numCurrentDivisions; i++)
            activeParams[j][i] = 0;
    
    if (numCurrentDivisions > 1)
        {
        paramCount = 0;
        for (j=0; j<NUM_LINKED; j++) /* loop over parameters */
            {
            isFirst = YES;
            for (i=0; i<numCurrentDivisions; i++)
                modelId[i] = 0;     
            for (i=0; i<numCurrentDivisions-1; i++) /* loop over partitions */
                {
                for (k=i+1; k<numCurrentDivisions; k++)
                    {
                    if (IsModelSame (j, i, k, &isApplicable1, &isApplicable2) == NO || linkTable[j][i] != linkTable[j][k])
                        {
                        /* we cannot link the parameters */
                        if (isApplicable1 == NO)
                            modelId[i] = -1;
                        if (isApplicable2 == NO)
                            modelId[k] = -1;
                        if (isApplicable1 == YES)
                            {
                            if (isFirst == YES && modelId[i] == 0)
                                {
                                modelId[i] = ++paramCount;
                                isFirst = NO;
                                }
                            else
                                {
                                if (modelId[i] == 0)
                                    modelId[i] = ++paramCount;
                                }
                            }
                        if (modelId[k] == 0 && isApplicable2 == YES)
                            modelId[k] = ++paramCount;
                        }
                    else
                        {
                        /* we can link the parameters */
                        if (isFirst == YES)
                            {
                            if (modelId[i] == 0)
                                modelId[i] = ++paramCount;
                            isFirst = NO;
                            }
                        else
                            {
                            if (modelId[i] == 0)
                                modelId[i] = ++paramCount;
                            }
                        modelId[k] = modelId[i];
                        }
                    }
                }
            for (i=0; i<numCurrentDivisions; i++)
                activeParams[j][i] = modelId[i];
            }
        }
    else
        {
        /* if we have only one partition, then we do things a bit differently */
        paramCount = 0;
        for (j=0; j<NUM_LINKED; j++) /* loop over parameters */
            {
            IsModelSame (j, 0, 0, &isApplicable1, &isApplicable2);
            if (isApplicable1 == YES)
                activeParams[j][0] = ++paramCount;
            else
                activeParams[j][0] = -1;
            }
        }
        
    /* Check that the same report format is specified for all partitions with the same rate multiplier */
    for (i=0; i<numCurrentDivisions; i++)
        check[i] = NO;
    for (i=0; i<numCurrentDivisions; i++)
        {
        m = activeParams[P_RATEMULT][i];
        if (m == -1 || check[i] == YES)
            continue;
        isSame = YES;
        for (j=i+1; j<numCurrentDivisions; j++)
            {
            if (activeParams[P_RATEMULT][j] == m)
                {
                check[i] = YES;
                if (strcmp(modelParams[i].ratemultFormat,modelParams[j].ratemultFormat)!= 0)
                    {
                    isSame = NO;
                    strcpy (modelParams[j].ratemultFormat, modelParams[i].ratemultFormat);
                    }
                }
            }
        if (isSame == NO)
            {
            YvyraPrint ("%s   WARNING: Report format for ratemult (parameter %d) varies across partitions.\n", spacer);
            YvyraPrint ("%s      yvyra will use the format for the first partition, which is %s.\n", spacer, modelParams[i].ratemultFormat);
            }
        }
       
    /* probably a good idea to clean up link table here */
    paramCount = 0;
    for (j=0; j<NUM_LINKED; j++)
        {
        for (i=0; i<numCurrentDivisions; i++)
            check[i] = NO;
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (check[i] == NO && activeParams[j][i] > 0)
                {
                m = activeParams[j][i];
                paramCount++;
                for (k=i; k<numCurrentDivisions; k++)
                    {
                    if (check[k] == NO && activeParams[j][k] == m)
                        {
                        activeParams[j][k] = paramCount;
                        check[k] = YES;
                        }
                    }
                }
            }
        }

    free (check);

    return (NO_ERROR);
}


/*------------------------------------------------------------------------
|
|   SetUpMoveTypes: Set up structs holding info on each move type
|
------------------------------------------------------------------------*/
void SetUpMoveTypes (void)
{
    /* Register the move type here when new move functions are added 
       Remember to check that the number of move types does not exceed NUM_MOVE_TYPES
       defined in bayes.h.         */
    int         i;
    MoveType    *mt;

    /* reset move types */
    for (i=0; i<NUM_MOVE_TYPES; i++)
        {
        mt = &moveTypes[i];
        mt->level = DEVELOPER;
        mt->numTuningParams = 0;
        mt->minimum[0] = mt->minimum[1] = -1000000000.0;
        mt->maximum[0] = mt->maximum[1] =  1000000000.0;
        mt->tuningParam[0] = mt->tuningParam[1] = 0.0;
        mt->nApplicable = 0;
        mt->name = mt->tuningName[0] = mt->tuningName[1] = "";
        mt->paramName = "";
        mt->subParams = NO;
        mt->relProposalProb = 0.0;
        mt->parsimonyBased = NO;
        mt->isApplicable = &IsApplicable;
        mt->Autotune = NULL;
        mt->targetRate = -1.0;
        }

    /* Moves are in alphabetic order after parameter name, which matches the name of a move function if
       there is a separate move function for the parameter. See proposal.h for declaration of move functions.
       Since 2010-10-04, some parameters use generalized move functions and do not have their own. */
    
    i = 0;
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = SYMPI_UNI;
    mt->applicableTo[1] = SYMPI_EXP;
    mt->applicableTo[2] = SYMPI_UNI_MS;
    mt->applicableTo[3] = SYMPI_EXP_MS;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_Beta;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* so-called lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_BrLen */
    mt = &moveTypes[i++];
    mt->name = "Random brlen hit with multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";  
    mt->applicableTo[0] = BRLENS_UNI;
    mt->applicableTo[1] = BRLENS_EXP;
    mt->applicableTo[2] = BRLENS_GamDir;
    mt->applicableTo[3] = BRLENS_iGmDir;
    mt->applicableTo[4] = BRLENS_twoExp;
    mt->nApplicable = 5;  // was 2
    mt->moveFxn = &Move_BrLen;
    mt->relProposalProb = 20.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (2.0);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;  /* change it smaller to avoid overflow in the exp function, same for following "smaller" */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_BMcorr */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = BMCORR_UNI;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_BMcorr;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* window size */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_BMsigma */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = BMSIGMA_UNI;
    mt->applicableTo[1] = BMSIGMA_GAMMA;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_BMsigma;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_ClockRate_M */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = CLOCKRATE_NORMAL;
    mt->applicableTo[1] = CLOCKRATE_LOGNORMAL;
    mt->applicableTo[2] = CLOCKRATE_GAMMA;
    mt->applicableTo[3] = CLOCKRATE_EXP;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_ClockRate_M;
    mt->relProposalProb = 4.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_Extinction */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = EXTRATE_BETA;
    mt->applicableTo[1] = EXTRATE_EXP;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_Extinction;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* window size */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_Fossilization */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = FOSLRATE_BETA;
    mt->applicableTo[1] = FOSLRATE_EXP;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_Fossilization;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* window size */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_AddBranch, for fossilization prior */
    mt = &moveTypes[i++];
    mt->name = "Add branch for FossilizedBD";
    mt->shortName = "AddBranch";
    // mt->subParams = YES;
    mt->applicableTo[0] = BRLENS_CLOCK_FOSSIL;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_AddBranch;
    mt->relProposalProb = 10.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level =STANDARD_USER;
    mt->isApplicable = &IsApplicable_AncestralFossil;

    /* Move_DelBranch, for fossilization prior */
    mt = &moveTypes[i++];
    mt->name = "Delete branch for FossilizedBD";
    mt->shortName = "DelBranch";
    // mt->subParams = YES;
    mt->applicableTo[0] = BRLENS_CLOCK_FOSSIL;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_DelBranch;
    mt->relProposalProb = 10.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level =STANDARD_USER;
    mt->isApplicable = &IsApplicable_AncestralFossil;

    /* Move_ExtSPR */
    mt = &moveTypes[i++];
    mt->name = "Extending SPR";
    mt->shortName = "ExtSPR";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->tuningName[1] = "Multiplier tuning parameter";
    mt->shortTuningName[1] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->applicableTo[2] = TOPOLOGY_RNCL_UNIFORM;
    mt->applicableTo[3] = TOPOLOGY_RNCL_CONSTRAINED;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_ExtSPR;
    mt->relProposalProb = 5.0;
    mt->numTuningParams = 2;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->tuningParam[1] = 2.0 * log (1.05);   /* lambda */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->minimum[1] = 0.00001;
    mt->maximum[1] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_ExtSPR1 */
    mt = &moveTypes[i++];
    mt->name = "Extending SPR variant 1";
    mt->shortName = "ExtSPR1";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->tuningName[1] = "Multiplier tuning parameter";
    mt->shortTuningName[1] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ExtSPR1;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 2;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->tuningParam[1] = 2.0 * log (1.05);   /* lambda */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->minimum[1] = 0.00001;
    mt->maximum[1] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_FiveTaxaOrMore;

    /* Move_ExtSPRClock */
    mt = &moveTypes[i++];
    mt->name = "Extending SPR for clock trees";
    mt->shortName = "ExtSPRClock";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->applicableTo[4] = TOPOLOGY_RCL_UNIFORM;
    mt->applicableTo[5] = TOPOLOGY_RCL_CONSTRAINED;
    mt->applicableTo[6] = TOPOLOGY_RCCL_UNIFORM;
    mt->applicableTo[7] = TOPOLOGY_RCCL_CONSTRAINED;
    mt->nApplicable = 8;
    mt->moveFxn = &Move_ExtSPRClock;
    mt->relProposalProb = 10.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;

    /* Move_ExtSPRClock_Fossil */
    mt = &moveTypes[i++];
    mt->name = "Extending fossil SPR for clock trees";
    mt->shortName = "ExtSPRClockFossil";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->applicableTo[4] = TOPOLOGY_RCL_UNIFORM;
    mt->applicableTo[5] = TOPOLOGY_RCL_CONSTRAINED;
    mt->applicableTo[6] = TOPOLOGY_RCCL_UNIFORM;
    mt->applicableTo[7] = TOPOLOGY_RCCL_CONSTRAINED;
    mt->nApplicable = 8;
    mt->moveFxn = &Move_ExtSPRClock_Fossil;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->parsimonyBased = NO;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;
    
    /* Move_ExtSS */
    mt = &moveTypes[i++];
    mt->name = "Extending subtree swapper";
    mt->shortName = "ExtSS";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->tuningName[1] = "Multiplier tuning parameter";
    mt->shortTuningName[1] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ExtSS;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 2;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->tuningParam[1] = 2.0 * log (1.05); /* lambda */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->minimum[1] = 0.00001;
    mt->maximum[1] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_ExtSSClock */
    mt = &moveTypes[i++];
    mt->name = "Extending subtree swapper";
    mt->shortName = "ExtSsClock";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->applicableTo[4] = TOPOLOGY_RCL_UNIFORM;
    mt->applicableTo[5] = TOPOLOGY_RCL_CONSTRAINED;
    mt->applicableTo[6] = TOPOLOGY_RCCL_UNIFORM;
    mt->applicableTo[7] = TOPOLOGY_RCCL_CONSTRAINED;
    mt->nApplicable = 8;
    mt->moveFxn = &Move_ExtSSClock;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->parsimonyBased = NO;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;

    /* Move_ExtTBR */
    mt = &moveTypes[i++];
    mt->name = "Extending TBR";
    mt->shortName = "ExtTBR";
    mt->subParams = YES;
    mt->tuningName[0] = "Extension probability";
    mt->shortTuningName[0] = "p_ext";
    mt->tuningName[1] = "Multiplier tuning parameter";
    mt->shortTuningName[1] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->applicableTo[2] = TOPOLOGY_RNCL_UNIFORM;
    mt->applicableTo[3] = TOPOLOGY_RNCL_CONSTRAINED;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_ExtTBR;
    mt->relProposalProb = 5.0;
    mt->numTuningParams = 2;
    mt->tuningParam[0] = 0.5; /* extension probability */
    mt->tuningParam[1] = 2.0 * log (1.05);   /* lambda */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 0.99999;
    mt->minimum[1] = 0.00001;
    mt->maximum[1] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FiveTaxaOrMore;

    /* Move_GeneRate_Dir */
    mt = &moveTypes[i++];
    mt->name = "Dirichlet proposal";
    mt->shortName = "Dirichlet";
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->applicableTo[0] = GENETREERATEMULT_DIR;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_GeneRate_Dir;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 1000.0; /* alphaPi */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 1000000.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25;

    /* Move_GeneTree1/2/3 removed: species tree not supported in yvyra */

    /* Move_Growth_M */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = GROWTH_UNI;
    mt->applicableTo[1] = GROWTH_EXP;
    mt->applicableTo[2] = GROWTH_NORMAL;
    mt->nApplicable = 3;
    mt->moveFxn = &Move_Growth_M;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log(1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_Local */
    mt = &moveTypes[i++];
    mt->name = "BAMBE's LOCAL";
    mt->shortName = "Local";
    mt->subParams = YES;
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_Local;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_LocalClock */
    mt = &moveTypes[i++];
    mt->name = "Modified LOCAL for clock trees";
    mt->shortName = "LocalClock";
    mt->subParams = YES;
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_LocalClock;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (2.0);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;

    /* Move_MixtureRates */
    mt = &moveTypes[i++];
    mt->name = "Dirichlet proposal";
    mt->shortName = "Dirichlet";
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->applicableTo[0] = MIXTURE_RATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_MixtureRates;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 100.0; /* alphaPi per rate */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 10000.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25;
    
    /* Move_MixtureRates_Slider */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = MIXTURE_RATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_MixtureRates_Slider;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.2;  /* window size (change in proportions) */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 1.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;
    
    /* Move_NNI */
    mt = &moveTypes[i++];
    mt->name = "NNI move for parsimony trees";
    mt->shortName = "ParsNNI";
    mt->applicableTo[0] = TOPOLOGY_PARSIMONY_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_PARSIMONY_CONSTRAINED;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_NNI;
    mt->relProposalProb = 10.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;    /* no extra parsimony scores are needed */
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_NNIClock */
    mt = &moveTypes[i++];
    mt->name = "NNI move for clock trees";
    mt->shortName = "NNIClock";
    mt->subParams = YES;
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->applicableTo[4] = TOPOLOGY_RCL_UNIFORM;
    mt->applicableTo[5] = TOPOLOGY_RCL_CONSTRAINED;
    mt->applicableTo[6] = TOPOLOGY_RCCL_UNIFORM;
    mt->applicableTo[7] = TOPOLOGY_RCCL_CONSTRAINED;
    mt->nApplicable = 8;
    mt->moveFxn = &Move_NNIClock;
    mt->relProposalProb = 20.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;

    /* Move_NNI */
    mt = &moveTypes[i++];
    mt->name = "NNI move";
    mt->shortName = "NNI";
    mt->subParams = YES;
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->applicableTo[2] = TOPOLOGY_RNCL_UNIFORM;
    mt->applicableTo[3] = TOPOLOGY_RNCL_CONSTRAINED;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_NNI;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_NNI_Hetero */
    mt = &moveTypes[i++];
    mt->name = "NNI move for trees with independent brlens";
    mt->shortName = "MultNNI";
    mt->subParams = YES;
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HETERO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HETERO;
    mt->nApplicable = 2; /* 3; */
    mt->moveFxn = &Move_NNI_Hetero;
    mt->relProposalProb = 15.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_NodeSlider */
    mt = &moveTypes[i++];
    mt->name = "Node slider (uniform on possible positions)";
    mt->shortName = "Nodeslider";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = BRLENS_UNI;
    mt->applicableTo[1] = BRLENS_EXP;
    mt->applicableTo[2] = BRLENS_GamDir;
    mt->applicableTo[3] = BRLENS_iGmDir;
    mt->applicableTo[4] = BRLENS_twoExp;
    mt->nApplicable = 5;
    mt->moveFxn = &Move_NodeSlider;
    mt->relProposalProb = 7.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;

    /* Move_NodeSliderClock */
    mt = &moveTypes[i++];
    mt->name = "Node depth window slider (clock-constrained)";
    mt->shortName = "NodesliderClock";
    mt->tuningName[0] = "Window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = BRLENS_CLOCK_UNI;
    mt->applicableTo[1] = BRLENS_CLOCK_COAL;
    mt->applicableTo[2] = BRLENS_CLOCK_BD;
    mt->applicableTo[3] = BRLENS_CLOCK_FOSSIL;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_NodeSliderClock;
    mt->relProposalProb = 30.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1; /* window size */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased eraser version 1";
    mt->shortName = "pEraser1";
    mt->subParams = YES;
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->tuningName[1] = "parsimony warp factor";
    mt->shortTuningName[1] = "warp";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_ParsEraser1;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 2;
    mt->tuningParam[0] = 0.5; /* alphaPi */
    mt->tuningParam[1] = 0.1; /* warp */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 10000.0;
    mt->minimum[1] = 0.00001;
    mt->maximum[1] = 0.99999;
    mt->parsimonyBased = YES;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_FiveTaxaOrMore;

    /* Move_ParsSPR asym */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased SPR";
    mt->shortName = "ParsSPR";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->tuningName[1] = "reweighting probability";
    mt->shortTuningName[1] = "r";
    mt->tuningName[2] = "typical branch length";
    mt->shortTuningName[2] = "v_t";
    mt->tuningName[3] = "multiplier tuning parameter";
    mt->shortTuningName[3] = "lambda";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ParsSPR;
    mt->relProposalProb = 0.0;      /* other variants appear to be more efficient */
    mt->numTuningParams = 4;
    mt->tuningParam[0] = 0.1;              /* warp */
    mt->tuningParam[1] = 0.05;             /* upweight and downweight probability */
    mt->tuningParam[2] = 0.03;             /* typical branch length */
    mt->tuningParam[3] = 2.0 * log (1.05); /* multiplier tuning parameter lambda */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 1.0;
    mt->minimum[1] = 0.0;
    mt->maximum[1] = 0.3;
    mt->minimum[2] = 0.0001;
    mt->maximum[2] = 0.5;
    mt->minimum[3] = 0.0001;
    mt->maximum[3] = 20.0;
    mt->parsimonyBased = YES;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_ParsSPR1 e^{-S} */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased SPR variant 1";
    mt->shortName = "ParsSPR";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->tuningName[1] = "reweighting probability";
    mt->shortTuningName[1] = "r";
    mt->tuningName[2] = "multiplier tuning parameter";
    mt->shortTuningName[2] = "lambda";
    mt->tuningName[3] = "moving distance";
    mt->shortTuningName[3] = "d";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ParsSPR1;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 4;
    mt->tuningParam[0] = 0.5;              /* warp */
    mt->tuningParam[1] = 0.05;             /* upweight and downweight probability */
    mt->tuningParam[2] = 2.0 * log (1.05); /* multiplier tuning parameter lambda */
    mt->tuningParam[3] = 10.0;             /* distance to move picked branch */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 5.0;
    mt->minimum[1] = 0.0;
    mt->maximum[1] = 0.3;
    mt->minimum[2] = 0.0001;
    mt->maximum[2] = 20.0;
    mt->minimum[3] = 2.0;
    mt->maximum[3] = 1000.0;
    mt->parsimonyBased = YES;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_ParsSPR2 S/N */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased SPR variant 2";
    mt->shortName = "ParsSPR2";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->tuningName[1] = "reweighting probability";
    mt->shortTuningName[1] = "r";
    mt->tuningName[2] = "typical branch length";
    mt->shortTuningName[2] = "v_t";
    mt->tuningName[3] = "multiplier tuning parameter";
    mt->shortTuningName[3] = "lambda";
    mt->tuningName[4] = "moving distance";
    mt->shortTuningName[4] = "d";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ParsSPR2;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 5;
    mt->tuningParam[0] = 0.1;              /* warp */
    mt->tuningParam[1] = 0.05;             /* upweight and downweight probability */
    mt->tuningParam[2] = 0.03;             /* typical branch length */
    mt->tuningParam[3] = 2.0 * log (1.05); /* multiplier tuning parameter lambda */
    mt->tuningParam[4] = 10.0;             /* distance to move picked branch */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 1.0;
    mt->minimum[1] = 0.0;
    mt->maximum[1] = 0.3;
    mt->minimum[2] = 0.0001;
    mt->maximum[2] = 0.5;
    mt->minimum[3] = 0.0001;
    mt->maximum[3] = 20.0;
    mt->minimum[4] = 2.0;
    mt->maximum[4] = 1000.0;
    mt->parsimonyBased = YES;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_FourTaxaOrMore;

    /* Move_ParsSPRClock */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased SPR for clock trees";
    mt->shortName = "ParsSPRClock";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->applicableTo[4] = TOPOLOGY_RCL_UNIFORM;
    mt->applicableTo[5] = TOPOLOGY_RCL_CONSTRAINED;
    mt->applicableTo[6] = TOPOLOGY_RCCL_UNIFORM;
    mt->applicableTo[7] = TOPOLOGY_RCCL_CONSTRAINED;
    mt->nApplicable = 8;
    mt->moveFxn = &Move_ParsSPRClock;
    mt->relProposalProb = 10.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* warp */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 1.0;
    mt->parsimonyBased = YES;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;

    /* Move_ParsSPRClock_Fossil */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased fossil SPR for clock trees";
    mt->shortName = "ParsSPRClockFossil";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->applicableTo[0] = TOPOLOGY_CL_UNIFORM;
    mt->applicableTo[1] = TOPOLOGY_CCL_UNIFORM;
    mt->applicableTo[2] = TOPOLOGY_CL_CONSTRAINED;
    mt->applicableTo[3] = TOPOLOGY_CCL_CONSTRAINED;
    mt->applicableTo[4] = TOPOLOGY_RCL_UNIFORM;
    mt->applicableTo[5] = TOPOLOGY_RCL_CONSTRAINED;
    mt->applicableTo[6] = TOPOLOGY_RCCL_UNIFORM;
    mt->applicableTo[7] = TOPOLOGY_RCCL_CONSTRAINED;
    mt->nApplicable = 8;
    mt->moveFxn = &Move_ParsSPRClock_Fossil;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* warp */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 1.0;
    mt->parsimonyBased = YES;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_ThreeTaxaOrMore;
    
    /* Move_ParsTBR1 e^{-S} */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased TBR variant 1";
    mt->shortName = "ParsTBR";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->tuningName[1] = "reweighting probability";
    mt->shortTuningName[1] = "r";
    mt->tuningName[2] = "multiplier tuning parameter";
    mt->shortTuningName[2] = "lambda";
    mt->tuningName[3] = "moving distance";
    mt->shortTuningName[3] = "d";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ParsTBR1;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 4;
    mt->tuningParam[0] = 0.5;              /* warp */
    mt->tuningParam[1] = 0.05;             /* upweight and downweight probability */
    mt->tuningParam[2] = 2.0 * log (1.05); /* multiplier tuning parameter lambda */
    mt->tuningParam[3] = 5.0;              /* distance to move picked branch */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 5.0;
    mt->minimum[1] = 0.0;
    mt->maximum[1] = 0.3;
    mt->minimum[2] = 0.0001;
    mt->maximum[2] = 20.0;
    mt->minimum[3] = 2.0;
    mt->maximum[3] = 1000.0;
    mt->parsimonyBased = YES;
    mt->level = STANDARD_USER;
    mt->isApplicable = &IsApplicable_FiveTaxaOrMore;

    /* Move_ParsTBR2 S/N */
    mt = &moveTypes[i++];
    mt->name = "Parsimony-biased TBR variant 2";
    mt->shortName = "ParsTBR2";
    mt->subParams = YES;
    mt->tuningName[0] = "parsimony warp factor";
    mt->shortTuningName[0] = "warp";
    mt->tuningName[1] = "reweighting probability";
    mt->shortTuningName[1] = "r";
    mt->tuningName[2] = "typical branch length";
    mt->shortTuningName[2] = "v_t";
    mt->tuningName[3] = "multiplier tuning parameter";
    mt->shortTuningName[3] = "lambda";
    mt->tuningName[4] = "moving distance";
    mt->shortTuningName[4] = "d";
    mt->applicableTo[0] = TOPOLOGY_NCL_UNIFORM_HOMO;
    mt->applicableTo[1] = TOPOLOGY_NCL_CONSTRAINED_HOMO;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_ParsTBR2;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 5;
    mt->tuningParam[0] = 0.1;              /* warp */
    mt->tuningParam[1] = 0.05;             /* upweight and downweight probability */
    mt->tuningParam[2] = 0.05;             /* typical branch length */
    mt->tuningParam[3] = 2.0 * log (1.05); /* multiplier tuning parameter lambda */
    mt->tuningParam[4] = 5.0;              /* distance to move picked branch */
    mt->minimum[0] = 0.0;
    mt->maximum[0] = 1.0;
    mt->minimum[1] = 0.0;
    mt->maximum[1] = 0.3;
    mt->minimum[2] = 0.0001;
    mt->maximum[2] = 0.5;
    mt->minimum[3] = 0.0001;
    mt->maximum[3] = 20.0;
    mt->minimum[4] = 2.0;
    mt->maximum[4] = 1000.0;
    mt->parsimonyBased = YES;
    mt->level = DEVELOPER;
    mt->isApplicable = &IsApplicable_FiveTaxaOrMore;

    /* Move_Pinvar */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = PINVAR_UNI;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_Pinvar;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* window size */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 0.999;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_Popsize_M */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = POPSIZE_UNI;
    mt->applicableTo[1] = POPSIZE_LOGNORMAL;
    mt->applicableTo[2] = POPSIZE_NORMAL;
    mt->applicableTo[3] = POPSIZE_GAMMA;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_PopSize_M;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log(1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_RateMult_Dir */
    mt = &moveTypes[i++];
    mt->name = "Dirichlet proposal";
    mt->shortName = "Dirichlet";
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->applicableTo[0] = RATEMULT_DIR;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_RateMult_Dir;
    mt->relProposalProb = 0.75;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 50.0; /* alphaPi per site */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 10000.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25;

    /* Move_RateMult_Slider */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = RATEMULT_DIR;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_RateMult_Slider;
    mt->relProposalProb = 0.75;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.05;  /* window size */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_RateShape_M */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = SHAPE_UNI;
    mt->applicableTo[1] = SHAPE_EXP;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_RateShape_M;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 10.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = SPECRATE_UNI;
    mt->applicableTo[1] = SPECRATE_EXP;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_Speciation;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.1;  /* window size */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_Speciation_M */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = SPECRATE_UNI;
    mt->applicableTo[1] = SPECRATE_EXP;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_Speciation_M;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_SpeciesTree removed: species tree not supported in yvyra */

    /* Move_Statefreqs */
    mt = &moveTypes[i++];
    mt->name = "Dirichlet proposal";
    mt->shortName = "Dirichlet";
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->applicableTo[0] = PI_DIR;
    mt->applicableTo[1] = DIRPI_DIRxDIR;
    mt->applicableTo[2] = DIRPI_DIRxFIXED;
    mt->applicableTo[3] = DIRPI_MIX;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_Statefreqs;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 100.0; /* alphaPi per state */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 10000.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25;

    /* Move_Statefreqs_Slider */
    mt = &moveTypes[i++];
    mt->name = "Sliding window";
    mt->shortName = "Slider";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = PI_DIR;
    mt->applicableTo[1] = DIRPI_DIRxDIR;
    mt->applicableTo[2] = DIRPI_DIRxFIXED;
    mt->applicableTo[3] = DIRPI_MIX;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_Statefreqs_Slider;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.2;  /* window size (change in proportions) */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 1.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_StatefreqsRoot */
    mt = &moveTypes[i++];
    mt->name = "Dirichlet proposal (root)";
    mt->shortName = "Dirichlet_root";
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->applicableTo[0] = DIRPI_DIRxDIR;
    mt->applicableTo[1] = DIRPI_FIXEDxDIR;
    mt->applicableTo[2] = DIRPI_MIX;
    mt->nApplicable = 3;
    mt->moveFxn = &Move_StatefreqsRoot;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 100.0; /* alphaPi per state */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 10000.0;
    mt->parsimonyBased = NO; 
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25; 

    /* Move_StatefreqsRoot_Slider */
    mt = &moveTypes[i++];
    mt->name = "Sliding window (root)";
    mt->shortName = "Slider_root";
    mt->tuningName[0] = "Sliding window size";
    mt->shortTuningName[0] = "delta";
    mt->applicableTo[0] = DIRPI_DIRxDIR;
    mt->applicableTo[1] = DIRPI_FIXEDxDIR;
    mt->applicableTo[2] = DIRPI_MIX;
    mt->nApplicable = 3;
    mt->moveFxn = &Move_StatefreqsRoot_Slider;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 0.2;  /* window size (change in proportions) */
    mt->minimum[0] = 0.00001;
    mt->maximum[0] = 1.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneSlider;
    mt->targetRate = 0.25;

    /* Move_Statefreqs_SplitMerge */
    mt = &moveTypes[i++];
    mt->name = "RJ between stationary and directional model";
    mt->shortName = "RJ_Stat-Dir";
    mt->paramName = "Statefrmod";
    mt->applicableTo[0] = DIRPI_MIX;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_Statefreqs_SplitMerge;
    mt->relProposalProb = 0.5;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 100.0; /* alphaPi per state */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 10000.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25;

    /* Move_StatefreqsSymDirMultistate */
    mt = &moveTypes[i++];
    mt->name = "Dirichlet proposal";
    mt->shortName = "Dirichlet";
    mt->paramName = "Pi_symdir";
    mt->tuningName[0] = "Dirichlet parameter";
    mt->shortTuningName[0] = "alpha";
    mt->applicableTo[0] = SYMPI_FIX_MS;
    mt->applicableTo[1] = SYMPI_UNI_MS;
    mt->applicableTo[2] = SYMPI_EXP_MS;
    mt->nApplicable = 3;
    mt->moveFxn = &Move_StatefreqsSymDirMultistate;
    mt->relProposalProb = 5.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 50.0; /* alphaPi */
    mt->minimum[0] = 0.001;
    mt->maximum[0] = 10000.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneDirichlet;
    mt->targetRate = 0.25;
    mt = &moveTypes[i++];
    mt->name = "Tree stretch";
    mt->shortName = "TreeStretch";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = BRLENS_CLOCK_UNI;
    mt->applicableTo[1] = BRLENS_CLOCK_BD;
    mt->applicableTo[2] = BRLENS_CLOCK_COAL;
    mt->applicableTo[3] = BRLENS_CLOCK_FOSSIL;
    mt->nApplicable = 4;
    mt->moveFxn = &Move_TreeStretch;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log(1.01); /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 2.0 * log(2.0);
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_TreeLen, by Jeremy Brown */
    mt = &moveTypes[i++];
    mt->name = "Whole treelength hit with multiplier";
    mt->shortName = "TLMultiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";  
    mt->applicableTo[0] = BRLENS_UNI;
    mt->applicableTo[1] = BRLENS_EXP;
    mt->applicableTo[2] = BRLENS_GamDir;
    mt->applicableTo[3] = BRLENS_iGmDir;
    mt->applicableTo[4] = BRLENS_twoExp;
    mt->nApplicable = 5;  // was 2
    mt->moveFxn = &Move_TreeLen;
    mt->relProposalProb = 3.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (2.0);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_AddDeleteCPPEvent */
    mt = &moveTypes[i++];
    mt->name = "Random addition/deletion of CPP event";
    mt->shortName = "Add_delete";
    mt->applicableTo[0] = CPPEVENTS;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_AddDeleteCPPEvent;
    mt->relProposalProb = 1.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;

    /* Move_CPPEventPosition */
    mt = &moveTypes[i++];
    mt->name = "Random draw of CPP event position from prior";
    mt->shortName = "Prior_draw_pos";
    mt->applicableTo[0] = CPPEVENTS;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_CPPEventPosition;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;

    /* Move_CPPRate */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = CPPRATE_EXP;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_CPPRate;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_CPPRateMultiplier_M */
    mt = &moveTypes[i++];
    mt->name = "Random CPP rate multiplier hit with multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = CPPEVENTS;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_CPPRateMultiplier_M;
    mt->relProposalProb = 0.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_CPPRateMultiplierRnd */
    mt = &moveTypes[i++];
    mt->name = "Random draw of CPP rate multiplier from prior";
    mt->shortName = "Prior_draw_mult";
    mt->applicableTo[0] = CPPEVENTS;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_CPPRateMultiplierRnd;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;

    /* Move_Nu */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = TK02VAR_EXP;
    mt->applicableTo[1] = TK02VAR_UNI;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_Nu;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_TK02BranchRate */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = TK02BRANCHRATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_TK02BranchRate;
    mt->relProposalProb = 25.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_WNVar */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = WNVAR_EXP;
    mt->applicableTo[1] = WNVAR_UNI;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_WNVar;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_WNBranchRate */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = WNBRANCHRATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_WNBranchRate;
    mt->relProposalProb = 25.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_IgrVar */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = IGRVAR_EXP;
    mt->applicableTo[1] = IGRVAR_UNI;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_IgrVar;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_IgrBranchRate */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = IGRBRANCHRATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_IgrBranchRate;
    mt->relProposalProb = 25.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_IlnVar */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = ILNVAR_EXP;
    mt->applicableTo[1] = ILNVAR_UNI;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_IlnVar;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;
    
    /* Move_IlnBranchRate */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = ILNBRANCHRATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_IlnBranchRate;
    mt->relProposalProb = 25.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

   /* Move_MixedVar */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = MIXEDVAR_EXP;
    mt->applicableTo[1] = MIXEDVAR_UNI;
    mt->nApplicable = 2;
    mt->moveFxn = &Move_MixedVar;
    mt->relProposalProb = 2.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.5);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;                 /* smaller */
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;
    
    /* Move_MixedBranchRate */
    mt = &moveTypes[i++];
    mt->name = "Multiplier";
    mt->shortName = "Multiplier";
    mt->tuningName[0] = "Multiplier tuning parameter";
    mt->shortTuningName[0] = "lambda";
    mt->applicableTo[0] = MIXEDBRCHRATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_MixedBranchRate;
    mt->relProposalProb = 25.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0 * log (1.1);  /* lambda */
    mt->minimum[0] = 0.0001;
    mt->maximum[0] = 20.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneMultiplier;
    mt->targetRate = 0.25;

    /* Move_RelaxedClockModel */
    mt = &moveTypes[i++];
    mt->name = "rjMCMC among relaxed clock models";
    mt->shortName = "RJ_Clocks";
    mt->tuningName[0] = "Ratio between ILN and IGR variances";
    mt->shortTuningName[0] = "ratio";  // w = sigma_L/sigma_G
    mt->applicableTo[0] = MIXEDBRCHRATES;
    mt->nApplicable = 1;
    mt->moveFxn = &Move_RelaxedClockModel;
    mt->relProposalProb = 30.0;
    mt->numTuningParams = 1;
    mt->tuningParam[0] = 2.0;
    mt->minimum[0] = 0.01;
    mt->maximum[0] = 100.0;
    mt->parsimonyBased = NO;
    mt->level = STANDARD_USER;
    mt->Autotune = &AutotuneRJClocks;
    mt->targetRate = 1.0;
    
    numMoveTypes = i;
    
    assert(numMoveTypes < NUM_MOVE_TYPES);
}


/* ShowModel: Display model on screen */
int ShowModel (void)
{
    int         i, j, ns;

    YvyraPrint ("%s   Model settings:\n\n", spacer);
    for (i=0; i<numCurrentDivisions; i++)
        {
        ns = 0;

        if (numCurrentDivisions > 1)
            YvyraPrint ("%s      Settings for partition %d --\n", spacer, i+1);
        else
            YvyraPrint ("%s      Data not partitioned --\n", spacer);
        
        if (modelParams[i].dataType == DNA)
            {
            YvyraPrint ("%s         Datatype  = DNA\n", spacer);
            ns = 4;
            }
        else if (modelParams[i].dataType == RNA)
            {
            YvyraPrint ("%s         Datatype  = RNA\n", spacer);
            ns = 4;
            }
        else if (modelParams[i].dataType == PROTEIN)
            {
            YvyraPrint ("%s         Datatype  = Protein\n", spacer);
            ns = 20;
            }
        else if (modelParams[i].dataType == RESTRICTION)
            {
            YvyraPrint ("%s         Datatype  = Restriction\n", spacer);
            ns = 2;
            }
        else if (modelParams[i].dataType == STANDARD)
            {
            YvyraPrint ("%s         Datatype  = Standard\n", spacer);
            ns = MAX_STD_STATES;
            }
        else if (modelParams[i].dataType == CONTINUOUS)
            {
            YvyraPrint ("%s         Datatype  = Continuous\n", spacer);
            }
            
        if (modelSettings[i].dataType == CONTINUOUS)
            {
            /* begin description of continuous models */
              if (!strcmp(modelParams[i].brownCorrPr, "Fixed") && AreDoublesEqual(modelParams[i].brownCorrFix, 0.0, ETA)==YES)
                YvyraPrint ("%s         Model     = Independent Brownian motion\n", spacer);
            else
                YvyraPrint ("%s         Model     = Correlated Brownian motion\n", spacer);
            /* end description of continuous models */
            }
        else
            {
            /* begin description of discrete models */
            if (!strcmp(modelParams[i].parsModel, "Yes"))
                {
                YvyraPrint ("%s         Parsmodel = %s\n", spacer, modelParams[i].parsModel);
                }
            else
                {
                /* dna characters in this partition */
                if (modelSettings[i].dataType == DNA || modelSettings[i].dataType == RNA)
                    {
                    /* general form of the rate matrix */ 
                    YvyraPrint ("%s         Nucmodel  = %s\n", spacer, modelParams[i].nucModel);
                
                    /* constraints on rates of substitution */
                    YvyraPrint ("%s         Nst       = %s\n", spacer, modelParams[i].nst);
                    if (!strcmp(modelParams[i].nst, "2"))
                        {
                        if (!strcmp(modelParams[i].tRatioPr,"Beta"))
                            {
                            YvyraPrint ("%s                     Transition and transversion  rates, expressed\n", spacer);
                            YvyraPrint ("%s                     as proportions of the rate sum, have a\n", spacer);
                            YvyraPrint ("%s                     Beta(%1.2lf,%1.2lf) prior\n", spacer, modelParams[i].tRatioDir[0], modelParams[i].tRatioDir[1]);
                            }
                        else
                            {
                            YvyraPrint ("%s                     Transition/transversion rate ratio is fixed to %1.2lf.\n", spacer, modelParams[i].tRatioFix);
                            }
                        }
                    else if (!strcmp(modelParams[i].nst, "6"))
                        {
                        if (!strcmp(modelParams[i].revMatPr,"Dirichlet"))
                            {
                            YvyraPrint ("%s                     Substitution rates, expressed as proportions\n", spacer);
                            YvyraPrint ("%s                     of the rate sum, have a Dirichlet prior\n", spacer);
                            YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
                                modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
                                modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
                            }
                        else
                            {
                            YvyraPrint ("%s                     Substitution rates are fixed to be \n", spacer);
                            YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf).\n", spacer,
                                modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2],
                                modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);
                            }
                        }
                    else if (!strcmp(modelParams[i].nst, "Mixed"))
                        {
                        if (!strcmp(modelParams[i].revMatPr,"Dirichlet"))
                            {
                            YvyraPrint ("%s                     Substitution rates, expressed as proportions\n", spacer);
                            YvyraPrint ("%s                     of the rate sum, have a Dirichlet prior\n", spacer);
                            YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
                                modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
                                modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
                            }
                        else
                            {
                            YvyraPrint ("%s                     Substitution rates are fixed to be \n", spacer);
                            YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf).\n", spacer,
                                modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2],
                                modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);
                            }
                        }
                    
                    if (!strcmp(modelParams[i].nucModel,"Codon"))
                        {
                        /* what is the distribution on the nonsyn./syn. rate ratio */
                        if (!strcmp(modelParams[i].omegaVar, "Equal"))
                            {
                            if (!strcmp(modelParams[i].omegaPr,"Dirichlet"))
                                {
                                YvyraPrint ("%s                     Nonsynonymous and synonymous rates, expressed\n", spacer);
                                YvyraPrint ("%s                     as proportions of the rate sum, have a\n", spacer);
                                YvyraPrint ("%s                     Dirichlet(%1.2lf,%1.2lf) prior\n", spacer, modelParams[i].omegaDir[0], modelParams[i].omegaDir[1]);
                                }
                            else
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio is fixed to %1.2lf.\n", spacer, modelParams[i].omegaFix);
                                }
                            }
                        else if (!strcmp(modelParams[i].omegaVar, "Ny98"))
                            {
                            if (!strcmp(modelParams[i].ny98omega1pr, "Beta"))
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for purifying selection\n", spacer);
                                YvyraPrint ("%s                     (class 1) has a Beta(%1.2lf,%1.2lf) on the interval (0,1).\n", spacer, modelParams[i].ny98omega1Beta[0], modelParams[i].ny98omega1Beta[1]);
                                }
                            else
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for purifying selection\n", spacer);
                                YvyraPrint ("%s                     (class 1) is fixed to %1.2lf.\n", spacer, modelParams[i].ny98omega1Fixed);
                                }
                            YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for neutral selection\n", spacer);
                            YvyraPrint ("%s                     (class 2) is fixed to 1.0.\n", spacer);
                            if (!strcmp(modelParams[i].ny98omega3pr, "Uniform"))
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for positive selection\n", spacer);
                                YvyraPrint ("%s                     is uniformly distributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].ny98omega3Uni[0], modelParams[i].ny98omega3Uni[1]);
                                }
                            else if (!strcmp(modelParams[i].ny98omega3pr, "Exponential"))
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for positive selection\n", spacer);
                                YvyraPrint ("%s                     is exponentially distributed with parameter (%1.2lf).\n", spacer, modelParams[i].ny98omega3Exp);
                                }
                            else
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for positive \n", spacer);
                                YvyraPrint ("%s                     selection is fixed to %1.2lf.\n", spacer, modelParams[i].ny98omega3Fixed);
                                }
                            }
                        else if (!strcmp(modelParams[i].omegaVar, "M3"))
                            {
                            if (!strcmp(modelParams[i].m3omegapr, "Exponential"))
                                {
                                YvyraPrint ("%s                     Nonsynonymous and synonymous rates for the tree classes of\n", spacer);
                                YvyraPrint ("%s                     omega are exponentially distributed random variables.\n", spacer);
                                }
                            else if (!strcmp(modelParams[i].m3omegapr, "Fixed"))
                                {
                                YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for the three omega\n", spacer);
                                YvyraPrint ("%s                     are fixed to %1.2lf, %1.2lf, and %1.2lf.\n", spacer, modelParams[i].m3omegaFixed[0], modelParams[i].m3omegaFixed[1], modelParams[i].m3omegaFixed[2]);
                                }
                            }
                        else if (!strcmp(modelParams[i].omegaVar, "M10"))
                            {
                            YvyraPrint ("%s                     Nonsynonymous/synonymous rate ratio for purifying \n", spacer);
                            YvyraPrint ("%s                     selection (class 1) has a Beta(alpha1,beta1) on the \n", spacer);
                            YvyraPrint ("%s                     interval (0,1). Nonsynonymous/synonymous rate ratio \n", spacer);
                            YvyraPrint ("%s                     for positive selection (class 2) has an offset \n", spacer);
                            YvyraPrint ("%s                     Gamma(alpha2,beta2) on the interval (1,Infinity).\n", spacer);
                            }
                            
                        /* genetic code that is used (if nucmodel=codon) */
                        YvyraPrint ("%s         Code      = %s\n", spacer, modelParams[i].geneticCode);
                        }
                    }
                /* amino acid characters in this partition */
                else if (modelSettings[i].dataType == PROTEIN)
                    {
                    if (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA)
                        YvyraPrint ("%s         Nucmodel  = %s\n", spacer, modelParams[i].nucModel);
                    /* constraints on rates of substitution in 20 X 20 matrix */
                    if (!strcmp(modelParams[i].aaModelPr, "Mixed"))
                        YvyraPrint ("%s         Aamodel   = Mixture of models with fixed rate matrices\n", spacer);
                    else
                        YvyraPrint ("%s         Aamodel   = %s\n", spacer, modelParams[i].aaModel);
                    /* revmat rates */
                    if (!strcmp(modelParams[i].aaModelPr, "Mixed"))
                        YvyraPrint ("%s                     Substitution rates come from the mixture of models\n", spacer);
                    else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && (!strcmp(modelParams[i].aaModel, "Poisson") ||
                                !strcmp(modelParams[i].aaModel, "Equalin")))
                        YvyraPrint ("%s                     Substitution rates are fixed to be equal\n", spacer);
                    else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && strcmp(modelParams[i].aaModel, "Gtr")!=0)
                        YvyraPrint ("%s                     Substitution rates are fixed to the %s rates\n", spacer, modelParams[i].aaModel);
                    else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && !strcmp(modelParams[i].aaModel, "Gtr"))
                        {
                        if (!strcmp(modelParams[i].aaRevMatPr,"Dirichlet"))
                            {
                            for (j=0; j<190; j++)
                                if (AreDoublesEqual(modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[j], 0.00001) == NO)
                                    break;
                            if (j == 190)
                                {
                                YvyraPrint ("%s                     Substitution rates have a Dirichlet(%1.2lf,%1.2lf,...) prior\n",
                                    spacer, modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[0]);
                                }
                            else
                                {
                                YvyraPrint ("%s                     Substitution rates have a Dirichlet(\n", spacer);
                                for (j=0; j<190; j++)
                                    {
                                    if (j % 10 == 0)
                                        YvyraPrint ("%s                        ", spacer);
                                    YvyraPrint ("%1.2lf", modelParams[i].aaRevMatDir[j]);
                                    if (j == 189)
                                        YvyraPrint (") prior\n");
                                    else if ((j+1) % 10 == 0)
                                        YvyraPrint (",\n");
                                    else
                                        YvyraPrint (",");
                                    }
                                }
                            }
                        else /* if (!strcmp(modelParams[i].aaRevMatPr,"Fixed")) */
                            {
                            for (j=0; j<190; j++)
                                if (AreDoublesEqual(modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[j], 0.00001) == NO)
                                    break;
                            if (j == 190)
                                {
                                YvyraPrint ("%s                     Substitution rates are fixed to (%1.1le,%1.1le,...)\n",
                                    spacer, modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[0]);
                                }
                            else
                                {
                                YvyraPrint ("%s                     Substitution rates are fixed to (\n", spacer);
                                for (j=0; j<190; j++)
                                    {
                                    if (j % 10 == 0)
                                        YvyraPrint ("%s                        ", spacer);
                                    YvyraPrint ("%1.1le", modelParams[i].aaRevMatFix[j]);
                                    if (j == 189)
                                        YvyraPrint (")\n");
                                    else if ((j+1) % 10 == 0)
                                        YvyraPrint (",\n");
                                    else
                                        YvyraPrint (",");
                                    }
                                }
                            }
                        }
                    }
                /* restriction site or morphological characters in this partition */
                else if (modelSettings[i].dataType == RESTRICTION || modelSettings[i].dataType == STANDARD)
                    {
                    /* what type of characters are sampled? */
                    YvyraPrint ("%s         Coding    = %s\n", spacer, modelParams[i].codingString);
                    }
                    
                /* is there rate variation in a single site across the tree? */
                if (((modelSettings[i].dataType == DNA || modelSettings[i].dataType == RNA) && !strcmp(modelParams[i].nucModel, "4by4")) || modelSettings[i].dataType == PROTEIN)
                    {
                    /* do rates change on tree accoding to covarion model? */
                    YvyraPrint ("%s         Covarion  = %s\n", spacer, modelParams[i].covarionModel);
                    if (!strcmp(modelParams[i].covarionModel, "Yes"))
                        {
                        /* distribution on switching parameters, if appropriate */
                        if (!strcmp(modelParams[i].covSwitchPr,"Uniform"))
                            {
                            YvyraPrint ("%s                     Switching rates have independent uniform dist-\n", spacer);
                            YvyraPrint ("%s                     ributions on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1]);
                            }
                        else if (!strcmp(modelParams[i].covSwitchPr,"Exponential"))
                            {
                            YvyraPrint ("%s                     Switching rates have independent exponential\n", spacer);
                            YvyraPrint ("%s                     distributions with parameters (%1.2lf).\n", spacer, modelParams[i].covswitchExp);
                            }
                        else
                            {
                            YvyraPrint ("%s                     Switching rates are fixed to %1.2lf and %1.2lf.\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[0]);
                            }
                        ns *= 2;
                        }
                    }

                /* now, let's deal with variation in omega */
                if ((modelParams[i].dataType == DNA || modelParams[i].dataType == RNA) && !strcmp(modelParams[i].nucModel,"Codon"))
                    {
                    YvyraPrint ("%s         Omegavar  = %s\n", spacer, modelParams[i].omegaVar);
                    if (!strcmp(modelParams[i].geneticCode, "Universal"))
                        ns = 61;
                    else if (!strcmp(modelParams[i].geneticCode, "Vertmt"))
                        ns = 60;
                    else if (!strcmp(modelParams[i].geneticCode, "Invermt"))
                        ns = 62;
                    else if (!strcmp(modelParams[i].geneticCode, "Mycoplasma"))
                        ns = 62;
                    else if (!strcmp(modelParams[i].geneticCode, "Yeast"))
                        ns = 62;
                    else if (!strcmp(modelParams[i].geneticCode, "Ciliate"))
                        ns = 63;
                    else if (!strcmp(modelParams[i].geneticCode, "Echinoderm"))
                        ns = 62;
                    else if (!strcmp(modelParams[i].geneticCode, "Euplotid"))
                        ns = 62;
                    else if (!strcmp(modelParams[i].geneticCode, "Metmt"))
                        ns = 62;
                    }

                /* what assumptions are made about the state frequencies? */
                if (modelParams[i].dataType != CONTINUOUS)
                    {
                    if (modelParams[i].dataType == STANDARD)
                        YvyraPrint ("%s         # States  = Variable, up to %d\n", spacer, MAX_STD_STATES);
                    else if (modelSettings[i].numStates != modelSettings[i].numModelStates)
                        YvyraPrint ("%s         # States  = %d (in the model)\n", spacer, modelSettings[i].numModelStates);
                    else
                        YvyraPrint ("%s         # States  = %d\n", spacer, ns);
                    if (modelSettings[i].dataType == STANDARD)
                        {
                        if (!strcmp(modelParams[i].symPiPr,"Fixed"))
                            {
                            if (AreDoublesEqual(modelParams[i].symBetaFix, -1.0, ETA)==YES)
                                YvyraPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
                            else
                                YvyraPrint ("%s                     Symmetric Dirichlet alpha is fixed to %1.2lf\n", spacer, modelParams[i].symBetaFix);
                            }
                        else if (!strcmp(modelParams[i].symPiPr,"Uniform"))
                            {
                            YvyraPrint ("%s                     Symmetric Dirichlet alpha has a Uniform(%1.2lf,%1.2lf) prior\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1]);
                            }
                        else
                            {
                            YvyraPrint ("%s                     Symmetric Dirichlet alpha has a Exponential(%1.2lf) prior\n", spacer, modelParams[i].symBetaExp);
                            }
                        }
                    else if (modelSettings[i].dataType == RESTRICTION)
                        {
                        /* distribution on state frequencies for restriction site model */
                        if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
                            {
                            YvyraPrint ("%s                     State frequencies have a Dirichlet (%1.2lf,%1.2lf) prior\n", spacer,
                                modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1]);
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
                            {
                            YvyraPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
                            {
                            YvyraPrint ("%s                     State frequencies are fixed(%1.2lf,%1.2lf)\n", spacer,
                                modelParams[i].stateFreqsFix[0], modelParams[i].stateFreqsFix[1]);
                            YvyraPrint ("%s                     State frequencies have been fixed by the user\n", spacer);
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
                            {
                            YvyraPrint ("%s                     State frequencies are fixed(%1.2lf,%1.2lf)\n", spacer,
                                modelParams[i].stateFreqsFix[0], modelParams[i].stateFreqsFix[1]);
                            YvyraPrint ("%s                     State frequencies have been fixed to the empirical frequencies in the data\n", spacer);
                            }
                        if (!strcmp(modelParams[i].statefreqModel,"Directional") || !strcmp(modelParams[i].statefreqModel,"Mixed"))
                            {
                            YvyraPrint ("%s                     State frequencies are potentially different for the root (directional model)\n", spacer);
                            if (!strcmp(modelParams[i].statefreqModel,"Directional") && !strcmp(modelParams[i].rootFreqPr,"Fixed"))
                                YvyraPrint ("%s                     Root state frequencies are fixed(%1.2lf,%1.2lf)\n", spacer,
                                    modelParams[i].rootFreqsFix[0], modelParams[i].rootFreqsFix[1]);
                            else if (!strcmp(modelParams[i].statefreqModel,"Directional") && !strcmp(modelParams[i].rootFreqPr,"Dirichlet"))
                                YvyraPrint ("%s                     Root state frequencies have a Dirichlet (%1.2lf,%1.2lf) prior\n", spacer,
                                    modelParams[i].rootFreqsDir[0], modelParams[i].rootFreqsDir[1]);
                            }
                        else if (!strcmp(modelParams[i].statefreqModel,"Mixed"))
                            {
                            YvyraPrint ("%s                     State frequencies are potentially different for the root\n", spacer);
                            }
                        }
                    else if (modelSettings[i].dataType == PROTEIN)
                        {
                        /* distribution on state frequencies for aminoacid model */
                        if (!strcmp(modelParams[i].aaModelPr, "Fixed") && (strcmp(modelParams[i].aaModel, "Equalin")==0 ||
                            strcmp(modelParams[i].aaModel, "Gtr")==0))
                            {
                            if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
                                {
                                YvyraPrint ("%s                     State frequencies have a Dirichlet prior\n", spacer);
                                YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,", spacer,
                                    modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1], modelParams[i].stateFreqsDir[2],
                                    modelParams[i].stateFreqsDir[3], modelParams[i].stateFreqsDir[4]);
                                YvyraPrint ("%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,\n",
                                    modelParams[i].stateFreqsDir[5], modelParams[i].stateFreqsDir[6], modelParams[i].stateFreqsDir[7],
                                    modelParams[i].stateFreqsDir[8], modelParams[i].stateFreqsDir[9]);
                                YvyraPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,", spacer,
                                    modelParams[i].stateFreqsDir[10], modelParams[i].stateFreqsDir[11], modelParams[i].stateFreqsDir[12],
                                    modelParams[i].stateFreqsDir[13], modelParams[i].stateFreqsDir[14]);
                                YvyraPrint ("%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n",
                                    modelParams[i].stateFreqsDir[15], modelParams[i].stateFreqsDir[16], modelParams[i].stateFreqsDir[17],
                                    modelParams[i].stateFreqsDir[18], modelParams[i].stateFreqsDir[19]);
                                }
                            else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
                                {
                                YvyraPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
                                }
                            else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
                                {
                                YvyraPrint ("%s                     State frequencies have been fixed by the user\n", spacer);
                                }
                            else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
                                {
                                YvyraPrint ("%s                     State frequencies have been fixed to the empirical frequencies in the data\n", spacer);
                                }
                            }
                        else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && !strcmp(modelParams[i].aaModel, "Poisson"))
                            {
                            YvyraPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
                            }
                        else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && strcmp(modelParams[i].aaModel, "Equalin") && strcmp(modelParams[i].aaModel, "Poisson"))
                            {
                            YvyraPrint ("%s                     State frequencies are fixed to the %s frequencies\n", spacer, modelParams[i].aaModel);
                            }
                        else
                            {
                            YvyraPrint ("%s                     State frequencies come from the mixture of models\n", spacer);
                            }
                        }
                    else
                        {
                        /* distribution on state frequencies for all other models */
                        if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
                            {
                            YvyraPrint ("%s                     State frequencies have a Dirichlet prior\n", spacer);
                            if (!strcmp(modelParams[i].nucModel, "Doublet"))
                                {
                                YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,\n", spacer,
                                    modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1], modelParams[i].stateFreqsDir[2],
                                    modelParams[i].stateFreqsDir[3]);
                                YvyraPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf,\n", spacer,
                                    modelParams[i].stateFreqsDir[4], modelParams[i].stateFreqsDir[5], modelParams[i].stateFreqsDir[6],
                                    modelParams[i].stateFreqsDir[7]);
                                YvyraPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf,\n", spacer,
                                    modelParams[i].stateFreqsDir[8], modelParams[i].stateFreqsDir[9], modelParams[i].stateFreqsDir[10],
                                    modelParams[i].stateFreqsDir[11]);
                                YvyraPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
                                    modelParams[i].stateFreqsDir[12], modelParams[i].stateFreqsDir[13], modelParams[i].stateFreqsDir[14],
                                    modelParams[i].stateFreqsDir[15]);
                                }
                            else if (!strcmp(modelParams[i].nucModel, "4by4"))
                                {
                                YvyraPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
                                    modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1], modelParams[i].stateFreqsDir[2],
                                    modelParams[i].stateFreqsDir[3]);
                                }
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
                            {
                            YvyraPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
                            {
                            YvyraPrint ("%s                     State frequencies have been fixed by the user\n", spacer);
                            }
                        else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
                            {
                            YvyraPrint ("%s                     State frequencies have been fixed to the empirical frequencies in the data\n", spacer);
                            }
                        }
                    }
                else
                    YvyraPrint ("%s         # States  = Infinity\n", spacer);

                /* now, let's deal with rate variation across sites */
                if (modelSettings[i].dataType != CONTINUOUS)
                    {
                    if (((modelSettings[i].dataType == DNA || modelSettings[i].dataType == RNA) && strcmp(modelParams[i].nucModel,"Codon")!=0) ||
                          modelSettings[i].dataType == PROTEIN || modelSettings[i].dataType == RESTRICTION || modelSettings[i].dataType == STANDARD)
                        {
                        if (!strcmp(modelParams[i].covarionModel, "No"))
                            YvyraPrint ("%s         Rates     = %s\n", spacer, modelParams[i].ratesModel);
                        else
                            {
                            if (!strcmp(modelParams[i].ratesModel, "Propinv"))
                                YvyraPrint ("%s         Rates     = Equal ", spacer);
                            else if (!strcmp(modelParams[i].ratesModel, "Invgamma"))
                                YvyraPrint ("%s         Rates     = Gamma ", spacer);
                            else
                                YvyraPrint ("%s         Rates     = %s ", spacer, modelParams[i].ratesModel);
                            YvyraPrint ("(+ Propinv induced by covarion model)\n");
                            }
                        
                        if ((modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD) && !strcmp(modelParams[i].ratesModel, "Adgamma"))
                            {
                            
                            }
                        else
                            {
                            if (!strcmp(modelParams[i].ratesModel, "Gamma") || !strcmp(modelParams[i].ratesModel, "Invgamma") ||
                                !strcmp(modelParams[i].ratesModel, "LNorm") || !strcmp(modelParams[i].ratesModel, "Adgamma") ||
                                !strcmp(modelParams[i].ratesModel, "Kmixture"))
                                {
                                /* how many categories is the continuous gamma/lnorm approximated by? or how many components are there in the mixture? */
                                if (!strcmp(modelParams[i].ratesModel, "Kmixture"))
                                    YvyraPrint ("%s                     There are %d components in the mixture.\n", spacer, modelParams[i].numMixtCats);
                                else if (!strcmp(modelParams[i].ratesModel, "Lnorm"))
                                    YvyraPrint ("%s                     The distribution is approximated using %d categories.\n", spacer, modelParams[i].numLnormCats);
                                else
                                    YvyraPrint ("%s                     The distribution is approximated using %d categories.\n", spacer, modelParams[i].numGammaCats);
                                /* distribution on shape parameter, if appropriate */
                                if (!strcmp(modelParams[i].shapePr,"Uniform"))
                                    {
                                    YvyraPrint ("%s                     Shape parameter is uniformly distributed\n", spacer);
                                    YvyraPrint ("%s                     on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1]);
                                    }
                                else if (!strcmp(modelParams[i].shapePr,"Exponential"))
                                    {
                                    YvyraPrint ("%s                     Shape parameter is exponentially\n", spacer);
                                    YvyraPrint ("%s                     distributed with parameter (%1.2lf).\n", spacer, modelParams[i].shapeExp);
                                    }
                                else
                                    {
                                    YvyraPrint ("%s                     Shape parameter is fixed to %1.2lf.\n", spacer, modelParams[i].shapeFix);
                                    }
                                }

                            if ((!strcmp(modelParams[i].ratesModel, "Propinv") || !strcmp(modelParams[i].ratesModel, "Invgamma")) && !strcmp(modelParams[i].covarionModel, "No"))
                                {
                                /* distribution on pInvar parameter, if appropriate */
                                if (!strcmp(modelParams[i].pInvarPr,"Uniform"))
                                    {
                                    YvyraPrint ("%s                     Proportion of invariable sites is uniformly dist-\n", spacer);
                                    YvyraPrint ("%s                     ributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1]);
                                    }
                                else
                                    {
                                    YvyraPrint ("%s                     Proportion of invariable sites is fixed to %1.2lf.\n", spacer, modelParams[i].pInvarFix);
                                    }
                                }
                            if (!strcmp(modelParams[i].ratesModel, "Adgamma"))
                                {
                                /* distribution on correlation parameter, if appropriate */
                                if (!strcmp(modelParams[i].adGammaCorPr,"Uniform"))
                                    {
                                    YvyraPrint ("%s                     Rate correlation parameter is uniformly dist-\n", spacer);
                                    YvyraPrint ("%s                     ributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].adgCorrUni[0], modelParams[i].adgCorrUni[1]);
                                    }
                                else
                                    {
                                    YvyraPrint ("%s                     Rate correlation parameter is fixed to %1.2lf.\n", spacer, modelParams[i].adgCorrFix);
                                    }
                                }
                            }
                        }
                    }
                }
            /* end description of discrete models */
            }

        if (i != numCurrentDivisions - 1)
            YvyraPrint ("\n");
        
        }

    YvyraPrint ("\n");
    ShowParameters (NO, NO, NO);
    
    return (NO_ERROR);
}


/*------------------------------------------------------------------------------
|
|   ShowMoves: Show applicable moves
|
------------------------------------------------------------------------------*/
int ShowMoves (int used)
{
    int             i, k, run, chain, chainIndex, areRunsSame, areChainsSame, numPrintedMoves;
    MCMCMove        *mv;

    numPrintedMoves = 0;
    for (i=0; i<numApplicableMoves; i++)
        {
        mv = moves[i];
        
        for (k=0; k<numGlobalChains; k++)
            {
            if (mv->relProposalProb[k] > 0.000001)
                break;
            }

        if (k == numGlobalChains && used == YES)
            continue;

        if (k < numGlobalChains && used == NO)
            continue;

        numPrintedMoves++;
        
        /* print move number and name */
        YvyraPrint ("%s   %4d -- Move        = %s\n", spacer, numPrintedMoves, mv->name);
        
        /* print move type */
        YvyraPrint ("%s           Type        = %s\n", spacer, mv->moveType->name);

        /* print parameter */
        if (mv->parm->nSubParams > 0)
            YvyraPrint ("%s           Parameters  = %s [param. %d] (%s)\n", spacer, mv->parm->name,
                mv->parm->index+1, mv->parm->paramTypeName);
        else
            YvyraPrint ("%s           Parameter   = %s [param. %d] (%s)\n", spacer, mv->parm->name,
                mv->parm->index+1, mv->parm->paramTypeName);
        for (k=0; k<mv->parm->nSubParams; k++)
            YvyraPrint ("%s                         %s [param. %d] (%s)\n", spacer, mv->parm->subParams[k]->name,
                mv->parm->subParams[k]->index+1, mv->parm->subParams[k]->paramTypeName);

        /* print tuning parameters */
        for (k=0; k<mv->moveType->numTuningParams; k++)
            {
            if (k==0)
                YvyraPrint ("%s           Tuningparam = %s (%s)\n", spacer, mv->moveType->shortTuningName[k], mv->moveType->tuningName[k]);
            else
                YvyraPrint ("%s                         %s (%s)\n", spacer, mv->moveType->shortTuningName[k], mv->moveType->tuningName[k]);
            }
        
        /* loop over tuning parameters */
        for (k=0; k<mv->moveType->numTuningParams; k++)
            {
            /* find if tuning parameters are different for different runs */
            areRunsSame = YES;
            for (run=1; run<chainParams.numRuns; run++)
                {
                for (chain=0; chain<chainParams.numChains; chain++)
                    {
                    chainIndex = run*chainParams.numChains + chain;
                    if (AreDoublesEqual (mv->tuningParam[chainIndex][k], mv->tuningParam[chain][k], 0.000001) == NO)
                        {
                        areRunsSame = NO;
                        break;
                        }
                    }
                if (areRunsSame == NO)
                    break;
                }
        
            /* now print values */
            for (run=0; run<chainParams.numRuns; run++)
                {
                if (areRunsSame == YES && run >= 1)
                    break;

                /* find out if chains are different within this run */
                areChainsSame = YES;
                for (chain=1; chain<chainParams.numChains; chain++)
                    {
                    chainIndex = run*chainParams.numChains + chain;
                    if (AreDoublesEqual (mv->tuningParam[chainIndex][k], mv->tuningParam[chainIndex-chain][k],0.000001) == NO)
                        {
                        areChainsSame = NO;
                        break;
                        }
                    }
                /* now we can print the values */
                for (chain=0; chain<chainParams.numChains; chain++)
                    {
                    chainIndex = run*chainParams.numChains + chain;
                    if (areChainsSame == YES && chain >= 1)
                        break;
                    
                    if (run == 0 && chain == 0)
                        YvyraPrint ("%s%22s = %1.3lf", spacer, mv->moveType->shortTuningName[k], mv->tuningParam[chainIndex][k]);
                    else
                        YvyraPrint ("%s                         %1.3lf", spacer, mv->tuningParam[chainIndex][k]);

                    if (areChainsSame == NO && areRunsSame == YES)
                        YvyraPrint ("  [chain %d]\n", chain+1);
                    else if (areChainsSame == YES && areRunsSame == NO)
                        YvyraPrint ("  [run %d]\n", run+1);
                    else if (areChainsSame == NO && areRunsSame == NO)
                        YvyraPrint ("  [run %d, chain %d]\n", run+1, chain+1);
                    else
                        YvyraPrint ("\n");
                    }
                }
            }   /* next tuning parameter */

        /* print target acceptance rate for autotuning */
        if (mv->moveType->targetRate > 0.0 && mv->moveType->targetRate < 1.0)
            {

            /* first find out if the targets are different in different runs */         
            areRunsSame = YES;
            for (run=1; run<chainParams.numRuns; run++)
                {
                for (chain=0; chain<chainParams.numChains; chain++)
                    {
                    chainIndex = run*chainParams.numChains + chain;
                    if (AreDoublesEqual (mv->targetRate[chainIndex], mv->targetRate[chain], 0.000001) == NO)
                        {
                        areRunsSame = NO;
                        break;
                        }
                    }
                if (areRunsSame == NO)
                    break;
                }
        
            /* now print values */
            for (run=0; run<chainParams.numRuns; run++)
                {
                if (areRunsSame == YES && run >= 1)
                    break;

                /* find out if chains are different within this run */
                areChainsSame = YES;
                for (chain=1; chain<chainParams.numChains; chain++)
                    {
                    chainIndex = run*chainParams.numChains + chain;
                    if (AreDoublesEqual (mv->targetRate[chainIndex], mv->targetRate[chainIndex-chain], 0.000001) == NO)
                        {
                        areChainsSame = NO;
                        break;
                        }
                    }
                /* now we can print the values */
                for (chain=0; chain<chainParams.numChains; chain++)
                    {
                    chainIndex = run*chainParams.numChains + chain;
                    if (areChainsSame == YES && chain >= 1)
                        break;
                    
                    if (run == 0 && chain == 0)
                        YvyraPrint ("%s           Targetrate  = %1.3lf", spacer, mv->targetRate[chainIndex]);
                    else
                        YvyraPrint ("%s                         %1.3lf", spacer, mv->targetRate[chainIndex]);

                    if (areChainsSame == NO && areRunsSame == YES)
                        YvyraPrint ("  [chain %d]\n", chain+1);
                    else if (areChainsSame == YES && areRunsSame == NO)
                        YvyraPrint ("  [run %d]\n", run+1);
                    else if (areChainsSame == NO && areRunsSame == NO)
                        YvyraPrint ("  [run %d, chain %d]\n", run+1, chain+1);
                    else
                        YvyraPrint ("\n");
                    }
                }
            }

        
        /* finally print the relative proposal probability */
        
        /* first find out if the probabilities are different in different runs */           
        areRunsSame = YES;
        for (run=1; run<chainParams.numRuns; run++)
            {
            for (chain=0; chain<chainParams.numChains; chain++)
                {
                chainIndex = run*chainParams.numChains + chain;
                if (AreDoublesEqual (mv->relProposalProb[chainIndex], mv->relProposalProb[chain], 0.000001) == NO)
                    {
                    areRunsSame = NO;
                    break;
                    }
                }
            if (areRunsSame == NO)
                break;
            }
    
        /* now print values */
        for (run=0; run<chainParams.numRuns; run++)
            {
            if (areRunsSame == YES && run >= 1)
                break;

            /* find out if chains are different within this run */
            areChainsSame = YES;
            for (chain=1; chain<chainParams.numChains; chain++)
                {
                chainIndex = run*chainParams.numChains + chain;
                if (AreDoublesEqual (mv->relProposalProb[chainIndex], mv->relProposalProb[chainIndex-chain], 0.000001) == NO)
                    {
                    areChainsSame = NO;
                    break;
                    }
                }
            /* now we can print the values */
            for (chain=0; chain<chainParams.numChains; chain++)
                {
                chainIndex = run*chainParams.numChains + chain;
                if (areChainsSame == YES && chain >= 1)
                    break;
                
                if (run == 0 && chain == 0)
                    YvyraPrint ("%s           Rel. prob.  = %1.1lf", spacer, mv->relProposalProb[chainIndex]);
                else
                    YvyraPrint ("%s                         %1.1lf", spacer, mv->relProposalProb[chainIndex]);

                if (areChainsSame == NO && areRunsSame == YES)
                    YvyraPrint ("  [chain %d]\n", chain+1);
                else if (areChainsSame == YES && areRunsSame == NO)
                    YvyraPrint ("  [run %d]\n", run+1);
                else if (areChainsSame == NO && areRunsSame == NO)
                    YvyraPrint ("  [run %d, chain %d]\n", run+1, chain+1);
                else
                    YvyraPrint ("\n");
                }
            }
        YvyraPrint ("\n");
        }   /* next move */
        
    if (numPrintedMoves == 0)
        {
        if (used == YES)
            {
            YvyraPrint ("%s      No moves currently used for this analysis. All parameters\n", spacer);
            YvyraPrint ("%s      will be fixed to their starting values.\n\n", spacer);
            }
        else
            {
            YvyraPrint ("%s      No additional moves available for this model.\n\n", spacer);
            }
        }

    return (NO_ERROR);
}


/*------------------------------------------------------------------------------
|
|   ShowParameters: Show parameter table and parameter info
|
------------------------------------------------------------------------------*/
int ShowParameters (int showStartVals, int showMoves, int showAllAvailable)
{
    int             a, b, d, i, j, k, m, n, run, chain, shouldPrint, isSame, areRunsSame, areChainsSame, nValues,
                    chainIndex, refIndex, numPrinted, numMovedChains, printedCol, screenWidth=100;
    Param           *p;
    Model           *mp;
    ModelInfo       *ms;
    YFlt          *value, *refValue, *subValue;
    MCMCMove        *mv;
    
    YvyraPrint ("%s   Active parameters: \n\n", spacer);
    if (numCurrentDivisions > 1)
        { 
        YvyraPrint ("%s                          Partition(s)\n", spacer);
        YvyraPrint ("%s      Parameters        ", spacer);
        for (i=0; i<numCurrentDivisions; i++)
            YvyraPrint (" %2d", i+1);
        YvyraPrint ("\n");
        YvyraPrint ("%s      ------------------", spacer);
        for (i=0; i<numCurrentDivisions; i++)
            YvyraPrint ("---");
        YvyraPrint ("\n");
        }
    else
        {
        YvyraPrint ("%s      Parameters\n", spacer);
        YvyraPrint ("%s      ---------------------\n", spacer);
        }
    for (j=0; j<NUM_LINKED; j++)
        {
        shouldPrint = NO;
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (activeParams[j][i] == -1)
                {}
            else
                shouldPrint = YES;
            }
        if (shouldPrint == NO)
            continue;
        
        else if (j == P_PI)
            {
            YvyraPrint ("%s      Statefreq         ", spacer);
            }
        else if (j == P_MIXTURE_RATES)
            {
            YvyraPrint ("%s      Mixturerates      ", spacer);
            }
        else if (j == P_SHAPE)
            {
            YvyraPrint ("%s      Shape             ", spacer);
            }
        else if (j == P_PINVAR)
            {
            YvyraPrint ("%s      Pinvar            ", spacer);
            }
        else if (j == P_CORREL)
            {
            YvyraPrint ("%s      Correlation       ", spacer);
            }
        else if (j == P_SWITCH)
            {
            YvyraPrint ("%s      Switchrates       ", spacer);
            }
        else if (j == P_RATEMULT)
            {
            YvyraPrint ("%s      Ratemultiplier    ", spacer);
            }
        else if (j == P_GENETREERATE)
            {
            YvyraPrint ("%s      Generatemult      ", spacer);
            }
        else if (j == P_TOPOLOGY)
            {
            YvyraPrint ("%s      Topology          ", spacer);
            }
        else if (j == P_BRLENS)
            {
            YvyraPrint ("%s      Brlens            ", spacer);
            }
        else if (j == P_SPECRATE)
            {
            YvyraPrint ("%s      Speciationrate    ", spacer);
            }
        else if (j == P_EXTRATE)
            {
            YvyraPrint ("%s      Extinctionrate    ", spacer);
            }
        else if (j == P_FOSLRATE)
            {
            YvyraPrint ("%s      Fossilizationrate ", spacer);
            }
        else if (j == P_POPSIZE)
            {
            YvyraPrint ("%s      Popsize           ", spacer);
            }
        else if (j == P_GROWTH)
            {
            YvyraPrint ("%s      Growthrate        ", spacer);
            } 
        else if (j == P_BMCORR)
            {
            YvyraPrint ("%s      Brownian corr.    ", spacer);
            }
        else if (j == P_BMSIGMA)
            {
            YvyraPrint ("%s      Brownian sigma    ", spacer);
            }
        else if (j == P_CPPRATE)
            {
            YvyraPrint ("%s      Cpprate           ", spacer);
            }
        else if (j == P_CPPMULTDEV)
            {
            YvyraPrint ("%s      Cppmultdev        ", spacer);
            }
        else if (j == P_CPPEVENTS)
            {
            YvyraPrint ("%s      Cppevents         ", spacer);
            }
        else if (j == P_TK02VAR)
            {
            YvyraPrint ("%s      TK02var           ", spacer);
            }
        else if (j == P_TK02BRANCHRATES)
            {
            YvyraPrint ("%s      TK02branchrates   ", spacer);
            }
        else if (j == P_WNVAR)
            {
            YvyraPrint ("%s      WNvar             ", spacer);
            }
        else if (j == P_WNBRANCHRATES)
            {
            YvyraPrint ("%s      WNbranchrates     ", spacer);
            }
        else if (j == P_IGRVAR)
            {
            YvyraPrint ("%s      IGRvar            ", spacer);
            }
        else if (j == P_IGRBRANCHRATES)
            {
            YvyraPrint ("%s      IGRbranchrates    ", spacer);
            }
        else if (j == P_ILNVAR)
            {
            YvyraPrint ("%s      ILNvar            ", spacer);
            }
        else if (j == P_ILNBRANCHRATES)
            {
            YvyraPrint ("%s      ILNbranchrates    ", spacer);
            }
        else if (j == P_MIXEDVAR)
            {
            YvyraPrint ("%s      Mixedvar          ", spacer);
            }
        else if (j == P_MIXEDBRCHRATES)
            {
            YvyraPrint ("%s      Mixedbranchrates  ", spacer);
            }
        else if (j == P_CLOCKRATE)
            {
            YvyraPrint ("%s      Clockrate         ", spacer);
            }
        else if (j == P_SPECIESTREE)
            {
            YvyraPrint ("%s      Speciestree       ", spacer);
            }
        else
            {
            YvyraPrint ("%s      ERROR: Someone forgot to name parameter type %d", spacer, j);
            return (ERROR);
            }
        
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (activeParams[j][i] == -1)
                YvyraPrint ("  .");
            else
                YvyraPrint (" %2d", activeParams[j][i]);
            }
        YvyraPrint ("\n");
        }
    if (numCurrentDivisions > 1)
        { 
        YvyraPrint ("%s      ------------------", spacer);
        for (i=0; i<numCurrentDivisions; i++)
            YvyraPrint ("---");
        YvyraPrint ("\n");
        }
    else
        {
        YvyraPrint ("%s      ---------------------\n", spacer);
        }
    
    YvyraPrint ("\n");
    
    if (numCurrentDivisions > 1)
        YvyraPrint ("%s      Parameters can be linked or unlinked across partitions using 'link' and 'unlink'\n\n", spacer);
    
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        j = p->paramType;
        
        mp = &modelParams[p->relParts[0]];
        ms = &modelSettings[p->relParts[0]];
        
        /* print parameter number and name */
        YvyraPrint ("%s   %4d --  Parameter  = %s\n", spacer, i+1, p->name);
        YvyraPrint ("%s            Type       = %s\n", spacer, p->paramTypeName);
        /* print prior */
        if (j == P_PI)
            {
            if (ms->dataType == STANDARD)
                {
                if (!strcmp(mp->symPiPr, "Uniform"))
                    YvyraPrint ("%s            Prior      = Symmetric dirichlet with uniform(%1.2lf,%1.2lf) variance parameter\n", spacer, mp->symBetaUni[0], mp->symBetaUni[1]);
                else if (!strcmp(mp->symPiPr, "Exponential"))
                    YvyraPrint ("%s            Prior      = Symmetric dirichlet with exponential(%1.2lf) variance parameter\n", spacer, mp->symBetaExp);
                else
                    { /* mp->symBetaFix == -1 */
                    if (AreDoublesEqual(mp->symBetaFix, -1.0, ETA)==YES)
                        YvyraPrint ("%s            Prior      = Symmetric dirichlet with all parameters fixed to infinity\n", spacer);
                    else
                        YvyraPrint ("%s            Prior      = Symmetric dirichlet with all parameters fixed to %1.2lf\n", spacer, mp->symBetaFix);
                    }
                }
            else if (ms->dataType == PROTEIN)
                {
                if (!strcmp(mp->aaModelPr, "Fixed") && (!strcmp(mp->aaModel, "Equalin") || !strcmp(mp->aaModel, "Gtr")))
                    {
                    if (!strcmp(mp->stateFreqPr,"Dirichlet"))
                        {
                        YvyraPrint ("%s            Prior      = Dirichlet\n", spacer);
                        }
                    else if (!strcmp(mp->stateFreqPr,"Fixed") && !strcmp(mp->stateFreqsFixType,"Equal"))
                        {
                        YvyraPrint ("%s            Prior      = Fixed (equal frequencies)\n", spacer);
                        }
                    else if (!strcmp(mp->stateFreqPr,"Fixed") && !strcmp(mp->stateFreqsFixType,"User"))
                        {
                        YvyraPrint ("%s            Prior      = Fixed (user-specified)\n", spacer);
                        }
                    else if (!strcmp(mp->stateFreqPr,"Fixed") && !strcmp(mp->stateFreqsFixType,"Empirical"))
                        {
                        YvyraPrint ("%s            Prior      = Fixed (empirical frequencies)\n", spacer);
                        }
                    }
                else if (!strcmp(mp->aaModelPr, "Fixed") && !strcmp(mp->aaModel, "Poisson"))
                    {
                    YvyraPrint ("%s            Prior      = Fixed (equal frequencies)\n", spacer);
                    }
                else if (!strcmp(mp->aaModelPr, "Fixed") && strcmp(mp->aaModel, "Equalin") && strcmp(mp->aaModel, "Poisson"))
                    {
                    YvyraPrint ("%s            Prior      = Fixed (%s frequencies)\n", spacer, mp->aaModel);
                    }
                else
                    {
                    YvyraPrint ("%s            Prior      = Fixed (from mixture of models)\n", spacer);
                    }
                }
            else
                {
                if (!strcmp(mp->stateFreqPr,"Dirichlet"))
                    YvyraPrint ("%s            Prior      = Dirichlet\n", spacer);
                else
                    YvyraPrint ("%s            Prior      = Fixed\n", spacer);
                }
            }
        else if (j == P_MIXTURE_RATES)
            {
            YvyraPrint ("%s            Prior      = Ordered flat Dirichlet distribution\n", spacer);
            }
        else if (j == P_SHAPE)
            {
            if (!strcmp(mp->shapePr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->shapeUni[0], mp->shapeUni[1]);
            else if (!strcmp(mp->shapePr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->shapeExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->shapeFix);
            }
        else if (j == P_PINVAR)
            {
            if (!strcmp(mp->pInvarPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->pInvarUni[0], mp->pInvarUni[1]);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->pInvarFix);
            }
        else if (j == P_CORREL)
            {
            if (!strcmp(mp->adGammaCorPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->adgCorrUni[0], mp->adgCorrUni[1]);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->adgCorrFix);
            }
        else if (j == P_SWITCH)
            {
            if (!strcmp(mp->covSwitchPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->covswitchUni[0], mp->covswitchUni[1]);
            else if (!strcmp(mp->covSwitchPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->covswitchExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf,%1.2lf)\n", spacer, mp->covswitchFix[0], mp->covswitchFix[1]);
            }
        else if (j == P_RATEMULT)
            {
            if (!strcmp(mp->ratePr,"Fixed"))
                YvyraPrint ("%s            Prior      = Fixed(1.0)\n", spacer);
            else
                {
                YvyraPrint ("%s            Prior      = Dirichlet(", spacer);
                for (d=n=0; d<numCurrentDivisions; d++)
                    {
                    if (activeParams[j][d] == i+1)
                        n++;
                    }
                for (d=m=0; d<numCurrentDivisions; d++)
                    {
                    if (activeParams[j][d] == i+1)
                        {
                        m++;
                        if (m < n)
                            YvyraPrint ("%1.2lf,", modelParams[d].ratePrDir);
                        else
                            YvyraPrint ("%1.2lf)\n", modelParams[d].ratePrDir);
                        }
                    }
                }
            }
        else if (j == P_GENETREERATE)
            {
            if (!strcmp(mp->generatePr,"Fixed"))
                YvyraPrint ("%s            Prior      = Fixed(1.0)\n", spacer);
            else
                {
                YvyraPrint ("%s            Prior      = Dirichlet(", spacer);
                printedCol = (int)(strlen(spacer)) + 25 + 10;
                for (n=0; n<numTrees-1; n++)
                    {
                    if (printedCol + 5 > screenWidth)
                        {
                        YvyraPrint ("\n%s                                   ", spacer);
                        printedCol = (int)(strlen(spacer)) + 25 + 10;
                        }
                    if (n == numTrees-2)
                        YvyraPrint ("1.00)\n");
                    else
                        YvyraPrint ("1.00,");
                    printedCol += 5;
                    }
                }
            }
        else if (j == P_TOPOLOGY)
            {
            if (!strcmp(mp->topologyPr,"Uniform"))
                YvyraPrint ("%s            Prior      = All topologies equally probable a priori\n", spacer);
            else if (!strcmp(mp->topologyPr,"Speciestree"))
                YvyraPrint ("%s            Prior      = Topology constrained to fold within species tree\n", spacer);
            else if (!strcmp(mp->topologyPr,"Constraints"))
                YvyraPrint ("%s            Prior      = Prior on topology obeys the following constraints:\n", spacer);
            else
                YvyraPrint ("%s            Prior      = Prior is fixed to the topology of user tree '%s'\n", spacer,
                                                    userTree[mp->topologyFix]->name);
            if (!strcmp(mp->topologyPr,"Constraints"))
                {
                for (a=0; a<numDefinedConstraints; ++a)
                    {
                    if (mp->activeConstraints[a] == NO)
                        continue;
                    if (definedConstraintsType[a] == HARD)
                        YvyraPrint ("%s                         -- Hard constraint \"%s\"\n", spacer, constraintNames[a]);
                    else if (definedConstraintsType[a] == PARTIAL)
                        YvyraPrint ("%s                         -- Partial constraint \"%s\"\n", spacer, constraintNames[a]);
                    else /* if (true) */
                        YvyraPrint ("%s                         -- Negative constraint \"%s\"\n", spacer, constraintNames[a]);
                    }
                }
            }
        else if (j == P_BRLENS)
            {
            if (!strcmp(mp->parsModel, "Yes"))
                YvyraPrint ("%s            Prior      = Reconstructed using parsimony\n", spacer);
            else
                {
                if (!strcmp(mp->brlensPr,"Unconstrained"))
                    {
                    YvyraPrint ("%s            Prior      = Unconstrained:%s", spacer, mp->unconstrainedPr);
                    if (!strcmp(mp->unconstrainedPr, "Uniform"))
                        YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->brlensUni[0], mp->brlensUni[1]);
                    else if (!strcmp(mp->unconstrainedPr, "Exponential"))
                        YvyraPrint ("(%1.1lf)\n", mp->brlensExp);
                    else if (!strcmp(mp->unconstrainedPr, "twoExp"))
                        YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->brlens2Exp[0], mp->brlens2Exp[1]);
                    else
                        YvyraPrint ("(%1.1lf,%1.4lf,%1.1lf,%1.1lf)\n", mp->brlensDir[0], mp->brlensDir[1], mp->brlensDir[2], mp->brlensDir[3]);
                    }
                else if (!strcmp(mp->brlensPr, "Clock"))
                    {
                    YvyraPrint ("%s            Prior      = Clock:%s\n", spacer, mp->clockPr);
                    if ((!strcmp(mp->clockPr, "Uniform") ||
                         !strcmp(mp->clockPr, "Birthdeath") ||
                         !strcmp(mp->clockPr, "Fossilization")) && !strcmp(mp->nodeAgePr, "Unconstrained"))
                        {
                        YvyraPrint ("%s                         Tree age has a %s distribution\n", spacer, mp->treeAgePr.name);
                        }
                    if (!strcmp(mp->nodeAgePr,"Calibrated"))
                        {
                        b = 0;
                        for (a=0; a<numTaxa; a++)
                            {
                            if (taxaInfo[a].isDeleted == NO && tipCalibration[a].prior != unconstrained)
                                b++;
                            }
                        for (a=0; a<numDefinedConstraints; a++)
                            {
                            if (mp->activeConstraints[a] == YES && nodeCalibration[a].prior != unconstrained)
                                b++;
                            }
                        if (b > 1)
                            YvyraPrint ("%s                         Node depths are constrained by the following age constraints:\n", spacer);
                        else
                            YvyraPrint ("%s                         Node depths are calibrated by the following age constraint:\n", spacer);
                        for (a=0; a<numTaxa; a++)
                            {
                            if (taxaInfo[a].isDeleted == NO && tipCalibration[a].prior != unconstrained)
                                {
                                YvyraPrint ("%s                         -- The age of terminal \"%s\" is %s\n", spacer, taxaNames[a],
                                    tipCalibration[a].name);
                                }
                            }
                        for (a=b=0; a<numDefinedConstraints; a++)
                            {
                            if (mp->activeConstraints[a] == YES && nodeCalibration[a].prior != unconstrained)
                                {
                                YvyraPrint ("%s                         -- The age of node '%s' is %s\n", spacer,
                                    constraintNames[a], nodeCalibration[a].name);
                                for (k=0; k<numTaxa; k++)
                                    if (taxaInfo[k].isDeleted == NO && IsBitSet(k,definedConstraint[a]) == NO)
                                        break;
                                if (k == numTaxa)
                                    b = 1;          /* root is calibrated */
                                }
                            }
                        if (b == 0) /* we need to use default calibration for root for uniform and birthdeath */
                            {
                            if (!strcmp(mp->clockPr,"Uniform") || !strcmp(mp->clockPr, "Birthdeath") || !strcmp(mp->clockPr,"Fossilization"))
                                {
                                YvyraPrint ("%s                         -- Tree age has a %s distribution\n", spacer, mp->treeAgePr.name);
                                }
                            }
                        }
                    else
                        YvyraPrint ("%s                         Node ages are not constrained\n", spacer);
                    }
                else
                    {
                    assert (!strcmp(mp->brlensPr, "Fixed"));
                    YvyraPrint ("%s            Prior      = Fixed, branch lengths are fixed to the ones of the user tree '%s'\n", spacer,
                                                    userTree[mp->topologyFix]->name);
                    }
                }
            }
        else if (j == P_SPECRATE)
            {
            if (!strcmp(mp->speciationPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->speciationUni[0], mp->speciationUni[1]);
            else if (!strcmp(mp->speciationPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->speciationExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->speciationFix);
            }
        else if (j == P_EXTRATE)
            {
            if (!strcmp(mp->extinctionPr,"Beta"))
                YvyraPrint ("%s            Prior      = Beta(%1.2lf,%1.2lf)\n", spacer, mp->extinctionBeta[0], mp->extinctionBeta[1]);
            else if (!strcmp(mp->extinctionPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->extinctionExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->extinctionFix);
            }
        else if (j == P_FOSLRATE)
            {
            if (!strcmp(mp->fossilizationPr,"Beta"))
                YvyraPrint ("%s            Prior      = Beta(%1.2lf,%1.2lf)\n", spacer, mp->fossilizationBeta[0], mp->fossilizationBeta[1]);
            else if (!strcmp(mp->fossilizationPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->fossilizationExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->fossilizationFix);
            }
        else if (j == P_POPSIZE)
            {
            if (!strcmp(mp->popSizePr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->popSizeUni[0], mp->popSizeUni[1]);
            else if (!strcmp(mp->popSizePr,"Lognormal"))
                YvyraPrint ("%s            Prior      = Lognormal(%1.2lf,%1.2lf)\n", spacer, mp->popSizeLognormal[0], mp->popSizeLognormal[1]);
            else if (!strcmp(mp->popSizePr,"Normal"))
                YvyraPrint ("%s            Prior      = Truncated Normal(%1.2lf,%1.2lf)\n", spacer, mp->popSizeNormal[0], mp->popSizeNormal[1]);
            else if (!strcmp(mp->popSizePr,"Gamma"))
                YvyraPrint ("%s            Prior      = Gamma(%1.2lf,%1.2lf)\n", spacer, mp->popSizeGamma[0], mp->popSizeGamma[1]);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.5lf)\n", spacer, mp->popSizeFix);
            if (!strcmp(mp->topologyPr,"Speciestree"))
                {
                if (!strcmp(mp->popVarPr,"Equal") || !strcmp(mp->popSizePr,"Fixed"))
                    YvyraPrint ("%s                         Population size the same across species tree\n", spacer);
                else
                    YvyraPrint ("%s                         Population size varies across branches in species tree\n", spacer);
                }
            }
        else if (j == P_GROWTH)
            {
            if (!strcmp(mp->growthPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->growthUni[0], mp->growthUni[1]);
            else if (!strcmp(mp->growthPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->growthExp);
            else if (!strcmp(mp->growthPr,"Normal"))
                YvyraPrint ("%s            Prior      = Normal(%1.2lf,%1.2lf)\n", spacer, mp->growthNorm[0], mp->growthNorm[1]);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->growthFix);
            } 
        else if (j == P_BMCORR)
            {
            if (!strcmp(mp->brownCorrPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->brownCorrUni[0], mp->brownCorrUni[1]);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->brownCorrFix);
            }
        else if (j == P_BMSIGMA)
            {
            if (!strcmp(mp->brownScalePr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->brownScaleUni[0], mp->brownScaleUni[1]);
            else if (!strcmp(mp->brownScalePr,"Gamma"))
                YvyraPrint ("%s            Prior      = Gamma(%1.2lf,%1.2lf)\n", spacer, mp->brownScaleGamma[0], mp->brownScaleGamma[1]);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->brownScaleFix);
            }
        else if (j == P_CPPRATE)
            {
            if (!strcmp(mp->cppRatePr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->cppRateExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->cppRateFix);
            }
        else if (j == P_CPPMULTDEV)
            {
            if (!strcmp(mp->cppMultDevPr,"Fixed"))
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->cppMultDevFix);
            }
        else if (j == P_CPPEVENTS)
            {
            YvyraPrint ("%s            Prior      = Poisson (%s) [Events]\n", spacer, modelSettings[p->relParts[0]].cppRate->name);
            YvyraPrint ("%s                         Lognormal (0.00,%1.2lf) [Rate multipliers]\n", spacer, mp->cppMultDevFix);
            }
        else if (j == P_TK02VAR)
            {
            if (!strcmp(mp->tk02varPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->tk02varUni[0], mp->tk02varUni[1]);
            else if (!strcmp(mp->tk02varPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->tk02varExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->tk02varFix);
            }
        else if (j == P_TK02BRANCHRATES)
            {
            YvyraPrint ("%s            Prior      = LogNormal (expectation = r_0, variance = %s * v) \n", spacer, modelSettings[p->relParts[0]].tk02var->name);
            YvyraPrint ("%s                            [r_0 is beginning rate of branch, v is branch length]\n", spacer);
            }
        else if (j == P_WNVAR)
            {
            if (!strcmp(mp->wnvarPr, "Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->wnvarUni[0], mp->wnvarUni[1]);
            else if (!strcmp(mp->wnvarPr, "Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->wnvarExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->wnvarFix);
            }
        else if (j == P_WNBRANCHRATES)
            {
            YvyraPrint ("%s            Prior      = Gamma (expectation = 1.0, variance = %s / v) \n", spacer, modelSettings[p->relParts[0]].wnvar->name);
            YvyraPrint ("%s                            [v is branch length (t * c)]\n", spacer);
            }
        else if (j == P_IGRVAR)
            {
            if (!strcmp(mp->igrvarPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->igrvarUni[0], mp->igrvarUni[1]);
            else if (!strcmp(mp->igrvarPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->igrvarExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->igrvarFix);
            }
        else if (j == P_IGRBRANCHRATES)
            {
            YvyraPrint ("%s            Prior      = Gamma (expectation = 1.0, variance = %s) \n", spacer, modelSettings[p->relParts[0]].igrvar->name);
            }
        else if (j == P_ILNVAR)
            {
            if (!strcmp(mp->ilnvarPr, "Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->ilnvarUni[0], mp->ilnvarUni[1]);
            else if (!strcmp(mp->ilnvarPr, "Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->ilnvarExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->ilnvarFix);
            }
        else if (j == P_ILNBRANCHRATES)
            {
            YvyraPrint ("%s            Prior      = LogNormal (expectation = 1.0, variance = %s) \n", spacer, modelSettings[p->relParts[0]].ilnvar->name);
            }
        else if (j == P_MIXEDVAR)
            {
            if (!strcmp(mp->mixedvarPr,"Uniform"))
                YvyraPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, mp->mixedvarUni[0], mp->mixedvarUni[1]);
            else if (!strcmp(mp->mixedvarPr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->mixedvarExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, mp->mixedvarFix);
            }
        else if (j == P_MIXEDBRCHRATES)
            {
            YvyraPrint ("%s            Prior      = Mixed IGR and ILN (mean = 1, variance = %s) \n", spacer, modelSettings[p->relParts[0]].mixedvar->name);
            YvyraPrint ("%s                         Uniform prior on relaxed clock models [Pr(IGR)=Pr(ILN)=0.5]\n", spacer);
            }
        else if (j == P_CLOCKRATE)
            {
            if (!strcmp(mp->clockRatePr,"Normal"))
                YvyraPrint ("%s            Prior      = Normal(%1.6lf,%1.6lf)\n", spacer, mp->clockRateNormal[0], mp->clockRateNormal[1]);
            else if (!strcmp(mp->clockRatePr,"Lognormal"))
                YvyraPrint ("%s            Prior      = Lognormal(%1.2lf,%1.2lf)\n", spacer, mp->clockRateLognormal[0], mp->clockRateLognormal[1]);
            else if (!strcmp(mp->clockRatePr,"Gamma"))
                YvyraPrint ("%s            Prior      = Gamma(%1.2lf,%1.2lf)\n", spacer, mp->clockRateGamma[0], mp->clockRateGamma[1]);
            else if (!strcmp(mp->clockRatePr,"Exponential"))
                YvyraPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, mp->clockRateExp);
            else
                YvyraPrint ("%s            Prior      = Fixed(%1.6lf)\n", spacer, mp->clockRateFix);
            if (!strcmp(mp->clockVarPr,"Strict"))
                YvyraPrint ("%s                         The clock rate is constant throughout the tree (strict clock)\n", spacer);
            else if (!strcmp(mp->clockVarPr,"Cpp"))
                YvyraPrint ("%s                         The clock rate varies according to a CPP model\n", spacer);
            else if (!strcmp(mp->clockVarPr,"TK02"))
                YvyraPrint ("%s                         The clock rate varies according to a autocorrelated lognormal model\n", spacer);
            else if (!strcmp(mp->clockVarPr,"WN"))
                YvyraPrint ("%s                         The clock rate varies according to a white noise model\n", spacer);
            else if (!strcmp(mp->clockVarPr,"IGR"))
                YvyraPrint ("%s                         The clock rate varies according to an independent gamma model\n", spacer);
            else if (!strcmp(mp->clockVarPr,"ILN"))
                YvyraPrint ("%s                         The clock rate varies according to an independent lognormal model\n", spacer);
            else /* if (!strcmp(mp->clockVarPr,"Mixed")) */
                YvyraPrint ("%s                         The clock rate varies according to mixed IGR and ILN models\n", spacer);
            }
        else if (j == P_SPECIESTREE)
            {
            YvyraPrint ("%s            Prior      = Uniform on topologies and branch lengths\n", spacer);
            }
                
        /* print partitions */
        if (numCurrentDivisions > 1)
            {
            if (p->nRelParts == 1)
                YvyraPrint ("%s            Partition  = %d\n", spacer, p->relParts[0]+1);
            else if (p->nRelParts == 2)
                {
                YvyraPrint ("%s            Partitions = %d and %d\n", spacer, p->relParts[0]+1, p->relParts[1]+1);                    
                }
            else if (p->nRelParts == numCurrentDivisions)
                {
                YvyraPrint ("%s            Partitions = All\n", spacer);                  
                }           
            else /* if (p->nRelParts > 2) */
                {
                YvyraPrint ("%s            Partitions = ", spacer);
                for (j=0; j<p->nRelParts; j++)
                    {
                    if (j == p->nRelParts - 2)
                        YvyraPrint ("%d, and ", p->relParts[j]+1);
                    else if (j == p->nRelParts - 1)
                        YvyraPrint ("%d\n", p->relParts[j]+1);
                    else
                        YvyraPrint ("%d, ", p->relParts[j]+1);
                    }               
                }
            }

        /* show subparams */
        if (p->nSubParams > 0)
            {
            if (p->nSubParams == 1)
                YvyraPrint ("%s            Subparam.  = %s\n", spacer, p->subParams[0]->name);
            else
                {
                printedCol = 0;
                for (k=0; k<p->nSubParams; k++)
                    {
                    if (k == 0)
                        {
                        YvyraPrint ("%s            Subparams  = %s", spacer, p->subParams[k]->name);
                        printedCol = (int)(strlen(spacer)) + 25 + (int)(strlen(p->subParams[k]->name));
                        }
                    else if (k == p->nSubParams - 1)
                        {
                        if (printedCol + 5 > screenWidth)
                            YvyraPrint ("\n%s                         and ", spacer);
                        else if (printedCol + (int)(strlen(p->subParams[k]->name)) + 5 > screenWidth)
                            YvyraPrint (" and \n%s                         ", spacer);
                        else
                            YvyraPrint (" and ");
                        YvyraPrint ("%s\n", p->subParams[k]->name);
                        }
                    else
                        {
                        if (printedCol + (int)(strlen(p->subParams[k]->name)) + 2 > screenWidth)
                            {
                            YvyraPrint (", \n%s                         ", spacer);
                            printedCol = (int)(strlen(spacer)) + 25;
                            }
                        else
                            {
                            YvyraPrint (", ");
                            printedCol += 2;
                            }
                        YvyraPrint ("%s", p->subParams[k]->name);
                        printedCol += (int)strlen(p->subParams[k]->name);
                        }
                    }
                }
            }

        /* show used moves */
        if (showMoves == YES && (p->printParam == YES || p->paramType == P_TOPOLOGY || p->paramType == P_BRLENS))
            {
            /* check if runs are same */
            areRunsSame = YES;
            for (run=1; run<chainParams.numRuns; run++)
                {
                for (chain=0; chain<chainParams.numChains; chain++)
                    {
                    for (k=0; k<numApplicableMoves; k++)
                        {
                        mv = moves[k];
                        if (mv->parm != p)
                            continue;
                        chainIndex = run*chainParams.numChains + chain;
                        refIndex = chain;
                        if (AreDoublesEqual (mv->relProposalProb[chainIndex], mv->relProposalProb[refIndex], 0.000001) == NO)
                            areRunsSame = NO;
                        }
                    }
                }

            for (run=0; run<chainParams.numRuns; run++)
                {
                if (run > 0 && areRunsSame == YES)
                    continue;
                
                /* check if chains are same */
                areChainsSame = YES;
                for (chain=1; chain<chainParams.numChains; chain++)
                    {
                    for (k=0; k<numApplicableMoves; k++)
                        {
                        mv = moves[k];
                        if (mv->parm != p)
                            continue;
                        chainIndex = run*chainParams.numChains + chain;
                        refIndex = run*chainParams.numChains;
                        if (AreDoublesEqual (mv->relProposalProb[chainIndex], mv->relProposalProb[refIndex], 0.000001) == NO)
                            areChainsSame = NO;
                        }
                    }

                /* now print moves */
                for (chain=0; chain<chainParams.numChains; chain++)
                    {
                    if (chain > 0 && areChainsSame == YES)
                        continue;
                    if (run==0 && chain==0)
                        YvyraPrint ("%s            Moves      = ", spacer);
                    else
                        YvyraPrint ("%s                         ", spacer);
                    numPrinted = 0;
                    printedCol = (int)(strlen(spacer)) + 25;
                    for (k=0; k<numApplicableMoves; k++)
                        {
                        mv = moves[k];
                        if (mv->parm != p)
                            continue;
                        chainIndex = run*chainParams.numChains + chain;
                        if (mv->relProposalProb[chainIndex] <= 0.000001)
                            continue;
                        if (numPrinted == 0)
                            {
                            YvyraPrint ("%s <prob %.1f>", mv->name, mv->relProposalProb[chainIndex]);
                            printedCol += 9 + (int)strlen(mv->name) + (int)(log10(mv->relProposalProb[chainIndex])) + 3;
                            }
                        else
                            {
                            if (printedCol + 11 + (int)(strlen(mv->name)) + (int)(log10(mv->relProposalProb[chainIndex])) + 3 > screenWidth)
                                {
                                YvyraPrint (", \n%s                         ", spacer);
                                printedCol = 25 + (int)(strlen(spacer));
                                }
                            else
                                {
                                YvyraPrint (", ");
                                printedCol += 2;
                                }
                            YvyraPrint ("%s <prob %.1f>", mv->name, mv->relProposalProb[chainIndex]);
                            printedCol += (9 + (int)(strlen(mv->name)) + (int)(log10(mv->relProposalProb[chainIndex])) + 3);
                            }
                        numPrinted++;
                        }

                    if (numPrinted == 0 && p->paramType == P_BRLENS)
                        {
                        for (k=0; k<numParams; k++)
                            if (params[k].tree == p->tree)
                                break;
                        YvyraPrint ("For moves, see param. '%s'", params[k].name);
                        }
                    else if (numPrinted == 0)
                        YvyraPrint ("WARNING! No moves -- param. fixed to startvals");

                    if (areRunsSame == YES && areChainsSame == YES)
                        YvyraPrint ("\n");
                    else if (areRunsSame == YES && areChainsSame == NO)
                        YvyraPrint (" [chain %d]\n", chain+1);
                    else if (areRunsSame == NO && areChainsSame == YES)
                        YvyraPrint (" [run %d]\n", run+1);
                    else /* if (areRunsSame == NO && areChainsSame == NO) */
                            YvyraPrint (" [run %d, chain %d]\n", run+1, chain+1);
                    }
                }
            }
            
        /* show available moves */
        if (showAllAvailable == YES && (p->printParam == YES || p->paramType == P_TOPOLOGY || p->paramType == P_BRLENS))
            {
            /* loop over moves */
            numPrinted = 0;
            printedCol = 0;
            for (k=0; k<numApplicableMoves; k++)
                {
                mv = moves[k];
                if (mv->parm != p)
                    continue;
                numMovedChains = 0;
                for (run=0; run<chainParams.numRuns; run++)
                    {
                    for (chain=0; chain<chainParams.numChains; chain++)
                        {
                        chainIndex = run*chainParams.numChains + chain;
                        if (mv->relProposalProb[chainIndex] > 0.000001)
                            numMovedChains++;
                        }
                    }
                if (numMovedChains == 0)
                    {
                    if (numPrinted == 0)
                        {
                        YvyraPrint ("%s            Not used   = ", spacer);
                        printedCol = (int)(strlen(spacer)) + 25;
                        }
                    else if (printedCol + 2 + (int)(strlen(mv->moveType->shortName)) > screenWidth)
                        {
                        YvyraPrint (", \n%s                         ", spacer);
                        printedCol = (int)(strlen(spacer)) + 25;
                        }
                    else
                        {
                        YvyraPrint (", ");
                        printedCol += 2;
                        }
                    YvyraPrint ("%s", mv->moveType->shortName);
                    printedCol += (int)strlen(mv->moveType->shortName);
                    numPrinted++;
                    }
                }
                if (numPrinted > 0)
                    YvyraPrint ("\n");
            }

        /* show startvals */
        if (showStartVals == YES && (p->printParam == YES || p->paramType == P_TOPOLOGY || p->paramType == P_BRLENS || p->paramType == P_SPECIESTREE || p->paramType == P_POPSIZE))
            {                   
            if (p->paramType == P_TOPOLOGY || p->paramType == P_BRLENS || p->paramType == P_SPECIESTREE)
                {
                /* check if they are the same */
                areRunsSame = YES;
                if (p->paramType == P_TOPOLOGY)
                    {
                    for (run=1; run<chainParams.numRuns; run++)
                        {
                        for (chain=0; chain<chainParams.numChains; chain++)
                            {
                            if (AreTopologiesSame (GetTree (p, run*chainParams.numChains + chain, 0), GetTree (p, chain, 0)) == NO)
                                {
                                areRunsSame = NO;
                                break;
                                }
                            }
                        }
                    }
                else if (p->paramType == P_BRLENS)
                    {
                    for (run=1; run<chainParams.numRuns; run++)
                        {
                        for (chain=0; chain<chainParams.numChains; chain++)
                            {
                            if (AreTreesSame (GetTree (p, run*chainParams.numChains + chain, 0), GetTree (p, chain, 0)) == NO)
                                {
                                areRunsSame = NO;
                                break;
                                }
                            }
                        }
                    }
                else if (p->paramType == P_SPECIESTREE)
                    {
                    for (run=1; run<chainParams.numRuns; run++)
                        {
                        for (chain=0; chain<chainParams.numChains; chain++)
                            {
                            if (AreTreesSame (GetTree (p, run*chainParams.numChains + chain, 0), GetTree (p, chain, 0)) == NO)
                                {
                                areRunsSame = NO;
                                break;
                                }
                            }
                        }
                    }
                
                /* print trees */
                for (run=0; run<chainParams.numRuns; run++)
                    {
                    if (run > 0 && areRunsSame == YES)
                        break;
                    for (chain=0; chain<chainParams.numChains; chain++)
                        {
                        if (run == 0 && chain == 0)
                            YvyraPrint ("%s            Startvals  = tree '%s'", spacer, GetTree (p, run*chainParams.numChains+chain, 0)->name);
                        else
                            YvyraPrint ("%s                         tree '%s'", spacer, GetTree (p, run*chainParams.numChains+chain, 0)->name);
                        if (chainParams.numChains > 1 && areRunsSame == YES)
                            YvyraPrint ("  [chain %d]", chain+1);
                        else if (chainParams.numChains == 1 && areRunsSame == NO)
                            YvyraPrint ("  [run %d]", run+1);
                        else if (areRunsSame == NO)
                            YvyraPrint ("  [run %d, chain %d]", run+1, chain+1);
                        YvyraPrint ("\n");
                        }
                    }               
                }   /* end topology and brlens parameters */

            else
                {
                /* run of the mill parameter */
                if (p->paramType == P_CLOCKRATE)
                    {
                     for (j=0; j<numGlobalChains; j++)
                        {
                        if (UpdateClockRate(-1.0, j) == ERROR)
                            {
                            YvyraPrint ("%s            Warning: There is no appropriate clock rate that would satisfy all calibrated trees for run:%d chain%d. Some of calibration, trees or clockprior needs to be changed. ", spacer, j/chainParams.numChains, j%chainParams.numChains);
                            }
                        }
                    }
                areRunsSame = YES;
                for (run=1; run<chainParams.numRuns; run++)
                    {
                    for (chain=0; chain<chainParams.numChains; chain++)
                        {
                        if ((p->paramType == P_PI && modelParams[p->relParts[0]].dataType != STANDARD))
                            {
                            nValues = p->nSubValues;
                            value = GetParamSubVals (p, run*chainParams.numChains+chain, 0);
                            refValue = GetParamSubVals (p, chain, 0);
                            }
                        else
                            {
                            nValues = p->nValues;
                            value = GetParamVals (p, run*chainParams.numChains+chain, 0);
                            refValue = GetParamVals (p, chain, 0);
                            }
                        for (k=0; k<nValues; k++)
                            {
                            if (AreDoublesEqual(value[k], refValue[k], 0.000001) == NO)
                                {
                                areRunsSame = NO;
                                break;
                                }                           
                            }
                        if (areRunsSame == NO)
                            break;                          
                        }
                    if (areRunsSame == NO)
                        break;
                    }

                /* print values */
                for (run=0; run<chainParams.numRuns; run++)
                    {
                    if (run > 0 && areRunsSame == YES)
                        break;
                    areChainsSame = YES;
                    for (chain=1; chain<chainParams.numChains; chain++)
                        {
                        if ((p->paramType == P_PI && modelParams[p->relParts[0]].dataType != STANDARD))
                            {
                            nValues = p->nSubValues;
                            value = GetParamSubVals (p, run*chainParams.numChains+chain, 0);
                            refValue = GetParamSubVals (p, run*chainParams.numChains, 0);
                            }
                        else
                            {
                            nValues = p->nValues;
                            value = GetParamVals (p, run*chainParams.numChains+chain, 0);
                            refValue = GetParamVals (p, run*chainParams.numChains, 0);
                            }
                        for (k=0; k<nValues; k++)
                            {
                            if (AreDoublesEqual(value[k], refValue[k], 0.000001) == NO)
                                {
                                areChainsSame = NO;
                                break;
                                }                           
                            }
                        }
                    for (chain=0; chain<chainParams.numChains; chain++)
                        {
                        if (areChainsSame == YES && chain > 0)
                            continue;

                        if ((p->paramType == P_PI && modelParams[p->relParts[0]].dataType != STANDARD))
                            {
                            nValues = p->nSubValues;
                            value = GetParamSubVals (p, run*chainParams.numChains+chain, 0);
                            }
                        else
                            {
                            nValues = p->nValues;
                            value = GetParamVals (p, run*chainParams.numChains+chain, 0);
                            }

                        if (run == 0 && chain == 0)
                            YvyraPrint ("%s            Startvals  = (%1.3lf", spacer, value[0]);
                        else
                            YvyraPrint ("%s                         (%1.3lf", spacer, value[0]);
                        
                        for (k=1; k<nValues; k++)
                            {
                            if (k%10==0)
                                YvyraPrint (",\n%s                          %1.3lf", spacer, value[k]);
                            else
                                YvyraPrint (",%1.3lf", value[k]);
                            }
                        YvyraPrint (")");
                        
                        if (areChainsSame == YES && areRunsSame == NO)
                            YvyraPrint ("  [run %d]", run+1);
                        else if (areChainsSame == NO && areRunsSame == YES)
                            YvyraPrint ("  [chain %d]", chain+1);
                        else if (areChainsSame == NO && areRunsSame == NO)
                            YvyraPrint ("  [run %d, chain %d]", run+1, chain+1);
                        YvyraPrint ("\n");
                        }   
                    }
                }
            }   /* end print start values */

        YvyraPrint ("\n");
        }   /* next parameter */

    return (NO_ERROR);
}


int Unlink (void)
{
    return (NO_ERROR);
}


/*-------------------------------------------------
|
|   UpdateClockRate:    Update clockRate of the given chain. Above all it will enforce fixed clockrate prior if it is set.
|                       Error will be returned if fixed clockrate prior may not be respected.
|   @param clockRate    is the new clockRate to setup. Clock rate value could be set as positive, 0.0 or negative value. 
|                       The function does the fallowing depending on one of this three values:
|                        positive    - check that this 'positive' value is suitable rate. At the end re-enforce (update) the 'positive' value as clock rate on all trees. 
|                        0.0         - check if current rate is suitable, if not update it with minimal suitable value. At the end re-enforce (update) the resulting clock rate on all trees. 
|                        negative    - check if current rate is suitable, if not update it with minimal suitable value. At the end re-enforce (update) the resulting clock rate ONLY if clock rate was changed 
|   @return             ERROR if clockRate can not be set up, NO_ERROR otherwise. 
|
--------------------------------------------------*/
int UpdateClockRate (YFlt clockRate, int chain)
{
    int         i, updateTrees;
    YFlt      *clockRatep;
    Tree        *t, *t_calibrated=NULL;
    YFlt      mintmp,maxtmp,minClockRate,maxClockRate;

    clockRatep=NULL;
    minClockRate = 0.0;
    maxClockRate = YFLT_MAX;

    for (i=0; i<numTrees; i++)
        {
        t = GetTreeFromIndex(i, chain, 0);
        if (t->isCalibrated == NO)
            continue;

        if (clockRatep == NULL)
            {
            clockRatep = GetParamVals(modelSettings[t->relParts[0]].clockRate, chain, 0);
            t_calibrated = t;
            assert (clockRatep);
            }

        findAllowedClockrate (t, &mintmp, &maxtmp);

        if (minClockRate < mintmp)
            minClockRate = mintmp;

        if (maxClockRate > maxtmp)
            maxClockRate = maxtmp;

        }
    /* clock rate is the same for all trees of a given chain */
    if (clockRatep != NULL)
        {
        if (minClockRate > maxClockRate)
            {
            YvyraPrint ("%s   ERROR: Calibrated trees require compatible clockrates but they are incompatible for run %d, chain %d.\n",
                          spacer, chain/chainParams.numChains + 1, chain%chainParams.numChains + 1);
            *clockRatep=0;
            return (ERROR);
            }

        if (!strcmp(modelParams[t_calibrated->relParts[0]].clockRatePr, "Fixed"))
            {
            if (clockRate < 0.0 && AreDoublesEqual (*clockRatep, modelParams[t_calibrated->relParts[0]].clockRateFix, 0.0001) == YES)
                {
                updateTrees = NO;
                }
            else
                {
                updateTrees = YES;
                }
            *clockRatep = modelParams[t_calibrated->relParts[0]].clockRateFix;
            if ((*clockRatep < minClockRate && AreDoublesEqual (*clockRatep, minClockRate, 0.0001) == NO) ||
                (*clockRatep > maxClockRate && AreDoublesEqual (*clockRatep, maxClockRate, 0.0001) == NO))
                {
                YvyraPrint ("%s   ERROR: Calibrated trees require clockrate in range from %f to %f, while clockrate prior is fixed to %f for run %d chain %d.\n",
                              spacer, minClockRate, maxClockRate, *clockRatep, chain/chainParams.numChains + 1, chain%chainParams.numChains + 1);
                *clockRatep=0;
                return (ERROR);
                }
            if (clockRate > 0.0)
                {
                if (AreDoublesEqual (*clockRatep, clockRate, 0.0001) == NO)
                    {
                    YvyraPrint ("%s   ERROR: Requested clockrate:%f does not match fixed clockrate prior:%f.\n", spacer, clockRate, *clockRatep);
                    *clockRatep=0;
                    return (ERROR);
                    }
                }
            }
        else {
            /* clock prior is not fixed */
            updateTrees = YES;
            if (clockRate > 0.0)
                {
                *clockRatep = clockRate;
                if ((*clockRatep < minClockRate && AreDoublesEqual (*clockRatep, minClockRate, 0.0001) == NO) ||
                    (*clockRatep > maxClockRate && AreDoublesEqual (*clockRatep, maxClockRate, 0.0001) == NO))
                    {
                    YvyraPrint ("%s   ERROR: Calibrated trees require clockrate in range from %f to %f, while requested clockrate is %f for run %d chain %d.\n",
                                  spacer, minClockRate, maxClockRate, clockRate, chain/chainParams.numChains + 1, chain%chainParams.numChains + 1);
                    *clockRatep=0;
                    return (ERROR);
                    }
                }
            else if (clockRate == 0.0) 
                {
                if ((*clockRatep < minClockRate && AreDoublesEqual (*clockRatep, minClockRate, 0.0001) == NO) ||
                    (*clockRatep > maxClockRate && AreDoublesEqual (*clockRatep, maxClockRate, 0.0001) == NO))
                    {
                    *clockRatep = minClockRate;
                    }
                }
            else // if (clockRate < 0.0)
                {
                if ((*clockRatep < minClockRate && AreDoublesEqual (*clockRatep, minClockRate, 0.0001) == NO) ||
                    (*clockRatep > maxClockRate && AreDoublesEqual (*clockRatep, maxClockRate, 0.0001) == NO))
                    {
                    *clockRatep = minClockRate;
                    }
                else
                    {
                    updateTrees = NO;
                    }
                }
            }

        
        if (updateTrees == YES)
            {
            for (i=0; i<numTrees; i++)
                {
                t = GetTreeFromIndex(i, chain, 0);
                if (t->isCalibrated == NO)
                    continue;
                UpdateTreeWithClockrate (t,*clockRatep);
                }
            }
        }

    return (NO_ERROR);
}


/*----------------------------------------------
|
|   UpdateCppEvolLength: Recursive function to
|      update evolLength of one node for Cpp
|      relaxed clock model
|
-----------------------------------------------*/
int UpdateCppEvolLength (int *nEvents, YFlt **pos, YFlt **rateMult, YFlt *evolLength, TreeNode *p, YFlt baseRate)
{
    int     i;
    YFlt  endRate;

    if (p != NULL)
        {
#   ifdef DEBUG_CPP
        if (baseRate < POS_MIN || baseRate > POS_INFINITY)
            {
            printf ("baseRate out of bounds (%.15e for node %d\n", baseRate, p->index);
            return (ERROR);
            }
#   endif
        p->upDateTi = YES;
        p->upDateCl = YES;
        if (nEvents[p->index] == 0)
            {
            evolLength[p->index] = p->length * baseRate;
            }
        else
            {
            /* note that event positions are from the top of the branch (more recent)
               to the bottom of the branch (older) */
            /* The algorithm below successively multiplies in the more basal rate multipliers,
               starting from the top of the branch. The length of the branch is first assumed
               to be 1.0; at the end we multiply the result with the true length of the branch. */
            evolLength[p->index] = pos[p->index][0] * rateMult[p->index][0];
            for (i=1; i<nEvents[p->index]; i++)
                {
                evolLength[p->index] += (pos[p->index][i] - pos[p->index][i-1]);
                evolLength[p->index] *= rateMult[p->index][i];
                }
            evolLength[p->index] += (1.0 - pos[p->index][nEvents[p->index]-1]);
            evolLength[p->index] *= baseRate;
            evolLength[p->index] *= p->length;
            }

        /* calculate end rate; we can do this in any order */
        endRate = baseRate;
        for (i=0; i<nEvents[p->index]; i++)
            endRate *= rateMult[p->index][i];

#   ifdef DEBUG_CPP
        if (endRate < POS_MIN || endRate > POS_INFINITY)
            {
            printf ("endRate out of bounds (%.15e for node %d)\n", endRate, p->index);
            return (ERROR);
            }
        if (p->anc != NULL && p->anc->anc != NULL && (evolLength[p->index] < POS_MIN || evolLength[p->index] > POS_INFINITY))
            {
            printf ("Effective branch length out of bounds (%.15e for node %d)\n", evolLength[p->index], p->index);
            return (ERROR);
            }
#   endif
        /* call left and right descendants */
        if (UpdateCppEvolLength (nEvents, pos, rateMult, evolLength, p->left, endRate)==ERROR)
            return (ERROR);
        if (UpdateCppEvolLength (nEvents, pos, rateMult, evolLength, p->right, endRate)==ERROR)
            return (ERROR);
        }

    return (NO_ERROR);
}


/*-------------------------------------------------
 |
 |  UpdateCppEvolLengths: Recalculate effective
 |      evolutionary lengths and set update flags
 |      for Cpp relaxed clock model
 |
 --------------------------------------------------*/
int UpdateCppEvolLengths (Param *param, TreeNode *p, int chain)
{
    int         i, *nEvents;
    TreeNode    *q;
    YFlt      baseRate = 1.0, **pos, **rateMult, *evolLength;
    
    i = 2*chain + state[chain];
    nEvents = param->nEvents[i];
    pos = param->position[i];
    rateMult = param->rateMult[i];
    evolLength = GetParamSubVals (param, chain, state[chain]);
    
    q = p->anc;
    while (q->anc != NULL)
        {
        for (i=0; i<nEvents[q->index]; i++)
            baseRate *= rateMult[q->index][i];
        q = q->anc;
        }
    
    if (UpdateCppEvolLength (nEvents, pos, rateMult, evolLength, p, baseRate)==ERROR)
        return (ERROR);
    
    return(NO_ERROR);
}


/* UpdateTK02EvolLengths: update branch lengths for tk02 model */
int UpdateTK02EvolLengths (Param *param, Tree *t, int chain)
{
    int         i;
    YFlt      *tk02Rate, *brlens;
    TreeNode    *p;
    
    tk02Rate = GetParamVals (param, chain, state[chain]);
    brlens = GetParamSubVals (param, chain, state[chain]);
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        brlens[p->index] = p->length * (tk02Rate[p->index] + tk02Rate[p->anc->index]) / 2.0;
        }
    
    return (NO_ERROR);
}


/* UpdateIndBranchLengths: update branch lengths for independent rates model */
int UpdateIndBrachLengths (Param *param, Tree *t, int chain)
{
    int         i;
    YFlt      *indRate, *brlens;
    TreeNode    *p;
    
    indRate = GetParamVals (param, chain, state[chain]);
    brlens = GetParamSubVals (param, chain, state[chain]);
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        brlens[p->index] = p->length * indRate[p->index];
        }

    return (NO_ERROR);
}

#ifndef MCMC_H_
#define MCMC_H_

int     AddToPrintString (char *tempStr);
void    AutotuneDirichlet (YFlt acceptanceRate, YFlt targetRate, int batch, YFlt *alphaPi, YFlt minTuning, YFlt maxTuning);
void    AutotuneMultiplier (YFlt acceptanceRate, YFlt targetRate, int batch, YFlt *lambda, YFlt minTuning, YFlt maxTuning);
void    AutotuneSlider (YFlt acceptanceRate, YFlt targetRate, int batch, YFlt *width, YFlt minTuning, YFlt maxTuning);
void    AutotuneRJClocks (YFlt acceptanceRate, YFlt targetRate, int batch, YFlt *w, YFlt minTuning, YFlt maxTuning);
int     DoMcmc (void);
int     DoMcmcp (void);
int     DoMcmcParm (char *parmName, char *tkn);
int     DoSs (void);
int     DoSsp (void);
int     DoSsParm (char *parmName, char *tkn);
int     ExhaustiveParsimonySearch (Tree *t, int chain, TreeInfo *tInfo);
YFlt  GetParsDP (Tree *t, TreeNode *p, int chain);
void    GetParsFP (Tree *t, TreeNode *p, int chain);
int     GetParsimonyBrlens (Tree *t, int chain, YFlt *brlens);
YFlt  GetParsimonyLength (Tree *t, int chain);
void    GetParsimonySubtreeRootstate (Tree *t, TreeNode *root, int chain);
YFlt  GetRate (int division, int chain);
int     LnBirthDeathPriorPr (Tree *t, YFlt clockRate, YFlt *prob, YFlt sR, YFlt eR, char *sS, YFlt sF);
int     LnCoalescencePriorPr (Tree *t, YFlt *prob, YFlt theta, YFlt growth);
YFlt  LnUniformPriorPr (Tree *t, YFlt clockRate);
int     LnFossilizationPriorPr (Tree *t, YFlt clockRate, YFlt *prob, YFlt *sR, YFlt *eR, YFlt *fR, YFlt sF, char *sS);
int     LogClockTreePriorRatio (Param *param, int chain, YFlt *lnPriorRatio);
YFlt  LogDirPrior (Tree *t, ModelParams *mp, int PV);
YFlt  LogOmegaPrior (YFlt w1, YFlt w2, YFlt w3);
FILE    *OpenNewMBPrintFile (char *fileName);
int     ResetScalersPartition (int *isScalerNode, Tree* t, unsigned rescaleFreq);
int     SafeSprintf (char **target, int *targetLen, char *fmt, ...);
int     SetFilePositions (int samplePos);
YFlt  TreeLength (Param *param, int chain);

#endif  /* MCMC_H_ */

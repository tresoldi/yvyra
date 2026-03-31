#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

//#define TIMING_ANALIZ
#if defined (TIMING_ANALIZ)
    static clock_t         CPUCondLikeDown;
    static clock_t         CPUScalers;
    static clock_t         CPUScalersRemove;
    static clock_t         CPUCondLikeRoot;
    static clock_t         CPULilklihood;

    #define TIME(X1,CPUtime)\
        {CPUTimeStart = clock();\
        X1;\
        CPUtime += (clock()-CPUTimeStart);}
#else
    #define TIME(X1,CPUtime)\
        X1;
#endif

void      CopySiteScalers (ModelInfo *m, int chain);
void      ResetSiteScalers (ModelInfo *m, int chain);
void      FlipCijkSpace (ModelInfo *m, int chain);
void      FlipCondLikeSpace (ModelInfo *m, int chain, int nodeIndex);
void      FlipNodeScalerSpace (ModelInfo *m, int chain, int nodeIndex);
void      FlipSiteScalerSpace (ModelInfo *m, int chain);
void      FlipTiProbsSpace (ModelInfo *m, int chain, int nodeIndex);

int       CondLikeDown_Std (TreeNode *p, int division, int chain);
int       CondLikeRoot_Std (TreeNode *p, int division, int chain);
int       CondLikeScaler_Std (TreeNode *p, int division, int chain);
int       CondLikeUp_Std (TreeNode *p, int division, int chain);
void      LaunchLogLikeForDivision (int chain, int d, YFlt* lnL);
int       Likelihood_Pars (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats);
int       Likelihood_ParsStd (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats);
int       Likelihood_Std (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats);
int       Likelihood_Cont (TreeNode *p, int division, int chain, YFlt *lnL, int whichSitePats);
int       TiProbs_Std (TreeNode *p, int division, int chain);

#endif  /* LIKELIHOOD_H_ */

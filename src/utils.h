#ifndef UTILS_H_
#define UTILS_H_

/* M_PI and M_PI_2 not part of standard C */
#ifndef M_PI
#    define M_PI    3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#    define M_PI_2  1.57079632679489661923132169163975144
#endif

struct YComplex
{
    YFlt re;
    YFlt im;
};
typedef struct YComplex YComplex;

typedef struct
    {
    YFlt     mean;
    YFlt     median;
    YFlt     lower;
    YFlt     upper;
    YFlt     var;
    YFlt     PSRF;
    YFlt     avrESS;
    YFlt     minESS;
    }
    Stat;


/* For explanation why the following two macros exists, see
 * http://stackoverflow.com/questions/38569628/calling-a-free-wrapper-dereferencing-type-punned-pointer-will-break-strict-al
 */

#define SAFEFREE(ptr) (ptr = SafeFree(ptr))
#define ALIGNEDSAFEFREE(ptr) (ptr = AlignedSafeFree(ptr))

int      AddBitfield (BitsLong ***list, int listLen, int *set, int setLen);
#if defined (SSE_ENABLED)
void    *AlignedMalloc (size_t size, size_t alignment);
void    *AlignedSafeFree (void *ptr);
#endif
int      AreBitfieldsEqual (BitsLong *p, BitsLong *q, int length);
int      Bit (int n, BitsLong *p);
void     ClearBit (int i, BitsLong *bits);
void     ClearBits (BitsLong *bits, int nLongs);
void     CopyBits (BitsLong *dest, BitsLong *source, int nLongs);
int      CopyResults (FILE *toFile, char *fromFileName, long long lastGen);
int      CopyProcessSsFile (FILE *toFile, char *fromFileName, int lastStep, YFlt *marginalLnLSS, YFlt *splitfreqSS);
int      CopyTreeResults (FILE *toFile, char *fromFileName, long long lastGen, int *treeNum);
int      FirstTaxonInPartition (BitsLong *partition, int length);
long     FirstTree (FILE *fp, char *lineBuf, int longestLine);
int      Flip01 (int x);
void     FlipBits (BitsLong *partition, int length, BitsLong *mask);
void     FlipOneBit (int n, BitsLong *p);
int      FromGrowthFxnToIndex (int *growthFxn);
void     FromIndexToGrowthFxn (int index, int *growthFxn);
void     GetIntSummary (int **vals, int nRows, int *rowCount, Stat *theStats, int HPD);
int      GetKFromGrowthFxn (int *growthFxn);
void     GetSummary (YFlt **vals, int nRows, int *rowCount, Stat *theStats, int HPD);
int      HarmonicArithmeticMeanOnLogs (YFlt *vals, int nVals, YFlt *mean, YFlt *harm_mean);
int      IsBitSet (int i, BitsLong *bits);
int      IsConsistentWith (const char *token, const char *expected);
int      IsPartNested (BitsLong *smaller, BitsLong *larger, int length);
int      IsPartCompatible (BitsLong *smaller, BitsLong *larger, int length);
int      IsSectionEmpty (BitsLong *bitField1, BitsLong *bitField2, int length);
int      IsUnionEqThird (BitsLong *bitField1, BitsLong *bitField2, BitsLong *bitField3, int length);
long     LastBlock (FILE *fp, char *lineBuf, int longestLine);
int      LineTermType (FILE *fp);
YFlt   LnDirichlet (YFlt *alphai, YFlt *xi, int lengthi);
int      LongestLine (FILE *fp);
void     LowerUpperMedian (YFlt *vals, int nVals, YFlt *lower, YFlt *upper, YFlt *median);
void     LowerUpperMedianHPD (YFlt *vals, int nVals, YFlt *lower, YFlt *upper, YFlt *median);
void     MeanVariance (YFlt *vals, int nVals, YFlt *mean, YFlt *var);
void     MeanVarianceLog (YFlt *vals, int nVals, YFlt *mean, YFlt *var, YFlt *varEst);
int      NextTaxonInPartition (int currentTaxon, BitsLong *partition, int length);
int      NBits (int x);
int      NumBits (BitsLong *x, int len);
char    *PrintNum (YFlt num);
void     YvyraPrint (char *format, ...);
void     YvyraPrintf (FILE *f, char *format, ...);
FILE    *OpenBinaryFileR (char *name);
FILE    *OpenTextFileA (char *name);
FILE    *OpenTextFileR (char *name);
FILE    *OpenTextFileRQuait (char *name);
FILE    *OpenTextFileW (char *name);
YFlt   PotentialScaleReduction (YFlt **vals, int nRows, int *count);
void     EstimatedSampleSize (YFlt **vals, int nRuns, int *count, YFlt *returnESS);
void    *SafeCalloc (size_t n, size_t s);
int      SafeFclose (FILE **fp);
void    *SafeFree (void *ptr);
void    *SafeMalloc (size_t s);
void    *SafeRealloc (void *ptr, size_t s);
char    *SafeStrcat (char **target, const char *source);
char    *SafeStrcpy (char **target, const char *source);
void     SetBit (int i, BitsLong *bits);
void     SortInts (int *item, int *assoc, int count, int descendingOrder);
void     SortInts2 (int *item, int *assoc, int left, int right, int descendingOrder);
void     SortYFlt_Asc (YFlt *item, int left, int right);
void     SortYFlt_Des (YFlt *item, int left, int right);
int      StrCmpCaseInsensitiveLen (const char *s, const char *t, size_t len);
int      StrCmpCaseInsensitive (char *s, char *t);
void     StripComments (char *s);
FILE    *TestOpenTextFileR (char *name);
void     UpdateGrowthFxn (int *growthFxn);
int      UpperTriangIndex (int i, int j, int size);
int      WantTo (const char *msg);

/* tree utility functions */
int       AddToTreeList (TreeList *treeList, Tree *tree);
Tree     *AllocateTree (int numTaxa);
Tree     *AllocateFixedTree (int numTaxa, int isRooted);
int       AllocateTreePartitions (Tree *t);
PolyTree *AllocatePolyTree (int numTaxa);
int       AllocatePolyTreePartitions (PolyTree *pt);
int       AllocatePolyTreeRelClockParams (PolyTree *pt, int nBSets, int nESets);
int       AreTopologiesSame (Tree *t1, Tree *t2);
int       AreTreesSame (Tree *t1, Tree *t2);
int       BuildConstraintTree (Tree *t, PolyTree *pt, char **localTaxonNames);
int       BuildRandomRTopology (Tree *t, RandLong *seed);
int       BuildRandomUTopology (Tree *t, RandLong *seed);
int       CheckConstraints (Tree *t);
int       CheckSetConstraints (Tree *t);
void      ColorClusters (TreeNode *p, int *index);
void      CopySubtreeToTree (Tree *subtree, Tree *t);
int       CopyToPolyTreeFromPolyTree (PolyTree *to, PolyTree *from);
int       CopyToSpeciesTreeFromPolyTree (Tree *to, PolyTree *from);
int       CopyToTreeFromPolyTree (Tree *to, PolyTree *from);
void      CopyPolyNodes (PolyNode *p, PolyNode *q, int nLongsNeeded);
int       CopyToTreeFromTree (Tree *to, Tree *from);
void      CopyTreeNodes (TreeNode *p, TreeNode *q, int nLongsNeeded);
void      CopyTreeToSubtree (Tree *t, Tree *subtree);
int       Deroot (PolyTree *pt);
void      EraseTreeList (TreeList *treeList);
void      findAllowedClockrate (Tree *t, YFlt *minClockRate, YFlt *maxClockRate);
void      FreePolyTree (PolyTree *pt);
void      FreePolyTreePartitions (PolyTree *pt);
void      FreePolyTreePopSizeParams (PolyTree *pt);
void      FreePolyTreeRelClockParams (PolyTree *pt);
void      FreeTree (Tree *t);
void      FreeTreePartitions (Tree *pt);
void      GetDatedNodeDepths (TreeNode *p, YFlt *nodeDepths);
void      GetDatedNodes (TreeNode *p, TreeNode **datedNodes);
void      GetDownPass (Tree *t);
void      GetNodeDownPass (Tree *t, TreeNode *p, int *i, int *j);
void      GetPolyAges (PolyTree *t);
void      GetPolyDepths (PolyTree *t);
void      GetPolyDownPass (PolyTree *t);
void      GetPolyNodeDownPass (PolyTree *t, PolyNode *p, int *i, int *j);
int       GetRandomEmbeddedSubtree (Tree *t, int nTerminals, RandLong *seed, int *nEmbeddedTrees);
int       GetFromTreeList (TreeList *treeList, Tree *tree);
int       InitBrlens (Tree *t, YFlt v);
int       InitCalibratedBrlens (Tree *t, YFlt minLength, RandLong *seed);
int       InitClockBrlens (Tree *t);
int       IsCalibratedClockSatisfied (Tree *t,YFlt *minClockRate,YFlt *maxClockRate , YFlt tol);
int       IsClockSatisfied (Tree *t, YFlt tol);
int       IsTreeConsistent (Param *param, int chain, int state);
int       LabelTree (Tree *t, char **taxonNames);
void      Mark (TreeNode *p);
void      MarkDistance (TreeNode *p, int YESorNO, int dist, int *n);
void      MarkUnconstrained (TreeNode *p);
int       MoveCalculationRoot (Tree *t, int outgroup);
int       MovePolyCalculationRoot (PolyTree *t, int outgroup);
int       NumConstrainedTips (TreeNode *p);
int       NumDatedTips (TreeNode *p);
void      OrderTips (PolyTree *t);
void      PrintNewick (char **s, int *len, Tree *t);
void      PrintNodes (Tree *t);
void      PrintPolyNodes (PolyTree *pt);
void      PrintTranslateBlock (FILE *fp, Tree *t);
int       PrunePolyTree (PolyTree *pt);
int       RandPerturb (Tree *t, int nPert, RandLong *seed);
int       RandResolve (Tree *tt, PolyTree *t, RandLong *seed, int destinationIsRooted);
int       ResetBrlensFromTree (Tree *tree, Tree *vTree);
void      ResetIntNodeIndices (PolyTree *t);
void      ResetPolyTree (PolyTree *t);
void      ResetPolyTreePartitions (PolyTree *pt);
void      ResetPolyTreeRelClockParams (PolyTree *pt);
int       ResetRootHeight (Tree *t, YFlt rootHeight);
void      ResetTipIndices (PolyTree *pt);
int       ResetTopology (Tree *t, char *s);
int       ResetTopologyFromTree (Tree *tree, Tree *top);
int       ResetTopologyFromPolyTree (Tree *tree, PolyTree *top);
void      ResetTreePartitions (Tree *t);
int       RetrieveRTopology (Tree *t, int *order);
int       RetrieveRTree (Tree *t, int *order, YFlt *brlens);
int       RetrieveRTreeWithIndices (Tree *t, int *order, YFlt *brlens);
int       RetrieveUTopology (Tree *t, int *order);
int       RetrieveUTree (Tree *t, int *order, YFlt *brlens);
void      SetDatedNodeAges (Param* param, int chain, int state);
void      SetNodeDepths (Tree *t);
int       SetTreeNodeAges (Param *param, int chain, int state);
int       ShowPolyNodes (PolyTree *pt);
int       ShowTree (Tree *t);
int       StoreRPolyTopology (PolyTree *t, int *order);
int       StoreRPolyTree (PolyTree *t, int *order, YFlt *brlens);
int       StoreRTopology (Tree *t, int *order);
int       StoreRTree (Tree *t, int *order, YFlt *brlens);
int       StoreRTreeWithIndices (Tree *t, int *order, YFlt *brlens);
int       StoreUPolyTopology (PolyTree *t, int *order);
int       StoreUPolyTree (PolyTree *t, int *order, YFlt *brlens);
int       StoreUTopology (Tree *t, int *order);
int       StoreUTree (Tree *t, int *order, YFlt *brlens);
YFlt    TreeLen (Tree *t);
void      Unmark (TreeNode *p);
void      UpdateTreeWithClockrate (Tree *t, YFlt clockRate);
void      WriteEventTree (TreeNode *p, int chain, Param *param);
void      WriteEventTreeToPrintString (TreeNode *p, int chain, Param *param, int printAll);
void      WriteNoEvtTreeToPrintString (TreeNode *p, int chain, Param *param, int showBrlens, int isRooted);
void      WriteEvolTree (TreeNode *p, int chain, Param *param);
void      WriteTopologyToFile (FILE *fp, TreeNode *p, int isRooted);

/* math utility functions */
YComplex **AllocateSquareComplexMatrix (int dim);
YFlt  **AllocateSquareDoubleMatrix (int dim);
int     **AllocateSquareIntegerMatrix (int dim);
int       AutodGamma (YFlt *M, YFlt rho, int K);
void      BetaBreaks (YFlt alpha, YFlt beta, YFlt *values, int K);
YFlt    BetaQuantile (YFlt alpha, YFlt beta, YFlt x);
void      CalcCijk (int dim, YFlt *c_ijk, YFlt **u, YFlt **v);
void      CopyComplexMatrices (int dim, YComplex **from, YComplex **to);
void      CopyDoubleMatrices (int dim, YFlt **from, YFlt **to);
void      DirichletRandomVariable (YFlt *alp, YFlt *z, int n, RandLong *seed);
int       DiscreteGamma (YFlt *rK, YFlt alfa, YFlt beta, int K, int median);
void      FreeSquareComplexMatrix (YComplex **m);
void      FreeSquareDoubleMatrix (YFlt **m);
void      FreeSquareIntegerMatrix (int **m);
int       GetEigens (int dim, YFlt **q, YFlt *eigenValues, YFlt *eigvalsImag, YFlt **eigvecs, YFlt **inverseEigvecs, YComplex **Ceigvecs, YComplex **CinverseEigvecs);
YFlt    LnFactorial (int value);
YFlt    LnGamma (YFlt alp);
YFlt    LnPriorProbExponential (YFlt val, YFlt *params);
YFlt    LnPriorProbExponential_Param_Mean (YFlt val, YFlt *params);
YFlt    LnPriorProbFix (YFlt val, YFlt *params);
YFlt    LnPriorProbGamma (YFlt val, YFlt *params);
YFlt    LnPriorProbGamma_Param_Mean_Sd (YFlt val, YFlt *params);
YFlt    LnPriorProbLognormal (YFlt val, YFlt *params);
YFlt    LnPriorProbLognormal_Param_Mean_Sd (YFlt val, YFlt *params);
YFlt    LnPriorProbNormal (YFlt val, YFlt *params);
YFlt    LnPriorProbOffsetExponential (YFlt val, YFlt *params);
YFlt    LnPriorProbOffsetExponential_Param_Offset_Mean (YFlt val, YFlt *params);
YFlt    LnPriorProbOffsetGamma (YFlt val, YFlt *params);
YFlt    LnPriorProbOffsetGamma_Param_Offset_Mean_Sd (YFlt val, YFlt *params);
YFlt    LnPriorProbOffsetLognormal (YFlt val, YFlt *params);
YFlt    LnPriorProbOffsetLognormal_Param_Offset_Mean_Sd (YFlt val, YFlt *params);
YFlt    LnPriorProbTruncatedNormal (YFlt val, YFlt *params);
YFlt    LnPriorProbTruncatedNormal_Param_Trunc_Mean_Sd (YFlt val, YFlt *params);
YFlt    LnPriorProbUniform (YFlt val, YFlt *params);
YFlt    LnProbRatioExponential (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioExponential_Param_Mean (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioFix (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioGamma (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioGamma_Param_Mean_Sd (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioLognormal (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioLognormal_Param_Mean_Sd (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioNormal (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioOffsetExponential (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioOffsetExponential_Param_Offset_Mean (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioOffsetGamma (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioOffsetGamma_Param_Offset_Mean_Sd (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioOffsetLognormal (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioOffsetLognormal_Param_Offset_Mean_Sd (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioTruncatedNormal (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioTruncatedNormal_Param_Trunc_Mean_Sd (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbRatioUniform (YFlt newX, YFlt oldX, YFlt *params);
YFlt    LnProbGamma (YFlt alpha, YFlt beta, YFlt x);
YFlt    LnProbTruncGamma (YFlt alpha, YFlt beta, YFlt x, YFlt min, YFlt max);
YFlt    LnProbLogNormal (YFlt mu, YFlt sigma, YFlt x);
YFlt    LnProbLogNormal_Mean_LogVar (YFlt mean, YFlt sigma2, YFlt x);
YFlt    LnProbLogNormal_Mean_Var (YFlt mean, YFlt var, YFlt x);
YFlt    LnProbNormal (YFlt mu, YFlt sigma, YFlt x);
YFlt    LnRatioLogNormal (YFlt mu, YFlt sigma, YFlt xNew, YFlt xOld);
YFlt    LogNormalRandomVariable (YFlt mean, YFlt var, RandLong *seed);
YFlt    MaximumValue (YFlt x, YFlt y);
YFlt    MinimumValue (YFlt x, YFlt y);
void      MultiplyMatrices (int dim, YFlt **a, YFlt **b, YFlt **result);
int       MultiplyMatrixNTimes (int dim, YFlt **Mat, int power, YFlt **Result);
YFlt    PointNormal (YFlt prob);
YFlt    PsiGammaLnProb (YFlt alpha, YFlt value);
YFlt    PsiGammaLnRatio (YFlt alpha, YFlt numerator, YFlt denominator);
YFlt    PsiGammaRandomVariable (YFlt alpha, RandLong *seed);
YFlt    QuantileGamma (YFlt x, YFlt alfa, YFlt beta);
YFlt    RandomNumber (RandLong *seed);
YFlt    QuantileLogNormal (YFlt prob, YFlt mu, YFlt sigma);
int       DiscreteLogNormal (YFlt *rK, YFlt sigma, int K, int median);
YFlt    LogNormalPoint (YFlt x, YFlt mu, YFlt sigma);

/* qsort utility function */
int       cmpYFlt(const void *a, const void *b);

#endif  /* UTILS_H_ */

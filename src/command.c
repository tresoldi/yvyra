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
#include "command.h"
#include "model.h"
#include "mcmc.h"
#include "sumpt.h"
#include "diagnostics.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

#define NUMCOMMANDS                     71    /* The total number of commands in the program  */
#define NUMPARAMS                       292   /* The total number of parameters  */
#define PARAM(i, s, f, l)               p->string = s;    \
                                        p->fp = f;        \
                                        p->valueList = l; \
                                        p++;
#define HIDE                            0
#define SHOW                            1

/* Debugging options */
#undef SHOW_TOKENS
#undef ECHO_PROCESSED_COMMANDS

/* Local function prototypes */
int      AddNameSet(NameSet **nameSetList, int numNameSets, char **nameSet, int numNames);
int      AddToSet (int i, int j, int k, int id);
int      AllocCharacters (void);
int      AllocMatrix (void);
int      AllocTaxa (void);
char     ChangeCase (char c);
int      CharacterCode (char ch, int *charCode, int chType);
int      CharacterNumber (int charCode, int chType);
int      CheckInitialPartitions (void);
int      Dex (TreeNode *p);
int      DoAbout (void);
int      DoAcknowledgments (void);
int      DoBeginParm (char *parmName, char *tkn);
int      DoBreaks (void);
int      DoBreaksParm (char *parmName, char *tkn);
int      DoCalibrate (void);
int      DoCalibrateParm (char *parmName, char *tkn);
int      DoCharset (void);
int      DoCharsetParm (char *parmName, char *tkn);
int      DoCharStat (void);
int      DoCitations (void);
int      DoConstraint (void);
int      DoConstraintParm (char *parmName, char *tkn);
int      DoCtype (void);
int      DoCtypeParm (char *parmName, char *tkn);
int      DoUsertype (void);
int      DoUsertypeParm (char *parmName, char *tkn);
int      DoUsercost (void);
int      DoWtset (void);
int      DoWtsetParm (char *parmName, char *tkn);
int      DoTipweights (void);
int      DoTipweightsParm (char *parmName, char *tkn);
int      DoCharlabels (void);
int      DoCharlabelsParm (char *parmName, char *tkn);
int      DoAsrentropy (void);
int      DoConvergence (void);
int      DoDelete (void);
int      DoDeleteParm (char *parmName, char *tkn);
int      DoDimensions (void);
int      DoDimensionsParm (char *parmName, char *tkn);
int      DoDisclaimer (void);
int      DoEndBlock (void);
int      DoExecuteParm (char *parmName, char *tkn);
int      DoExclude (void);
int      DoExcludeParm (char *parmName, char *tkn);
int      DoFormat (void);
int      DoFormatParm (char *parmName, char *tkn);
int      DoHelp (void);
int      DoHelpParm (char *parmName, char *tkn);
int      DoInclude (void);
int      DoIncludeParm (char *parmName, char *tkn);
int      DoLog (void);
int      DoLogParm (char *parmName, char *tkn);
int      DoManual (void);
int      DoManualParm (char *parmName, char *tkn);
int      DoMatrix (void);
int      DoMatrixParm (char *parmName, char *tkn);
int      DoNexusParm (char *parmName, char *tkn);
int      DoOutgroup (void);
int      DoOutgroupParm (char *parmName, char *tkn);
int      DoPairs (void);
int      DoPairsParm (char *parmName, char *tkn);
int      DoPartition (void);
int      DoPartitionParm (char *parmName, char *tkn);
int      DoRestore (void);
int      DoRestoreParm (char *parmName, char *tkn);
int      DoSet (void);
int      DoSetParm (char *parmName, char *tkn);
int      DoShowMatrix (void);
int      DoShowUserTrees (void);
int      DoSpeciespartition (void);
int      DoSpeciespartitionParm (char *parmName, char *tkn);
int      DoTaxaset (void);
int      DoTaxasetParm (char *parmName, char *tkn);
int      DoTaxaStat (void);
int      DoTaxlabels (void);
int      DoTaxlabelsParm (char *parmName, char *tkn);
int      DoTranslate (void);
int      DoTranslateParm (char *parmName, char *tkn);
int      DoTree (void);
int      DoTreeParm (char *parmName, char *tkn);
int      DoUserTree (void);
int      DoUserTreeParm (char *parmName, char *tkn);
int      DoVersion (void);
int      FindValidParam (char *tk, int *numMatches);
int      FreeCharacters (void);
int      FreeMatrix (void);
int      FreeTaxa (void);
int      GetNumPartDivisions (int n);
int      GetUserHelp (char *helpTkn);
int      IsAmbig (int charCode, int dType);
int      IsMissing (int charCode, int dType);
int      MBResID (char nuc);
int      NucID (char nuc);
void     PrintSettings (char *command);
void     PrintYesNo (int yn, char s[4]);
int      ProtID (char aa);
int      SetPartition (int part);
int      SetSpeciespartition (int part);
int      SetTaxaFromTranslateTable (void);
int      StandID (char nuc);
void     WhatVariableExp (BitsLong exp, char *st);
YFlt   WhichCont (int x);

/* globals */
int             autoClose;             /* autoclose                                     */
int             autoOverwrite;         /* Overwrite or append outputfiles when nowarnings=yes */
Calibration     *calibrationPtr;       /* ptr to calibration being set                  */
CharInformation *charInfo;             /* holds critical information about characters   */
UserType        userTypes[MAX_NUM_USERTYPES]; /* user-defined rate matrices               */
int             numUserTypes = 0;      /* number of defined user types                  */
int             currentUserTypeIndex = -1; /* usertype index being applied by ctype     */
YFlt          *tipWeights = NULL;    /* confidence weights for tips                  */
int             hasTipWeights = NO;    /* any confidence-weighted tips?                */
char            **charLabels = NULL;   /* character names from charlabels command      */
char            ***charStateNames = NULL; /* per-char state names from YAML          */
int             *charNStateNames = NULL;  /* per-char count of named states          */
BitsLong        **charSet;             /* holds information about defined charsets      */
char            **charSetNames;        /* holds names of character sets                 */
Comptree        comptreeParams;        /* holds parameters for comparetree command      */
char            **constraintNames;     /* holds names of constraints                    */
int             dataType;              /* type of data                                  */
Calibration     defaultCalibration;    /* default calibration                           */
BitsLong        **definedConstraint;          /* bitfields representing taxa sets of defined constraints                                             */
BitsLong        **definedConstraintTwo;       /* bitfields representing second taxa sets of defined constraints (used for PARTIAL constraints)       */
BitsLong        **definedConstraintPruned;    /* bitfields representing taxa sets of defined constraints after deleted taxa are removed              */
BitsLong        **definedConstraintTwoPruned; /* bitfields representing second taxa sets of defined constraints for PARTIAL constraints after deleted*/
                                              /* taxa are removed and for NEGATIVE constraint it contains complements of definedConstraintPruned     */
int             echoMB;                /* flag used by Manual to prevent echoing        */
BitsLong        expecting;             /* variable denoting expected token type         */
int             foundNewLine;          /* whether a new line has been found             */
int             inComment;             /* flag for whether input stream is commented    */
int             inComparetreeCommand;  /* flag set whenever you enter comparetree cmd   */
int             inferAncStates;        /* should ancestral states be inferred (y/n)     */
int             inferSiteOmegas;       /* should site omegas be inferred (y/n)          */
int             inferSiteRates;        /* should site rates be inferred (y/n)           */
int             inferSiteLikes;        /* should site lnLs be output (y/n)              */
int             inYvyraBlock;        /* flag for whether we are in a yvyra block    */
int             inSumtCommand;         /* flag set whenever you enter sumt cmd          */
int             inTreesBlock;          /* flag for whether we are in a trees block      */
int             inValidCommand;        /* a useful flag set whenever you enter a cmd    */
int             isInAmbig, isInPoly;   /* flags whether we are within () or {}          */
int             isTaxsetDef;           /* is a taxon set defined                        */
int             isTranslateDef;        /* is a translation block defined                */
int             isTranslateDiff;       /* is translate different from current taxaset?  */
char            logFileName[100];      /* name of the log file                          */
int             logToFile;             /* should screen output be logged to a file      */
FILE            *logFileFp;            /* file pointer to log file                      */
int             longIntegerSize;       /* size of an unsigned integer                   */
char            manFileName[100];      /* name of the file for the command help info    */
int             *matrix;               /* matrix containing original data               */
int             matrixHasPoly;         /* flag for whether matrix has polymorphisms     */
int             memAllocs[NUM_ALLOCS]; /* allocated memory flags                        */
int             mode;                  /* mode of program (interactive/noninteractive)  */
Calibration     *nodeCalibration;      /* holds information about node calibrations     */
int             noWarn;                /* no warnings on overwriting files              */
int             numChar;               /* number of characters in character matrix      */
int             numCharSets;           /* number of character sets                      */
int             numComments;           /* counts how deeply nested a comment is         */
int             numDefinedConstraints; /* number of constraints defined                 */
int             numDefinedPartitions;  /* number of partitions defined                  */
int             numDefinedSpeciespartitions;  /* number of speciespartitions defined    */
int             numNamedTaxa;          /* number of named taxa during parsing of cmd    */
int             numOpenExeFiles;       /* number of execute files open                  */
int             numSpecies;            /* number of species in current speciespartition */
int             numTaxa;               /* number of taxa in character matrix            */
int             numTaxaSets;           /* number of taxa sets                           */
int             numTranslates;         /* number of taxa in active translate block      */
int             outGroupNum;           /* number of outgroup taxon                      */
ParmInfo        paramTable[NUMPARAMS]; /* information on parameters                     */
char            **partitionNames;      /* hold names of partitions (first is "default") */
int             **partitionId;         /* holds information about defined partitions    */
int             partitionNum;          /* index of current partition                    */
Plot            plotParams;            /* holds parameters for plot command             */
int             precision;             /* precision of samples and summary stats        */
int             quitOnError;           /* quit on error?                                */
int             replaceLogFile;        /* should logfile be replace/appended to         */
int             scientific;            /* use scientific format for samples ?           */
char            spacer[10];            /* holds blanks for printing indentations        */
NameSet         *speciesNameSets;      /* hold species name sets, one for each speciespartition     */
int             **speciespartitionId;  /* holds info about defined speciespartitions    */
char            **speciespartitionNames;    /* hold names of speciespartitions (first is "default") */
int             speciespartitionNum;   /* index of current speciespartition             */
Sump            sumpParams;            /* holds parameters for sump command             */
Sumt            sumtParams;            /* holds parameters for sumt command             */
Sumss           sumssParams;           /* holds parameters for sumss command            */
TaxaInformation *taxaInfo;             /* holds critical information about taxa         */
char            **taxaNames;           /* holds name of taxa                            */
BitsLong        **taxaSet;             /* holds information about defined taxasets      */
char            **taxaSetNames;        /* holds names of taxa sets                      */
int             *tempActiveConstraints;/* temporarily holds active constraints size allocated        */
enum ConstraintType  *definedConstraintsType;  /* Store type of constraint              */
int             *tempSet;              /* temporarily holds defined set                 */
int             *tempSetNeg;           /* holds bitset of negative set of taxa for partial constraint*/
int             theAmbigChar;          /* int containing ambiguous character            */
Calibration     *tipCalibration;       /* holds information about node calibrations     */
char            **transFrom;           /* translation block information                 */
char            **transTo;             /* translation block information                 */
int             userBrlensDef;         /* are the branch lengths on user tree defined   */


/* local (to this file) */
char            *tokenP, token[CMD_STRING_LENGTH], *cmdStr=NULL;
Calibration     defaultCalibration = {
                    "Unconstrained",      /* name */
                    unconstrained,        /* prior */
                    { -1.0, -1.0, -1.0 }, /* priorParams */
                    NULL,                 /* LnPriorProb */
                    NULL,                 /* LnPriorRatio */
                    -1.0,                 /* min */
                    -1.0                  /* max */
                };

CmdType     commands[] =
            {
            /*  Information on commands initialization:
             
                    1 = Command number (cmdNumber)
                    2 = Command name (string)
                    3 = Special command (YES/NO) (specialCmd) 
                    4 = Pointer to finishing function (fp)
                    5 = Number of valid parameters (numParms)
                    6 = List of valid parameters (parmList) 
                    7 = Expecting (2^TokenType) (expect) (PARAMETER = 4; SEMICOLON = 32; ALPHA = 16384; 
                        ALPHA | QUESTIONMARK | DASH | NUMBER | ASTERISK | EXCLAMATIONMARK | PERCENT | WEIRD | SEMICOLON = 11715360;
                        ALPHA | QUESTIONMARK | DASH | NUMBER | ASTERISK | EXCLAMATIONMARK | PERCENT | WEIRD | VERTICALBAR | SEMICOLON | LEFTPAR | RIGHTPAR | LEFTCURL | RIGHTCURL = 649252640;
                        PARAMETER | SEMICOLON = 36; NUMBER | ALPHA = 49152; ALPHA | SEMICOLON = 16416; EQUALSIGN = 8; NUMBER = 32768)
                    8 = Description of the command (cmdDescription)
                    9 = Where should the command be used (cmdUse) (IN_CMD = used from command line or yvyra block; IN_FILE = used in data block or in tree block)
                   10 = Should the command be shown when "help" is typed (hiding).
             
              #1                 #2   #3                 #4  #5                                                                                                #6        #7                                                            #8       #9   #10
             -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
            {  0,               "#",  NO,              NULL,  1,                                                                                              {0},        4,                                                           "", IN_FILE, HIDE },
            {  1,           "About",  NO,           DoAbout,  0,                                                                                             {-1},       32,                                      "Describes the program",  IN_CMD, SHOW },
            {  2, "Acknowledgments",  NO, DoAcknowledgments,  0,                                                                                             {-1},       32,                              "Shows program acknowledgments",  IN_CMD, SHOW },
            {  3,           "Begin",  NO,              NULL,  6,                                                                              {1,2,3,201,226,227},        4,                         "Denotes beginning of block in file", IN_FILE, SHOW },
            {  4,       "Calibrate",  NO,       DoCalibrate,  1,                                                                                            {119},        4,               "Assigns dates to terminals or interior nodes",  IN_CMD, SHOW },
            {  5,         "Charset",  NO,         DoCharset,  1,                                                                                             {15},        4,                          "Assigns a group of sites to a set",  IN_CMD, SHOW },
            {  6,        "Charstat",  NO,        DoCharStat,  0,                                                                                             {-1},       32,                                 "Shows status of characters",  IN_CMD, SHOW },
            {  7,       "Citations",  NO,       DoCitations,  0,                                                                                             {-1},       32,                   "Citation of program, models, and methods",  IN_CMD, SHOW },
            {  8,     "Comparetree",  NO,     DoCompareTree,  7,                                                                    {127,128,129,130,221,222,223},       36,                     "Compares the trees from two tree files",  IN_CMD, SHOW },
            {  9,      "Constraint",  NO,      DoConstraint,  1,                                                                                             {66},        4,                      "Defines a constraint on tree topology",  IN_CMD, SHOW },
            { 10,           "Ctype",  NO,           DoCtype,  1,                                                                                             {65},        4,                        "Assigns ordering for the characters",  IN_CMD, SHOW },
            { 11,      "Databreaks", YES,          DoBreaks,  1,                                                                                             {93},    32768,           "Defines data breaks for autodiscrete gamma model",  IN_CMD, SHOW },
            { 12,          "Delete", YES,          DoDelete,  1,                                                                                             {47},    49152,                             "Deletes taxa from the analysis",  IN_CMD, SHOW },
            { 13,      "Dimensions",  NO,      DoDimensions,  2,                                                                                            {4,5},        4,                           "Defines size of character matrix", IN_FILE, SHOW },
            { 14,      "Disclaimer",  NO,      DoDisclaimer,  0,                                                                                             {-1},       32,                               "Describes program disclaimer",  IN_CMD, SHOW },
            { 15,             "End",  NO,        DoEndBlock,  0,                                                                                             {-1},       32,                             "Denotes end of a block in file", IN_FILE, SHOW },
            { 16,        "Endblock",  NO,        DoEndBlock,  0,                                                                                             {-1},       32,                 "Alternative way of denoting end of a block", IN_FILE, SHOW },
            { 17,         "Exclude", YES,         DoExclude,  1,                                                                                             {45},    49152,                           "Excludes sites from the analysis",  IN_CMD, SHOW },
            { 18,         "Execute", YES,         DoExecute,  1,                                                                                             {12},    16384,                                            "Executes a file",  IN_CMD, SHOW },
            { 19,          "Format",  NO,          DoFormat,  7,                                                                             {6,7,8,9,10,219,220},        4,                     "Defines character format in data block", IN_FILE, SHOW },
            { 20,            "Help", YES,            DoHelp,  1,                                                                                             {50},    16416,                  "Provides detailed description of commands",  IN_CMD, SHOW },
            { 21,         "Include", YES,         DoInclude,  1,                                                                                             {46},    49152,                                             "Includes sites",  IN_CMD, SHOW },
            { 22,            "Link",  NO,            DoLink, 30,  {55,56,57,58,59,60,61,62,63,72,73,74,75,76,105,118,193,194,195,196,197,242,243,252,253,255,256,
                                                                                                                                                     270,273,274},        4,               "Links parameters across character partitions",  IN_CMD, SHOW },
            { 23,             "Log",  NO,             DoLog,  5,                                                                                 {85,86,87,88,89},        4,                               "Logs screen output to a file",  IN_CMD, SHOW },
            { 24,            "Lset",  NO,            DoLset, 20,                                     {28,29,30,31,32,33,34,40,51,52,53,90,91,131,188,189,276,277,280,282},4,                "Sets the parameters of the likelihood model",  IN_CMD, SHOW },
            { 25,          "Manual",  NO,          DoManual,  1,                                                                                            {126},       36,                  "Prints a command reference to a text file",  IN_CMD, SHOW },
            { 26,          "Matrix", YES,          DoMatrix,  1,                                                                                             {11},649252640,                 "Defines matrix of characters in data block", IN_FILE, SHOW },
            { 27,            "Mcmc",  NO,            DoMcmc, 46,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215},       36,                   "Starts Markov chain Monte Carlo analysis",  IN_CMD, SHOW },
            { 28,           "Mcmcp",  NO,           DoMcmcp, 46,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215},        4,     "Sets parameters of a chain (without starting analysis)",  IN_CMD, SHOW },
            { 29,        "Outgroup", YES,        DoOutgroup,  1,                                                                                             {78},    49152,                                     "Changes outgroup taxon",  IN_CMD, SHOW },
            { 30,           "Pairs", YES,           DoPairs,  1,                                                                                             {92},    32768,        "Defines nucleotide pairs (doublets) for stem models",  IN_CMD, SHOW },
            { 31,       "Partition",  NO,       DoPartition,  1,                                                                                             {16},        4,                              "Assigns a character partition",  IN_CMD, SHOW },
            { 32,            "Plot",  NO,            DoPlot,  6,                                                                        {106,107,108,109,224,225},       36,                        "Plots parameters from MCMC analysis",  IN_CMD, SHOW },
            { 33,           "Prset",  NO,           DoPrset, 44,  {35,36,37,38,39,41,42,43,44,54,64,67,68,69,70,71,77,100,101,102,103,104,110,111,117,120,121,133,
                                                                                                 168,172,173,174,183,184,185,218,241,246,247,251,254,269,271,272},        4,                         "Sets the priors for the parameters",  IN_CMD, SHOW },
            { 34,         "Propset",  NO,         DoPropset,  1,                                                                                            {186},        4,          "Sets proposal probabilities and tuning parameters",  IN_CMD, SHOW },
            { 35,            "Quit",  NO,            DoQuit,  0,                                                                                             {-1},       32,                                          "Quits the program",  IN_CMD, SHOW },
            { 36,          "Report",  NO,          DoReport, 10,                                                        {122,123,124,125,134,135,136,192,217,283},        4,                 "Controls how model parameters are reported",  IN_CMD, SHOW },
            { 37,         "Restore", YES,         DoRestore,  1,                                                                                             {48},    49152,                                              "Restores taxa",  IN_CMD, SHOW },
            { 38,             "Set",  NO,             DoSet, 23,   {13,14,94,145,170,171,179,181,182,216,229,233,234,235,236,238,239,240,245,268,275,278,279},        4,      "Sets run conditions and defines active data partition",  IN_CMD, SHOW },
            { 40,      "Showmatrix",  NO,      DoShowMatrix,  0,                                                                                             {-1},       32,                             "Shows current character matrix",  IN_CMD, SHOW },
            { 41,   "Showmcmctrees",  NO,   DoShowMcmcTrees,  0,                                                                                             {-1},       32,                          "Shows trees used in MCMC analysis",  IN_CMD, SHOW },
            { 42,       "Showmodel",  NO,       DoShowModel,  0,                                                                                             {-1},       32,                                       "Shows model settings",  IN_CMD, SHOW },
            { 43,       "Showmoves",  NO,       DoShowMoves,  1,                                                                                            {180},       36,                              "Shows moves for current model",  IN_CMD, SHOW },
            { 44,      "Showparams",  NO,      DoShowParams,  0,                                                                                             {-1},       32,                          "Shows parameters in current model",  IN_CMD, SHOW },
            { 45,   "Showusertrees",  NO,   DoShowUserTrees,  0,                                                                                             {-1},       32,                                   "Shows user-defined trees",  IN_CMD, SHOW },
            { 46,"Speciespartition",  NO,DoSpeciespartition,  1,                                                                                            {244},        4,                   "Defines a partition of tips into species",  IN_CMD, SHOW },
            { 47,              "Ss",  NO,              DoSs, 50,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215,248,249,250,257},       36,                             "Starts stepping-stone sampling",  IN_CMD, SHOW },
            { 48,             "Ssp",  NO,             DoSsp, 50,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215,248,249,250,257},       36,"Sets parameters of stepping-stone analysis (without starting)",IN_CMD, SHOW },
            { 49,       "Startvals",  NO,       DoStartvals,  1,                                                                                            {187},        4,                         "Sets starting values of parameters",  IN_CMD, SHOW },
            { 50,            "Sump",  NO,            DoSump, 14,                                          {96,97,137,138,139,140,141,161,162,176,178,211,212,231},       36,                   "Summarizes parameters from MCMC analysis",  IN_CMD, SHOW },
            { 51,           "Sumss",  NO,           DoSumSs, 10,                                                        {258,259,260,261,262,263,264,265,266,267},       36,         "Summarizes parameters from stepping-stone analysis",  IN_CMD, SHOW },
            { 52,            "Sumt",  NO,            DoSumt, 21,                {80,81,82,95,146,147,163,164,165,167,175,177,204,205,206,207,208,209,210,230,232},       36,                        "Summarizes trees from MCMC analysis",  IN_CMD, SHOW },
            { 53,        "Taxastat",  NO,        DoTaxaStat,  0,                                                                                             {-1},       32,                                       "Shows status of taxa",  IN_CMD, SHOW },
            { 54,          "Taxset",  NO,         DoTaxaset,  1,                                                                                             {49},        4,                           "Assigns a group of taxa to a set",  IN_CMD, SHOW },
            { 55,       "Taxlabels", YES,       DoTaxlabels,  1,                                                                                            {228},    49152,                                       "Defines taxon labels", IN_FILE, SHOW },
            { 56,       "Translate", YES,       DoTranslate,  1,                                                                                             {83},    49152,                         "Defines alternative names for taxa", IN_FILE, SHOW },
            { 57,            "Tree",  NO,            DoTree,  1,                                                                                             {79},        4,                                             "Defines a tree", IN_FILE, SHOW },
            { 58,          "Unlink",  NO,          DoUnlink, 30,  {55,56,57,58,59,60,61,62,63,72,73,74,75,76,105,118,193,194,195,196,197,242,243,252,253,255,256,
                                                                                                                                                     270,273,274},        4,             "Unlinks parameters across character partitions",  IN_CMD, SHOW },
            { 59,        "Usertree", YES,        DoUserTree,  1,                                                                                            {203},        8,                                 "Defines a single user tree",  IN_CMD, HIDE },
            { 60,         "Version",  NO,         DoVersion,  0,                                                                                             {-1},       32,                                      "Shows program version",  IN_CMD, SHOW },
            { 61,      "Compareref",  NO,     DoCompRefTree,  7,                                                                    {127,128,129,130,221,222,223},       36,                     "Compares the tree to the reference trees",  IN_CMD, HIDE },
            { 62,        "Usertype",  NO,        DoUsertype,  1,                                                                                            {284},        4,                "Defines a custom rate matrix for characters",  IN_CMD, SHOW },
            { 63,        "Usercost",  NO,       DoUsercost,  1,                                                                                            {284},        4,             "Defines a custom cost matrix for characters",  IN_CMD, SHOW },
            { 64,           "Wtset",  NO,          DoWtset,  1,                                                                                            {286},        4,                          "Sets character weights",  IN_CMD, SHOW },
            { 65,      "Tipweights",  NO,    DoTipweights,  1,                                                                                            {287},        4,              "Sets confidence weights for tip states",  IN_CMD, SHOW },
            { 66,        "Chardiag",  NO,      DoChardiag,  1,                                                                                            {288},        4,                   "Character contribution diagnostics",  IN_CMD, SHOW },
            { 67,     "Sensitivity",  NO,   DoSensitivity,  1,                                                                                            {289},        4,              "Weight sensitivity analysis for clades",  IN_CMD, SHOW },
            { 68,      "Charlabels", YES,   DoCharlabels,  1,                                                                                            {290},    32800,                         "Sets character names",  IN_CMD, SHOW },
            { 69,     "Asrentropy",  NO,   DoAsrentropy,  0,                                                                                             {-1},       32,   "Ancestral state entropy from posterior trees",  IN_CMD, SHOW },
            { 70,   "Convergence",  NO, DoConvergence,  0,                                                                                             {-1},       32,   "Comprehensive convergence diagnostics",  IN_CMD, SHOW },
            /* NOTE: If you add a command here, make certain to change NUMCOMMANDS (above, in this file) appropriately! */
            { 999,             NULL,  NO,              NULL,  0,                                                                                             {-1},       32,                                                           "",  IN_CMD, HIDE }  
            };
int                 inDataBlock, inForeignBlock, isInterleaved, isFirstMatrixRead, isFirstInterleavedBlock, 
                    taxonCount, fromI, toJ, everyK, foundDash, foundSlash, foundFirst, isMixed, whichPartition,
                    isNegative, numDivisions, charOrdering, foundExp, foundColon, isFirstNode, nextAvailableNode,
                    pairId, firstPair, inTaxaBlock, inCharactersBlock, foundEqual;
char                gapId, missingId, matchId, tempSetName[100], **tempNames;
CmdType             *commandPtr; /* Points to the commands array entry which corresponds to currently processed command */
ParmInfoPtr         paramPtr;    /* Points to paramTable table array entry which corresponds to currently processed parameter of current command */
TreeNode            *pPtr, *qPtr;

enum ConstraintType constraintType; /* Used only in processing of constraint command to indicate the type of constraint */


int AddToGivenSet (int i, int j, int k, int id, int *Set)
{
    int     m, n;
    
    if (id <= 0)
        {
        YvyraPrint ("%s   The id for a temporary set should be greater than 0\n", spacer);
        return (ERROR);
        }
    
    if (i < 0 && j < 0)
        return (ERROR);
    else if (i < 0 && j >= 0)
        return (ERROR);
    else if (i >= 0 && j < 0)
        {
        if (k >= 0)
            return (ERROR);
        else
            {
            if (Set[i] != 0)
                {
                YvyraPrint ("%s   Character %d defined more than once\n", spacer, i+1);
                return (ERROR);
                }
            Set[i] = id;
            }
        }
    else if (i >= 0 && j >= 0)
        {
        if (k < 0)
            {
            for (m=i; m<=j; m++)
                {
                if (Set[m] != 0)
                    {
                    YvyraPrint ("%s   Character %d defined more than once\n", spacer, m+1);
                    return (ERROR);
                    }
                Set[m] = id;
                }
            }
        else
            {
            n = k;
            for (m=i; m<=j; m++)    
                {
                if (n % k == 0)
                    {
                    if (Set[m] != 0)
                        {
                        YvyraPrint ("%s   Character %d defined more than once\n", spacer, m+1);
                        return (ERROR);
                        }
                    Set[m] = id;
                    }
                n++;
                }
            }
        }

    return (NO_ERROR);
    
}


int AddToSet (int i, int j, int k, int id)
{
    return AddToGivenSet (i, j, k,id, tempSet);
}


/* AddNameSet: Push a name set onto the end of a list of name sets, with reallocation
      of list to hold the extra element. The calling function needs to keep track of
      the counter holding the length of the list. */
int AddNameSet (NameSet **nameSetList, int numNameSets, char **nameSet, int numNames)
{
    int     i;

    (*nameSetList) = (NameSet*) SafeRealloc ((void*)(*nameSetList), ((size_t)numNameSets+1)*sizeof(NameSet));

    (*nameSetList)[numNameSets].names    = NULL;
    (*nameSetList)[numNameSets].numNames = numNames;

    for (i=0; i<numNames; i++)
        AddString(&((*nameSetList)[numNameSets].names), i, nameSet[i]);
    
    return NO_ERROR;
}


/* AddString: Push a string onto the end of a list, with reallocation of list
      to hold the extra element. The calling function needs to keep track of
      the counter holding the length of the list. */
int AddString (char ***list, int len, char *token)
{
    (*list) = (char **) SafeRealloc ((void *)(*list), ((size_t)len+1)*sizeof(char*));
    if (!(*list))
        return ERROR;

    (*list)[len] = (char *) SafeCalloc ((strlen(token)+1), sizeof(char));
    if (!(*list)[len])
        return ERROR;

    strcpy ((*list)[len], token);
    
    return NO_ERROR;
}


int AllocCharacters (void)
{
    int     i, tempSetSize;

    if (memAllocs[ALLOC_MATRIX] == YES)
        goto errorExit;
    matrix = (int *) SafeMalloc((size_t)numTaxa * (size_t)numChar * sizeof(int));
    if (!matrix)
        {
        YvyraPrint ("%s   Problem allocating matrix (%d)\n", spacer, numTaxa * numChar * sizeof(int));
        goto errorExit;
        }
    for (i=0; i<numTaxa * numChar; i++)
        matrix[i] = 0;
    memAllocs[ALLOC_MATRIX] = YES;

    if (memAllocs[ALLOC_CHARINFO] == YES)
        goto errorExit;
    charInfo = (CharInformation *) SafeMalloc ((size_t)numChar * sizeof(CharInformation));
    if (!charInfo)
        {
        YvyraPrint ("%s   Problem allocating charInfo (%d)\n", spacer, numChar * sizeof(CharInformation));
        goto errorExit;
        }
    for (i=0; i<numChar; i++)
        {
        charInfo[i].isExcluded = NO;
        charInfo[i].numStates = 0;
        charInfo[i].charType = 0;
        charInfo[i].isMissAmbig = NO;
        charInfo[i].ctype = UNORD;
        charInfo[i].userTypeIndex = -1;
        charInfo[i].charId = 0;
        charInfo[i].pairsId = 0;
        charInfo[i].bigBreakAfter = NO;
        charInfo[i].weight = 1.0;
        }
    memAllocs[ALLOC_CHARINFO] = YES;

    if (memAllocs[ALLOC_CHARSETS] == YES)
        goto errorExit;
    charSetNames = NULL;
    charSet = NULL;
    numCharSets = 0;
    memAllocs[ALLOC_CHARSETS] = YES;    /* safe to do free */

    if (memAllocs[ALLOC_PARTITIONS] == YES)
        goto errorExit;
    partitionNames = NULL;
    partitionId = (int**) SafeMalloc ((size_t)numChar * sizeof(int*));
    for (i=0; i<numChar; i++)
        partitionId[i] = (int *) SafeMalloc (sizeof(int));
    numDefinedPartitions = 0;   /* number of defined partitions */
    memAllocs[ALLOC_PARTITIONS] = YES;  /* safe to do free */

    if (memAllocs[ALLOC_PARTITIONVARS] == YES)
        goto errorExit;
    numVars           = NULL;
    tempLinkUnlinkVec = NULL;
    activeParts       = NULL;
    tempLinkUnlinkVec = NULL;
    tempNum           = NULL;
    linkTable[0]      = NULL;
    tempLinkUnlink[0] = NULL;
    for (i=0; i<NUM_LINKED; i++)
        {
        linkTable[i]      = NULL;
        tempLinkUnlink[i] = NULL;
        activeParams[i]   = NULL;
        }
    memAllocs[ALLOC_PARTITIONVARS] = YES;

    if (memAllocs[ALLOC_TMPSET] == NO)
        goto errorExit;
    if (numChar > numTaxa)
        tempSetSize = numChar;
    else
        tempSetSize = numTaxa;
    tempSet = (int *) SafeRealloc ((void *)tempSet, (size_t)tempSetSize * sizeof(int));
    tempSetNeg = (int *) SafeRealloc ((void *)tempSetNeg, (size_t)tempSetSize * sizeof(int));
    if (!tempSet || !tempSetNeg)
        {
        YvyraPrint ("%s   Problem reallocating tempSet (%d)\n", spacer, tempSetSize * sizeof(int));
        goto errorExit;
        }

    YvyraPrint ("%s   Allocated matrix\n", spacer);
    return (NO_ERROR);

    errorExit:
        YvyraPrint ("%s   Problem allocating matrix\n", spacer);
        FreeMatrix();
        return (ERROR);
}


int AllocMatrix (void)
{
    if (memAllocs[ALLOC_TAXA] == NO && AllocTaxa() == ERROR)
        return ERROR;
    else
        return (AllocCharacters());
}


int AllocTaxa (void)
{
    int             i;

    if (defTaxa == NO)
        {
        YvyraPrint ("%s   Number of taxa not defined\n", spacer);
        return (ERROR);
        }

    if (numTaxa == 0)
        {
        YvyraPrint ("%s   Number of taxa is 0\n", spacer);
        return (ERROR);
        }

    /* allocate space for taxa */
    if (memAllocs[ALLOC_TAXA] == YES)
        goto errorExit;

    taxaNames = NULL;           /* This variable is allocated in AddString */
    taxaInfo =
        (TaxaInformation *) SafeMalloc ((size_t) numTaxa *
                                        sizeof (TaxaInformation));

    if (!taxaInfo)
        goto errorExit;

    tipCalibration =
        (Calibration *) SafeMalloc ((size_t) numTaxa * sizeof (Calibration));

    if (!tipCalibration)
        {
        free (taxaInfo);
        taxaInfo = NULL;
        goto errorExit;
        }

    for (i = 0; i < numTaxa; i++)
        {
        taxaInfo[i].isDeleted = NO;
        taxaInfo[i].charCount = 0;
        }

    memAllocs[ALLOC_TAXA] = YES;

    /* taxa sets */
    if (memAllocs[ALLOC_TAXASETS] == YES)
        goto errorExit;

    taxaSetNames = NULL;
    taxaSet = NULL;
    numTaxaSets = 0;
    memAllocs[ALLOC_TAXASETS] = YES;    /* safe to free */

    /* species partitions; allocate space and set default species partition */
    if (memAllocs[ALLOC_SPECIESPARTITIONS] == YES)
        goto errorExit;

    speciespartitionNames = NULL;
    speciesNameSets = NULL;
    speciespartitionId =
        (int **) SafeMalloc ((size_t) numTaxa * sizeof (int *));

    for (i = 0; i < numTaxa; i++)
        {
        speciespartitionId[i] = (int *) SafeMalloc (sizeof (int));
        speciespartitionId[i][0] = i + 1;   /* 1-based taxon index, do not ask me why */
        }

    numDefinedSpeciespartitions = 0;    /* number of defined species partitions */
    memAllocs[ALLOC_SPECIESPARTITIONS] = YES;   /* safe to do free */

    /* constraints */
    if (memAllocs[ALLOC_CONSTRAINTS] == YES)
        goto errorExit;

    constraintNames = NULL;
    definedConstraintsType = NULL;
    definedConstraint = NULL;
    definedConstraintTwo = NULL;
    definedConstraintPruned = NULL;
    definedConstraintTwoPruned = NULL;
    numDefinedConstraints = 0;
    tempActiveConstraints = NULL;
    memAllocs[ALLOC_CONSTRAINTS] = YES; /* safe to free */

    /* translate table */
    transFrom = NULL;
    transTo = NULL;
    numTranslates = 0;

    /* tempSet */
    if (memAllocs[ALLOC_TMPSET] == YES)
        goto errorExit;

    tempSet = (int *) SafeMalloc ((size_t) numTaxa * sizeof (int));
    tempSetNeg = (int *) SafeMalloc ((size_t) numTaxa * sizeof (int));

    if (!tempSet || !tempSetNeg)
        goto errorExit;

    memAllocs[ALLOC_TMPSET] = YES;

    /* make sure previous user trees are freed */
    if (numUserTrees > 0)
        {
        YvyraPrint ("%s   Previous user trees not freed\n", spacer);
        goto errorExit;
        }

    YvyraPrint ("%s   Allocated taxon set\n", spacer);
    return NO_ERROR;

errorExit:
    YvyraPrint ("%s   Problem allocating taxon set\n", spacer);
    FreeTaxa();
    return ERROR;
}


char ChangeCase (char c)
{
    int     x;
    
    x = tolower(c);
    return (x);
}


int CharacterCode (char ch, int *charCode, int chType)
{
    if (chType == DNA || chType == RNA)
        {
        if ((*charCode = NucID (ch)) == -1)
            {
            YvyraPrint ("%s   Unrecognized DNA/RNA character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == PROTEIN)
        {
        if ((*charCode = ProtID (ch)) == -1)
            {
            YvyraPrint ("%s   Unrecognized Protein character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == RESTRICTION)
        {
        if ((*charCode = MBResID (ch)) == -1)
            {
            YvyraPrint ("%s   Unrecognized Restriction character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == STANDARD)
        {
        if ((*charCode = StandID (ch)) == -1)
            {
            YvyraPrint ("%s   Unrecognized Standard character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == CONTINUOUS)
        {
        YvyraPrint ("%s   CharacterCode function cannot check continuous characters\n", spacer);
        }
    else
        {
        YvyraPrint ("%s   Unrecognized character type (%d)\n", spacer, chType);
        return (ERROR);
        }
        
    return (NO_ERROR);
}


int CharacterNumber (int charCode, int chType)
{
    int i, x = charCode;
    
    if (chType == CONTINUOUS)
        return 0;

    for (i=0; x!=0; i++)
        x >>= 1;

    return (i);
}


int CheckInitialPartitions (void)
{
    int     i;
    
    for (i=0; i<numChar; i++)
        {
        if (partitionId[i][0] <= 0 || partitionId[i][0] > numDivisions)
            {
            YvyraPrint ("%s   The partition for site %d is incorrect\n", spacer, i+1); 
            return (ERROR);
            }
        }
        
    return (NO_ERROR);
}


int CheckStringValidity (char *s)
{
    int         i, numUnknownChars, tempNumComments, tempInComment;
    char        temp[100];

    i = 0;
    numUnknownChars = 0;
    tempNumComments = numComments;
    tempInComment = inComment;

    while (s[i] != '\0')
        {
        if (tempInComment == NO)
            {
            if (!IsIn(s[i],"=abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.;:,#()[]?-*/'\\'!%\"&~+^$@|{}`>< "))
                {
                if (IsWhite(s[i]) == 1 || IsWhite(s[i]) == 2)
                    {
                    
                    }
                else
                    {
                    if (commandPtr == NULL) 
                        return (ERROR);
                    YvyraPrint ("%s   Unknown character \"%c\" (ASCII code %d)\n", spacer, s[i], s[i]);
                    if (!strcmp(commandPtr->string,"Matrix"))
                        {
                        if (foundNewLine == NO)
                            {
                            YvyraPrint ("%s   The error is in character %d for taxon %s\n", spacer, taxaInfo[taxonCount-1].charCount+i+1, "???"); /* bug? */
                            }
                        else
                            {
                            if (taxonCount == 0)
                                YvyraPrint ("%s   The error is in the first taxon name\n", spacer);
                            else
                                {
                                strcpy(temp, taxaNames[taxonCount]);
                                if (isInterleaved == NO)
                                    YvyraPrint ("%s   The error is in the name of the taxon following taxon %s\n", spacer, temp);
                                else
                                    {
                                    YvyraPrint ("%s   The error is in the name of the taxon following taxon %s\n", spacer, temp);
                                    YvyraPrint ("%s   in one of the interleaved data blocks\n", spacer);
                                    }
                                }
                            }
                        }
                    else if (!strcmp(commandPtr->string,"Execute"))
                        {
                        YvyraPrint ("%s   Assuming irrelevant characters at beginning of file; processing continues\n", spacer);
                        return (NO_ERROR);
                        }
                    return (ERROR);
                    }
                }
            if (s[i]=='[')
                {
                tempInComment = YES;
                tempNumComments++;
                }
            }
        else if (tempInComment == YES)
            {
            if (s[i]==']')
                {
                tempNumComments--;
                if (tempNumComments == 0)
                    tempInComment = NO;
                }
            }
        i++;
        }
        
    if (numUnknownChars > 0)
        return (ERROR);
    else
        return (NO_ERROR);
}


/* CheckString: This function simply checks a vector of strings for a match against token.
          Upon return, matchIndex contains the index of the matched string. An
          ERROR is returned if there are no matches.  */
int CheckString (char **list, int len, char *token, int *matchIndex)
{
    int         i;      
        
    *matchIndex = -1;
    for (i=0; i<len; i++)
        {
        if (StrCmpCaseInsensitive(token,list[i]) == 0)
            {
            *matchIndex = i;
            return (NO_ERROR);
            }
        }

    return (ERROR); 
}


int Dex (TreeNode *p)
{
    return (p == NULL) ? -1 : p->index;
}


int DoAbout (void)
{
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");
    YvyraPrint ("   About yvyra                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra is a program for Bayesian phylogenetic inference of linguistic and      \n");
    YvyraPrint ("   cultural data. It is a fork of MrBayes 3.2.7a (Ronquist et al. 2012),        \n");
    YvyraPrint ("   stripped of molecular models and extended with features for discrete          \n");
    YvyraPrint ("   morphological and linguistic characters.                                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   Extensions over MrBayes:                                                      \n");
    YvyraPrint ("     - Asymmetric rate matrices via usercost/usertype commands                   \n");
    YvyraPrint ("     - Character weights (wtset) and confidence-weighted tips (tipweights)       \n");
    YvyraPrint ("     - Per-model-group ascertainment correction for heterogeneous models         \n");
    YvyraPrint ("     - Post-processing diagnostics: chardiag, sensitivity, asrentropy            \n");
    YvyraPrint ("     - Per-site log-likelihoods in dedicated .slk file with pattern weights      \n");
    YvyraPrint ("     - Character labels (charlabels) for diagnostic output                       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra uses Markov chain Monte Carlo (MCMC) to approximate the posterior      \n");
    YvyraPrint ("   probability distribution of trees. Trees sampled by MCMC are saved in a      \n");
    YvyraPrint ("   \".t\" file; model parameters are saved in a \".p\" file; per-site log-      \n");
    YvyraPrint ("   likelihoods (if requested) go to a \".slk\" file. Summarize results with     \n");
    YvyraPrint ("   the \"sumt\" and \"sump\" commands.                                           \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra was created by Tiago Tresoldi in 2026, building on the MrBayes         \n");
    YvyraPrint ("   codebase originally written by John P. Huelsenbeck and Fredrik Ronquist.     \n");
    YvyraPrint ("   Type 'citations' for papers to cite, and 'acknowledgments' for the full      \n");
    YvyraPrint ("   list of MrBayes contributors.                                                \n");
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoAcknowledgments (void)
{
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");
    YvyraPrint ("   Acknowledgments                                                               \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra is built on MrBayes, and we gratefully acknowledge all who contributed  \n");
    YvyraPrint ("   to that project over its 20+ year history.                                    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   MrBayes was originally written by John P. Huelsenbeck and Fredrik Ronquist.   \n");
    YvyraPrint ("   JPH and FR thank Gautam Altekar, Andrea Betancourt, Jon Bollback, Barry       \n");
    YvyraPrint ("   Hall, Jimmy McGuire, Rasmus Nielsen, David Swofford, Johan Nylander,          \n");
    YvyraPrint ("   Mikael Thollesson, and Derrick Zwickl for help during initial development.    \n");
    YvyraPrint ("   Important contributions came from Clemens Lakner, Sebastian Hoehna, Paul      \n");
    YvyraPrint ("   Lewis, Mark Holder, Julian Catchen, Bret Larget, Marc Suchard, Daniel         \n");
    YvyraPrint ("   Ayres, and Aaron Darling. Bug fixes and support by Paul van der Mark          \n");
    YvyraPrint ("   (2005-2007), Maxim Teslenko (2010-2012), and Chi Zhang (2012-2015).           \n");
    YvyraPrint ("   From 2015, Andreas Kahari and Johan Nylander at NBIS (https://nbis.se)       \n");
    YvyraPrint ("   maintained the MrBayes codebase.                                              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   JPH was supported by NSF grants DEB-007540 and MCB-0075404 and a Wenner-      \n");
    YvyraPrint ("   Gren scholarship. FR was supported by the Swedish Research Council.            \n");
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoBeginParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        /* set Data (inDataBlock) *************************************************************/
        if (!strcmp(parmName, "Data"))
            {
            if (FreeModel () == ERROR)
                return (ERROR);
            if (FreeMatrix () == ERROR)
                return (ERROR);
            YvyraPrint ("   Reading data block\n");
            inDataBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Characters (inCharactersBlock) *************************************************************/
        else if (!strcmp(parmName, "Characters"))
            {
            if (FreeModel () == ERROR)
                return (ERROR);
            if (FreeCharacters () == ERROR)
                return (ERROR);
            YvyraPrint ("   Reading characters block\n");
            inCharactersBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Taxa (inTaxaBlock) *************************************************************/
        else if (!strcmp(parmName, "Taxa"))
            {
            if (FreeModel () == ERROR)
                return (ERROR);
            if (FreeMatrix () == ERROR)
                return (ERROR);
            YvyraPrint ("   Reading taxa block\n");
            inTaxaBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Yvyra block (inYvyraBlock) — accept both block names ********************/
        else if (!strcmp(parmName, "Mrbayes") || IsSame("Yvyra", tkn) != DIFFERENT)
            {
            YvyraPrint ("   Reading yvyra block\n");
            inYvyraBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Trees (inTreesBlock) *******************************************************/
        else if (!strcmp(parmName, "Trees"))
            {
            YvyraPrint ("   Reading trees block\n");
            inTreesBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Foreign (inForeignBlock) *******************************************************/
        else
            {
            YvyraPrint ("   Skipping \"%s\" block\n", tkn);
            inForeignBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "");
            }
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoBreaks (void)
{
    int         i, numBreaks;
    
    numBreaks = 0;
    for (i=0; i<numChar; i++)
        {
        if (charInfo[i].bigBreakAfter == YES)
            {
            numBreaks++;
            }
        }
        
    if (numBreaks > 0)
        {
        if (numBreaks == 1)
            YvyraPrint ("%s   One data break found after character ", spacer, numBreaks);
        else
            YvyraPrint ("%s   %d data breaks found after characters: ", spacer, numBreaks);
        for (i=0; i<numChar; i++)
            {
            if (charInfo[i].bigBreakAfter == YES)
                {
                YvyraPrint ("%d ", i+1);
                }
            }
        YvyraPrint ("\n");

        if (numBreaks == 1)
            YvyraPrint ("%s   Successfully defined one break in data\n", spacer);
        else
            YvyraPrint ("%s   Successfully defined %d breaks in data\n", spacer, numBreaks);
        }
    else
        {
        YvyraPrint ("%s   No breaks in data found\n", spacer);
        }
        
    return (NO_ERROR);
}


int DoBreaksParm (char *parmName, char *tkn)
{
    int     i, tempInt;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can define breaks in the data\n", spacer);
        return (ERROR);
        }
            
    if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            for (i=0; i<numChar; i++)
                charInfo[i].bigBreakAfter = NO;
            return (ERROR);
            }
        if (tempInt == numChar)
            {
            YvyraPrint ("%s   Character number %d is the last character. yvyra will define the\n", spacer, tempInt);
            YvyraPrint ("%s   break, even though it doesn't make too much sense.\n", spacer);
            }
        tempInt--;
                    
        charInfo[tempInt].bigBreakAfter = YES;
        
        expecting  = (Expecting(NUMBER) | Expecting(SEMICOLON));
        }
    else
        {
        for (i=0; i<numChar; i++)
            charInfo[i].bigBreakAfter = NO;
        return (ERROR);
        }

    return (NO_ERROR);
}


int DoCalibrate (void)
{
    int         i;

    /* show calibration times (for debugging) */
#   if 0
    YvyraPrint ("Taxon ages\n");
    for (i=0; i<numTaxa; i++)
        YvyraPrint ("%4d  --  %s\n", i+1, tipCalibration[i].name);
    YvyraPrint ("Constraint ages\n");
    for (i=0; i<numDefinedConstraints; i++)
        {
        if (definedConstraintsType[i] != HARD)
            continue;
        YvyraPrint ("%4d  --  %s\n", i+1, nodeCalibration[i].name);
        }
#   endif

    /* Update model if calibrations enforced */
    for (i=0; i<numCurrentDivisions; i++)
        {
        if (!strcmp(modelParams[i].nodeAgePr,"Calibrated"))
            {
            if (SetUpAnalysis (&globalSeed) == ERROR)
                return (ERROR);
            break;
            }
        }

    return (NO_ERROR);
}


int DoCalibrateParm (char *parmName, char *tkn)
{
    static int              isTaxon, paramIndex;
    static char             nodeName[100], calName[100];
    static YFlt           priorParams[3];
    static enum CALPRIOR    calPrior;
    int                     howMany, index;
    char                    s[20], tempStr[100];
    YFlt                  tempD;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can calibrate nodes\n", spacer);
        return (ERROR);
        }
        
    if (expecting == Expecting(PARAMETER))
        {
        if (strcmp(parmName, "Xxxxxxxxxx") != 0)
            {
            YvyraPrint ("%s   Unexpected error - Wrong parmName in DoCalibrateParm\n", spacer);
            return (ERROR);
            }

        /* find taxon with this name */
        calibrationPtr = NULL;
        howMany = 0;

        /* first look in constraint names */
        if (CheckString (constraintNames, numDefinedConstraints, tkn, &index) != ERROR && definedConstraintsType[index] == HARD)
            {
            calibrationPtr = &nodeCalibration[index];
            howMany++;
            isTaxon = NO;
            strcpy (nodeName, tkn);
            }
        
        /* then look in terminal taxon names */
        if (CheckString (taxaNames, numTaxa, tkn, &index) != ERROR)
            {
            calibrationPtr = &tipCalibration[index];
            howMany++;
            isTaxon = YES;
            strcpy (nodeName, tkn);
            }

        /* return error if not found or ambiguous */
        if (howMany == 0)
            {
            YvyraPrint ("%s   No taxon or hard constraint named ""%s"" found. Note that only hard constraint can be calibrated.\n", spacer, tkn);
            return (ERROR);
            }
        else if (howMany > 1)
            {
            YvyraPrint ("%s   Both a taxon and a constraint named ""%s"" encountered -- please rename one\n", spacer, tkn);
            return (ERROR);
            }

        /* get ready to find the equal sign */
        expecting = Expecting(EQUALSIGN);
        }

    else if (expecting == Expecting(EQUALSIGN))
        {
        /* get ready to find the calibration prior */
        expecting = Expecting(ALPHA);
        }

    else if (expecting == Expecting(ALPHA))
        {
        /* set the calibration prior type */
        if (IsArgValid(tkn,tempStr) == NO_ERROR)
            {
            if (!strcmp (tempStr, "Unconstrained"))
                calPrior = unconstrained;
            else if (!strcmp (tempStr, "Fixed"))
                calPrior = fixed;
            else if (!strcmp (tempStr, "Uniform"))
                calPrior = uniform;
            else if (!strcmp (tempStr, "Offsetexponential"))
                calPrior = offsetExponential;
            else if (!strcmp (tempStr, "Truncatednormal"))
                calPrior = truncatedNormal;
            else if (!strcmp (tempStr, "Lognormal"))
                calPrior = logNormal;
            else if (!strcmp (tempStr, "Offsetlognormal"))
                calPrior = offsetLogNormal;
            else if (!strcmp (tempStr, "Gamma"))
                calPrior = standardGamma;
            else if (!strcmp (tempStr, "Offsetgamma"))
                calPrior = offsetGamma;

            if (calPrior == unconstrained)
                {
                /* reset the values of the calibration */
                YvyraPrint ("%s   Resetting previous calibration for ""%s""\n", spacer, nodeName);

                calibrationPtr->prior           = defaultCalibration.prior;
                calibrationPtr->priorParams[0]  = defaultCalibration.priorParams[0];
                calibrationPtr->priorParams[1]  = defaultCalibration.priorParams[1];
                calibrationPtr->priorParams[2]  = defaultCalibration.priorParams[2];
                calibrationPtr->LnPriorProb     = defaultCalibration.LnPriorProb;
                calibrationPtr->LnPriorRatio    = defaultCalibration.LnPriorRatio;
                calibrationPtr->min             = defaultCalibration.min;
                calibrationPtr->max             = defaultCalibration.max;
                strcpy(calibrationPtr->name, defaultCalibration.name);
            
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                strcpy (calName, tempStr);
                paramIndex = 0;
                priorParams[0] = priorParams[1] = priorParams[2] =  -1.0;
                expecting = Expecting(LEFTPAR);
                }
            }
        else
            {
            YvyraPrint ("%s   Invalid calibration prior argument \n", spacer);
            return (ERROR);
            }
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        strcat (calName, "(");
        expecting  = Expecting(NUMBER);
        }
    else if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%lf", &tempD);
        if (paramIndex == 0)
            {
            if (calPrior == logNormal)
                {
                if (tempD < 0.0)
                    {
                    YvyraPrint ("%s   Mean age must be nonnegative\n", spacer);
                    YvyraPrint ("%s   Parameters of the lognormal distribution used for dating are mean age\n", spacer);
                    YvyraPrint ("%s   and standard deviation, both specified on the linear scale, not as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == standardGamma)
                {
                if (tempD <= 0.0)
                    {
                    YvyraPrint ("%s   Mean parameter must be positive\n", spacer);
                    YvyraPrint ("%s   Parameters of the gamma distribution used for dating are mean age and\n", spacer);
                    YvyraPrint ("%s   standard deviation. In terms of the common shape (alpha) and rate (beta)\n", spacer);
                    YvyraPrint ("%s   parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    YvyraPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }
            else if (tempD < 0.0)
                {
                if (calPrior == fixed)
                    YvyraPrint ("%s   Fixed age must be nonnegative\n", spacer);
                else if (calPrior == uniform)
                    {
                    YvyraPrint ("%s   Minimum age must be nonnegative\n", spacer);
                    YvyraPrint ("%s   Parameters of the uniform are minimum age and maximum age.\n", spacer);
                    }
                else if (calPrior == truncatedNormal)
                    {
                    YvyraPrint ("%s   Offset (minimum or truncation) age must be nonnegative.\n", spacer);
                    YvyraPrint ("%s   Parameters of the truncated normal distribution are offset (minimum\n", spacer);
                    YvyraPrint ("%s   or truncation) age, mean age and standard deviation.\n", spacer);
                    }
                else if (calPrior == offsetGamma)
                    {
                    YvyraPrint ("%s   Offset age must be nonnegative\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset gamma distribution used for dating are offset age,\n", spacer);
                    YvyraPrint ("%s   mean age, and standard deviation. In terms of the common shape (alpha) and\n", spacer);
                    YvyraPrint ("%s   rate (beta) parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    YvyraPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    }
                else if (calPrior == offsetExponential)
                    {
                    YvyraPrint ("%s   Offset age must be nonnegative\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset exponential are offset age and mean age.\n", spacer);
                    }
                else if (calPrior == offsetLogNormal)
                    {
                    YvyraPrint ("%s   Offset age must be nonnegative\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset lognormal distribution are offset age, mean age,\n", spacer);
                    YvyraPrint ("%s   and standard deviation. All values are specified on the linear scale, not\n", spacer);
                    YvyraPrint ("%s   as log values.\n", spacer);
                    }
                return (ERROR);
                }
            priorParams[0] = tempD;
            if (calPrior == fixed)
                expecting = Expecting(RIGHTPAR);
            else
                expecting = Expecting(COMMA);
            }
        else if (paramIndex == 1)
            {
            if (calPrior == uniform)
                {
                if (tempD <= priorParams[0])
                    {
                    YvyraPrint ("%s   Maximum age of uniform distribution must be larger than minimum age\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetExponential)
                {
                if (tempD <= priorParams[0])
                    {
                    YvyraPrint ("%s   Mean age must be larger than offset age.\n", spacer);
                    YvyraPrint ("%s   yvyra now uses offset and mean rather than offset and rate\n", spacer);
                    YvyraPrint ("%s   as the parameters for the offset exponential distribution.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == truncatedNormal)
                {
                if (tempD <= priorParams[0])
                    {
                    YvyraPrint ("%s   Mean age must be larger than offset (truncation) age.\n", spacer);
                    YvyraPrint ("%s   Parameters of the truncated normal distribution are offset (minimum\n", spacer);
                    YvyraPrint ("%s   or truncation) age, mean age and standard deviation\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == logNormal)
                {
                if (tempD <= 0.0)
                    {
                    YvyraPrint ("%s   Standard deviation must be positive.\n", spacer);
                    YvyraPrint ("%s   Parameters of the lognormal distribution used for dating are mean age\n", spacer);
                    YvyraPrint ("%s   and standard deviation, both specified on the linear scale, not as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetLogNormal)
                {
                if (tempD <= priorParams[0])
                    {
                    YvyraPrint ("%s   Mean age must be larger than offset age.\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset lognormal distribution are offset age, mean age,\n", spacer);
                    YvyraPrint ("%s   and standard deviation. All values are specified on the linear scale, not\n", spacer);
                    YvyraPrint ("%s   as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == standardGamma)
                {
                if (tempD <= 0.0)
                    {
                    YvyraPrint ("%s   Standard deviation must be positive.\n", spacer);
                    YvyraPrint ("%s   Parameters of the gamma distribution used for dating are mean age and\n", spacer);
                    YvyraPrint ("%s   standard deviation. In terms of the common shape (alpha) and rate (beta)\n", spacer);
                    YvyraPrint ("%s   parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    YvyraPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetGamma)
                {
                if (tempD <= 0.0)
                    {
                    YvyraPrint ("%s   Mean age must be positive.\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset gamma distribution used for dating are offset age,\n", spacer);
                    YvyraPrint ("%s   mean age, and standard deviation. In terms of the common shape (alpha) and\n", spacer);
                    YvyraPrint ("%s   rate (beta) parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    YvyraPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }

            priorParams[1] = tempD;
            if (calPrior == uniform || calPrior == standardGamma || calPrior == logNormal || calPrior == offsetExponential)
                expecting = Expecting(RIGHTPAR);
            else
                expecting = Expecting(COMMA);
            }
        else /* if (paramIndex == 2) */
            {
            if (calPrior == offsetGamma)
                {
                if (tempD <= 0.0)
                    {
                    YvyraPrint ("%s   Standard deviation must be positive.\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset gamma distribution used for dating are offset age,\n", spacer);
                    YvyraPrint ("%s   mean age, and standard deviation. In terms of the common shape (alpha) and\n", spacer);
                    YvyraPrint ("%s   rate (beta) parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    YvyraPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetLogNormal)
                {
                if (tempD <= 0.0)
                    {
                    YvyraPrint ("%s   Standard deviation must be positive.\n", spacer);
                    YvyraPrint ("%s   Parameters of the offset lognormal distribution are offset age, mean age,\n", spacer);
                    YvyraPrint ("%s   and standard deviation. All values are specified on the linear scale, not\n", spacer);
                    YvyraPrint ("%s   as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            priorParams[2] = tempD;
            expecting = Expecting(RIGHTPAR);
            }
        sprintf (s, "%1.2lf", tempD);
        strcat (calName, s);
        }
    else if (expecting == Expecting(COMMA))
        {
        strcat (calName, ",");
        paramIndex++;
        expecting  = Expecting(NUMBER);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        strcat (calName, ")");
        if (isTaxon == YES)
            YvyraPrint ("%s   Setting age of taxon '%s' to %s\n", spacer, nodeName, calName);
        else
            YvyraPrint ("%s   Setting age of constraint node '%s' to %s\n", spacer, nodeName, calName);

        /* set calibration based on collected values and settings */
        strcpy(calibrationPtr->name, calName);
        calibrationPtr->priorParams[0]  = priorParams[0];
        calibrationPtr->priorParams[1]  = priorParams[1];
        calibrationPtr->priorParams[2]  = priorParams[2];
        calibrationPtr->prior           = calPrior;
        if (calPrior == fixed)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbFix;
            calibrationPtr->LnPriorRatio    = &LnProbRatioFix;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = priorParams[0];
            }
        else if (calPrior == uniform)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbUniform;
            calibrationPtr->LnPriorRatio    = &LnProbRatioUniform;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = priorParams[1];
            }
        else if (calPrior == offsetExponential)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbOffsetExponential_Param_Offset_Mean;
            calibrationPtr->LnPriorRatio    = &LnProbRatioOffsetExponential_Param_Offset_Mean;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == truncatedNormal)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbTruncatedNormal_Param_Trunc_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioTruncatedNormal_Param_Trunc_Mean_Sd;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == logNormal)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbLognormal_Param_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioLognormal_Param_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN;
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == offsetLogNormal)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbOffsetLognormal_Param_Offset_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioOffsetLognormal_Param_Offset_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN + priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == standardGamma)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbGamma_Param_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioGamma_Param_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN;
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == offsetGamma)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbOffsetGamma_Param_Offset_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioOffsetGamma_Param_Offset_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN + priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }

        /* get ready to find more calibrated nodes or taxa, if present */
        expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoCharset (void)
{
    /* first add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* add name to charSetNames */
    if (AddString (&charSetNames, numCharSets, tempSetName) == ERROR)
        {
        YvyraPrint ("%s   Problem adding charset %s to list\n", spacer, tempSetName);
        return (ERROR);
        }

    /* store charset */
    AddBitfield (&charSet, numCharSets, tempSet, numChar);

    /* increment number of char sets */
    numCharSets++;

    return (NO_ERROR);
}


int DoCharsetParm (char *parmName, char *tkn)
{
    int     i, index, tempInt, allDigit;
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before charsets can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check that the name of the charset is not a number */
            allDigit = YES;
            for (i=0; i<(int)strlen(tkn); i++)
                {
                if (tkn[i] == '0' || tkn[i] == '1' || tkn[i] == '2' || tkn[i] == '3' || tkn[i] == '4' || 
                    tkn[i] == '5' || tkn[i] == '6' || tkn[i] == '7' || tkn[i] == '8' || tkn[i] == '9' || tkn[i] == '.')
                    {}
                else
                    allDigit = NO;
                }
            if (allDigit == YES)
                {
                YvyraPrint ("%s   Charset name may not be a number\n", spacer);
                return (ERROR);
                }
            
            /* check size of charset name */
            if (strlen(tkn) > 99)
                {
                YvyraPrint ("%s   Charset name is too long\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a charset */
            if (numCharSets > 1)
                {
                if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                    {
                    /* if the charset name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    YvyraPrint ("%s   Charset name has been used previously\n", spacer);
                    return (ERROR);
                    }
                }
                
            /* add the name to the character set */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numChar; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            YvyraPrint ("%s   Defining charset called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* We are defining a character set in terms of another (called tkn, here). We should be able
           to find tkn in the list of character set names. If we cannot, then we have a problem and
           return an error. */
        if (numCharSets < 1)
            {
            YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
            {
            YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        /* add characters from charset "tkn" to new tempset */
        for (i=0; i<numChar; i++)
            {
            if (IsBitSet(i,charSet[index]) == YES)
                tempSet[i] = 1;
            }       
        fromI = toJ = everyK = -1;

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                    {
                    return (ERROR);
                    }
                }
                
            }

        
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoCharStat (void)
{
    int         i, j, numDivs;
    char        tempName[100];
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }
            
    if (numDefinedPartitions == 1)
        YvyraPrint ("%s   1 character partition defined:\n", spacer, numDefinedPartitions);
    else
        YvyraPrint ("%s   %d character partitions defined:\n", spacer, numDefinedPartitions);
    for (i=0; i<numDefinedPartitions; i++)
        {
        numDivs = GetNumPartDivisions (i);
        if (numDivs == 1)
            YvyraPrint ("%s      Partition %d (\"%s\") does not divide the characters\n", spacer, i+1, partitionNames[i]);
        else
            YvyraPrint ("%s      Partition %d (\"%s\") divides the characters into %d parts\n", spacer, i+1, partitionNames[i], numDivs);
        }
    YvyraPrint ("%s      Current partition is \"%s\"\n", spacer, partitionNames[partitionNum]);
    YvyraPrint ("\n");

    /* print out list of characters with information about each */
    YvyraPrint ("%s   Showing character status:\n\n", spacer);
    YvyraPrint ("%s                                                    Partition(s)\n", spacer);
    YvyraPrint ("%s      #      Type      In/Out    Ambiguity Order  ", spacer);
    for (i=0; i<numDefinedPartitions; i++)
        YvyraPrint (" %2d", i+1);
    YvyraPrint ("\n");
    YvyraPrint ("%s   -----------------------------------------------", spacer);
    for (i=0; i<numDefinedPartitions; i++)
        YvyraPrint ("---");
    YvyraPrint ("\n");
    for (i=0; i<numChar; i++)
        {
        YvyraPrint ("%s   %4d -- ", spacer, i+1);
                
        if (charInfo[i].charType == DNA)
            YvyraPrint ("   DNA");
        else if (charInfo[i].charType == RNA)
            YvyraPrint ("   RNA");
        else if (charInfo[i].charType == PROTEIN)
            YvyraPrint ("  Prot");
        else if (charInfo[i].charType == RESTRICTION)
            YvyraPrint ("  Rest");
        else if (charInfo[i].charType == STANDARD)
            YvyraPrint (" Stand");
        else if (charInfo[i].charType == CONTINUOUS)
            YvyraPrint ("  Cont");
            
        if (charInfo[i].charType == DNA)
            YvyraPrint ("   4");
        else if (charInfo[i].charType == RNA)
            YvyraPrint ("   4");
        else if (charInfo[i].charType == PROTEIN)
            YvyraPrint ("  20");
        else if (charInfo[i].charType == RESTRICTION)
            YvyraPrint ("   2");
        else if (charInfo[i].charType == STANDARD)
            YvyraPrint ("  %2d", charInfo[i].numStates);
        else if (charInfo[i].charType == CONTINUOUS)
            YvyraPrint (" Inf");
            
        if (charInfo[i].isExcluded == NO)
            YvyraPrint ("  Included");
        else
            YvyraPrint ("  Excluded");
            
        if (charInfo[i].isMissAmbig == YES)
            YvyraPrint ("  MissAmbig");
        else
            YvyraPrint ("       None");
            
        if (charInfo[i].ctype == UNORD)
            YvyraPrint (" Unord");
        else if (charInfo[i].ctype == ORD)
            YvyraPrint ("   Ord");
        else if (charInfo[i].ctype == DOLLO)
            YvyraPrint (" Dollo");
        else if (charInfo[i].ctype == IRREV)
            YvyraPrint (" Irrev");

        YvyraPrint ("  ");
            
        for (j=0; j<numDefinedPartitions; j++)
            YvyraPrint (" %2d", partitionId[i][j]);

        /* YvyraPrint ("%4d   ", charSet[i]);*/
        
        if (charInfo[i].pairsId > 0)
            {
            /* find paired character */
            for (j=0; j<numChar; j++)
                {
                if (i != j && charInfo[j].pairsId == charInfo[i].pairsId)
                    {
                    YvyraPrint (" (coupled with %d)", j+1);
                    break;
                    }
                }
            }
                    
        YvyraPrint ("\n");
        
        if (charInfo[i].bigBreakAfter == YES)
            {
            YvyraPrint ("%s   ", spacer);
            YvyraPrint ("     - - - - - - - - - - - - - - - - - - - -  \n");
            }
        
        /* we may want to pause */
        if (autoClose == NO)
            {
            if ((i+1) % 100 == 0)
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

    return (NO_ERROR);
}


int DoCitations (void)
{
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");
    YvyraPrint ("   Citations                                                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra is a fork of MrBayes. If you publish results obtained using yvyra,     \n");
    YvyraPrint ("   please cite MrBayes:                                                          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Ronquist, F. et al. 2012. MrBayes 3.2: Efficient Bayesian phylogenetic    \n");
    YvyraPrint ("         inference and model selection across a large model space.               \n");
    YvyraPrint ("         Syst. Biol. 61:539-542.                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   For the Mk model used for discrete morphological/linguistic characters:       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Lewis, P. O. 2001. A likelihood approach to estimating phylogeny from      \n");
    YvyraPrint ("         discrete morphological character data. Syst. Biol. 50:913-925.          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   If you use the parallel abilities of the program, also cite:                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Altekar, G., S. Dwarkadas, J. P. Huelsenbeck, and F. Ronquist. 2004.       \n");
    YvyraPrint ("         Parallel Metropolis-coupled Markov chain Monte Carlo for Bayesian       \n");
    YvyraPrint ("         phylogenetic inference. Bioinformatics 20:407-415.                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   You should also cite other papers for different ideas that are implemented    \n");
    YvyraPrint ("   in the program.  For example, the program performs Bayesian inference of      \n");
    YvyraPrint ("   phylogeny, an idea that was first proposed in the following papers:           \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Larget, B., and D. Simon. 1999. Markov chain Monte Carlo                   \n");
    YvyraPrint ("         algorithms for the Bayesian analysis of phylogenetic trees.             \n");
    YvyraPrint ("         Mol. Biol. Evol. 16:750-759.                                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Li, S. 1996. Phylogenetic tree construction using Markov chain             \n");
    YvyraPrint ("         Monte carlo. Ph. D. dissertation, Ohio State University, Columbus.      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Mau, B. 1996. Bayesian phylogenetic inference via Markov chain             \n");
    YvyraPrint ("         Monte carlo methods. Ph. D. dissertation, University of                 \n");
    YvyraPrint ("         Wisconsin, Madison.                                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Mau, B., and M. Newton. 1997. Phylogenetic inference for binary            \n");
    YvyraPrint ("         data on dendrograms using Markov chain Monte Carlo. Journal of          \n");
    YvyraPrint ("         Computational and Graphical Statistics 6:122-131.                       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Mau, B., M. Newton, and B. Larget. 1999. Bayesian phylogenetic             \n");
    YvyraPrint ("         inference via Markov chain Monte carlo methods. Biometrics. 55:1-12.    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Newton, M., B. Mau, and B. Larget. 1999. Markov chain Monte Carlo          \n");
    YvyraPrint ("         for the Bayesian analysis of evolutionary trees from aligned            \n");
    YvyraPrint ("         molecular sequences. In Statistics in molecular biology (F. Seillier-   \n");
    YvyraPrint ("         Moseiwitch, T. P. Speed, and M. Waterman, eds.). Monograph Series       \n");
    YvyraPrint ("         of the Institute of Mathematical Statistics.                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Rannala, B., and Z. Yang. 1996. Probability distribution of                \n");
    YvyraPrint ("         molecular evolutionary trees: a new method of phylogenetic              \n");
    YvyraPrint ("         inference. J. Mol. Evol. 43:304-311.                                    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z., and B. Rannala. 1997. Bayesian phylogenetic inference            \n");
    YvyraPrint ("         using DNA sequences: a Markov chain Monte Carlo method. Molecular       \n");
    YvyraPrint ("         Biology and Evolution. 14:717-724.                                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra uses Markov chain Monte Carlo (MCMC) to approximate the posterior     \n");
    YvyraPrint ("   probability of trees.  MCMC was developed in the following papers:            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Metropolis, N., A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller,          \n");
    YvyraPrint ("         and E. Teller. 1953. Equations of state calculations by fast            \n");
    YvyraPrint ("         computing machines. J. Chem. Phys. 21:1087-1091.                        \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Hastings, W. K. 1970. Monte Carlo sampling methods using Markov            \n");
    YvyraPrint ("         chains and their applications. Biometrika 57:97-109.                    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   In particular, yvyra implements a variant of MCMC that was described by     \n");
    YvyraPrint ("   Charles Geyer:                                                                \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Geyer, C. J. 1991. Markov chain Monte Carlo maximum likelihood.            \n");
    YvyraPrint ("         Pages 156-163 in Computing Science and Statistics: Proceed-             \n");
    YvyraPrint ("         ings of the 23rd Symposium on the Interface. (E. M. Keramidas,          \n");
    YvyraPrint ("         ed.). Fairfax Station: Interface Foundation.                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra implements a large number of DNA substitution models. These models    \n");
    YvyraPrint ("   are of three different structures.  The \"4by4\" models are the usual flavor  \n");
    YvyraPrint ("   of phylogenetic models.  The \"Doublet\" model was first proposed by          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Schoniger, M., and A. von Haeseler. 1994. A stochastic model and the       \n");
    YvyraPrint ("         evolution of autocorrelated DNA sequences. Molecular Phylogenetics      \n");
    YvyraPrint ("         and Evolution 3:240-247.                                                \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The program also implements codon models.  Two papers, published back-to-back \n");
    YvyraPrint ("   were the first to implement a codon model of DNA substitution in which the    \n");
    YvyraPrint ("   substitution process is modelled on the codon, not on a site-by-site basis:   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Goldman, N., and Z. Yang. 1994. A codon-based model of nucleotide          \n");
    YvyraPrint ("         substitution for protein coding DNA sequences. Molecular Biology        \n");
    YvyraPrint ("         and Evolution. 11:725-736.                                              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Muse, S., and B. Gaut. 1994. A likelihood approach for comparing           \n");
    YvyraPrint ("         synonymous and non-synonymous substitution rates, with application      \n");
    YvyraPrint ("         to the chloroplast genome. Molecular Biology and Evolution.             \n");
    YvyraPrint ("         11:715-724.                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The program can be used to detect positively selected amino-acid sites using  \n");
    YvyraPrint ("   a full hierarchical Bayes analysis.  The method is based on the excellent     \n");
    YvyraPrint ("   paper by Nielsen and Yang:                                                    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Nielsen, R., and Z. Yang. 1998. Likelihood models for detecting            \n");
    YvyraPrint ("         positively selected amino acid sites and applications to the HIV-1      \n");
    YvyraPrint ("         envelope gene. Genetics. 148:929-936.                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The previous four papers describe three different structures for the          \n");
    YvyraPrint ("   nucleotide models implemented in yvyra--the four-by-four models, the        \n");
    YvyraPrint ("   16-by-16 (doublet) models and the 64-by-64 (codon) models.  The program       \n");
    YvyraPrint ("   implements three different substitution models within each model structure.   \n");
    YvyraPrint ("   These include the nst=1 models:                                               \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Jukes, T., and C. Cantor. 1969. Evolution of protein molecules.            \n");
    YvyraPrint ("         Pages 21-132 in Mammalian Protein Metabolism. (H. Munro, ed.).          \n");
    YvyraPrint ("         Academic Press, New York.                                               \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Felsenstein, J. 1981. Evolutionary trees from DNA sequences: A             \n");
    YvyraPrint ("         maximum likelihood approach. Journal of Molecular Evolution             \n");
    YvyraPrint ("         17:368-376.                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   the nst=2 models:                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Kimura, M. 1980. A simple method for estimating evolutionary rates         \n");
    YvyraPrint ("         of base substitutions through comparative studies of nucleotide         \n");
    YvyraPrint ("         sequences. Journal of Molecular Evolution. 16:111-120.                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Hasegawa, M., T. Yano, and H. Kishino. 1984. A new molecular clock         \n");
    YvyraPrint ("         of mitochondrial DNA and the evolution of Hominoids. Proc.              \n");
    YvyraPrint ("         Japan Acad. Ser. B 60:95-98.                                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Hasegawa, M., H. Kishino, and T. Yano. 1985. Dating the human-ape          \n");
    YvyraPrint ("         split by a molecular clock of mitochondrial DNA. Journal of             \n");
    YvyraPrint ("         Molecular Evolution 22:160-174.                                         \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   and the the nst=6 models:                                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Tavare, S. 1986. Some probabilistic and statistical problems on the        \n");
    YvyraPrint ("         analysis of DNA sequences. Lect. Math. Life Sci. 17:57-86.              \n");
    YvyraPrint ("         17:368-376.                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra implements a large number of amino-acid models.  These include:       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Poisson --                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Bishop, M.J., and A.E. Friday. 1987. Tetrapod relationships: the           \n");
    YvyraPrint ("         molecular evidence. Pp. 123-139 in Molecules and morphology in          \n");
    YvyraPrint ("         evolution: conflict or compromise? (C. Patterson, ed.). Cambridge       \n");
    YvyraPrint ("         University Press, Cambridge, England.                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Jones --                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Jones, D.T., W. R. Taylor, and J. M. Thornton. 1992. The rapid generation  \n");
    YvyraPrint ("         of mutation data matrices from protein sequences. Comput. Appl.         \n");
    YvyraPrint ("         Biosci. 8:275-282.                                                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Dayhoff --                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Dayhoff, M.O., R.M. Schwartz, and B.C. Orcutt. 1978. A model of evol-      \n");
    YvyraPrint ("         utionary change in proteins. Pp. 345-352 in Atlas of protein sequence   \n");
    YvyraPrint ("         and structure. Vol. 5, Suppl. 3. National Biomedical Research           \n");
    YvyraPrint ("          Foundation, Washington, D.C.                                           \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Mtrev --                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Adachi, J. and M. Hasegawa. 1996. MOLPHY version 2.3: programs for         \n");
    YvyraPrint ("         molecular phylogenetics based on maximum likelihood.  Computer Science  \n");
    YvyraPrint ("         Monographs of Institute of Statistical Mathematics 28:1-150.            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Mtmam --                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Cao, Y., A. Janke, P.J. Waddell, M. Westerman, O. Takenaka, S. Murata,     \n");
    YvyraPrint ("         N. Okada, S. Paabo, and M. Hasegawa. 1998. Conflict amongst individual  \n");
    YvyraPrint ("         mitochondrial proteins in resolving the phylogeny of eutherian orders.  \n");
    YvyraPrint ("         Journal of Molecular Evolution                                          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z., R. Nielsen, and M. Hasegawa. 1998.  Models of amino acid         \n");
    YvyraPrint ("         substitution and applications to mitochondrial protein evolution        \n");
    YvyraPrint ("         Molecular Biology and Evolution 15:1600-1611.                           \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      WAG --                                                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Whelan, S. and Goldman, N. 2001. A general empirical model of protein      \n");
    YvyraPrint ("         evolution derived from multiple protein families using a maximum-       \n");
    YvyraPrint ("         likelihood approach. Molecular Biology and Evolution 18:691-699.        \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Rtrev --                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Dimmic M.W., J.S. Rest, D.P. Mindell, and D. Goldstein. 2002. RArtREV:     \n");
    YvyraPrint ("         An amino acid substitution matrix for inference of retrovirus and       \n");
    YvyraPrint ("         reverse transcriptase phylogeny. Journal of Molecular Evolution         \n");
    YvyraPrint ("         55: 65-73.                                                              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Cprev --                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Adachi, J., P. Waddell, W. Martin, and M. Hasegawa. 2000. Plastid          \n");
    YvyraPrint ("         genome phylogeny and a model of amino acid substitution for proteins    \n");
    YvyraPrint ("         encoded by chloroplast DNA. Journal of Molecular Evolution              \n");
    YvyraPrint ("         50:348-358.                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Blosum --                                                                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Henikoff, S., and J. G. Henikoff. 1992. Amino acid substitution            \n");
    YvyraPrint ("         matrices from protein blocks. Proc. Natl. Acad. Sci., U.S.A.            \n");
    YvyraPrint ("         89:10915-10919. The matrix implemented in yvyra is Blosum62.          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Vt --                                                                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Muller, T., and M. Vingron. 2000. Modeling amino acid replacement.         \n");
    YvyraPrint ("         Journal of Computational Biology 7:761-776.                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      LG --                                                                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Le, Si Q. & Gascuel, O. 2008 An improved general amino- acid replacement   \n");
    YvyraPrint ("         matrix. Mol. Biol. Evol. 25, 1307-1320.                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra implements a simple Jukes-Cantor-like model for restriction sites     \n");
    YvyraPrint ("   and other binary data.  A problem with some of these data is that there is a  \n");
    YvyraPrint ("   coding bias, such that certain characters are missing from any observable data\n");
    YvyraPrint ("   matrix.  It is impossible, for instance, to observe restriction sites that are\n");
    YvyraPrint ("   absent in all the studied taxa.  However, yvyra corrects for this coding    \n");
    YvyraPrint ("   bias according to an idea described in                                        \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Felsenstein, J. 1992. Phylogenies from restriction sites: A maximum-       \n");
    YvyraPrint ("         likelihood approach. Evolution 46:159-173.                              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The model used by yvyra for 'standard' or morphological data is based on    \n");
    YvyraPrint ("   the ideas originally presented by                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Lewis, P. O. 2001. A likelihood approach to estimating phylogeny from      \n");
    YvyraPrint ("         discrete morphological character data. Systematic Biology 50:913-925.   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   For both DNA sequence and amino-acid data, the program allows rates to change \n");
    YvyraPrint ("   under a covarion-like model, first described by Tuffley and Steel             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Tuffley, C., and M. Steel. 1998. Modeling the covarion hypothesis          \n");
    YvyraPrint ("         of nucleotide substitution. Mathematical Biosciences 147:63-91.         \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   and implemented by Huelsenbeck (2002)                                         \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Huelsenbeck, J. P. 2002. Testing a covariotide model of DNA sub-           \n");
    YvyraPrint ("         stitution. Molecular Biology and Evolution 19(5):698-707.               \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   Galtier (2001) implements a different variant of the covarion model in a      \n");
    YvyraPrint ("   paper that is worth reading:                                                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Galtier, N. 2001. Maximum-likelihood phylogenetic analysis under a         \n");
    YvyraPrint ("         covarion-like model. Mol. Biol. Evol. 18:866-873.                       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   A number of models are available that allow rates to vary across the          \n");
    YvyraPrint ("   characters.  The program implements the proportion of invariable sites model  \n");
    YvyraPrint ("   and two variants of gamma distributed rate variation.  Yang\'s (1993) paper   \n");
    YvyraPrint ("   is a good one to cite for implementing a gamma-distributed rates model.       \n");
    YvyraPrint ("   In the 1994 paper he provides a way to approximate the continuous gamma       \n");
    YvyraPrint ("   distribution:                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z. 1993. Maximum likelihood estimation of phylogeny from DNA         \n");
    YvyraPrint ("         sequences when substitution rates differ over sites. Molecular          \n");
    YvyraPrint ("         Biology and Evolution 10:1396-1401.                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z. 1994. Maximum likelihood phylogenetic estimation from DNA         \n");
    YvyraPrint ("         sequences with variable rates over sites: Approximate methods.          \n");
    YvyraPrint ("         Journal of Molecular Evolution 39:306-314.                              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The program also implements Yang\'s autocorrelated gamma model.  In this      \n");
    YvyraPrint ("   model, the rate at one site depends to some extent on the rate at an adjacent \n");
    YvyraPrint ("   site.  The appropriate citation for this model is:                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z. 1995. A space-time process model for the evolution of             \n");
    YvyraPrint ("         DNA sequences. Genetics 139:993-1005.                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The following two papers show how ancestral states on a tree can be           \n");
    YvyraPrint ("   reconstructed.  The Yang et al. paper implements an empirical Bayes approach  \n");
    YvyraPrint ("   while Huelsenbeck and Bollback use a pure, hierarchical Bayes approach. The   \n");
    YvyraPrint ("   method used in yvyra is the latter, since it integrates over uncertainty in \n");
    YvyraPrint ("   model parameters.                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z., S. Kumar, and M. Nei. 1995. A new method of inference of         \n");
    YvyraPrint ("         ancestral nucleotide and amino acid sequences. Genetics 141:1641        \n");
    YvyraPrint ("         1650.                                                                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Huelsenbeck, J. P., and J. P. Bollback. 2001. Empirical and hier-          \n");
    YvyraPrint ("         archical Bayesian estimation of ancestral states. Systematic            \n");
    YvyraPrint ("         Biology 50:351-366.                                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   You may also want to consult a more recent review of Bayesian reconstruction  \n");
    YvyraPrint ("   of ancestral states and character evolution:                                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Ronquist, F. 2004. Bayesian inference of character evolution. Trends in    \n");
    YvyraPrint ("         Ecology and Evolution 19: 475-481.                                      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra allows you to analyze gene tree - species tree problems using the     \n");
    YvyraPrint ("   multi-species coalescent approach originally proposed by Edwards et al:       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Edwards, S., L. Liu, and D. Pearl. 2007. High-resolution species trees     \n");
    YvyraPrint ("         without concatenation. Proc. Natl. Acad. Sci. USA 104: 5936-5941.       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The program implements an incredibly parameter rich model, first described by \n");
    YvyraPrint ("   Tuffley and Steel (1997), that orders trees in the same way as the so-called  \n");
    YvyraPrint ("   parsimony method of phylogenetic inference.  The appropriate citation is:     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Tuffley, C., and M. Steel. 1997. Links between maximum likelihood          \n");
    YvyraPrint ("         and maximum parsimony under a simple model of site substitution.        \n");
    YvyraPrint ("         Bull. Math. Bio. 59:581-607.                                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   yvyra implements five relaxed clock models: the Compound Poisson Process    \n");
    YvyraPrint ("   (CPP), the Autocorrelated Lognormal (TK02), the White Noise (WN), the         \n");
    YvyraPrint ("   Independent Gamma Rates (IGR), and the Independent Lognormal (ILN) models.    \n");
    YvyraPrint ("   The CPP model was first described by Huelsenbeck et al. (2000).  It is an     \n");
    YvyraPrint ("   autocorrelated discrete model of rate variation over time.  Instead of        \n");
    YvyraPrint ("   the modified gamma distribution originally proposed for the rate multipliers, \n");
    YvyraPrint ("   yvyra uses a lognormal distribution.  The extensions necessary to sample    \n");
    YvyraPrint ("   over tree space under this model are original to yvyra; the original paper  \n");
    YvyraPrint ("   only considered fixed trees.                                                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The TK02 model was first described by Thorne and Kishino (2002), and is a     \n");
    YvyraPrint ("   variant of a model presented by them earlier (Thorne et al., 1998).  It is an \n");
    YvyraPrint ("   autocorrelated continuous model, in which rates vary according to a lognormal \n");
    YvyraPrint ("   distribution.  Specifically, the rate of a descendant node is assumed to      \n");
    YvyraPrint ("   be drawn from a lognormal distribution with the mean being the rate of the    \n");
    YvyraPrint ("   ancestral node (and not the log of the rate of the ancestral node), and the   \n");
    YvyraPrint ("   variance (on the log scale) being proportional to the length of the branch    \n");
    YvyraPrint ("   separating the nodes (measured in terms of expected substitutions per site at \n");
    YvyraPrint ("   the base rate of the clock). Note that in MrBayes up to version 3.2.7,      \n");
    YvyraPrint ("   the variance was measured on the linear scale and not on the log scale, which \n");
    YvyraPrint ("   a slight difference from the original formulation of the model. Note also that\n");
    YvyraPrint ("   we factor out the base rate of the clock in yvyra, so the model applies to  \n");
    YvyraPrint ("   the deviations from this rate.                                                \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   In the WN model, the branch rates are modeled as being drawn independently    \n");
    YvyraPrint ("   from gamma distributions. The distributions are not identical as the variance \n");
    YvyraPrint ("   is proportional to the branch length. See Lepage et al. (2007) for details.   \n");
    YvyraPrint ("   Note that the WN model was named 'IGR' in previous versions of MrBayes (up  \n");
    YvyraPrint ("   version 3.2.7), but now 'IGR' refers to a slightly different model, in which  \n");
    YvyraPrint ("   the branch rates are drawn from independent and identically distributed       \n");
    YvyraPrint ("   (i.i.d.) gamma distributions. This naming change better reflects common usage \n");
    YvyraPrint ("   of the term 'independent' in the relaxed clock model literature.              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The ILN model is analogous to IGR but differs in that the branch rates follow \n");
    YvyraPrint ("   i.i.d. lognormal (instead of gamma) distributions (Drummond et al. 2006).     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Huelsenbeck, J. P., B. Larget, and D. Swofford. 2000. A compound Poisson   \n");
    YvyraPrint ("         process for relaxing the molecular clock. Genetics 154: 1879-1892.      \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Thorne, J. L., H. Kishino, and I. S. Painter. 1998. Estimating the rate    \n");
    YvyraPrint ("         of evolution of the rate of molecular evolution. Mol. Biol. Evol.       \n");
    YvyraPrint ("         15: 1647-1657.                                                          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Thorne, J. L., and H. Kishino. 2002. Divergence time and evolutionary      \n");
    YvyraPrint ("         rate estimation with multilocus data. Syst. Biol. 51: 689-702.          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Drummond, A. J., S. Y. W. Ho, M. J. Phillips, and A. Rambaut. 2006.        \n");
    YvyraPrint ("         Relaxed phylogenetics and dating with confidence. PLoS Biology          \n");
    YvyraPrint ("         4: 699-710.                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Lepage, T., D. Bryant, H. Philippe, and N. Lartillot. 2007. A general      \n");
    YvyraPrint ("         comparison of relaxed molecular clock models. Mol. Biol. Evol.          \n");
    YvyraPrint ("         24: 2669-2680.                                                          \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The standard tree proposals used by yvyra are described by Lakner et al.    \n");
    YvyraPrint ("   (2008).  The parsimony-biased tree proposals are still undescribed, although  \n");
    YvyraPrint ("   a rough outline of the idea is presented in the same paper.                   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Lakner, C., P. van der Mark, J. P. Huelsenbeck, B. Larget, and F. Ronquist.\n");
    YvyraPrint ("         2008. Efficiency of Markov chain Monte Carlo tree proposals in Bayesian \n");
    YvyraPrint ("         phylogenetics. Syst. Biol. 57: 86-103.                                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The topology convergence diagnostic used by yvyra, the average standard     \n");
    YvyraPrint ("   deviation of split frequencies, is described by Lakner et al. (2008).         \n");
    YvyraPrint ("   The potential scale reduction factor, the diagnostic used by yvyra for      \n");
    YvyraPrint ("   continous parameters, was first proposed by Gelman and Rubin (1992).  The     \n");
    YvyraPrint ("   auto-tuning mechanism used in yvyra is based on a paper by Roberts and      \n");
    YvyraPrint ("   Rosenthal (2009).                                                             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Gelman, A., and D. B. Rubin. 1992. Inference from iterative simulation     \n");
    YvyraPrint ("         using multiple sequences. Statistical Science 7: 457-472.               \n");
    YvyraPrint ("         Bull. Math. Bio. 59:581-607.                                            \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Lakner, C., P. van der Mark, J. P. Huelsenbeck, B. Larget, and F. Ronquist.\n");
    YvyraPrint ("         2008. Efficiency of Markov chain Monte Carlo tree proposals in Bayesian \n");
    YvyraPrint ("         phylogenetics. Syst. Biol. 57: 86-103.                                  \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Roberts, G. O., and J. S. Rosenthal. 2009. Examples of adaptive            \n");
    YvyraPrint ("         MCMC. Journal of Compuational and Graphical Statistics 18: 349-367.     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   The harmonic mean estimator of model likelihoods, used for Bayes factor       \n");
    YvyraPrint ("   testing, was discussed by Newton and Raftery (1996).  The more accurate       \n");
    YvyraPrint ("   stepping-stone algorithm was first proposed by Xie et al. (2011).  The paper  \n");
    YvyraPrint ("   by Lartillot and Philippe (2006) presents an interesting discussion of the    \n");
    YvyraPrint ("   shortcomings of the harmonic mean estimator and describes thermodynamic       \n");
    YvyraPrint ("   integration, a technique that is similar to the stepping-stone algorithm.     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Newton, M. A., and A. E. Raftery. 1994. Approximate Bayesian inference     \n");
    YvyraPrint ("         with the weighted likelihood bootstrap. J. R. Stat. Soc. B. 56. 3-48.   \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Lartillot, N., and H. Philippe. 2006. Computing Bayes factors using        \n");
    YvyraPrint ("         thermodynamic integration. Syst. Biol. 55: 195-207.                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Xie, W., P. O. Lewis, Y. Fan, L. Kuo, and M.-H. Chen. 2011. Improving      \n");
    YvyraPrint ("         marginal likelihood estimation for Bayesian phylogenetic model          \n");
    YvyraPrint ("         selection. Syst. Biol. 60: 150-160.                                     \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   For unconstrained branch lengths, yvyra implements the compound Dirichlet   \n");
    YvyraPrint ("   priors for branch lengths described by Rannala et al. (2012) and Zhang et al. \n");
    YvyraPrint ("   (2012).  Compared with the i.i.d. exponential and uniform priors for branch   \n");
    YvyraPrint ("   lengths in the previous versions of MrBayes, the Dirichlet priors appear more \n");
    YvyraPrint ("   reasonable and may avoid the problem of extremely long trees, as discussed    \n");
    YvyraPrint ("   by Brown et al. (2010) and Marshall (2010).  The two-exponential prior on     \n");
    YvyraPrint ("   internal and external branch lengths described by Yang & Rannala (2005) and   \n");
    YvyraPrint ("   Yang (2007) is also implemented in this version.                              \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Brown, J. M., S. M. Hedtke, A. R. Lemmon, and E. M. Lemmon. 2010. When     \n");
    YvyraPrint ("         trees  grow too long: investigating the causes of highly inaccurate     \n");
    YvyraPrint ("         Bayesian branch-length estimates. Syst. Biol. 59:145-161.               \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Marshall, D. C. 2010. Cryptic failure of partitioned Bayesian phylogenetic \n");
    YvyraPrint ("         analyses: lost in the land of long trees. Syst. Biol. 59:108-117.       \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Rannala, B., T. Zhu, and Z. Yang. 2012. Tail paradox, partial              \n");
    YvyraPrint ("         identifiability and influential priors in Bayesian branch length        \n");
    YvyraPrint ("         inference. Mol. Biol. Evol. 29:325-335.                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Zhang, C., B. Rannala, and Z. Yang. 2012. Robustness of compound Dirichlet \n");
    YvyraPrint ("         priors for Bayesian inference of branch lengths. Syst. Biol. 61:779-784.\n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z. 2007. Fair-balance paradox, star-tree paradox and Bayesian        \n");
    YvyraPrint ("         phylogenetics. Mol. Biol. Evol. 24:1639-1655.                           \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      Yang, Z., and B. Rannala. 2005. Branch-length prior influences Bayesian    \n");
    YvyraPrint ("         posterior probability of phylogeny. Syst. Biol. 54:455-470.             \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoConstraint (void)
{
    int         i, howMany;
    int         *tset;

    if (constraintType == PARTIAL)
        tset=tempSetNeg;
    else
        tset=tempSet;

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToGivenSet (fromI, toJ, everyK, 1, tset) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToGivenSet (fromI, toJ, everyK, 1, tset) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToGivenSet (fromI, toJ, everyK, 1, tset) == ERROR)
            return (ERROR);
        }
            
    /* check that this is not a stupid constraint */
    howMany = 0;
    for (i=0; i<numTaxa; i++)
        if (tempSet[i] != 0)
            howMany++;

    if (howMany == 0)
        {
        YvyraPrint ("%s   This constraint does not include any taxa and will not be defined\n", spacer);
        return (ERROR);
        }

    if (constraintType == HARD)
        {
        if (howMany == numTaxa)
            {
            /* We allow this so we can report states from and calibrate root */
            }
        
        } /*end constraintType == HARD */
    else if (constraintType == PARTIAL)
        {
        if (howMany == 1)
            {
            YvyraPrint ("%s   This partial constraint includes only one taxon. It is always satisfied and will not be defined.\n", spacer);
            return (ERROR);
            }

        howMany = 0;
        for (i=0; i<numTaxa; i++)
            {
            if (tempSetNeg[i] != 0)
                {
                howMany++;
                if (tempSetNeg[i] == tempSet[i])
                    {
                    YvyraPrint ("%s   Two sets of taxa in partial constraint are not allowed to intersect. Constraint will not be defined\n", spacer);
                    return (ERROR);
                    }
                }
            }
        if (howMany == 0)
            {
            YvyraPrint ("%s   This partial constraint does not include any taxa in the second set and will not be defined\n", spacer);
            return (ERROR);
            }
        }
    else if (constraintType == NEGATIVE)
        {
        if (howMany == 1)
            {
            YvyraPrint ("%s   Negative constraint should include more than one taxon. Constraint will not be defined\n", spacer);
            return (ERROR);
            }
        }

    /* add name to constraintNames */
    if (AddString (&constraintNames, numDefinedConstraints, tempSetName) == ERROR)
        {
        YvyraPrint ("%s   Problem adding constraint %s to list\n", spacer, tempSetName);
        return (ERROR);
        }

    /* store tempSet */
    AddBitfield (&definedConstraint, numDefinedConstraints, tempSet, numTaxa);
    if (constraintType == PARTIAL)
        {
        AddBitfield (&definedConstraintTwo, numDefinedConstraints, tempSetNeg, numTaxa);
        }
    else
        {
        definedConstraintTwo = (BitsLong **) SafeRealloc ((void *)(definedConstraintTwo), ((size_t)numDefinedConstraints+1)*sizeof(BitsLong *));
        if (definedConstraintTwo==NULL)
            return ERROR;
        definedConstraintTwo[numDefinedConstraints]=NULL;
        }
    
    /* add a default node calibration */
    nodeCalibration = (Calibration *) SafeRealloc ((void *)nodeCalibration, ((size_t)numDefinedConstraints+1)*sizeof(Calibration));
    nodeCalibration[numDefinedConstraints].prior            = defaultCalibration.prior;
    nodeCalibration[numDefinedConstraints].priorParams[0]   = defaultCalibration.priorParams[0];
    nodeCalibration[numDefinedConstraints].priorParams[1]   = defaultCalibration.priorParams[1];
    nodeCalibration[numDefinedConstraints].priorParams[2]   = defaultCalibration.priorParams[2];
    nodeCalibration[numDefinedConstraints].min              = defaultCalibration.min;
    nodeCalibration[numDefinedConstraints].max              = defaultCalibration.max;
    nodeCalibration[numDefinedConstraints].LnPriorProb      = defaultCalibration.LnPriorProb;
    nodeCalibration[numDefinedConstraints].LnPriorRatio     = defaultCalibration.LnPriorRatio;
    strcpy(nodeCalibration[numDefinedConstraints].name, defaultCalibration.name);

    /* increment number of defined constraints */
    numDefinedConstraints++;

    /* reallocate and initialize space for activeConstraints */
    for (i=0; i<numCurrentDivisions; i++)
        {
        modelParams[i].activeConstraints = (int *) SafeRealloc((void *)(modelParams[i].activeConstraints), (size_t)numDefinedConstraints*sizeof(int));
        modelParams[i].activeConstraints[numDefinedConstraints-1] = NO;
        }

    /* reallocate and initialize space for tempActiveConstraints */
    tempActiveConstraints = (int *) SafeRealloc((void *)(tempActiveConstraints), (size_t)numDefinedConstraints*sizeof(int));
    tempActiveConstraints[numDefinedConstraints-1] = NO;

    definedConstraintsType = (enum ConstraintType *) SafeRealloc((void *)(definedConstraintsType), (size_t)numDefinedConstraints*sizeof(enum ConstraintType));
    if (definedConstraintsType==NULL)
        return ERROR;
    definedConstraintsType[numDefinedConstraints-1] = constraintType;

    definedConstraintPruned = (BitsLong **) SafeRealloc ((void *)(definedConstraintPruned), (size_t)numDefinedConstraints*sizeof(BitsLong *));
    if (definedConstraintPruned==NULL)
        return ERROR;
    definedConstraintPruned[numDefinedConstraints-1]=NULL;

    definedConstraintTwoPruned = (BitsLong **) SafeRealloc ((void *)(definedConstraintTwoPruned), (size_t)numDefinedConstraints*sizeof(BitsLong *));
    if (definedConstraintTwoPruned==NULL)
        return ERROR;
    definedConstraintTwoPruned[numDefinedConstraints-1]=NULL;

    /* show taxset (for debugging) */
    // for (i=0; i<numTaxa; i++)
    //     YvyraPrint ("%4d  %4d\n", i+1, taxaInfo[i].constraints[numDefinedConstraints-1]);

    return (NO_ERROR);
}


int DoConstraintParm (char *parmName, char *tkn)
{
    int         i, index, tempInt;
    YFlt      tempD;
    static int  *tempSetCurrent;
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before constraints can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of constraint name */
            if (strlen(tkn) > 99)
                {
                YvyraPrint ("%s   Constraint name is too long\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a constraint */
            if (numDefinedConstraints > 0)
                {
                if (CheckString (constraintNames, numDefinedConstraints, tkn, &index) == ERROR)
                    {
                    /* an ERROR returned if the constraint name has not been used. we _want_ to be here */
                    }
                else
                    {
                    YvyraPrint ("%s   Constraint name '%s' has been used previously\n", spacer, tkn);
                    return (ERROR);
                    }
                }
                
            /* copy the name to the temporary constraint names string */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 0;

            constraintType = HARD; /* set default constrain type */
            tempSetCurrent=tempSet;
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            YvyraPrint ("%s   Defining constraint called '%s'\n", spacer, tkn);
            foundExp = NO;
            foundFirst = YES;
            foundEqual = NO;
            isNegative = NO;
            foundColon = NO;
            expecting = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(DASH);
            expecting |= Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        foundEqual = YES;
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        isNegative = NO;
        expecting = Expecting(NUMBER);
        expecting |= Expecting(DASH);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        isNegative = NO;
        foundExp = NO;
        expecting = Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(DASH))
        {
        if (foundExp == YES)
            isNegative = YES;
        else
            foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (foundFirst == YES && foundEqual == NO)
            {
            /* We are filling in the probability for the constraint. Specifically, we expect exp(number). */
            if (IsSame ("Partial", tkn) == SAME)
                {
                for (i=0; i<numTaxa; i++)
                    tempSetNeg[i] = 0;

                constraintType = PARTIAL;
                expecting = Expecting(EQUALSIGN);
                expecting |= Expecting(ALPHA);
                }
            else if (IsSame ("Hard", tkn) == SAME)
                {
                constraintType = HARD;
                expecting = Expecting(EQUALSIGN);
                expecting |= Expecting(ALPHA);
                }
            else if (IsSame ("Negative", tkn) == SAME)
                {
                constraintType = NEGATIVE;
                expecting = Expecting(EQUALSIGN);
                expecting |= Expecting(ALPHA);
                }
            else if (IsSame ("Exp", tkn) == SAME || IsSame ("Exp", tkn) == CONSISTENT_WITH)
                {
                foundExp = YES;
                foundDash = NO;
                isNegative = NO;
                expecting  = Expecting(LEFTPAR);
                }
            else
                {
                YvyraPrint ("%s   Do not understand %s\n", spacer, tkn);
                return (ERROR);
                }
            }
        else
            {
            /* We are defining a constraint in terms of a taxon set (called tkn, here) or we are referring to
               the taxon name. We should be able to find tkn in the list of taxon set names or in the list
               of taxon names. If we cannot, then we have a problem and return an error. */
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                if (numTaxaSets < 1)
                    {
                    YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                /* add taxa from taxset tkn to new tempSet */
                for (i=0; i<numTaxa; i++)
                    {
                    if (IsBitSet(i, taxaSet[index]) == YES)
                        {
                        tempSetCurrent[i] = 1;
                        }
                    }
                }
            else
                {
                tempSetCurrent[index] = 1;
                }
            fromI = toJ = everyK = -1;

            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            if (constraintType != PARTIAL || foundColon == YES)
                expecting |= Expecting(SEMICOLON);
            else
                expecting |= Expecting(COLON);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundFirst == YES && foundEqual == NO)
            {
            /* We are filling in the probability for the constraint. Specifically, we expect number. */
            sscanf (tkn, "%lf", &tempD);        
            if (foundExp == NO && tempD < 0.0)
                {
                YvyraPrint ("%s   The probability of a clade cannot be less than zero\n", spacer, tkn);
                return (ERROR);
                }
            if (isNegative == YES || foundDash == YES)
                tempD *= -1.0;
            if (foundExp == YES)
                {
                expecting  = Expecting(RIGHTPAR);
                }
            else
                {
                expecting  = Expecting(EQUALSIGN);
                }
            foundFirst = NO;
            foundDash = NO;
            }
        else
            {       
            if (strlen(tkn) == 1 && !strcmp(tkn, "."))
                {
                tempInt = numTaxa;
                }
            else
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt <= 0 || tempInt > numTaxa)
                    {
                    YvyraPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                    return (ERROR);
                    }
                }
            tempInt--;
            if (foundDash == YES)
                {
                if (fromI >= 0)
                    toJ = tempInt;
                else
                    {
                    YvyraPrint ("%s   Improperly formatted constraint\n", spacer);
                    return (ERROR);
                    }
                foundDash = NO;
                }
            else if (foundSlash == YES)
                {
                tempInt++;
                if (tempInt <= 1)
                    {
                    YvyraPrint ("%s   Improperly formatted constraint\n", spacer);
                    return (ERROR);
                    }
                if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                    everyK = tempInt;
                else
                    {
                    YvyraPrint ("%s   Improperly formatted constraint\n", spacer);
                    return (ERROR);
                    }
                foundSlash = NO;
                }
            else
                {
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToGivenSet (fromI, toJ, everyK, 1, tempSetCurrent) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    }
                else if (fromI < 0 && toJ < 0)
                    {
                    fromI = tempInt;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToGivenSet (fromI, toJ, everyK, 1, tempSetCurrent) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToGivenSet (fromI, toJ, everyK, 1, tempSetCurrent) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else
                    {
                    YvyraPrint ("%s   Improperly formatted constraint\n", spacer);
                        {
                        return (ERROR);
                        }
                    }
                }

            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(DASH);
            expecting |= Expecting(BACKSLASH);
            if (constraintType != PARTIAL || foundColon == YES)
                expecting |= Expecting(SEMICOLON);
            else
                expecting |= Expecting(COLON);
            }
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COLON))
        {
        if (foundColon == YES)
            {
            YvyraPrint ("%s   Improperly formatted constraint: two colon characters in constraint command.\n", spacer);
            return (ERROR);
            }

        /* add set to tempSet */
        if (fromI >= 0 && toJ < 0)
            {
            if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                return (ERROR);
            }
        else if (fromI >= 0 && toJ >= 0 && everyK < 0)
            {
            if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                return (ERROR);
            }
        else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
            {
            if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                return (ERROR);
            }
        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;

        foundColon = YES;
        tempSetCurrent = tempSetNeg;
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoCtype (void)
{
    int         i, foundIllegal, marks[5], numAppliedTo;

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with ctype */
    numAppliedTo = 0;
    for (i=0; i<5; i++)
        marks[i] = NO;
    for (i=0; i<numChar; i++)
        {
        if (tempSet[i] != 0)
            {
            foundIllegal = NO;
            if (charOrdering != UNORD)
                {
                if (charInfo[i].charType == DNA)
                    {
                    foundIllegal = YES;
                    if (marks[0] == NO)
                        YvyraPrint ("%s   Ctype not applied to DNA states which must be unordered\n", spacer);
                    marks[0] = YES;
                    }
                else if (charInfo[i].charType == RNA)
                    {
                    foundIllegal = YES;
                    if (marks[1] == NO)
                        YvyraPrint ("%s   Ctype not applied to RNA states which must be unordered\n", spacer);
                    marks[1] = YES;
                    }
                else if (charInfo[i].charType == PROTEIN)
                    {
                    foundIllegal = YES;
                    if (marks[2] == NO)
                        YvyraPrint ("%s   Ctype not applied to amino acid states which must be unordered\n", spacer);
                    marks[2] = YES;
                    }
                else if (charInfo[i].charType == RESTRICTION)
                    {
                    foundIllegal = YES;
                    if (marks[3] == NO)
                        YvyraPrint ("%s   Ctype not applied to restriction site states which must be unordered\n", spacer);
                    marks[3] = YES;
                    }
                else if (charInfo[i].charType == CONTINUOUS)
                    {
                    foundIllegal = YES;
                    if (marks[4] == NO)
                        YvyraPrint ("%s   Ctype not applied to continuous characters\n", spacer);
                    marks[4] = YES;
                    }
                }
            if (foundIllegal == NO)
                {
                charInfo[i].ctype = charOrdering;
                if (charOrdering == USERTYPE)
                    charInfo[i].userTypeIndex = currentUserTypeIndex;
                numAppliedTo++;
                }
            }
        }
    if (numAppliedTo > 0)
        {
        YvyraPrint ("%s   Ctype was applied to %d standard characters\n", spacer, numAppliedTo);
        }
    else
        {
        YvyraPrint ("%s   No standard characters found to apply ctype to\n", spacer);
        }

    /* mark that analysis needs rebuilding (will be done when mcmc starts) */
    if (defMatrix == YES)
        setUpAnalysisSuccess = NO;

    return (NO_ERROR);
}


/* Variables for usertype parsing */
static char     utName[MAX_USERTYPE_NAME];
static int      utNStates;
static int      utRateCount;
static int      utExpectingRates;
static YFlt   utRates[MAX_STD_STATES * MAX_STD_STATES];

int DoUsertype (void)
{
    int         n, i, j;
    UserType    *ut;

    if (numUserTypes >= MAX_NUM_USERTYPES)
        {
        YvyraPrint ("%s   Maximum number of user types (%d) exceeded\n", spacer, MAX_NUM_USERTYPES);
        return (ERROR);
        }

    if (utRateCount != utNStates * utNStates)
        {
        YvyraPrint ("%s   Expected %d rate values but got %d\n", spacer,
                      utNStates * utNStates, utRateCount);
        return (ERROR);
        }

    /* validate: off-diagonal rates must be non-negative */
    n = utNStates;
    for (i=0; i<n; i++)
        {
        for (j=0; j<n; j++)
            {
            if (i != j && utRates[i*n + j] < 0.0)
                {
                YvyraPrint ("%s   Rate[%d][%d] = %f is negative\n", spacer, i, j, utRates[i*n + j]);
                return (ERROR);
                }
            }
        }

    /* store the usertype */
    ut = &userTypes[numUserTypes];
    strcpy(ut->name, utName);
    ut->nStates = utNStates;
    ut->rates = (YFlt *) SafeCalloc (n * n, sizeof(YFlt));
    if (ut->rates == NULL)
        {
        YvyraPrint ("%s   Problem allocating usertype rates\n", spacer);
        return (ERROR);
        }
    for (i=0; i<n*n; i++)
        ut->rates[i] = utRates[i];

    YvyraPrint ("%s   Defined usertype '%s' with %d states\n", spacer, ut->name, ut->nStates);
    numUserTypes++;

    return (NO_ERROR);
}

int DoUsercost (void)
{
    int         n, i, j;

    /* Convert costs to rates: rate = 1/cost for off-diagonal, 0 for diagonal */
    n = utNStates;
    for (i=0; i<n; i++)
        {
        for (j=0; j<n; j++)
            {
            if (i == j)
                {
                utRates[i*n + j] = 0.0;  /* diagonal ignored */
                }
            else
                {
                YFlt cost = utRates[i*n + j];
                if (cost < 0.0)
                    {
                    YvyraPrint ("%s   Cost[%d][%d] = %f is negative\n", spacer, i, j, cost);
                    return (ERROR);
                    }
                else if (cost == 0.0)
                    {
                    /* cost 0 means forbidden transition */
                    utRates[i*n + j] = 0.0;
                    }
                else
                    {
                    utRates[i*n + j] = 1.0 / cost;
                    }
                }
            }
        }

    YvyraPrint ("%s   Converting cost matrix to rates (rate = 1/cost)\n", spacer);

    /* Now use the same storage logic as DoUsertype */
    return DoUsertype();
}

int DoUsertypeParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        /* first token: the name of the usertype */
        strncpy(utName, tkn, MAX_USERTYPE_NAME - 1);
        utName[MAX_USERTYPE_NAME - 1] = '\0';
        utNStates = 0;
        utRateCount = 0;
        utExpectingRates = NO;
        expecting = Expecting(LEFTPAR) | Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        /* optional (Standard) qualifier */
        expecting = Expecting(ALPHA);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (IsSame("Standard", tkn) != SAME && IsSame("Standard", tkn) != CONSISTENT_WITH)
            {
            YvyraPrint ("%s   Only 'Standard' datatype is supported for usertypes\n", spacer);
            return (ERROR);
            }
        expecting = Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        /* next we expect the number of states */
        utExpectingRates = NO;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (utExpectingRates == NO && utNStates == 0)
            {
            /* this is the nStates value */
            int n;
            sscanf (tkn, "%d", &n);
            if (n < 2 || n > MAX_STD_STATES)
                {
                YvyraPrint ("%s   Number of states (%d) must be between 2 and %d\n", spacer, n, MAX_STD_STATES);
                return (ERROR);
                }
            utNStates = n;
            utExpectingRates = YES;
            expecting = Expecting(NUMBER) | Expecting(SEMICOLON);
            }
        else
            {
            /* this is a rate value */
            YFlt val;
            sscanf (tkn, "%lf", &val);
            if (utRateCount >= utNStates * utNStates)
                {
                YvyraPrint ("%s   Too many rate values (expected %d)\n", spacer, utNStates * utNStates);
                return (ERROR);
                }
            utRates[utRateCount++] = val;
            expecting = Expecting(NUMBER) | Expecting(SEMICOLON);
            }
        }
    else
        return (ERROR);

    return (NO_ERROR);
}

/* Variables for wtset parsing */
static int      wtCharNum;     /* current character number being parsed */
static YFlt   wtValue;       /* current weight value */
static int      wtExpectColon; /* expecting : between char and weight */

int DoWtset (void)
{
    YvyraPrint ("%s   Character weights set\n", spacer);

    /* mark that analysis needs rebuilding */
    if (defMatrix == YES)
        setUpAnalysisSuccess = NO;

    return (NO_ERROR);
}

int DoWtsetParm (char *parmName, char *tkn)
{
    int     tempInt;

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before weights can be set\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        /* first token: * or name */
        wtCharNum = -1;
        wtValue = 1.0;
        wtExpectColon = NO;
        expecting = Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%d", &tempInt);
        if (wtExpectColon == NO)
            {
            /* this is a character number (1-based) */
            if (tempInt < 1 || tempInt > numChar)
                {
                YvyraPrint ("%s   Character number %d out of range (1-%d)\n", spacer, tempInt, numChar);
                return (ERROR);
                }
            wtCharNum = tempInt - 1;  /* convert to 0-based */
            wtExpectColon = YES;
            expecting = Expecting(COLON);
            }
        else
            {
            /* this is a weight value (after colon) */
            YFlt w;
            sscanf (tkn, "%lf", &w);
            if (w < 0.0)
                {
                YvyraPrint ("%s   Weight must be non-negative (got %f)\n", spacer, w);
                return (ERROR);
                }
            charInfo[wtCharNum].weight = w;
            wtExpectColon = NO;
            expecting = Expecting(NUMBER) | Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(COLON))
        {
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}

/* Variables for charlabels parsing */
static int      clCharNum = -1;     /* character index being labeled (0-based) */
static int      clExpectEquals = NO;/* expecting = after character number */

int DoCharlabels (void)
{
    int i, count = 0;

    if (charLabels == NULL)
        {
        YvyraPrint ("%s   No character labels set\n", spacer);
        return (ERROR);
        }

    for (i=0; i<numChar; i++)
        if (charLabels[i] != NULL)
            count++;

    YvyraPrint ("%s   Character labels set for %d of %d characters\n", spacer, count, numChar);

    return (NO_ERROR);
}

static int StoreCharLabel (int charIdx, char *label)
{
    int i;

    /* allocate charLabels array if needed */
    if (charLabels == NULL)
        {
        charLabels = (char **) SafeCalloc ((size_t)numChar, sizeof(char *));
        if (charLabels == NULL)
            {
            YvyraPrint ("%s   Problem allocating charLabels\n", spacer);
            return (ERROR);
            }
        for (i=0; i<numChar; i++)
            charLabels[i] = NULL;
        }

    /* free old label if any */
    if (charLabels[charIdx] != NULL)
        free (charLabels[charIdx]);

    /* store the label */
    charLabels[charIdx] = (char *) SafeCalloc (strlen(label) + 1, sizeof(char));
    if (charLabels[charIdx] == NULL)
        {
        YvyraPrint ("%s   Problem allocating character label\n", spacer);
        return (ERROR);
        }
    strcpy (charLabels[charIdx], label);

    return (NO_ERROR);
}

int DoCharlabelsParm (char *parmName, char *tkn)
{
    int     tempInt;

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before character labels can be set\n", spacer);
        return (ERROR);
        }

    if ((expecting & Expecting(NUMBER)) == Expecting(NUMBER))
        {
        if (clExpectEquals == NO)
            {
            /* this is a character number (1-based) */
            sscanf (tkn, "%d", &tempInt);
            if (tempInt < 1 || tempInt > numChar)
                {
                YvyraPrint ("%s   Character number %d out of range (1-%d)\n", spacer, tempInt, numChar);
                return (ERROR);
                }
            clCharNum = tempInt - 1;
            clExpectEquals = YES;
            expecting = Expecting(EQUALSIGN);
            }
        else
            {
            /* a numeric label (after =) */
            if (clCharNum < 0)
                return (ERROR);
            if (StoreCharLabel(clCharNum, tkn) == ERROR)
                return (ERROR);
            clCharNum = -1;
            clExpectEquals = NO;
            expecting = Expecting(NUMBER) | Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        clExpectEquals = NO;
        expecting = Expecting(ALPHA) | Expecting(NUMBER);
        }
    else if ((expecting & Expecting(ALPHA)) == Expecting(ALPHA))
        {
        if (clCharNum < 0)
            {
            YvyraPrint ("%s   Internal error: no character index set\n", spacer);
            return (ERROR);
            }

        if (StoreCharLabel(clCharNum, tkn) == ERROR)
            return (ERROR);

        clCharNum = -1;
        expecting = Expecting(NUMBER) | Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}

/* Variables for tipweights parsing */
static int      twTaxon = -1;       /* taxon index for current tipweights command */
static int      twChar = -1;        /* character index (0-based) */
static int      twState = -1;       /* current state being weighted */
static int      twExpectColon = NO; /* expecting : after state number */

int DoTipweights (void)
{
    if (twTaxon < 0 || twChar < 0)
        {
        YvyraPrint ("%s   Tipweights command incomplete\n", spacer);
        return (ERROR);
        }

    hasTipWeights = YES;
    YvyraPrint ("%s   Tip weights set for taxon %d, character %d\n", spacer, twTaxon+1, twChar+1);

    /* mark that analysis needs rebuilding */
    if (defMatrix == YES)
        setUpAnalysisSuccess = NO;

    return (NO_ERROR);
}

int DoTipweightsParm (char *parmName, char *tkn)
{
    int     tempInt, i;
    YFlt  tempFloat;

    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before tip weights can be set\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        /* first token: taxon name */
        /* find the taxon */
        twTaxon = -1;
        for (i=0; i<numTaxa; i++)
            {
            if (IsSame(taxaNames[i], tkn) == SAME || IsSame(taxaNames[i], tkn) == CONSISTENT_WITH)
                {
                twTaxon = i;
                break;
                }
            }
        if (twTaxon < 0)
            {
            YvyraPrint ("%s   Could not find taxon '%s'\n", spacer, tkn);
            return (ERROR);
            }
        twChar = -1;
        twExpectColon = NO;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (twChar < 0)
            {
            /* this is the character number (1-based) */
            sscanf (tkn, "%d", &tempInt);
            if (tempInt < 1 || tempInt > numChar)
                {
                YvyraPrint ("%s   Character number %d out of range (1-%d)\n", spacer, tempInt, numChar);
                return (ERROR);
                }
            twChar = tempInt - 1;
            expecting = Expecting(EQUALSIGN);
            }
        else if (twExpectColon == NO)
            {
            /* this is a state number */
            sscanf (tkn, "%d", &tempInt);
            if (tempInt < 0 || tempInt >= MAX_STD_STATES)
                {
                YvyraPrint ("%s   State %d out of range (0-%d)\n", spacer, tempInt, MAX_STD_STATES-1);
                return (ERROR);
                }
            twState = tempInt;
            twExpectColon = YES;
            expecting = Expecting(COLON);
            }
        else
            {
            /* this is a weight value (after colon) */
            sscanf (tkn, "%lf", &tempFloat);
            if (tempFloat < 0.0 || tempFloat > 1.0)
                {
                YvyraPrint ("%s   Weight must be between 0.0 and 1.0 (got %f)\n", spacer, tempFloat);
                return (ERROR);
                }

            /* allocate tipWeights if needed */
            if (tipWeights == NULL)
                {
                tipWeights = (YFlt *) SafeCalloc ((size_t)numTaxa * numChar * MAX_STD_STATES, sizeof(YFlt));
                if (tipWeights == NULL)
                    {
                    YvyraPrint ("%s   Problem allocating tipWeights\n", spacer);
                    return (ERROR);
                    }
                }

            /* store the weight */
            tipWeights[(twTaxon * numChar + twChar) * MAX_STD_STATES + twState] = tempFloat;

            /* also set the bit in the matrix so this state is recognized */
            matrix[pos(twTaxon, twChar, numChar)] |= (1 << twState);

            twExpectColon = NO;
            expecting = Expecting(NUMBER) | Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}

int DoCtypeParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before typesets can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            if (IsSame ("Ordered", tkn) == SAME || IsSame ("Ordered", tkn) == CONSISTENT_WITH)
                charOrdering = ORD;
            else if (IsSame ("Unordered", tkn) == SAME || IsSame ("Unordered", tkn) == CONSISTENT_WITH)
                charOrdering = UNORD;
            else if (IsSame ("Dollo", tkn) == SAME || IsSame ("Dollo", tkn) == CONSISTENT_WITH)
                charOrdering = DOLLO;
            else if (IsSame ("Irreversible", tkn) == SAME || IsSame ("Irreversible", tkn) == CONSISTENT_WITH)
                charOrdering = IRREV;
            else
                {
                /* check if this is a user-defined type name */
                int utIdx;
                int foundUt = NO;
                for (utIdx = 0; utIdx < numUserTypes; utIdx++)
                    {
                    if (IsSame(userTypes[utIdx].name, tkn) == SAME ||
                        IsSame(userTypes[utIdx].name, tkn) == CONSISTENT_WITH)
                        {
                        charOrdering = USERTYPE;
                        currentUserTypeIndex = utIdx;
                        foundUt = YES;
                        break;
                        }
                    }
                if (foundUt == NO)
                    {
                    YvyraPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
                    return (ERROR);
                    }
                }
            
            /* clear tempSet */
            for (i=0; i<numChar; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            YvyraPrint ("%s   Setting characters to %s\n", spacer, tkn);
            expecting = Expecting(COLON);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting  = Expecting(ALPHA) | Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* first, check that we are not trying to put in another character ordering */
        if (IsSame ("Ordered", tkn) == SAME || IsSame ("Ordered", tkn) == CONSISTENT_WITH)
            {
            YvyraPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        else if (IsSame ("Unordered", tkn) == SAME || IsSame ("Unordered", tkn) == CONSISTENT_WITH)
            {
            YvyraPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        else if (IsSame ("Dollo", tkn) == SAME || IsSame ("Dollo", tkn) == CONSISTENT_WITH)
            {
            YvyraPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        else if (IsSame ("Irreversible", tkn) == SAME || IsSame ("Irreversible", tkn) == CONSISTENT_WITH)
            {
            YvyraPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        
        /* We are defining a type set in terms of another (called tkn, here). We should be able
           to find tkn in the list of character set names. If we cannot, then we have a problem and
           return an error. */
        if (IsSame ("All", tkn) == SAME)
            {
            for (i=0; i<numChar; i++)
                tempSet[i] = 1;
            fromI = toJ = everyK = -1;
            expecting = Expecting(SEMICOLON);
            }
        else
            {
            if (numCharSets < 1)
                {
                YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                {
                YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
                
            /* add characters from charset tkn to new tempset */
            for (i=0; i<numChar; i++)
                {
                if (IsBitSet(i, charSet[index]) == YES)
                    tempSet[i] = 1;
                }
            fromI = toJ = everyK = -1;
            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted ctype\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                YvyraPrint ("%s   Improperly formatted ctype\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted ctype\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted ctype\n", spacer);
                    {
                    return (ERROR);
                    }
                }
                
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoDelete (void)
{
    int         i, alreadyDone;

    YvyraPrint ("%s   Excluding taxa\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with taxaset */
    alreadyDone = NO;
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] == 1)
            {
            if (taxaInfo[i].isDeleted == YES && alreadyDone == NO)
                {
                YvyraPrint ("%s   Some taxa already excluded\n", spacer);
                alreadyDone = YES;
                }
            taxaInfo[i].isDeleted = YES;
            }
        }

    SetLocalTaxa ();
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    /* show tempSet (for debugging) */
#   if 0
    for (i=0; i<numTaxa; i++)
        YvyraPrint ("%4d  %4d\n", i+1, tempSet[i]);
#   endif

    return (NO_ERROR);
}


int DoDeleteParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can delete taxa\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = NO;
        for (i=0; i<numTaxa; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 1;
            }
        else
            {
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                /* we are using a pre-defined taxa set */
                if (numTaxaSets < 1)
                    {
                    YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                /* add taxa from taxset tkn to new tempSet */
                for (i=0; i<numTaxa; i++)
                    {
                    if (IsBitSet (i, taxaSet[index]) == YES)
                        tempSet[i] = 1;
                    }
                }
            else
                {
                /* we found the taxon name */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                    
                tempSet[index] = 1;
                }
            }
        foundDash = NO;
        fromI = toJ = everyK = -1;
        
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && !strcmp(tkn, "."))
            tempInt = numTaxa;
        else
            {
            sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numTaxa)
                {
                YvyraPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                return (ERROR);
                }
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted delete set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted delete set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
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


int DoDimensions (void)
{
    if (inDataBlock == NO && inTaxaBlock == NO && inCharactersBlock == NO)
        {
        YvyraPrint ("%s   Dimensions can only be defined in a data, characters or taxa block\n", spacer);
        return (ERROR);
        }

    /* other problems are detected already when reading in DoDimensionsParm */
    if (inDataBlock == YES && (defTaxa == NO || defChars == NO))
        {
        YvyraPrint ("%s   Expecting both Ntax and Nchar to be defined in a data block\n", spacer);
        return (ERROR);
        }

    /* allocate matrix */
    if (inTaxaBlock == YES)
        {
        if (AllocTaxa () == ERROR)
            return ERROR;
        YvyraPrint ("%s   Defining new set of %d taxa\n", spacer, numTaxa);
        }

    if (inCharactersBlock == YES)
        {
        if (AllocMatrix() == ERROR)
            return (ERROR);
        YvyraPrint ("%s   Defining new character matrix with %d characters\n", spacer, numChar);
        }

    if (inDataBlock == YES)
        {
        if (AllocMatrix() == ERROR)
            return (ERROR);
        YvyraPrint ("%s   Defining new matrix with %d taxa and %d characters\n", spacer, numTaxa, numChar);
        }

    return (NO_ERROR);
}


int DoDimensionsParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Ntax (numTaxa) *****************************************************************/
        if (!strcmp(parmName, "Ntax"))
            {
            if (inCharactersBlock == YES)
                {
                YvyraPrint ("%s   You cannot define ntax in a characters block\n");
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &numTaxa);
                defTaxa = YES;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nchar (numChar) ****************************************************************/
        else if (!strcmp(parmName, "Nchar"))
            {
            if (inTaxaBlock == YES)
                {
                YvyraPrint ("%s   You cannot define nchar in a taxa block\n");
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &numChar);
                defChars = YES;
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


int DoDisclaimer (void)
{
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");
    YvyraPrint ("   Disclaimer                                                                    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   Copyright 2003 by John P. Huelsenbeck and Fredrik Ronquist                    \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   This software package is provided \"as is\" and without a warranty of any     \n");
    YvyraPrint ("   kind. In no event shall the authors be held responsible for any damage        \n");
    YvyraPrint ("   resulting from the use of this software. The program--including source code,  \n");
    YvyraPrint ("   example data sets, and executables--is distributed free of charge for         \n");
    YvyraPrint ("   academic use only.                                                            \n");
    YvyraPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoEndBlock (void)
{
    if (inYvyraBlock == YES)
        {
        YvyraPrint ("   Exiting yvyra block\n");
        inYvyraBlock = NO;
        }
    else if (inDataBlock == YES)
        {
        YvyraPrint ("   Exiting data block\n");
        inDataBlock = NO;
        }
    else if (inCharactersBlock == YES)
        {
        YvyraPrint ("   Exiting characters block\n");
        inCharactersBlock = NO;
        }
    else if (inTaxaBlock == YES)
        {
        YvyraPrint ("   Exiting taxa block\n");
        if (numNamedTaxa < numTaxa)
            {
            YvyraPrint ("%s   Leaving taxa block without taxon labels being defined\n", spacer);
            FreeTaxa();
            }
        inTaxaBlock = NO;
        }
    else if (inTreesBlock == YES)
        {
        YvyraPrint ("   Exiting trees block\n");
        inTreesBlock = NO;
        ResetTranslateTable();
        }
    else if (inForeignBlock == YES)
        {
        YvyraPrint ("   Exiting foreign block\n");
        inForeignBlock = NO;
        }
    else
        {
        YvyraPrint ("   Unknown \"end\" statement\n");
        return (ERROR);
        }

    strcpy(spacer,"");  /* reset indentation */
    return (NO_ERROR);
}


int DoExecute (void)
{
    int         rc, cmdLine, lineTerm, longestLineLength, nErrors;
    char        *s;
    FILE        *fp;
    CmdType     *oldCommandPtr;
    char        *oldTokenP, oldToken[CMD_STRING_LENGTH];
#   if defined (MPI_ENABLED)
    int         sumErrors;
#   endif
        
    nErrors = 0;
    cmdLine = 0;
    numOpenExeFiles++;
    s = NULL;
    
    if (numOpenExeFiles > 1)
        YvyraPrint ("\n%s   Executing file \"%s\"...\n\n", spacer, inputFileName);
    else
        YvyraPrint ("%s   Executing file \"%s\"\n", spacer, inputFileName);

    /* Save old command ptr, token pointer and token */
    oldCommandPtr = commandPtr;
    oldTokenP     = tokenP;
    strcpy(oldToken, token);

    /* open binary file */
    if ((fp = OpenBinaryFileR(inputFileName)) == NULL)
        nErrors++;

    /* set indentation to 0 */
    strcpy (spacer, "");
    
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif

    /* find out what type of line termination is used */
    lineTerm = LineTermType (fp);
    if (lineTerm == LINETERM_MAC)
        YvyraPrint ("%s   Macintosh line termination\n", spacer);
    else if (lineTerm == LINETERM_DOS)
        YvyraPrint ("%s   DOS line termination\n", spacer);
    else if (lineTerm == LINETERM_UNIX)
        YvyraPrint ("%s   UNIX line termination\n", spacer);
    else
        {
        YvyraPrint ("%s   Unknown line termination\n", spacer);
        nErrors++;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif
            
    /* find length of longest line */
    longestLineLength = LongestLine (fp);
    YvyraPrint ("%s   Longest line length = %d\n", spacer, longestLineLength);
    longestLineLength += 10;
    
    /* allocate a string long enough to hold a line */
    s = (char *)SafeMalloc((size_t)longestLineLength * sizeof(char));
    if (!s)
        {
        YvyraPrint ("%s   Problem allocating string for reading file\n", spacer);
        nErrors++;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif

    /* close binary file */
    SafeFclose (&fp);
    
    /* open text file */
    if ((fp = OpenTextFileR(inputFileName)) == NULL)
        {
        nErrors++;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif
    
    /* parse file, reading each line in turn */
    YvyraPrint ("%s   Parsing file\n", spacer);

    inYvyraBlock = inDataBlock = inForeignBlock = inTreesBlock = NO;
    foundNewLine = NO;
    expecting = Expecting(COMMAND);
    cmdLine = 0;

    /* read lines into s until end of file */
    while ( fgets(s, longestLineLength, fp) != NULL)
        {
        foundNewLine = YES;
        cmdLine++;

        /* process string if not empty */
        if (strlen(s) > 1)
            {
            /* check that all characters in the string are valid */
            if (CheckStringValidity (s) == ERROR)
                {
                nErrors++;
                }
#           if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
                goto errorExit;
                }
#           else
            if (nErrors > 0)
                goto errorExit;
#           endif
                
            /* interpret commands on line */
            rc = ParseCommand (s);
            if (rc == ERROR)
                nErrors++;
#           if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
                goto errorExit;
                }
#           else
            if (nErrors > 0)
                goto errorExit;
#           endif
            if (rc == NO_ERROR_QUIT)
                nErrors++;
#           if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                goto quitExit;
#           else
            if (nErrors > 0)
                goto quitExit;
#           endif
            }
        }
    
    YvyraPrint ("%s   Reached end of file\n", spacer);

    if (inComment == YES)
        nErrors++;

#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        YvyraPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif

    if (s)
        free (s);
    SafeFclose (&fp);
    numOpenExeFiles--;

    if (numOpenExeFiles > 0)
        {
        inYvyraBlock = YES;
        YvyraPrint ("\n   Returning execution to calling file ...\n\n");
        strcpy (spacer, "   ");
        }
    else
        strcpy (spacer, "");

    commandPtr = oldCommandPtr;

    return (NO_ERROR);
    
    quitExit:
        if (s)
            free (s);
        SafeFclose (&fp);
        numOpenExeFiles--;
        if (numOpenExeFiles > 0)
            {
            inYvyraBlock = YES;
            strcpy (spacer, "   ");
            }
        else
            strcpy (spacer, "");

        commandPtr = oldCommandPtr;
        tokenP     = oldTokenP;
        strcpy(token, oldToken);

        return (NO_ERROR_QUIT);
            
    errorExit:
        if (inComment == YES)
            {
            YvyraPrint ("%s   ERROR: Reached end of file while in comment.\n", spacer);
            inComment = NO;
            numComments = 0;
            }
        if (fp)
            {
            YvyraPrint ("%s   The error occurred when reading char. %d-%d on line %d\n", spacer,
                (size_t)(tokenP-s)-strlen(token)+1, (size_t)(tokenP-s), cmdLine);
            YvyraPrint ("%s      in the file '%s'\n", spacer, inputFileName);
            }
        if (s)
            free (s);
        SafeFclose (&fp);

        /* make sure we exit the block we were reading from correctly */
        if (inYvyraBlock == YES)
            inYvyraBlock = NO;
        else if (inDataBlock == YES)
            inDataBlock = NO;
        else if (inTreesBlock == YES)
            {
            inTreesBlock = NO;
            ResetTranslateTable();
            }
        else if (inForeignBlock == YES)
            inForeignBlock = NO;

        /* make sure correct return if we came from yvyra block in another execute file */
        if (numOpenExeFiles > 1)
            {
            inYvyraBlock = YES;
            YvyraPrint ("\n   Returning execution to calling file ...\n\n");
            strcpy (spacer, "   ");
            }
        else
            {
            strcpy (spacer, "");
            YvyraPrint ("\n   Returning execution to command line ...\n\n");
            }

        numOpenExeFiles--;  /* we increase the value above even if no file is successfully opened */

        /* restore state of globals */
        commandPtr = oldCommandPtr;
        tokenP = oldTokenP;
        strcpy(token, oldToken);

        return (ERROR);
}


int DoExecuteParm (char *parmName, char *tkn)
{
    if (strlen(tkn)>99)
        {
        YvyraPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
        YvyraPrint ("%s      '%s'\n", spacer,tkn);
        YvyraPrint ("%s   has %d characters.\n", spacer, strlen(tkn));
        return (ERROR);
        }
    strcpy (inputFileName, tkn);
    
    expecting = Expecting (SEMICOLON);

    return (NO_ERROR);
}


int DoExclude (void)
{
    int         i, alreadyDone;

    YvyraPrint ("%s   Excluding character(s)\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with charset */
    alreadyDone = NO;
    for (i=0; i<numChar; i++)
        {
        if (tempSet[i] == 1)
            {
            if (charInfo[i].isExcluded == YES && alreadyDone == NO)
                {
                YvyraPrint ("%s   Some characters already excluded\n", spacer);
                alreadyDone = YES;
                }
            charInfo[i].isExcluded = YES;
            }
        }
        
    foundFirst = NO;

    /* reset analysis to recompress data */
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    return (NO_ERROR);
}


int DoExcludeParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can exclude characters\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;
        for (i=0; i<numChar; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                tempSet[i] = 1;
            }
        else if (IsSame ("Missambig", tkn) == SAME || IsSame ("Missambig", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                {
                if (charInfo[i].isMissAmbig == YES)
                    tempSet[i] = 1;
                }
            }
        else
            {
            /* we are using a pre-defined character set */
            if (numCharSets < 1)
                {
                YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                {
                YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add characters from charset tkn to new tempSet */
            for (i=0; i<numChar; i++)
                {
                if (IsBitSet(i, charSet[index]) == YES)
                    tempSet[i] = 1;
                }
            fromI = toJ = everyK = -1;
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted exclude set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                YvyraPrint ("%s   Improperly formatted exclude set\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted exclude set\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted exclude set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoFormat (void)
{
    if (inDataBlock == NO && inCharactersBlock == NO)
        {
        YvyraPrint ("%s   Formats can only be defined in a data or characters block\n", spacer);
        return (ERROR);
        }

    return CheckInitialPartitions();
}


int DoFormatParm (char *parmName, char *tkn)
{
    int         i, tempInt;
    char        tempStr[100];
    
    if (inDataBlock == NO && inCharactersBlock == NO)
        {
        YvyraPrint ("%s   Formats can only be defined in a data or characters block\n", spacer);
        return (ERROR);
        }
    if (defTaxa == NO || defChars == NO)
        {
        YvyraPrint ("%s   The dimensions of the matrix must be defined before the format\n", spacer);
        return (ERROR);
        }
    
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        if (!strcmp(parmName, "Interleave"))
            {
            expecting = Expecting(EQUALSIGN) | Expecting(PARAMETER) | Expecting(SEMICOLON);
            isInterleaved = YES;
            }
        }
    else
        {
        /* set Datatype (dataType) ************************************************************/
        if (!strcmp(parmName, "Datatype"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (isMixed == NO)
                        {
                        if (!strcmp(tempStr, "Standard"))
                            dataType = STANDARD;
                        else if (!strcmp(tempStr, "Dna") || !strcmp(tempStr, "Rna") || !strcmp(tempStr, "Protein") || !strcmp(tempStr, "Restriction"))
                            {
                            YvyraPrint ("%s   Molecular datatypes (DNA, RNA, Protein, Restriction) are not supported in yvyra.\n", spacer);
                            YvyraPrint ("%s   yvyra only supports 'datatype=Standard' for morphological/linguistic data.\n", spacer);
                            return ERROR;
                            }
                        else if (!strcmp(tempStr, "Continuous"))
                            {
                            YvyraPrint ("%s   Continuous datatype is not yet supported in yvyra.\n", spacer);
                            return ERROR;
                            }
                        else if (!strcmp(tempStr, "Mixed"))
                            {
                            dataType = MIXED;
                            isMixed = YES;
                            for (i=0; i<numChar; i++)
                                tempSet[i] = 0;
                            fromI = toJ = everyK = -1;
                            foundDash = foundSlash = NO;
                            numDivisions = 0;
                            YvyraPrint ("%s   Data is Mixed\n", spacer);
                            }
                        if (dataType == MIXED)
                            expecting = Expecting(LEFTPAR);
                        else
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                        }
                    else
                        {
                        if (!strcmp(tempStr, "Standard"))
                            dataType = STANDARD;
                        else if (!strcmp(tempStr, "Dna") || !strcmp(tempStr, "Rna") || !strcmp(tempStr, "Protein") || !strcmp(tempStr, "Restriction"))
                            {
                            YvyraPrint ("%s   Molecular datatypes (DNA, RNA, Protein, Restriction) are not supported in yvyra.\n", spacer);
                            YvyraPrint ("%s   yvyra only supports 'datatype=Standard' for morphological/linguistic data.\n", spacer);
                            return ERROR;
                            }
                        else if (!strcmp(tempStr, "Continuous"))
                            {
                            YvyraPrint ("%s   Continuous datatype is not yet supported in yvyra.\n", spacer);
                            return ERROR;
                            }
                        else if (!strcmp(tempStr, "Mixed"))
                            {
                            YvyraPrint ("%s   Cannot have mixed datatype within a mixed datatype\n", spacer);
                            return (ERROR);
                            }
                        expecting = Expecting(COLON);
                        for (i=0; i<numChar; i++)
                            tempSet[i] = 0;
                        fromI = toJ = everyK = -1;
                        foundDash = foundSlash = NO;
                        }
                    if (isMixed == NO)
                        {
                        numDivisions = 1;
                        for (i=0; i<numChar; i++)
                            {
                            charInfo[i].charType = dataType;
                            partitionId[i][0] = numDivisions;
                            }
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid data type argument\n", spacer);
                    return (ERROR);
                    }
                if (isMixed == NO)
                    YvyraPrint ("%s   Data is %s\n", spacer, tempStr);
                else if (strcmp(tempStr, "Mixed"))
                    YvyraPrint ("%s      Data for partition %d is %s\n", spacer, numDivisions+1, tempStr);
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(COLON))
                {
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                if (strlen(tkn) == 1 && tkn[0] == '.')
                    tempInt = numChar;
                else
                    sscanf (tkn, "%d", &tempInt);
                if (tempInt <= 0 || tempInt > numChar)
                    {
                    YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
                    return (ERROR);
                    }
                tempInt--;
                if (foundDash == YES)
                    {
                    if (fromI >= 0)
                        toJ = tempInt;
                    else
                        {
                        YvyraPrint ("%s   Improperly formatted partition\n", spacer);
                        return (ERROR);
                        }
                    foundDash = NO;
                    }
                else if (foundSlash == YES)
                    {
                    tempInt++;
                    if (tempInt <= 1)
                        {
                        YvyraPrint ("%s   Improperly formatted partition\n", spacer);
                        return (ERROR);
                        }
                    if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                        everyK = tempInt;
                    else
                        {
                        YvyraPrint ("%s   Improperly formatted partition\n", spacer);
                        return (ERROR);
                        }
                    foundSlash = NO;
                    }
                else
                    {
                    if (fromI >= 0 && toJ < 0)
                        {
                        if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                            return (ERROR);
                        fromI = tempInt;
                        }
                    else if (fromI < 0 && toJ < 0)
                        {
                        fromI = tempInt;
                        }
                    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                        {
                        if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                            return (ERROR);
                        fromI = tempInt;
                        toJ = everyK = -1;
                        }
                    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                        {
                        if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                            return (ERROR);
                        fromI = tempInt;
                        toJ = everyK = -1;
                        }
                    else
                        {
                        YvyraPrint ("%s   Improperly formatted partition\n", spacer);
                            {
                            return (ERROR);
                            }
                        }
                        
                    }
                expecting  = Expecting(NUMBER);
                expecting |= Expecting(DASH);
                expecting |= Expecting(BACKSLASH);
                expecting |= Expecting(COMMA);
                expecting |= Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(BACKSLASH))
                {
                foundSlash = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(COMMA))
                {
                /* add set to tempSet */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                for (i=0; i<numChar; i++)
                    {
                    if (tempSet[i] == numDivisions)
                        charInfo[i].charType = dataType;
                    }

                /* merge tempSet */
                for (i=0; i<numChar; i++)
                    {
                    if (tempSet[i] != 0)
                        {
                        if (partitionId[i][0] == 0)
                            {
                            charInfo[i].charType = dataType;
                            partitionId[i][0] = numDivisions + 1;
                            }
                        else
                            {
                            YvyraPrint ("%s   Improperly formatted partition (same character found in multiple partitions)\n", spacer);
                            return (ERROR);
                            }
                        }
                    }
                
                /* increment number of partitions */
                numDivisions++;             
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                /* add set to tempSet */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                    
                /* merge tempSet */
                for (i=0; i<numChar; i++)
                    {
                    if (tempSet[i] != 0)
                        {
                        if (partitionId[i][0] == 0)
                            {
                            charInfo[i].charType = dataType;
                            partitionId[i][0] = numDivisions + 1;
                            }
                        else
                            {
                            YvyraPrint ("%s   Improperly formatted partition (same character found in multiple partitions)\n", spacer);
                            return (ERROR);
                            }
                        }
                    }
                
                /* increment number of partitions */
                numDivisions++;             
                if (isMixed == YES)
                    dataType = MIXED;
                    
                if (numDivisions > 1)
                    YvyraPrint ("%s   There are a total of %d default data divisions\n", spacer, numDivisions);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Interleave (isInterleaved) *****************************************************/
        else if (!strcmp(parmName, "Interleave"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        isInterleaved = YES;
                    else
                        isInterleaved = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for interleaved data\n", spacer);
                    return (ERROR);
                    }
                if (isInterleaved == YES)
                    YvyraPrint ("%s   Data matrix is interleaved\n", spacer);
                else
                    YvyraPrint ("%s   Data matrix is not interleaved\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Gap (gapId) ********************************************************************/
        else if (!strcmp(parmName, "Gap"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting  = Expecting(ALPHA);
                expecting |= Expecting(QUESTIONMARK);
                expecting |= Expecting(DASH);
                expecting |= Expecting(NUMBER);
                expecting |= Expecting(ASTERISK);
                expecting |= Expecting(EXCLAMATIONMARK);
                expecting |= Expecting(PERCENT);
                expecting |= Expecting(WEIRD);
                expecting |= Expecting(VERTICALBAR);
                }
            else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
                     ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
                     ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
                     ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
                     ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
                     ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
                     ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
                     ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)) ||
                     ((expecting & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR)))
                {
                if (strlen(tkn) == 1)
                    {
                    if (tkn[0] == matchId || tkn[0] == missingId)
                        {
                        YvyraPrint ("%s   Gap character matches matching or missing characters\n", spacer);
                        return (ERROR);
                        }
                    gapId = tkn[0];
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid gap argument %s\n", spacer, tkn);
                    return (ERROR);
                    }
                YvyraPrint ("%s   Gaps coded as %s\n", spacer, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Missing (missingId) ************************************************************/
        else if (!strcmp(parmName, "Missing"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting  = Expecting(ALPHA);
                expecting |= Expecting(QUESTIONMARK);
                expecting |= Expecting(DASH);
                expecting |= Expecting(NUMBER);
                expecting |= Expecting(ASTERISK);
                expecting |= Expecting(EXCLAMATIONMARK);
                expecting |= Expecting(PERCENT);
                expecting |= Expecting(WEIRD);
                expecting |= Expecting(VERTICALBAR);
                }
            else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
                     ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
                     ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
                     ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
                     ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
                     ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
                     ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
                     ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)) ||
                     ((expecting & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR)))
                {
                if (strlen(tkn) == 1)
                    {
                    if (tkn[0] == gapId || tkn[0] == matchId)
                        {
                        YvyraPrint ("%s   Missing character matches matching or gap characters\n", spacer);
                        return (ERROR);
                        }
                    missingId = tkn[0];
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid missing argument %s\n", spacer, tkn);
                    return (ERROR);
                    }
                YvyraPrint ("%s   Missing data coded as %s\n", spacer, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Matchchar (matchId) ************************************************************/
        else if (!strcmp(parmName, "Matchchar"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting  = Expecting(ALPHA);
                expecting |= Expecting(QUESTIONMARK);
                expecting |= Expecting(DASH);
                expecting |= Expecting(NUMBER);
                expecting |= Expecting(ASTERISK);
                expecting |= Expecting(EXCLAMATIONMARK);
                expecting |= Expecting(PERCENT);
                expecting |= Expecting(WEIRD);
                expecting |= Expecting(VERTICALBAR);
                }
            else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
                     ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
                     ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
                     ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
                     ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
                     ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
                     ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
                     ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)) ||
                     ((expecting & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR)))
                {
                if (strlen(tkn) == 1)
                    {
                    if (tkn[0] == gapId || tkn[0] == missingId)
                        {
                        YvyraPrint ("%s   Matching character matches gap or missing characters\n", spacer);
                        return (ERROR);
                        }
                    matchId = tkn[0];
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid matchchar argument %s\n", spacer, tkn);
                    return (ERROR);
                    }
                YvyraPrint ("%s   Matching characters coded as %s\n", spacer, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* skip Symbols ***************************************************************/
        else if (!strcmp(parmName, "Symbols"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                YvyraPrint ("%s   WARNING: yvyra does not support 'symbols' specification; default symbols assumed\n", spacer);
                readWord=YES;
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* on Equate return ERROR ***************************************************************/
        else if (!strcmp(parmName, "Equate"))
            {
            YvyraPrint ("%s   ERROR: yvyra does not support 'Equate' macros; please remove or comment out\n", spacer);
            return (ERROR);
            }
        else
            return (ERROR);
        }

    return (NO_ERROR);
}


int DoHelp (void)
{
    int         i, j, longestDescription;
    CmdType     *p;

    if (foundFirst == NO)
        {
        longestDescription = 0;
        for (i=1; i<NUMCOMMANDS; i++)
            {
            p = commands + i;
            if ((int)strlen(p->string) > longestDescription)
                longestDescription = (int) strlen(p->string);
            }
        
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Commands that are available from the command                                  \n");
        YvyraPrint ("   line or from a yvyra block include:                                         \n");
        YvyraPrint ("                                                                                 \n");
        for (i=1; i<NUMCOMMANDS; i++)
            {
            p = commands + i;
            if (p->cmdUse == IN_CMD && p->hiding == SHOW)
                {
                YvyraPrint ("   %s", p->string);
                for (j=0; j<longestDescription - (int) strlen(p->string); j++)
                YvyraPrint (" ");
                YvyraPrint (" -- %s\n", p->cmdDescription);
                }
            }
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Commands that should be in a NEXUS file (data                                 \n");
        YvyraPrint ("   block, trees block or taxa block) include:                                                \n");
        YvyraPrint ("                                                                                 \n");
        for (i=1; i<NUMCOMMANDS; i++)
            {
            p = commands + i;
            if (p->cmdUse == IN_FILE && p->hiding == SHOW)
                {
                YvyraPrint ("   %s", p->string);
                for (j=0; j<longestDescription - (int) strlen(p->string); j++)
                YvyraPrint (" ");
                YvyraPrint (" -- %s\n", p->cmdDescription);
                }
            }
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that this program supports the use of the shortest unambiguous           \n"); 
        YvyraPrint ("   spelling of the above commands (e.g., \"exe\" instead of \"execute\").        \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    foundFirst = NO;

    return (NO_ERROR);
}


int DoHelpParm (char *parmName, char *tkn)
{
    int         i, j, tkLen, targetLen, numDiff, numMatches;
    CmdType     *p, *q=NULL;

    if (expecting == Expecting(ALPHA))
        {
        p = commands + 0;
        tkLen = (int) strlen(tkn);
        numMatches = 0;
        for (i=0; i<NUMCOMMANDS; i++)
            {
            targetLen = (int) strlen(p->string);
            if (tkLen <= targetLen)
                {
                for (j=0, numDiff=0; j<tkLen; j++)
                    {
                    if (ChangeCase(tkn[j]) != ChangeCase(p->string[j]))
                        numDiff++;
                    }
                if (numDiff == 0)
                    {
                    numMatches++;
                    q = p;
                    if (tkLen == targetLen)
                        break;
                    }
                }       
            p++;
            }
        if (numMatches == 0)
            {
            YvyraPrint ("%s   Could not find command \"%s\"\n", spacer, tkn);
            return (ERROR);
            }
        else if (numMatches == 1)
            {
            if (GetUserHelp (q->string) == ERROR)
                {
                YvyraPrint ("%s   Problem getting help for command \"%s\"\n", spacer, q->string);
                }
            }
        else 
            {
            YvyraPrint ("%s   Ambiguous command \"%s\"\n", spacer, tkn);
            return (ERROR);
            }
            
        expecting = Expecting(SEMICOLON);
        foundFirst = YES;
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoInclude (void)
{
    int         i, alreadyDone;

    YvyraPrint ("%s   Including character(s)\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with excludedChars */
    alreadyDone = NO;
    for (i=0; i<numChar; i++)
        {
        if (tempSet[i] == 1)
            {
            if (charInfo[i].isExcluded == NO && alreadyDone == NO)  
                {
                YvyraPrint ("%s   Some characters already included\n", spacer);
                alreadyDone = YES;
                }
            charInfo[i].isExcluded = NO;
            }
        }

    /* reset analysis to recompress data */
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    return (NO_ERROR);
}


int DoIncludeParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can include characters\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;
        for (i=0; i<numChar; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                tempSet[i] = 1;
            }
        else if (IsSame ("Missambig", tkn) == SAME || IsSame ("Missambig", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                {
                if (charInfo[i].isMissAmbig == YES)
                    tempSet[i] = 1;
                }
            }
        else
            {
            /* we are using a pre-defined character set */
            if (numCharSets < 1)
                {
                YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                {
                YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add characters from charset tkn to new tempSet */
            for (i=0; i<numChar; i++)
                {
                if (IsBitSet(i, charSet[index]) == YES)
                    tempSet[i] = 1;
                }
            fromI = toJ = everyK = -1;
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted include set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                YvyraPrint ("%s   Improperly formatted include set\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted include set\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted include set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoLog (void)
{
    if (logToFile == YES)
        {
        SafeFclose (&logFileFp);
        if (replaceLogFile == YES)
            {
            if ((logFileFp = OpenTextFileW (logFileName)) == NULL)  
                {
                logToFile = NO;
                return (ERROR);
                }
            }
        else
            {
            if ((logFileFp = OpenTextFileA (logFileName)) == NULL)  
                {
                logToFile = NO;
                return (ERROR);
                }
            }
        YvyraPrint ("%s   Logging screen output to file \"%s\"\n", spacer, logFileName);
        }
    else
        {
        SafeFclose (&logFileFp);
        YvyraPrint ("%s   Terminating log output\n", spacer);
        }

    return (NO_ERROR);
}


int DoLogParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Start"))
            {
            if (logToFile == YES)
                YvyraPrint ("%s   Logging to file is already on\n", spacer, logFileName);
            else
                logToFile = YES;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else if (!strcmp(parmName, "Stop"))
            {
            if (logToFile == NO)
                YvyraPrint ("%s   Logging to file is already off\n", spacer, logFileName);
            else
                logToFile = NO;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else if (!strcmp(parmName, "Replace"))
            {
            replaceLogFile = YES;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else if (!strcmp(parmName, "Append"))
            {
            replaceLogFile = NO;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else
            expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (logFileName, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else
            {
            YvyraPrint ("%s   Unknown parameter in Log\n", spacer);
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


int DoManual (void)
{
    int     i, j, logSetting;
    char    title[100];
    FILE    *fp, *logfp;
    CmdType *p;
    
    /* try to open file, return error if present */
    if ((fp = OpenTextFileRQuait(manFileName)) != NULL)
        {
        YvyraPrint ("%s   File \"%s\" already exists \n", spacer, manFileName);
        SafeFclose(&fp);
        return (ERROR);
        }

    /* try to open file for writing, return error if not possible */
    if ((fp = OpenTextFileW(manFileName)) == NULL)
        return (ERROR);

    /* print message */
    YvyraPrint ("%s   Producing command reference file \"%s\"\n", spacer, manFileName);

    /* temporarily disable normal logging and switch echoing off */
    logSetting = logToFile;
    logfp = logFileFp;
    echoMB = NO;
    logToFile = YES;
    logFileFp = fp;
    
    /* produce command reference file */
    /* header */
    strcpy (title, "Command Reference for yvyra ver. ");
    strcat (title, VERSION_NUMBER);

    i = (70 - (int) strlen (title)) / 2;
    j = 70 - i - (int) strlen(title);

    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("      %*c%s%*c      \n", i, ' ', title, j, ' ');
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("                   (c) John P. Huelsenbeck, Fredrik Ronquist                     \n");
    YvyraPrint ("                               and Maxim Teslenko                                \n");
    YvyraPrint ("                                                                                 \n");

    /* summary */
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   ***************************************************************************   \n");
    YvyraPrint ("   *                                                                         *   \n");
    YvyraPrint ("   *  1. Command summary                                                     *   \n");
    YvyraPrint ("   *                                                                         *   \n");
    YvyraPrint ("   ***************************************************************************   \n");
    YvyraPrint ("                                                                                 \n");
    foundFirst = NO;
    if (DoHelp() == ERROR)
        {
        YvyraPrint ("%s   Could not produce command reference summary\n", spacer);
        goto errorExit;
        }
    
    /* list of yvyra commands */
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   ***************************************************************************   \n");
    YvyraPrint ("   *                                                                         *   \n");
    YvyraPrint ("   *  2. yvyra commands                                                    *   \n");
    YvyraPrint ("   *                                                                         *   \n");
    YvyraPrint ("   ***************************************************************************   \n");
    YvyraPrint ("                                                                                 \n");
    for (i=1; i<NUMCOMMANDS; i++)
        {
        p = commands + i;
        if (p->cmdUse == IN_CMD && p->hiding == SHOW)
            {
            if (GetUserHelp(p->string)==ERROR)
                goto errorExit;
            }
        }

    /* list of data or tree block commands */
    YvyraPrint ("                                                                                 \n");
    YvyraPrint ("   ***************************************************************************   \n");
    YvyraPrint ("   *                                                                         *   \n");
    YvyraPrint ("   *  3. 'Data' or 'tree' block commands (in #NEXUS file)                    *   \n");
    YvyraPrint ("   *                                                                         *   \n");
    YvyraPrint ("   ***************************************************************************   \n");
    YvyraPrint ("                                                                                 \n");
    for (i=1; i<NUMCOMMANDS; i++)
        {
        p = commands + i;
        if (p->cmdUse == IN_FILE && p->hiding == SHOW)
            {
            if (GetUserHelp(p->string) == ERROR)
                goto errorExit;
            }
        }

    /* return logging to previous settings and switch echoing on */
    SafeFclose (&fp);
    logToFile = logSetting;
    logFileFp = logfp;
    echoMB = YES;

    YvyraPrint ("%s   Successfully produced command reference file \"%s\"\n", spacer, manFileName);

    return (NO_ERROR);

    errorExit:
        SafeFclose (&fp);
        logToFile = logSetting;
        logFileFp = logfp;
        echoMB = YES;

        return (ERROR);
}


int DoManualParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (manFileName, tkn);
                expecting = Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else
            {
            YvyraPrint ("%s   Unknown parameter in Manual\n", spacer);
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


int DoMatrix (void)
{
    int         i, j, hasMissingAmbig;
    
    if (taxonCount != numTaxa)
        {
        YvyraPrint ("%s   Problem with number of taxa read in (%d taxa read in, while expecting %d)\n", spacer, taxonCount, numTaxa);
        FreeMatrix();
        return (ERROR);
        }
    for (i=0; i<numTaxa; i++)
        {
        if (taxaInfo[i].charCount != numChar)
            {
            YvyraPrint ("%s   Problem with number of characters read in (%d expected for taxon %d, %d read in)\n", spacer, numChar, i, taxaInfo[i].charCount);
            FreeMatrix();
            return (ERROR);
            }
        }
        
    /* find out which characters have missing or ambiguous states (one time only, so no special function) */
    for (i=0; i<numChar; i++)
        {
        hasMissingAmbig = NO;
        for (j=0; j<numTaxa; j++)
            {
            if (IsMissing (matrix[pos(j,i,numChar)], charInfo[i].charType) == YES)
                hasMissingAmbig = YES;
            if (IsAmbig (matrix[pos(j,i,numChar)], charInfo[i].charType) == YES)
                hasMissingAmbig = YES;
            }
        if (hasMissingAmbig == YES)
            charInfo[i].isMissAmbig = YES;
        }

    YvyraPrint ("%s   Successfully read matrix\n", spacer);
    if (matrixHasPoly == YES)
        YvyraPrint ("%s   Matrix contains polymorphisms, interpreted as ambiguity\n", spacer);
    defMatrix = YES;
    isTaxsetDef = YES;

    /* add name of default partition */
    if (AddString (&partitionNames, 0, "Default") == ERROR)
        {
        YvyraPrint ("%s   Problem adding Default name to partition list\n", spacer);
        return (ERROR);
        }
    numDefinedPartitions = 1;

    if (numDefinedSpeciespartitions == 0)   /* the default species partition could have been added already in DoTaxLabels */
        {
        /* add default speciespartition name to list of valid speciespartitions */
        if (AddString (&speciespartitionNames, 0, "Default") == ERROR)
            {
            YvyraPrint ("%s   Problem adding Default speciespartition to list\n", spacer);
            return (ERROR);
            }

        /* add default species name set */
        AddNameSet(&speciesNameSets, 0, taxaNames, numTaxa);

        /* set number of defined speciespartitions to 1 */
        numDefinedSpeciespartitions = 1;
        }
        
    if (SetPartition (0) == ERROR)
        return ERROR;

    if (SetSpeciespartition (0) == ERROR)
        return ERROR;

    if (numCurrentDivisions == 1)
        YvyraPrint ("%s   Setting default partition (does not divide up characters)\n", spacer);
    else
        YvyraPrint ("%s   Setting default partition, dividing characters into %d parts\n", spacer, numCurrentDivisions);
    
    if (SetModelDefaults () == ERROR)
        return (ERROR);

    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);

    /* set default names for some output file names based on processed file */
    strcpy (sumtParams.sumtFileName, inputFileName);
    strcpy (sumtParams.sumtOutfile, inputFileName);
    strcpy (sumpParams.sumpFileName, inputFileName);
    strcpy (sumpParams.sumpOutfile, inputFileName);
    strcpy (comptreeParams.comptOutfile, inputFileName);

    if (chainParams.numRuns == 1)
        {
        sprintf (comptreeParams.comptFileName1, "%s.t", inputFileName);
        sprintf (comptreeParams.comptFileName2, "%s.t", inputFileName);
        }
    else /* if (chainParams.numRuns > 1) */
        {
        sprintf (comptreeParams.comptFileName1, "%s.run1.t", inputFileName);
        sprintf (comptreeParams.comptFileName2, "%s.run2.t", inputFileName);
        }

    if (chainParams.numRuns == 1)
        sprintf (plotParams.plotFileName, "%s.p", inputFileName);
    else /* if (chainParams.numRuns > 1) */
        sprintf (plotParams.plotFileName, "%s.run1.p", inputFileName);

    strcpy (chainParams.chainFileName, inputFileName);

    if (chainParams.numRuns > 1)
        YvyraPrint ("%s   Setting output file names to \"%s.run<i>.<p|t>\"\n", spacer, chainParams.chainFileName);
    else
        YvyraPrint ("%s   Setting output file names to \"%s.<p|t>\"\n", spacer, chainParams.chainFileName);

#   if 0
    for (i=0; i<numChar; i++)
        {
        int     j;
        YvyraPrint ("%4d -- ", i+1);
        for (j=0; j<numTaxa; j++)
            YvyraPrint ("%2d ", matrix[pos(j,i,numChar)]);
        YvyraPrint ("\n");
        }
#   endif

    return (NO_ERROR);
}


int DoMatrixParm (char *parmName, char *tkn)
{
    int             i, j, charCode=0, index;
    YFlt          charValue;

    expecting  = Expecting(ALPHA);
    expecting |= Expecting(QUESTIONMARK);
    expecting |= Expecting(DASH);
    expecting |= Expecting(NUMBER);
    expecting |= Expecting(ASTERISK);
    expecting |= Expecting(EXCLAMATIONMARK);
    expecting |= Expecting(PERCENT);
    expecting |= Expecting(WEIRD);
    expecting |= Expecting(VERTICALBAR);
    expecting |= Expecting(SEMICOLON);
    expecting |= Expecting(LEFTPAR);
    expecting |= Expecting(RIGHTPAR);
    expecting |= Expecting(LEFTCURL);
    expecting |= Expecting(RIGHTCURL);

    if (defTaxa == NO || defChars == NO)
        {
        YvyraPrint ("%s   Number of taxa and characters needs to be defined before matrix is read\n", spacer);
        goto errorExit;
        }
    if (inDataBlock == NO && inCharactersBlock == NO)
        {
        YvyraPrint ("%s   Must be in data or characters block to read in character matrix\n", spacer);
        goto errorExit;
        }

    if (isFirstMatrixRead == YES)
        {
        foundNewLine = YES;
        isFirstInterleavedBlock = YES;
        taxonCount = 0;
        isNegative = NO;
        }
    isFirstMatrixRead = NO;
    
    /* allow line breaks in non-interleaved matrices */
    if (isInterleaved == NO)
        {
        if (foundNewLine == YES && taxonCount > 0)
            {
            if (taxaInfo[taxonCount-1].charCount < numChar)
                foundNewLine = NO;
            }
        }

    if (taxonCount >= numTaxa && foundNewLine == YES)
        {
        if (isInterleaved == YES)
            {
            taxonCount = 0;
            isFirstInterleavedBlock = NO;
            }
        else
            {
            YvyraPrint ("%s   Too many taxa in matrix\n", spacer);
            goto errorExit;
            }
        }
    
    if (taxaInfo[0].charCount > 4010)
        i = 1;  /* FIXME: Not used (from clang static analyzer) */

    if (foundNewLine == YES)
        {
        /* Should be a taxon. */
        if (isFirstInterleavedBlock == YES)
            {
            /* If this is the first interleaved block, then we need to add the taxon
               to the set of taxon names unless there is already a defined taxon set. */
            if (strlen(tkn)>99)
                {
                YvyraPrint ("%s   Taxon name %s is too long. Maximum 99 characters is allowed.\n", spacer, tkn);
                goto errorExit;
                }
            if (isTaxsetDef == NO && AddString (&taxaNames, taxonCount, tkn) == ERROR)
                {
                YvyraPrint ("%s   Problem adding taxon %s to taxon set\n", spacer, tkn);
                goto errorExit;
                }
            if (numTaxa < 10)
                YvyraPrint ("%s   Taxon %d -> %s\n", spacer, taxonCount+1, tkn);
            else if (numTaxa < 100 && numTaxa >= 10)
                YvyraPrint ("%s   Taxon %2d -> %s\n", spacer, taxonCount+1, tkn);
            else if (numTaxa < 1000 && numTaxa >= 100)
                YvyraPrint ("%s   Taxon %3d -> %s\n", spacer, taxonCount+1, tkn);
            else
                YvyraPrint ("%s   Taxon %4d -> %s\n", spacer, taxonCount+1, tkn);
            }
        else
            {
            /* If this is not the first interleaved block, then we need to
               check to see if taxon name is present and in correct place. */
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                YvyraPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
                goto errorExit;
                }
            if (index != taxonCount)
                {
                YvyraPrint ("%s   Could not find taxon %s in correct position in list of taxa\n", spacer, tkn);
                goto errorExit;
                }
            }
        foundNewLine = NO;
        isNegative = NO;
        taxonCount++;
        }
    else
        {
        /* Should be a character (either continuous or otherwise). */
        if (charInfo[taxaInfo[taxonCount-1].charCount].charType == CONTINUOUS)
            {
            /* If we have a CONTINUOUS character, then the entire token should either be
               a number or a dash (for a negative sign). */
            if (!strcmp(tkn, "-"))
                {
                /* Dealing with a negative number. We will multiply the next tkn, which
                   had better be a number, by -1. */
                isNegative = YES;
                }
            else
                {
                /* We have a number, we hope. */
                if (tkn[0] == matchId)
                    {
                    /* If the token is a matchchar, then things are simple. */
                    if (taxonCount == 1)
                        {
                        YvyraPrint ("%s   Matching characters cannot be in first taxon\n", spacer);
                        goto errorExit;
                        }
                    charCode = matrix[pos(0,taxaInfo[taxonCount-1].charCount,numChar)];
                    matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount,numChar)] = charCode;
                    }
                else
                    {
                    /* Otherwise, we have a number. Check that it is a valid number first... */
                    if (!IsIn(tkn[0],"0123456789."))
                        {
                        YvyraPrint ("%s   Expecting a number for the continuous character\n", spacer);
                        goto errorExit;
                        }
                    /* ... and then put the character into the matrix. Note that matrix
                       is defined as an integer, but we may have floating precision continuous
                       characters. To get around this, we multiply the value of the character
                       by 1000 before putting it into matrix. We will divide by 1000 later on
                       when/if we use the characters. */
                    sscanf (tkn, "%lf", &charValue);
                    charValue *= 1000.0;
                    if (isNegative == YES)
                        {
                        charValue *= -1.0;
                        isNegative = NO;
                        }
                    /*YvyraPrint ("%d \n", (int)charValue);*/
                    matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = (int)charValue;
                    }
                }
            }
        else
            {
            /* Otherwise, we are dealing with a run-of-the-mill character, and we
               cannot expect the entire token to contain only a single character. We
               must, therefore, go through the token character-by-character. */
            i = 0;
            while (tkn[i] != '\0')
                {
                /*YvyraPrint ("%c", tkn[i]);*/
                if (tkn[i] == matchId)
                    {
                    if (taxonCount == 1)
                        {
                        YvyraPrint ("%s   Matching characters cannot be in first taxon\n", spacer);
                        goto errorExit;
                        }
                    charCode = matrix[pos(0,taxaInfo[taxonCount-1].charCount,numChar)];
                    matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
                    }
                else
                    {
                    if ((tkn[i] == ')' && isInAmbig == YES) || (tkn[i] == '}' && isInPoly == YES))
                        {
                        isInAmbig = isInPoly = NO;
                        charCode = theAmbigChar;
                        j = CharacterNumber (charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType);
                        if (j > charInfo[taxaInfo[taxonCount-1].charCount].numStates)
                            charInfo[taxaInfo[taxonCount-1].charCount].numStates = j;
                        matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
                        theAmbigChar = 0;
                        }
                    else if ((tkn[i] == '(' && isInAmbig == YES) || (tkn[i] == '{' && isInPoly == YES))
                        {
                        if (isInAmbig == YES)
                            YvyraPrint ("%s   Found an inappropriate \"(\"\n", spacer);
                        else
                            YvyraPrint ("%s   Found an inappropriate \"{\"\n", spacer);
                        goto errorExit;
                        }
                    else if (isInAmbig == YES || isInPoly == YES)
                        {
                        if (tkn[i] == ',')
                            expecting |= Expecting (COMMA);
                        else 
                            {
                            if (CharacterCode(tkn[i], &charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType) == ERROR)
                                goto errorExit;
                            if (charCode == MISSING || charCode == GAP)
                                goto errorExit;
                            theAmbigChar |= charCode;
                            expecting ^= Expecting (COMMA);
                            }
                        }
                    else if (tkn[i] == '{' && isInPoly == NO && isInAmbig == NO)
                        {
                        isInPoly = YES;
                        matrixHasPoly = YES;
                        theAmbigChar = 0;
                        }   
                    else if (tkn[i] == '(' && isInPoly == NO && isInAmbig == NO)
                        {
                        isInAmbig = YES;
                        theAmbigChar = 0;
                        }
                    else if (tkn[i] == '(' && isInPoly == NO && isInAmbig == NO)
                        {
                        isInAmbig = YES;
                        theAmbigChar = 0;
                        }
                    else
                        {
                        if (CharacterCode(tkn[i], &charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType) == ERROR)
                            {
                            YvyraPrint ("%s   Error while reading character position %d (charCode %d)\n", spacer, taxaInfo[taxonCount-1].charCount+1, charCode);
                            goto errorExit;
                            }
                        if (charCode != MISSING && charCode != GAP)
                            {
                            j = CharacterNumber (charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType);
                            if (j > charInfo[taxaInfo[taxonCount-1].charCount].numStates)
                                charInfo[taxaInfo[taxonCount-1].charCount].numStates = j;
                            }
                        matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
                        }
                    }
                i++; 
                }
            }
        }

    return (NO_ERROR);
    errorExit:
        numTaxa=taxonCount;
        FreeMatrix();
        return (ERROR);
}


int DoNexusParm (char *parmName, char *tkn)
{
    if (!strcmp(parmName, "NEXUS"))
        {
        YvyraPrint ("%s   Expecting NEXUS formatted file\n", spacer);
        expecting = Expecting(COMMAND);
        }
    else
        {
        YvyraPrint ("%s   Found %s\n", spacer, tkn);
        return (ERROR);
        }
    
    return (NO_ERROR);
}


int DoOutgroup (void)
{
    YvyraPrint ("%s   Setting outgroup to taxon \"%s\"\n", spacer, taxaNames[outGroupNum]);
    return (NO_ERROR);
}


int DoOutgroupParm (char *parmName, char *tkn)
{
    int     index, tempInt;

    if (expecting == Expecting(ALPHA))
        {
        if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
            {
            YvyraPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
            return (ERROR);
            }
        outGroupNum = index;
        
        expecting = Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
            {
            /* OK, as we expect, the taxon is not a digit. So, now we assume that
               the user is assigning the outgroup by its number */
            sscanf (tkn, "%d", &tempInt);
            if (tempInt < 1 || tempInt > numTaxa)
                {
                YvyraPrint ("%s   Taxon number %d is out of range\n", spacer, tempInt);
                return (ERROR);
                }
            outGroupNum = tempInt - 1;
            }
        else
            {
            outGroupNum = index;
            }
        
        expecting = Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoPairs (void)
{
    YvyraPrint ("\n");
    YvyraPrint ("%s   Successfully defined character pairings\n", spacer);

    defPairs = YES;
    foundFirst = NO;
    
    return (NO_ERROR);
}


int DoPairsParm (char *parmName, char *tkn)
{
    int     i, tempInt;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can define pairs of characters\n", spacer);
        return (ERROR);
        }
    
    if (defPairs == YES)
        {
        YvyraPrint ("%s   Character pairs have been previously defined \n", spacer);
        YvyraPrint ("%s   Now overwriting old pairings\n", spacer);
        for (i=0; i<numChar; i++)
            charInfo[i].pairsId = 0;
        defPairs = NO;
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        pairId = 1;
        firstPair = YES;
        foundFirst = YES;
        YvyraPrint ("%s   Defining character pairings:\n\n", spacer);
        YvyraPrint ("%s      Pair --  First Second \n", spacer);
        }

    if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            YvyraPrint ("\n");
            YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            for (i=0; i<numChar; i++)
                charInfo[i].pairsId = 0;
            return (ERROR);
            }
        tempInt--;
        
        if (charInfo[tempInt].pairsId != 0)
            {
            YvyraPrint ("\n");
            YvyraPrint ("%s   Character number %d has already been included in a pairing\n", spacer, tempInt+1);
            for (i=0; i<numChar; i++)
                charInfo[i].pairsId = 0;
            return (ERROR);
            }
        if (charInfo[tempInt].charType != DNA && charInfo[tempInt].charType != RNA)
            {
            YvyraPrint ("\n");
            YvyraPrint ("%s  Pairings may only include nucleotide data\n", spacer);
            if (charInfo[tempInt].charType == PROTEIN)
                YvyraPrint ("%s  Character %d is an amino acid character\n", spacer, tempInt+1);
            else if (charInfo[tempInt].charType == RESTRICTION)
                YvyraPrint ("%s  Character %d is a restriction site character\n", spacer, tempInt+1);
            else if (charInfo[tempInt].charType == STANDARD)
                YvyraPrint ("%s  Character %d is a \"standard\" character\n", spacer, tempInt+1);
            else if (charInfo[tempInt].charType == CONTINUOUS)
                YvyraPrint ("%s  Character %d is a continuously varying character\n", spacer, tempInt+1);
            for (i=0; i<numChar; i++)
                charInfo[i].pairsId = 0;
            return (ERROR);
            }
            
        charInfo[tempInt].pairsId = pairId;
        
        if (firstPair == YES)
            {
            YvyraPrint ("%s      %4d --  %5d  ", spacer, pairId, tempInt+1);
            expecting  = Expecting(COLON);
            firstPair = NO;
            }
        else
            {
            YvyraPrint ("%5d\n", tempInt+1);
            expecting  = (Expecting(COMMA) | Expecting(SEMICOLON));
            firstPair = YES;
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        pairId++;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting = Expecting(NUMBER);
        }
    else
        {
        for (i=0; i<numChar; i++)
            charInfo[i].pairsId = 0;
        return (ERROR);
        }

    return (NO_ERROR);
}


int DoPartition (void)
{
    int     i, *partTypes;
        
    /* add set to tempSet */
    if (fromI >= 0)
        if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
            return (ERROR);

    /* check that all characters are included */
    for (i=0; i<numChar; i++)
        {
        /* YvyraPrint ("%4d %4d \n", i, tempSet[i]); */
        if (tempSet[i] == 0)
            {
            YvyraPrint ("%s   Character %d not included in partition\n", spacer, i+1);
            return (ERROR);
            }
        }

            
    /* check how many partitions were found against how many were expected */
    if (whichPartition != numDivisions - 1)
        {
        YvyraPrint ("%s   Did not find correct number of partitions (expecting %d, found %d)\n", spacer, numDivisions, whichPartition + 1);
        return (ERROR);
        }

    partTypes = (int *) SafeCalloc (numDivisions, sizeof(int));
    if (!partTypes)
        return ERROR;
    
    /* make certain that the partition labels go from 1 - numDivisions, inclusive */
    for (i=0; i<numChar; i++)
        partTypes[tempSet[i] - 1] = -1; //partTypes is temporary used here not as an indicator of partition type 
    for (i=0; i<numDivisions; i++)
        {
        if (partTypes[i] == 0)
            {
            YvyraPrint ("%s   Could not find a single character for division %d\n", spacer, i+1);
            return (ERROR);
            }
        }

    /* check if partition overruns data types */
    for (i=0; i<numChar; i++)
        {
        if (partTypes[ tempSet[i]-1 ] == -1)
            partTypes[ tempSet[i]-1 ] = charInfo[i].charType;
        else
            {
            if (partTypes[ tempSet[i]-1 ] != charInfo[i].charType)
                {
                YvyraPrint ("%s   There are two different data types for partition division %d\n", spacer, tempSet[i]);
                free (partTypes);
                return (ERROR);
                }
            }
        }
    free (partTypes);

    /* add name to list of valid partitions */
    if (AddString (&partitionNames, numDefinedPartitions, tempSetName) == ERROR)
        {
        YvyraPrint ("%s   Problem adding partition %s to list\n", spacer, tempSetName);
        return (ERROR);
        }
        
    /* add new partition */
    for (i=0; i<numChar; i++) {
        partitionId[i] = (int *) SafeRealloc ((void *)(partitionId[i]), ((size_t)numDefinedPartitions + 1) * sizeof(int));
        if (!partitionId[i])
            return ERROR;
    }

    /* set new partition */
    for (i=0; i<numChar; i++)
        partitionId[i][numDefinedPartitions] = tempSet[i];

    /* increment number of defined partitions */
    numDefinedPartitions++;
    
    return (NO_ERROR);
}


int DoPartitionParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before partitions can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        /* set Partition () ******************************************************************/
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of partition name */
            if (strlen(tkn) > 99)
                {
                YvyraPrint ("%s   Partition name is too long. Max 100 characters\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a partition */
            if (numDefinedPartitions > 1)
                {
                if (CheckString (partitionNames, numDefinedPartitions, tkn, &index) == ERROR)
                    {
                    /* if the partition name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    YvyraPrint ("%s   Partition name '%s' has been used previously\n", spacer, tkn);
                    return (ERROR);
                    }
                }
                
            /* add the name temporarily to tempSetName */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numChar; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            whichPartition = 0;
            foundFirst = NO;
            numDivisions = 0;
            YvyraPrint ("%s   Defining partition called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* We are defining a partition in terms of a character set (called tkn, here). We should be able
           to find tkn in the list of character set names. If we cannot, then we have a problem and
           return an error. */
        if (numCharSets < 1)
            {
            YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
            {
            YvyraPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        /* add characters from charset tkn to new tempSet */
        for (i=0; i<numChar; i++)
            {
            if (IsBitSet (i, charSet[index]) == YES)
                tempSet[i] = whichPartition + 1;
            }
        fromI = toJ = everyK = -1;

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(COMMA);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundFirst == NO)
            {
            sscanf (tkn, "%d", &tempInt);
            numDivisions = tempInt;
            expecting  = Expecting(COLON);
            foundFirst = YES;
            }
        else
            {
            if (strlen(tkn) == 1 && tkn[0] == '.')
                tempInt = numChar;
            else
                sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numChar)
                {
                YvyraPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
                return (ERROR);
                }
            tempInt--;
            if (foundDash == YES)
                {
                if (fromI >= 0)
                    toJ = tempInt;
                else
                    {
                    YvyraPrint ("%s   Improperly formatted partition\n", spacer);
                    return (ERROR);
                    }
                foundDash = NO;
                }
            else if (foundSlash == YES)
                {
                tempInt++;
                if (tempInt <= 1)
                    {
                    YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                    return (ERROR);
                    }
                if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                    everyK = tempInt;
                else
                    {
                    YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                    return (ERROR);
                    }
                foundSlash = NO;
                }
            else
                {
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    }
                else if (fromI < 0 && toJ < 0)
                    {
                    fromI = tempInt;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else
                    {
                    YvyraPrint ("%s   Improperly formatted charset\n", spacer);
                        {
                        return (ERROR);
                        }
                    }
                }
            
            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(SEMICOLON);
            expecting |= Expecting(DASH);
            expecting |= Expecting(BACKSLASH);
            expecting |= Expecting(COMMA);
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        /* add set to tempSet */
        if (fromI >= 0)
            if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                return (ERROR);

        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;
        whichPartition++;
        if (whichPartition > numDivisions)
            {
            YvyraPrint ("%s   Too many partitions of the data (expecting %d)\n", spacer, numDivisions);
            return (ERROR);
            }
        expecting  = Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting  = Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoRestore (void)
{
    int         i, alreadyDone;

    YvyraPrint ("%s   Restore taxa\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with excludedTaxa */
    alreadyDone = NO;
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] == 1)
            {
            if (taxaInfo[i].isDeleted == NO && alreadyDone == NO)
                {
                YvyraPrint ("%s   Some taxa already included\n", spacer);
                alreadyDone = YES;
                }
            taxaInfo[i].isDeleted = NO;
            }
        }

    SetLocalTaxa();
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    /* show tempSet (for debugging) */
#   if 0
    for (i=0; i<numTaxa; i++)
        YvyraPrint ("%4d  %4d\n", i+1, tempSet[i]);
#   endif
        
    return (NO_ERROR);
}


int DoRestoreParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before you can restore taxa\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = NO;
        for (i=0; i<numTaxa; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 1;
            }
        else
            {
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                /* we are using a pre-defined taxa set */
                if (numTaxaSets < 1)
                    {
                    YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                /* add taxa from taxset tkn to new tempSet */
                for (i=0; i<numTaxa; i++)
                    {
                    if (IsBitSet (i, taxaSet[index]) == YES)
                        tempSet[i] = 1;
                    }
                }
            else
                {
                /* we found the taxon name */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                tempSet[index] = 1;
                }
            fromI = toJ = everyK = -1;
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && !strcmp(tkn, "."))
            {
            tempInt = numTaxa;
            }
        else
            {
            sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numTaxa)
                {
                YvyraPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                return (ERROR);
                }
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted restore set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted restore set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
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


int DoSet (void)
{
    return (NO_ERROR);
}


int DoSetParm (char *parmName, char *tkn)
{
    int         index;
    char        tempStr[100];
    int         tempI;

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Autoclose (autoClose) **********************************************************/
        if (!strcmp(parmName, "Autoclose"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        autoClose = YES;
                    else
                        autoClose = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for autoclose\n", spacer);
                    return (ERROR);
                    }
                if (autoClose == YES)
                    YvyraPrint ("%s   Setting autoclose to yes\n", spacer);
                else
                    YvyraPrint ("%s   Setting autoclose to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nowarnings (noWarn) **********************************************************/
        else if (!strcmp(parmName, "Nowarnings"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        noWarn = YES;
                    else
                        noWarn = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for nowarnings\n", spacer);
                    return (ERROR);
                    }
                if (noWarn == YES)
                    YvyraPrint ("%s   Setting nowarnings to yes\n", spacer);
                else
                    YvyraPrint ("%s   Setting nowarnings to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Quitonerror (quitOnError) **************************************************/
        else if (!strcmp(parmName, "Quitonerror"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        quitOnError = YES;
                    else
                        quitOnError = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for quitonerror\n", spacer);
                    return (ERROR);
                    }
                if (quitOnError == YES)
                    YvyraPrint ("%s   Setting quitonerror to yes\n", spacer);
                else
                    YvyraPrint ("%s   Setting quitonerror to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Autoreplace (autoOverwrite) **************************************************/
        else if (!strcmp(parmName, "Autoreplace"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        {
                        autoOverwrite = YES;
                        YvyraPrint ("%s   Setting autoreplace to yes\n", spacer);
                        }
                    else
                        {
                        autoOverwrite = NO;
                        YvyraPrint ("%s   Setting autoreplace to no\n", spacer);
                        }
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for autoreplace\n", spacer);
                    return (ERROR);
                    }                   
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }           
        /* set Scientific (scientific) *********************************************/
        else if (!strcmp(parmName, "Scientific"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        scientific = YES;
                    else
                        scientific = NO;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for Scientific\n", spacer);
                    return (ERROR);
                    }
                if (scientific == YES)
                    YvyraPrint ("%s   Setting Scientific to Yes\n", spacer);
                else
                    YvyraPrint ("%s   Setting Scientific to No\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Userlevel (userLevel) **********************************************************/
        else if (!strcmp(parmName, "Userlevel"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Standard"))
                        userLevel = STANDARD_USER;
                    else if (!strcmp (tempStr,"Developer"))
                        userLevel = DEVELOPER;
                    }
                else
                    {
                    YvyraPrint ("%s   Invalid argument for userlevel\n", spacer);
                    return (ERROR);
                    }
                YvyraPrint ("%s   Setting userlevel to %s\n", spacer, tempStr);
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Npthreads (number of pthreads) ****************************************************/
        else if (!strcmp(parmName, "Npthreads"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                nPThreads = tempI;
                YvyraPrint ("%s   Setting Npthreads to %d\n", spacer, nPThreads);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }
        /* set Precision (number of decimals) ****************************************************/
        else if (!strcmp(parmName, "Precision"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 3 || tempI > 15)
                    {
                    YvyraPrint ("%s   Precision must be in the range 3 to 15\n", spacer);
                    return ERROR;
                    }
                precision = tempI;
                YvyraPrint ("%s   Setting Precision to %d\n", spacer, precision);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }
        /* set Partition (partitionNum) *******************************************************/
        else if (!strcmp(parmName, "Partition"))
            {
            if (defMatrix == NO)
                {
                YvyraPrint ("%s   A character matrix must be defined first\n", spacer);
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA) | Expecting(NUMBER);
            else if (expecting == Expecting(ALPHA))
                {
                /* first check to see if name is there */
                if (CheckString (partitionNames, numDefinedPartitions, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find \"%s\" as a defined partition\n", spacer, tkn);
                    return (ERROR);
                    }
                if (SetPartition (index) == ERROR)
                    return ERROR;
                if (numCurrentDivisions == 1)
                    YvyraPrint ("%s   Setting %s as the partition (does not divide up characters).\n", spacer, tkn); 
                else
                    YvyraPrint ("%s   Setting %s as the partition, dividing characters into %d parts.\n", spacer, tkn, numCurrentDivisions); 
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &index);
                if (index > numDefinedPartitions) 
                    {
                    YvyraPrint ("%s   Partition number %d is not a valid partition. Only %d partitions\n", spacer, index, numDefinedPartitions);
                    YvyraPrint ("%s   have been defined.\n", spacer);
                    return (ERROR);
                    }
                if (index < 1)
                    {
                    YvyraPrint ("%s   Partition number %d is not a valid partition. Must be between 1 and %d.\n", spacer, index+1, numDefinedPartitions);
                    return (ERROR);
                    }
                if (SetPartition (index) == ERROR)
                    return ERROR;
                if (numCurrentDivisions == 1)
                    YvyraPrint ("%s   Setting %s as the partition (does not divide up characters).\n", spacer, partitionNames[index]);
                else
                    YvyraPrint ("%s   Setting %s as the partition, dividing characters into %d parts.\n", spacer, partitionNames[index], numCurrentDivisions);
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Speciespartition (speciespartitionNum) *******************************************************/
        else if (!strcmp(parmName, "Speciespartition"))
            {
            if (defTaxa == NO)
                {
                YvyraPrint ("%s   A taxaset must be defined first\n", spacer);
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA) | Expecting(NUMBER);
            else if (expecting == Expecting(ALPHA))
                {
                /* first check to see if name is there */
                if (CheckString (speciespartitionNames, numDefinedSpeciespartitions, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find \"%s\" as a defined speciespartition\n", spacer, tkn);
                    return (ERROR);
                    }
                if (SetSpeciespartition (index) == ERROR)
                    return ERROR;
                YvyraPrint ("%s   Setting %s as the speciespartition, dividing taxa into %d species.\n", spacer, tkn, numSpecies);
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &index);
                if (index > numDefinedSpeciespartitions) 
                    {
                    YvyraPrint ("%s   Speciespartition number %d is not valid. Only %d speciespartitions\n", spacer, index, numDefinedSpeciespartitions);
                    YvyraPrint ("%s   have been defined.\n", spacer);
                    return (ERROR);
                    }
                if (index < 1)
                    {
                    YvyraPrint ("%s   Speciespartition number %d is not valid. Must be between 1 and %d.\n", spacer, index, numDefinedSpeciespartitions);
                    return (ERROR);
                    }
                if (SetSpeciespartition (index-1) == ERROR)
                    return ERROR;
                YvyraPrint ("%s   Setting %s as the speciespartition, dividing taxa into %d species.\n", spacer, speciespartitionNames[index-1], numSpecies);
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Seed (global variable globalSeed) ****************************************************/
        else if (!strcmp(parmName, "Seed"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI == 0 || tempI == 2147483647)
                    {
                    YvyraPrint ("%s   Error: Seed can be any natural number except 0 and 2147483647\n", spacer);
                    return (ERROR);
                    }
                globalSeed = tempI;
                YvyraPrint ("%s   Setting seed to %ld\n", spacer, globalSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }
        /* set Swapseed (global variable swapSeed) ***************************************************************/
        else if (!strcmp(parmName, "Swapseed"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI == 0 || tempI == 2147483647)
                    {
                    YvyraPrint ("%s   Error: Swapseed can be any natural number except 0 and 2147483647\n", spacer);
                    return (ERROR);
                    }
                swapSeed = tempI;
                YvyraPrint ("%s   Setting swapseed to %ld\n", spacer, swapSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Dir (global variable workingDir) ***************************************************************/
        else if (!strcmp(parmName, "Dir"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (strlen(tkn)>99)
                    {
                    YvyraPrint ("%s   Maximum allowed length of working directory name is 99 characters. The given name:\n", spacer);
                    YvyraPrint ("%s      '%s'\n", spacer,tkn);
                    YvyraPrint ("%s   has %d characters.\n", spacer,strlen(tkn));
                    return (ERROR);
                    }
                strcpy (workingDir, tkn);
#   if defined (WIN_VERSION)
                /* Reformat to Windows with trailing '\' */
                for (index=0; index<(int)strlen(workingDir); index++)
                    {
                    if (workingDir[index] == '/')
                        workingDir[index] = '\\';
                    }
                if (strlen(workingDir) > 0 && workingDir[strlen(workingDir)-1] != '\\')
                    strcat(workingDir,"\\");
#   else
                /* Reformat to Unix with trailing '/' */
                for (index=0; index<(int)strlen(workingDir); index++)
                    {
                    if (workingDir[index] == '\\')
                        workingDir[index] = '/';
                    }
                if (strlen(workingDir) > 0 && workingDir[strlen(workingDir)-1] != '/')
                    strcat(workingDir,"/");
#   endif
                YvyraPrint ("%s   Setting working directory to \"%s\"\n", spacer, workingDir);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Usebeagle (global variable BEAGLE usage) ***************************************************************/    
        else if (!strcmp(parmName, "Usebeagle"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Beagle resource number (global variable BEAGLE flag) ****************************************/
        else if (!strcmp(parmName, "Beagleresource"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting =  Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Beagle resources requirements (global variable BEAGLE flag) ****************************************/
        else if (!strcmp(parmName, "Beagledevice"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beagleprecision"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beagleopenmp"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beaglefreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }        
        else if (!strcmp(parmName, "Beaglesse"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beaglethreads"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beaglethreadcount"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beaglefloattips"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
#if 0
        else if (!strcmp(parmName, "Beaglevec"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
#endif
        else if (!strcmp(parmName, "Beaglescaling"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
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


int DoShowMatrix (void)
{
    int         i, j, nameLen, start, finish, ct, longestName;
    char        tempStr[100], stride;
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }
            
    longestName = 0;
    for (i=0; i<numTaxa; i++)
        {
        nameLen = (int) strlen(taxaNames[i]);
        if (nameLen > longestName)
            longestName = nameLen;
        }
            
    stride = 50;
    start = finish = 0;
    do
        {
        finish += stride;
        if (finish > numChar)
            finish = numChar;

        YvyraPrint ("%s   ", spacer);
        for (j=0; j<longestName; j++)
            YvyraPrint (" ");
        YvyraPrint ("  ");
        YvyraPrint ("%d\n", start+1);

        for (i=0; i<numTaxa; i++)
            {
            strcpy (tempStr, taxaNames[i]);
            nameLen = (int) strlen(tempStr);
            
            YvyraPrint ("%s   ", spacer);
            if (nameLen >= longestName)
                {
                for (j=0; j<longestName; j++)
                    YvyraPrint ("%c", tempStr[j]);
                }
            else
                {
                YvyraPrint ("%s", tempStr);
                for (j=0; j<longestName-nameLen; j++)
                    YvyraPrint (" ");
                }
            YvyraPrint ("  ");

            for (j=start; j<finish; j++)
                {
                ct = charInfo[j].charType;
                if (ct == STANDARD)
                    YvyraPrint ("%c", WhichStand(matrix[pos(i,j,numChar)]));
                else if (ct == CONTINUOUS)
                    {
                    if (WhichCont(matrix[pos(i,j,numChar)]) < 0.0)
                        YvyraPrint (" %2.2lf", WhichCont(matrix[pos(i,j,numChar)]));
                    else
                        YvyraPrint ("  %2.2lf", WhichCont(matrix[pos(i,j,numChar)]));
                    }
                else
                    {
                    YvyraPrint ("%s   Unknown data type\n", spacer);
                    return (ERROR);
                    }
                
                }
            YvyraPrint ("\n");
            }
        YvyraPrint ("\n");
        start = finish;
        } while (finish != numChar);

    return (NO_ERROR);
}


int DoShowUserTrees (void)
{
    int         i;

    if (numUserTrees == 0)
        {
        YvyraPrint ("%s   No user trees have been defined\n", spacer);
        }
    else
        {
        for (i=0; i<numUserTrees; i++)
            {
            YvyraPrint ("\n   Tree #%d -- '%s':\n\n", i+1, userTree[i]->name);
            ShowConTree (stdout, userTree[i], 70, NO);
            YvyraPrint ("\n");
            }
        }

    return (NO_ERROR);
}


int DoTaxlabels (void)
{
    isTaxsetDef = YES;

    /* add default speciespartition name to list of valid speciespartitions */
    if (AddString (&speciespartitionNames, 0, "Default") == ERROR)
        {
        YvyraPrint ("%s   Problem adding Default speciespartition to list\n", spacer);
        return (ERROR);
        }

    /* add default species name set */
    AddNameSet(&speciesNameSets, 0, taxaNames, numTaxa);

    /* set number of defined speciespartitions to 1 */
    numDefinedSpeciespartitions = 1;
        
    return (NO_ERROR);
}


int DoTaxlabelsParm (char *parmName, char *tkn)
{
    int         index;

    if (inTaxaBlock == NO)
        {
        YvyraPrint ("%s   You must be in a taxa block to read a taxlabels command\n", spacer);
        return (ERROR);
        }

    if (defTaxa == NO)
        {
        YvyraPrint ("%s   The number of taxa must be given before a set of taxon labels can be read\n", spacer);
        return ERROR;
        }

    if (isTaxsetDef == YES)
        {
        YvyraPrint ("%s   A set of taxon labels has already been defined\n", spacer);
        if (defMatrix == NO)
            if (WantTo ("Do you want to delete the current set of taxon labels") == NO)
                return (SKIP_COMMAND);
            else
                FreeTaxa();
        else
            if (WantTo ("Do you want to delete the current character matrix") == NO)
                return (SKIP_COMMAND);
            else
                FreeMatrix();
        }

    if (expecting == Expecting(ALPHA) ||
        expecting == Expecting(NUMBER))
        {
        if (CheckString (taxaNames, numNamedTaxa, tkn, &index) == ERROR)
            {
            if (strlen(tkn)>99)
                {
                YvyraPrint ("%s   Taxon name %s is too long. Maximum 99 characters is allowed.\n", spacer, tkn);
                return (ERROR);
                }
            if (AddString (&taxaNames, numNamedTaxa, tkn) == ERROR)
                {
                YvyraPrint ("%s   Problem adding label %s to list of taxon labels\n", spacer, tkn);
                return (ERROR);
                }
            numNamedTaxa++;
            }
        else
            {
            YvyraPrint ("%s   Taxon label '%s' is included twice in list of taxon labels\n", spacer, tkn);
            return (ERROR);
            }
        if (numNamedTaxa < numTaxa)
            {
            expecting = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            }
        else
            expecting |= Expecting(SEMICOLON);
        }

    return (NO_ERROR);
}


int DoSpeciespartition (void)
{
    int     i, *partCount;
        
    /* add set to tempSet */
    if (fromI >= 0)
        if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
            {
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }

    /* set numDivisions; not set while reading the speciespartition */
    numDivisions = whichPartition + 1;
    
    /* check that all species are included */
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] == 0)
            {
            YvyraPrint ("%s   Tip %d not included in speciespartition\n", spacer, i+1);
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        /*YvyraPrint ("%4d %4d \n", i, tempSet[i]);*/
        }

    partCount = (int *) SafeCalloc (numDivisions, sizeof(int));
    if (!partCount)
        {
        for (i=0; i<numDivisions; i++)
            free(tempNames[i]);
        free (tempNames);
        tempNames = NULL;
        return ERROR;
        }

    /* make certain that the partition labels go from 1 - numTaxa, inclusive */
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] < 1 || tempSet[i] > numTaxa)
            {
            YvyraPrint ("%s   Speciespartition index for tip %d out of bound (%d)\n", spacer, i+1, tempSet[i]);
            free (partCount);
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        partCount[tempSet[i] - 1]++;
        }
    for (i=0; i<numDivisions; i++)
        {
        if (partCount[i] == 0)
            {
            YvyraPrint ("%s   Could not find a single tip for species %d\n", spacer, i+1);
            free (partCount);
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        }
    free (partCount);

    /* add name to list of valid partitions */
    if (AddString (&speciespartitionNames, numDefinedSpeciespartitions, tempSetName) == ERROR)
        {
        YvyraPrint ("%s   Problem adding speciespartition %s to list\n", spacer, tempSetName);
        for (i=0; i<numDivisions; i++)
            free(tempNames[i]);
        free (tempNames);
        tempNames = NULL;
        return (ERROR);
        }

    /* add new partition */
    for (i=0; i<numTaxa; i++)
        {
        speciespartitionId[i] = (int *) SafeRealloc ((void *)(speciespartitionId[i]), ((size_t)numDefinedSpeciespartitions + 1) * sizeof(int));
        if (!speciespartitionId[i])
            {
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return ERROR;
            }
        }

    /* set new partition */
    for (i=0; i<numTaxa; i++)
        speciespartitionId[i][numDefinedSpeciespartitions] = tempSet[i];

    /* add new set of species names */
    AddNameSet(&speciesNameSets, numDefinedSpeciespartitions, tempNames, numDivisions);

    /* free species names */
    for (i=0; i<numDivisions; i++)
        free(tempNames[i]);
    free (tempNames);
    tempNames = NULL;

    /* increment number of defined partitions */
    numDefinedSpeciespartitions++;
    
    return (NO_ERROR);
}


int DoSpeciespartitionParm (char *parmName, char *tkn)
{
    int             i, index, tempInt;
    
    if (defTaxa == NO || numTaxa == 0)
        {
        YvyraPrint ("%s   A matrix or taxaset must be specified before partitions can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        /* set Speciespartition name ******************************************************************/
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of partition name */
            if (strlen(tkn) > 99)
                {
                YvyraPrint ("%s   Partition name is too long. Max 100 characters\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a partition */
            if (numDefinedSpeciespartitions > 0)
                {
                if (CheckString (speciespartitionNames, numDefinedSpeciespartitions, tkn, &index) == ERROR)
                    {
                    /* if the partition name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    YvyraPrint ("%s   Speciespartition name '%s' has been used previously\n", spacer, tkn);
                    return (ERROR);
                    }
                }
                
            /* add the name temporarily to tempSetName */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 0;
    
            /* make sure tempNames is NULL */
            assert (tempNames == NULL);

            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            whichPartition = 0;
            foundFirst = NO;
            numDivisions = 0;
            YvyraPrint ("%s   Defining speciespartition called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(ALPHA);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (foundFirst == NO)
            {
            AddString(&tempNames, whichPartition, tkn);
            foundFirst = YES;
            expecting = Expecting(COLON);
            }
        else
            {
            /* We are defining a species partition in terms of a tip name (called tkn, here). We should be able
               to find tkn in the list of tip names. If we cannot, then we have a problem and
               return an error. */
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                YvyraPrint ("%s   Could not find a tip called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add index of the tip named tkn to new tempSet */
            tempSet[index] = whichPartition + 1;
            fromI = toJ = everyK = -1;

            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(SEMICOLON);
            expecting |= Expecting(COMMA);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numTaxa;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numTaxa)
            {
            YvyraPrint ("%s   Tip number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
            for (i=0; i<whichPartition; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted speciespartition\n", spacer);
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                YvyraPrint ("%s   Improperly formatted speciespartition\n", spacer);
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted speciespartition\n", spacer);
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                    {
                    for (i=0; i<whichPartition; i++)
                        free(tempNames[i]);
                    free (tempNames);
                    tempNames = NULL;
                    return (ERROR);
                    }
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted speciespartition\n", spacer);
                    {
                    for (i=0; i<whichPartition; i++)
                        free(tempNames[i]);
                    free (tempNames);
                    tempNames = NULL;
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        expecting |= Expecting(COMMA);
        }
    else if (expecting == Expecting(COMMA))
        {
        /* add set to tempSet */
        if (fromI >= 0)
            if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                {
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }

        fromI = toJ = everyK = -1;
        foundDash = foundSlash = foundFirst = NO;
        whichPartition++;
        if (whichPartition > numTaxa)
            {
            YvyraPrint ("%s   Too many speciespartitions (expecting maximum %d speciespartitions)\n", spacer, numTaxa);
            for (i=0; i<whichPartition; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        expecting  = Expecting(ALPHA);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting  = Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        {
        for (i=0; i<whichPartition; i++)
            free(tempNames[i]);
        free (tempNames);
        tempNames = NULL;
        return (ERROR);
        }

    return (NO_ERROR);
}


int DoTaxaset (void)
{
    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* add name to taxaSetNames */
    if (AddString (&taxaSetNames, numTaxaSets, tempSetName) == ERROR)
        {
        YvyraPrint ("%s   Problem adding taxset %s to list\n", spacer, tempSetName);
        return (ERROR);
        }

    /* merge tempSet with taxaSet */
    AddBitfield (&taxaSet, numTaxaSets, tempSet, numTaxa);
    
    /* increment number of char sets */
    numTaxaSets++;
    
    /* show taxset (for debugging) */
#   if 0
    for (i=0; i<numTaxa; i++)
        YvyraPrint ("%4d  %4d\n", i+1, tempSet[i]);
#   endif

    return (NO_ERROR);
}


int DoTaxasetParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A matrix must be specified before taxsets can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of taxset name */
            if (strlen(tkn) > 99)
                {
                YvyraPrint ("%s   Taxset name is too long\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a taxset */
            if (numTaxaSets > 0)
                {
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    /* if the taxset name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    YvyraPrint ("%s   Taxset name has been used previously\n", spacer);
                    return (ERROR);
                    }
                }
            else if (numTaxaSets > 30)
                {
                YvyraPrint ("%s   You cannot define more than 30 taxsets\n", spacer);
                return (ERROR);
                }
                
            /* add the name to the taxa set */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            YvyraPrint ("%s   Defining taxset called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* We are defining a taxon set in terms of another (called tkn, here) or we are referring to
           the taxon name. We should be able to find tkn in the list of character set names or in the list
           of taxon names. If we cannot, then we have a problem and return an error. */
        if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
            {
            if (numTaxaSets < 1)
                {
                YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                {
                YvyraPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add taxa from taxset tkn to new tempSet */
            for (i=0; i<numTaxa; i++)
                {
                if (IsBitSet (i, taxaSet[index]) == YES)
                    tempSet[i] = 1;
                }
            }
        else
            {
            tempSet[index] = 1;
            }
        fromI = toJ = everyK = -1;

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && !strcmp(tkn, "."))
            {
            tempInt = numTaxa;
            }
        else
            {
            sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numTaxa)
                {
                YvyraPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                return (ERROR);
                }
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted taxset\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                YvyraPrint ("%s   Improperly formatted taxset\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                YvyraPrint ("%s   Improperly formatted taxset\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                YvyraPrint ("%s   Improperly formatted taxset\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoTaxaStat (void)
{
    int         i, j, maxLen, nameLen, nIncludedTaxa;
    char        tempName[100];
    
    if (defMatrix == NO)
        {
        YvyraPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }
        
    /* find maximum length of taxon name */
    maxLen = nIncludedTaxa = 0;
    for (i=0; i<numTaxa; i++)
        {
        strcpy (tempName, taxaNames[i]);
        if ((int)strlen(tempName) > maxLen)
            maxLen = (int) strlen(tempName);
        if (taxaInfo[i].isDeleted == NO)
            nIncludedTaxa++;
        }
            
    YvyraPrint ("%s   Showing taxon status:\n\n", spacer);
    if (nIncludedTaxa == numTaxa)
        YvyraPrint ("%s     Number of taxa        = %d (all of which are included)\n", spacer, numTaxa);
    else
        YvyraPrint ("%s     Number of taxa        = %d (of which %d are included)\n", spacer, numTaxa, nIncludedTaxa);
    YvyraPrint ("%s     Number of constraints = %d\n\n", spacer, numDefinedConstraints);
    
    if (numDefinedConstraints > 0)
        {
        for (j=0; j<numDefinedConstraints; j++)
            {
            strcpy (tempName, constraintNames[j]);

            /* for now, ignore the probability */
            if (definedConstraintsType[j] == HARD)
                YvyraPrint ("%s     %2d -- Trees with 'hard' constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
            else if (definedConstraintsType[j] == PARTIAL)
                YvyraPrint ("%s     %2d -- Trees with 'partial' constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
            else
                YvyraPrint ("%s     %2d -- Trees with 'negative' constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
            YvyraPrint ("%s           more probable than those without \n", spacer);
            }
        YvyraPrint ("\n");
        for (j=0; j<maxLen; j++)
            YvyraPrint (" ");
        YvyraPrint ("                             Constraints\n");
        }
    YvyraPrint ("%s     Taxon  ", spacer);
    for (j=0; j<maxLen; j++)
        YvyraPrint (" ");
    YvyraPrint ("   Inclusion");
    YvyraPrint ("   ");
    for (j=0; j<numDefinedConstraints; j++)
        YvyraPrint (" %2d", j+1);
    YvyraPrint ("\n");
    YvyraPrint ("%s   -------", spacer);
    for (j=0; j<maxLen; j++)
        YvyraPrint ("-");
    YvyraPrint ("--------------");
    
    if (numDefinedConstraints > 0)
        {
        YvyraPrint ("----");
        for (j=0; j<numDefinedConstraints; j++)
            YvyraPrint ("---");
        }
    YvyraPrint ("\n");
    for (i=0; i<numTaxa; i++)
        {
        strcpy (tempName, taxaNames[i]);
        nameLen = (int) strlen(tempName);
        
        if (i == outGroupNum)
            YvyraPrint ("%s ->%4d (%s) ", spacer, i+1, tempName);
        else
            YvyraPrint ("%s   %4d (%s) ", spacer, i+1, tempName);
        for (j=0; j<(maxLen-nameLen); j++)
            YvyraPrint (" ");
        YvyraPrint (" -- ");
        
        if (taxaInfo[i].isDeleted == YES)
            YvyraPrint ("Deleted ");
        else
            YvyraPrint ("Included");
            
        YvyraPrint ("    ");
            
        for (j=0; j<numDefinedConstraints; j++)
            {
            if (definedConstraintsType[j] == HARD)
                {
                if (IsBitSet(i, definedConstraint[j]) == NO)
                    YvyraPrint ("  .");
                else
                    YvyraPrint ("  *");
                }
            else if (definedConstraintsType[j] == PARTIAL)
                {
                if (IsBitSet(i, definedConstraint[j]) == YES)
                    YvyraPrint ("  +");
                else if (IsBitSet(i, definedConstraintTwo[j]) == YES)
                    YvyraPrint ("  -");
                else
                    YvyraPrint ("  .");
                }
            else if (definedConstraintsType[j] == NEGATIVE)
                {
                if (IsBitSet(i, definedConstraint[j]) == NO)
                    YvyraPrint ("  .");
                else
                    YvyraPrint ("  #");
                }
            }
        YvyraPrint ("\n");
        }
        
    YvyraPrint ("\n");
    YvyraPrint ("%s   '.' indicate that the taxon is not present in the constraint. \n", spacer);
    YvyraPrint ("%s   '*' indicate that the taxon is present in the 'hard' constraint. \n", spacer);
    YvyraPrint ("%s   '+' indicate that the taxon is present in the first group of 'partial' constraint. \n", spacer);
    YvyraPrint ("%s   '-' indicate that the taxon is present in the second group of 'partial' constraint. \n", spacer);
    YvyraPrint ("%s   '#' indicate that the taxon is present in the 'negative' constraint. \n", spacer);
    YvyraPrint ("%s   Arrow indicates current outgroup. \n", spacer);

    return (NO_ERROR);
}


int DoTranslate (void)
{
    int     i, j;

    if (inTreesBlock == NO)
        {
        YvyraPrint ("%s   You must be in a trees block to read a translate command\n", spacer);
        return (ERROR);
        }
    numTranslates++;    /* number of taxa in translate table */
    isTranslateDef = YES;

    isTranslateDiff = NO;
    if (isTaxsetDef == NO)
        SetTaxaFromTranslateTable();
    else
        {
        for (i=0; i<numTranslates; i++)
            {
            strcpy (token, transFrom[i]);
            if (CheckString (taxaNames, numTaxa, token, &j) == ERROR)
                {
                isTranslateDiff = YES;
                }
            }
        if (numTranslates != numTaxa)
            isTranslateDiff = YES;
        }

    return (NO_ERROR);
}


int DoTranslateParm (char *parmName, char *tkn)
{
    int         index;
    static int  whichTranslate;

    if (inTreesBlock == NO)
        {
        YvyraPrint ("%s   You must be in a trees block to read a translate command\n", spacer);
        return (ERROR);
        }

    if (isTranslateDef == YES)
        {
        YvyraPrint ("%s   A translation has already been defined for this tree block\n", spacer);
        return (ERROR);
        }
        
    if (expecting == Expecting(ALPHA) ||
        expecting == Expecting(NUMBER))
        {
        if (numTaxa == 0)
            {
            YvyraPrint ("%s   Data matrix should be defined before translation table could be set.\n", spacer);
            return (ERROR);
            }
        if (numTranslates == numTaxa)
            {
            YvyraPrint ("%s   Too many entries in translation table. Maximum number of taxon names to translate is %d\n", spacer,numTaxa);
            return (ERROR);
            }
        if (whichTranslate == 0)
            {
            if (CheckString (transTo, numTranslates, tkn, &index) == ERROR)
                {
                if (AddString (&transTo, numTranslates, tkn) == ERROR)
                    {
                    YvyraPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
                    return (ERROR);
                    }
                }
            else
                {
                YvyraPrint ("%s   Already found name (%s) in list\n", spacer, tkn);
                return (ERROR);
                }           
            whichTranslate++;
            expecting = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            }
        else 
            {
            if (CheckString (transFrom, numTranslates, tkn, &index) == ERROR)
                {
                if (AddString (&transFrom, numTranslates, tkn) == ERROR)
                    {
                    YvyraPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
                    return (ERROR);
                    }
                }
            else
                {
                YvyraPrint ("%s   Already found name (%s) in list\n", spacer, tkn);
                return (ERROR);
                }           
            whichTranslate = 0;
            expecting = Expecting(COMMA);
            expecting |= Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        numTranslates++;
        expecting = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }

    return (NO_ERROR);
}


int DoTree (void)
{
    readComment = NO;

    if (inSumtCommand == YES || inComparetreeCommand == YES)
        return (DoSumtTree ());

    return (NO_ERROR);
}


int DoTreeParm (char *parmName, char *tkn)
{
    int                 i, tempInt, index;
    YFlt              tempD;
    char                tempName[100];
    static BitsLong     lastExpecting; /* keep track of what we expected before a comment, in case we want to skip a comment */
    static char         *tempNameString=NULL; /* Contains multiple tokens which form name string of param set*/
    static int          foundAmpersand, foundColon, foundComment, foundE, foundB, foundN, foundFirst,
                        foundCurly, /* is set to YES when we are between two curly braces ONLY while processing CppEvent name */
                        foundClockrate, 
                        foundName, /*is set to YES when param set name token is found and set to NO once full param set name is processed*/
                        eSetIndex, /* is set in the beginning of reading CppEvent for a node/branch to the index of currently processed CppEvent set */
                        bSetIndex, eventIndex, treeIndex, nextIntNodeIndex;
    static PolyNode     *pp, *qq;
    static PolyTree     *t;
    
    /* This function will read in components of a tree description. We expect one of the following formats:
    
          tree <name> = [&R] <newick-description>;
          tree <name> = [&U] <newick-description>;
          tree <name> [&E CppEvents]  = [&R] [&clockrate = 1.23] ((1:0.021[&E CppEvents 2: (0.10 1.11,0.83 3.17)],...
          tree <name> [&B TK02Brlens] = [&R] [&clockrate = 1.23] ((1:0.021[&B TK02Brlens 0.019],...
          tree <name> [&B IgrBrlens]  = [&R] [&clockrate = 1.23] ((1:0.021[&B IgrBrlens 0.019],...
     
       Values will be stored in event sets that go with the tree and that are used to initialize the relaxed clock
       parameters before a run is started. Note that several sets of events can be stored with each tree.
    */

    if (isTaxsetDef == NO)
        {
        YvyraPrint ("%s   Taxon labels must be specified before a tree could be red in\n", spacer);
        return (ERROR);
        }
    if (inTreesBlock == NO)
        {
        YvyraPrint ("%s   You must be in a trees block to read a tree\n", spacer);
        return (ERROR);
        }
    
    if (expecting == Expecting(PARAMETER))
        {
        /* this is the name of the tree */
        if (inSumtCommand==YES || inComparetreeCommand == YES)
            {
            /* we are reading in a tree to sumt or comparetree counters */
            t = sumtParams.tree;
            ResetPolyTree (t);
            }
        else
            {
            /* we are reading in a user tree */
            /* check if the tree exists */
            treeIndex = 0;
            for (i=0; i<numUserTrees; i++)
                if (strcmp(tkn,userTree[i]->name) == 0)
                    break;
            treeIndex = i;
            if (treeIndex < numUserTrees)
                {
                YvyraPrint ("%s   Overwriting tree '%s'.\n", spacer, userTree[treeIndex]);
                FreePolyTree (userTree[treeIndex]);
                }
            if (treeIndex > MAX_NUM_USERTREES)
                {
                YvyraPrint ("%s   yvyra can only store %d user trees.\n", spacer, MAX_NUM_USERTREES);
                return (ERROR);
                }
            if ((userTree[treeIndex] = AllocatePolyTree (numTaxa)) == NULL)
                return (ERROR);
            t = userTree[treeIndex];
            }
        strncpy (t->name, tkn, 99);
        foundColon = foundAmpersand = foundEqual = foundComment = NO;
        foundE = foundB = foundN = foundFirst = foundClockrate = foundName = NO;
        eSetIndex = bSetIndex = eventIndex = 0;
        nextAvailableNode = 0;
        if (isTranslateDef == YES && isTranslateDiff == YES)
            nextIntNodeIndex = numTranslates;
        else
            nextIntNodeIndex = numTaxa;
        pp = &t->nodes[nextAvailableNode++];
        t->root = pp;
        t->isRooted = NO;  /* expect unrooted tree */
        t->isClock = NO;   /* expect nonclock tree */
        t->isCalibrated = NO;  /* expect uncalibrated tree */
        t->isRelaxed = NO;    /* expect strict clock if clock tree */
        t->clockRate = 0.0;     /* expect no clock rate */
        t->popSizeSet = NO;     
        readComment = YES;
        expecting = Expecting(EQUALSIGN) | Expecting(LEFTCOMMENT);
        lastExpecting = expecting;
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        if (foundClockrate == YES)
            expecting = Expecting(NUMBER);
        else
            {
            for (i=0; i<numTaxa; i++)
                tempSet[i] = NO;
            foundEqual = YES;
            expecting = Expecting(LEFTPAR) | Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        if (foundE == YES)
            {
            expecting = Expecting(NUMBER);
            }
        else
            {
            if (nextAvailableNode >= 2*numTaxa)
                {
                YvyraPrint ("%s   Too many nodes on tree '%s'\n", spacer, t->name);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            /* FIXME: t == NULL here (from clang static analyzer) */
            qq = &t->nodes[nextAvailableNode++];
            qq->anc = pp;
            pp->left = qq;
            pp->index = nextIntNodeIndex++;
            pp = qq;
            expecting = Expecting(LEFTPAR);
            expecting |= Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (foundAmpersand == YES)
            {
            if (strcmp(tkn,"E") == 0)
                {
                foundE = YES;
                expecting = Expecting(ALPHA);
                }
            else if (strcmp(tkn,"B") == 0)
                {
                foundB = YES;
                expecting = Expecting(ALPHA);
                }
            else if (strcmp(tkn,"N") == 0)
                {
                foundN = YES;
                expecting = Expecting(ALPHA);
                }
            else if (strcmp(tkn, "R") == 0)
                {
                t->isRooted = YES;
                t->isClock = YES;   /* assume clock if rooted */
                expecting = Expecting(RIGHTCOMMENT);
                }
            else if (strcmp(tkn, "D") == 0)
                {
                t->isRooted = YES;  /* we use 'D' for 'directed' to indicate a rooted non-clock tree */
                expecting = Expecting(RIGHTCOMMENT);
                }
            else if (strcmp(tkn, "U") == 0)
                {
                t->isRooted = NO;
                expecting = Expecting(RIGHTCOMMENT);
                }
            else if (strcmp(tkn, "clockrate") == 0)
                {
                t->isCalibrated = YES;
                foundClockrate = YES;
                expecting = Expecting(EQUALSIGN);
                }
            else
                {
                inComment = YES;
                numComments++;
                expecting = lastExpecting;
                }
            foundAmpersand = NO;
            }
        else if (foundName == YES && foundCurly == YES)
            {
            if (strcmp("all",tkn) == 0)
                {
                SafeStrcat (&tempNameString,tkn);
                expecting = Expecting(RIGHTCURL);
                }
            else
                {
                YvyraPrint ("%s   Unrecognized argument '%s'\n", spacer, tkn);
                return (ERROR);
                }
            }
        else if (foundE == YES) /* We have seen &E */
            {
            if (foundEqual == NO) /* We have not seen name before and we are in header */
                {
                t->nESets++;
                t->isRelaxed = YES;
                t->nEvents  = (int **) SafeRealloc ((void *)t->nEvents, t->nESets*sizeof(int *));
                t->position = (YFlt ***) SafeRealloc ((void *)t->position, t->nESets*sizeof(YFlt **));
                t->rateMult = (YFlt ***) SafeRealloc ((void *)t->rateMult, t->nESets*sizeof(YFlt **));
                t->nEvents[t->nESets-1]  = (int *) SafeCalloc (2*(size_t)numTaxa, sizeof(int));
                t->position[t->nESets-1] = (YFlt **) SafeCalloc (2*(size_t)numTaxa, sizeof(YFlt *));
                t->rateMult[t->nESets-1] = (YFlt **) SafeCalloc (2*(size_t)numTaxa, sizeof(YFlt *));
                t->eSetName = (char **) SafeRealloc ((void *)t->eSetName, t->nESets*sizeof(char **));
                }
            SafeStrcpy (&tempNameString,tkn);
            foundName = YES;
            expecting = Expecting(LEFTCURL);
            if (foundEqual == YES)
                expecting |= Expecting(NUMBER);
            else
                expecting |= Expecting(RIGHTCOMMENT);
            }
        else if (foundB == YES)
            {
            if (foundEqual == NO)
                {
                t->nBSets++;
                t->isRelaxed = YES;
                t->effectiveBrLen = (YFlt **) SafeRealloc ((void *)t->effectiveBrLen, (size_t)(t->nBSets)*sizeof(YFlt *));
                t->effectiveBrLen[t->nBSets-1] = (YFlt *) SafeCalloc (2*(size_t)numTaxa, sizeof(YFlt));
                for (i=0; i<2*numTaxa; i++)
                    t->effectiveBrLen[t->nBSets-1][i] = 1.0;
                t->bSetName = (char **) SafeRealloc ((void *)t->bSetName, (size_t)(t->nBSets)*sizeof(char *));
                t->bSetName[t->nBSets-1] = (char *) SafeCalloc (strlen(tkn)+1, sizeof(char));
                }
            SafeStrcpy (&tempNameString,tkn);
            foundName = YES;
            expecting = Expecting(LEFTCURL);
            if (foundEqual == YES)
                expecting |= Expecting(NUMBER);
            else
                expecting |= Expecting(RIGHTCOMMENT);
            }
        else if (foundN == YES)
            {
            if (foundEqual == NO)
                {
                if (t->popSizeSet == YES)
                    {
                    YvyraPrint ("%s   Cannot hold more than one population size set\n", spacer);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                t->popSizeSet = YES;
                if (isTranslateDef == YES && isTranslateDiff == YES)
                    t->popSize = (YFlt *) SafeCalloc (2*numTranslates, sizeof(YFlt));
                else
                    t->popSize = (YFlt *) SafeCalloc (2*numLocalTaxa, sizeof(YFlt));
                }
            SafeStrcpy (&tempNameString,tkn);
            foundName = YES;
            expecting = Expecting(LEFTCURL);
            if (foundEqual == YES)
                expecting |= Expecting(NUMBER);
            else
                expecting |= Expecting(RIGHTCOMMENT);
            }
        else   /* taxon name */
            {
            if (isTranslateDef == YES)
                {
                /* we are using the translation table */
                if (CheckString (transTo, numTranslates, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find token '%s' in taxon translation table\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                strcpy (tempName, transFrom[index]);
                if (isTranslateDiff == NO && CheckString (taxaNames, numTaxa, tempName, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find taxon '%s' in list of taxa\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                if (tempSet[index] == YES)
                    {
                    YvyraPrint ("%s   Taxon name '%s' already used in tree\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                tempSet[index] = YES;
                /* FIXME: pp is NULL here (from clang static analyzer) */
                strcpy (pp->label, tempName);
                pp->index = index;
                }
            else
                {
                /* Check to see if the name is in the list of taxon names. */
                if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find taxon '%s' in list of taxa\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                if (tempSet[index] == YES)
                    {
                    YvyraPrint ("%s   Taxon name '%s' already used in tree\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                tempSet[index] = YES;
                /* FIXME: pp is NULL here (from clang static analyzer) */
                strcpy (pp->label, tkn);
                pp->index = index;
                }
            expecting  = Expecting(COMMA);
            expecting |= Expecting(COLON);
            expecting |= Expecting(RIGHTPAR);
            }
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        if (foundE == YES)
            expecting = Expecting(RIGHTCOMMENT);
        else
            {
            if (pp->anc == NULL)
                {
                YvyraPrint ("%s   Incorrect tree format: cannot go down\n", spacer);//, tkn
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            if (pp->anc->left == pp)
                {
                YvyraPrint ("%s   Incorrect tree format: all nodes except tips should have more then one child. Either a single\n", spacer);
                YvyraPrint ("%s   taxon is surrounded with brackets or there is a clade surrounded by double brackets.\n", spacer);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            pp = pp->anc;
            if (pp->anc == NULL)
                {
                /* finish up tree */
                t->nNodes = nextAvailableNode;
                t->nIntNodes = t->nNodes;
                for (i=0; i<t->nNodes; i++)
                    {
                    if (t->nodes[i].left == NULL)
                        t->nIntNodes--;
                    }
                GetPolyDownPass(t);
                
                /* check that number of taxa is correct */
                if (t->isRooted == NO && t->nNodes-t->nIntNodes == t->nIntNodes + 1)
                    t->isRooted = YES;
                if ((t->isRooted == YES && t->nNodes-t->nIntNodes != t->nIntNodes + 1) ||
                    (t->isRooted == NO  && t->nNodes-t->nIntNodes != t->nIntNodes + 2))
                    {
                    /* we are protected from adding too many taxa by taxon-matching code above */
                    if (t->isRooted == YES && t->nNodes-t->nIntNodes == t->nIntNodes + 2)
                        {
                        YvyraPrint ("%s   The tree is declared as rooted (by comment [&R]) but\n", spacer);
                        YvyraPrint ("%s   the given tree has unrooted structure.\n", spacer);
                        }
                    else
                        YvyraPrint ("%s   Taxa missing in tree, or NOT a binary tree\n", spacer);

                    return (ERROR);
                    }

                /* check other properties */
                if (t->isClock == YES && t->isRooted == NO)
                    {
                    YvyraPrint ("%s   Tree has clock rate but is not rooted\n", spacer);
                    return (ERROR);
                    /* Note: any deviation from an ultrametric tree must be assumed to be due to dated
                       tips at this point */
                    }
                if (t->isRelaxed == YES && t->isClock == NO)
                    {
                    YvyraPrint ("%s   Tree has relaxed clock rates but is not a clock tree\n", spacer);
                    return (ERROR);
                    }
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    {
                    if (treeIndex == numUserTrees)
                        numUserTrees++;
                    YvyraPrint ("%s   Successfully read tree '%s'\n", spacer, userTree[treeIndex]->name);
                    }
                if (t->popSize == NULL)
                    {
                    readComment = NO;
                    expecting = Expecting(SEMICOLON);
                    }
                else
                    {
                    readComment = YES;
                    expecting = Expecting(LEFTCOMMENT);
                    lastExpecting = expecting;
                    }
                }
            else
                {
                expecting = Expecting(COMMA);
                expecting |= Expecting(COLON);
                expecting |= Expecting(RIGHTPAR);
                }
            }
        }
    else if (expecting == Expecting(COLON))
        {
        foundColon = YES;
        if (foundE == YES)
            expecting = Expecting(LEFTPAR);
        else
            expecting  = Expecting(NUMBER);
        expecting |= Expecting(LEFTCOMMENT);
        lastExpecting = expecting;
        }
    else if (expecting == Expecting(COMMA))
        {
        if (foundName == YES)
            {
            SafeStrcat (&tempNameString,",");
            expecting = Expecting(NUMBER);
            }
        else if (foundE == YES)
            {
            expecting = Expecting(NUMBER);
            }
        else
            {
            if (nextAvailableNode >= 2*numTaxa)
                {
                YvyraPrint ("%s   Too many nodes on tree '%s'\n", spacer, t->name);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            /* FIXME: t is NULL here (from clang static analyzer) */
            qq = &t->nodes[nextAvailableNode++];
            pp->sib = qq;
            qq->anc = pp->anc;
            pp = qq;
            expecting = Expecting(LEFTPAR);
            expecting |= Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundClockrate == YES)
            {
            sscanf (tkn, "%lf", &tempD);
            t->clockRate = tempD;
            foundClockrate = NO;
            expecting = Expecting(RIGHTCOMMENT);
            }
        else if (foundName == YES && foundCurly == YES)
            {
            /* still assembling name of a param set */
            SafeStrcat (&tempNameString,tkn);       
            expecting = Expecting(RIGHTCURL) | Expecting(COMMA);
            }
        else if (foundN == YES)
            {
            /* we only know now that name is complete if it does not have curlies in it */
            foundName = NO;

            if (strcmp(tempNameString,t->popSizeSetName) != 0)
                {
                YvyraPrint ("%s   Could not find population size set '%s'\n", spacer, tempNameString);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }

            sscanf (tkn, "%lf", &tempD);
            t->popSize[pp->index] = tempD;
            foundN = NO;
            expecting = Expecting(RIGHTCOMMENT);
            }
        else if (foundB == YES)
            {
            /* we only know now that name is complete if it does not have curlies in it */
            foundName = NO;

            /* find the right effective branch length set */
            for (i=0; i<t->nBSets; i++)
                if (strcmp(t->bSetName[i],tempNameString) == 0)
                    break;
            if (i == t->nBSets)
                {
                YvyraPrint ("%s   Could not find effective branch length set '%s'\n", spacer, tempNameString);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            bSetIndex = i;

            sscanf (tkn, "%lf", &tempD);
            t->effectiveBrLen[bSetIndex][pp->index] = tempD;
            foundB = NO;
            expecting = Expecting(RIGHTCOMMENT);
            }
        else if (foundE == YES)
            {
            if (foundColon == NO)
                {
                /* we only know now that name is complete if it does not have curlies in it */
                foundName = NO;

                /* find the right event set */
                for (i=0; i<t->nESets; i++)
                    if (strcmp(t->eSetName[i],tempNameString) == 0)
                        break;
                if (i == t->nESets)
                    {
                    YvyraPrint ("%s   Could not find event set '%s'\n", spacer, tempNameString);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                eSetIndex = i;

                sscanf (tkn, "%d", &tempInt);
                if (tempInt < 0)
                    {
                    YvyraPrint ("%s   Wrong number of events (%d) for event set '%s'\n", spacer, tempInt, t->eSetName[eSetIndex]);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                t->nEvents[eSetIndex][pp->index]  = tempInt;
                if (tempInt > 0)
                    {
                    t->position[eSetIndex][pp->index] = (YFlt *) SafeCalloc (tempInt, sizeof(YFlt));
                    t->rateMult[eSetIndex][pp->index] = (YFlt *) SafeCalloc (tempInt, sizeof(YFlt));
                    expecting = Expecting (COLON);
                    if (inSumtCommand == YES || inComparetreeCommand == YES)
                        expecting |= Expecting (RIGHTCOMMENT);  /* we allow empty event specifications in sumt and comparetree */
                    }
                else
                    expecting = Expecting (RIGHTCOMMENT);
                eventIndex = 0;
                }
            else if (foundFirst == NO)
                {
                /* processing the first number in the cpp event pair <position rate> */
                sscanf (tkn, "%lf", &tempD);
                t->position[eSetIndex][pp->index][eventIndex] = tempD;
                expecting = Expecting(NUMBER);
                foundFirst = YES;
                }
            else
                {
                /* processing the second number in the cpp event pair <position rate> */
                foundFirst = NO;
                sscanf (tkn, "%lf", &tempD);
                t->rateMult[eSetIndex][pp->index][eventIndex] = tempD;
                eventIndex++;
                if (eventIndex == t->nEvents[eSetIndex][pp->index])
                    {
                    expecting = Expecting(RIGHTPAR);
                    foundColon = NO;
                    }
                else
                    expecting = Expecting(COMMA);
                }
            }
        else if (foundColon == YES)
            {
            /* branch length */
            sscanf (tkn, "%lf", &tempD);
            pp->length = tempD;
            foundColon = NO;
            t->brlensDef = YES;
            expecting  = Expecting(COMMA);
            expecting |= Expecting(RIGHTPAR);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        else    /* taxon identifier */
            {
            if (isTranslateDef == YES)
                {
                /* we are using the translation table */
                if (CheckString (transTo, numTranslates, tkn, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find token '%s' in taxon translation table\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                strcpy (tempName, transFrom[index]);
                if (isTranslateDiff == NO && CheckString (taxaNames, numTaxa, tempName, &index) == ERROR)
                    {
                    YvyraPrint ("%s   Could not find taxon '%s' in list of taxa\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                if (tempSet[index] == YES)
                    {
                    YvyraPrint ("%s   Taxon name '%s' already used in tree\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                tempSet[index] = YES;
                /* FIXME: pp is NULL here (from clang static analyzer) */
                strcpy (pp->label, tempName);
                pp->index = index;
                }
            else
                {
                /* Simply use taxon number; first check to see if the name is in the list of taxon names. */
                if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                    {
                    /* The number could not be found as a taxon name in the list of taxon names. We will
                        assume that the user has then input taxa as numbers and not the names. */
                    sscanf (tkn, "%d", &index);
                    if (index < 1 || index > numTaxa)
                        {
                        YvyraPrint ("%s   Taxon number %d is out of range\n", spacer, index);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    index--;
                    if (tempSet[index] == YES)
                        {
                        YvyraPrint ("%s   Taxon name %d has already been used in tree '%s'\n", spacer, index+1, t->name);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    }
                else
                    {
                    /* The number is in the list of taxon names */
                    if (index < 0 || index >= numTaxa)
                        {
                        YvyraPrint ("%s   Taxon name %s could not be found\n", spacer, tkn);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    if (tempSet[index] == YES)
                        {
                        YvyraPrint ("%s   Taxon %d has already been used in tree '%s'\n", spacer, index+1, t->name);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    }
                tempSet[index] = YES;
                /* FIXME: pp is NULL here (from clang static analyzer) */
                strcpy (pp->label, taxaNames[index]);
                pp->index = index;
                }
            expecting  = Expecting(COMMA);
            expecting |= Expecting(COLON);
            expecting |= Expecting(RIGHTPAR);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(LEFTCOMMENT))
        {
        expecting = Expecting(AMPERSAND);
        foundComment = YES;
        }
    else if (expecting == Expecting(RIGHTCOMMENT))
        {
        if (foundEqual == NO)
            {
            /* We may have a complete name of a set of branch parameters, which needs to be recorded */
            if (foundName == YES)
                {
                if (foundE == YES)
                    {
                    t->eSetName[t->nESets-1] = (char *) SafeCalloc (strlen(tempNameString)+1,sizeof(char));
                    strcat(t->eSetName[t->nESets-1],tempNameString);
                    }
                else if (foundB == YES)
                    {
                    t->bSetName[t->nBSets-1] = (char *) SafeCalloc (strlen(tempNameString)+1,sizeof(char));
                    strcat(t->bSetName[t->nBSets-1],tempNameString);
                    }
                else if (foundN == YES)
                    {
                    t->popSizeSetName = (char *) SafeCalloc (strlen(tempNameString)+1,sizeof(char));
                    strcpy(t->popSizeSetName,tempNameString);
                    }
                foundName = NO;
                }
            expecting = Expecting(EQUALSIGN);
            }
        else
            {
            /* FIXME: pp is NULL here (from clang static analyzer) */
            if (pp->anc == NULL)
                {
                if (pp->left == NULL)
                    expecting = Expecting(LEFTPAR);
                else
                    expecting = Expecting(SEMICOLON);
                }
            else if (pp == pp->anc->left)
                expecting = Expecting(COMMA);
            else
                expecting = Expecting(RIGHTPAR);
            }
        foundE = foundB = foundN = NO;
        expecting |= Expecting(LEFTCOMMENT);
        }
    else if (expecting == Expecting(AMPERSAND))
        {
        foundAmpersand = YES;
        foundComment = NO;
        expecting = Expecting (ALPHA);
        }
    else if (foundComment == YES)
        {
        numComments++;
        foundComment = NO;
        }
    else if (expecting == Expecting(LEFTCURL))
        {
        if (foundName == YES)
            {
            foundCurly=YES;
            SafeStrcat (&tempNameString,"{");               
            expecting = Expecting(NUMBER) | Expecting(ALPHA);
            }
        else
            return(ERROR);
        }
    else if (expecting == Expecting(RIGHTCURL))
        {
        if (foundName == YES)
            {
            SafeStrcat (&tempNameString,"}");
            foundCurly=NO;
            if (foundEqual == NO)
                {
                /* We are processing a name of a set of branch params in the header of a tree.  */
                expecting = Expecting(RIGHTCOMMENT);
                }
            else
                {
                /* We are processing a param value of a branch param set  */
                expecting = Expecting(NUMBER);
                }
            }
        else
            return(ERROR);
        }

    return (NO_ERROR);
}


int DoUserTree (void)
{
    YvyraPrint ("%s   Usertree command deprecated. Define the tree in a treeblock and use 'Startvals' instead.\n", spacer);
    return (ERROR);
}


int DoUserTreeParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(LEFTPAR);
        expecting |= Expecting(RIGHTPAR);
        expecting |= Expecting(COLON);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        expecting = Expecting(LEFTPAR);
        expecting |= Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        expecting = Expecting(COLON);
        expecting |= Expecting(COMMA);
        expecting |= Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(NUMBER))
        {
        expecting = Expecting(COLON);
        expecting |= Expecting(COMMA);
        expecting |= Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COMMA))
        {
        expecting = Expecting(LEFTPAR);
        expecting |= Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        expecting = Expecting(RIGHTPAR);
        expecting |= Expecting(COMMA);
        expecting |= Expecting(COLON);
        expecting |= Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}

int DoVersion (void)
{
    YvyraPrint ("   ---------------------------------------------------------------------------\n");
    YvyraPrint ("   Version\n");
    YvyraPrint ("\n");
    YvyraPrint ("   yvyra %s\n", VERSION_NUMBER);
    YvyraPrint ("\n");
    YvyraPrint ("   Features: ");
#ifdef SSE_ENABLED
    YvyraPrint (" SSE");
#endif
#ifdef AVX_ENABLED
    YvyraPrint (" AVX");
#endif
#ifdef FMA_ENABLED
    YvyraPrint (" FMA");
#endif
#ifdef MPI_ENABLED
    YvyraPrint (" MPI");
#endif
#ifdef HAVE_LIBREADLINE
    YvyraPrint (" readline");
#endif
    YvyraPrint ("\n");
#if defined(HOST_TYPE) && defined(HOST_CPU)
    YvyraPrint ("   Host type: %s (CPU: %s)\n", HOST_TYPE, HOST_CPU);
#endif
#if defined(COMPILER_VENDOR) && defined(COMPILER_VERSION)
    YvyraPrint ("   Compiler:  %s %s\n", COMPILER_VENDOR, COMPILER_VERSION);
#endif
    YvyraPrint ("   ---------------------------------------------------------------------------\n");

    return (NO_ERROR);
}


BitsLong Expecting (int y)
{
    BitsLong x;
    
    x = (BitsLong) pow (2.0, (YFlt)y);
    
    return (x);
}


#ifdef HAVE_LIBREADLINE
/* This function is for commandline substitution: first word is always a command */
char *command_generator (const char *text, int state)
{
    static int      list_index, len;
    char           *command;
    char           *dupstring;

    if (state == 0)
        {
        list_index = 0;
        len = (int) strlen (text);
        }

    while ((command = commands[list_index].string) != NULL)
        {
        list_index++;

        if (StrCmpCaseInsensitiveLen (command, text, len) == 0)
            {
            /* memory is freed by the readline library so we need a strdup here */
            dupstring = SafeMalloc (strlen (command) + 1);
            strcpy (dupstring, command);
            return dupstring;
            }
        }

    return NULL;
}
#endif


int FindValidCommand (char *tk, int *numMatches)
{
    int             i, j, tkLen, targetLen, numDiff;
    CmdType         *p;

    p = commands + 0;
    tkLen = (int) strlen(tk);

    (*numMatches) = 0;
    for (i=0; i<NUMCOMMANDS; i++)
        {
        targetLen = (int) strlen(p->string);
        if (tkLen <= targetLen)
            {
            for (j=0, numDiff=0; j<tkLen; j++)
                {
                if (ChangeCase(tk[j]) != ChangeCase(p->string[j]))
                    numDiff++;
                }
            if (numDiff == 0)
                {
                (*numMatches)++;
                commandPtr = p;
                if (tkLen == targetLen)
                    break;
                }
            }
        p++;
        }

    inValidCommand = NO;
    if (*numMatches == 1)
        {
        inValidCommand = YES;
        return (NO_ERROR);
        }
    else
        return (ERROR);
}


int FindValidParam (char *tk, int *numMatches)
{
    int         i, j, tkLen, targetLen, numDiff;
    CmdType     *p;
    ParmInfoPtr q;

    if (commandPtr)
        p = commandPtr;
    else
        {
        YvyraPrint ("%s   Command pointer is NULL\n", spacer);
        return (ERROR);
        }
    tkLen = (int) strlen(tk);

    *numMatches = 0;
    for (i=0; i<p->numParms; i++)
        {
        q = paramTable + (p->parmList[i]);
        targetLen = (int) strlen(q->string);
        /* printf ("%s %d (%s %d)\n", q->string, targetLen, tk, p->numParms); */
        if (!strcmp(q->string, "Xxxxxxxxxx"))
            {
            (*numMatches)++;
            paramPtr = q;
            }
        else if (tkLen <= targetLen)
            {
            for (j=0, numDiff=0; j<tkLen; j++)
                {
                if (ChangeCase(tk[j]) != ChangeCase(q->string[j]))
                    numDiff++;
                }
            if (numDiff == 0)
                {
                (*numMatches)++;
                paramPtr = q;
                if (tkLen == targetLen)
                    break;
                }
            }   
        }
    
    if (*numMatches == 1)
        return (NO_ERROR);
    else
        return (ERROR);
}


int FreeCharacters (void)
{
    int     i, memoryLetFree;
    
    memoryLetFree = NO;

    if (memAllocs[ALLOC_TMPSET] == YES)
        {
        if (numChar > numTaxa)
            {
            tempSet = (int *) SafeRealloc ((void *)tempSet, (size_t)numTaxa*sizeof(int));
            tempSetNeg = (int *) SafeRealloc ((void *)tempSetNeg, (size_t)numTaxa*sizeof(int));
            }
        }
    if (memAllocs[ALLOC_MATRIX] == YES)
        {
        free (matrix);
        matrix = NULL;
        defMatrix = NO;
        memAllocs[ALLOC_MATRIX] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_CHARINFO] == YES)
        {
        free (charInfo);
        charInfo = NULL;
        memAllocs[ALLOC_CHARINFO] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_CHARSETS] == YES)
        {
        for (i=0; i<numCharSets; i++)
            {
            free (charSetNames[i]);
            free (charSet[i]);
            }
        free (charSetNames);
        free (charSet);
        charSetNames = NULL;
        charSet = NULL;
        numCharSets = 0;
        memAllocs[ALLOC_CHARSETS] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_PARTITIONS] == YES)
        {
        for (i=0; i<numDefinedPartitions; i++)
            free (partitionNames[i]);
        free (partitionNames);
        partitionNames = NULL;
        for (i=0; i<numChar; i++)
            free (partitionId[i]);
        free (partitionId);
        numDefinedPartitions = 0;
        memAllocs[ALLOC_PARTITIONS] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_PARTITIONVARS] == YES)
        {
        free (numVars);
        numVars = NULL;
        free (tempNum);
        tempNum = NULL;
        free (activeParams[0]);
        activeParams[0] = NULL;
        free (linkTable[0]);
        linkTable[0] = NULL;
        tempLinkUnlinkVec = NULL;
        activeParts = NULL;
        tempLinkUnlinkVec = NULL;
        for (i=0; i<NUM_LINKED; i++)
            {
            linkTable[i] = NULL;
            activeParams[i] = NULL;
            }
        memAllocs[ALLOC_PARTITIONVARS] = NO;
        memoryLetFree = YES;
        }

    /* Free charLabels */
    if (charLabels != NULL)
        {
        for (i=0; i<numChar; i++)
            if (charLabels[i] != NULL)
                free(charLabels[i]);
        free(charLabels);
        charLabels = NULL;
        }

    /* Free tipWeights */
    if (tipWeights != NULL)
        {
        free(tipWeights);
        tipWeights = NULL;
        hasTipWeights = NO;
        }

    ResetCharacterFlags();

    if (memoryLetFree == YES)
        YvyraPrint ("%s   Deleting previously defined characters\n", spacer);

    return (NO_ERROR);
}


int FreeMatrix (void)
{
    if (FreeCharacters() == ERROR)
        return ERROR;

    return (FreeTaxa());
}


int FreeTaxa (void)
{
    int             i, memoryLetFree;

    memoryLetFree = NO;

    if (memAllocs[ALLOC_SPECIESPARTITIONS] == YES)
        {
        for (i = 0; i < numDefinedSpeciespartitions; i++)
            SAFEFREE (speciespartitionNames[i]);

        SAFEFREE (speciespartitionNames);

        for (i = 0; i < numTaxa; i++)
            SAFEFREE (speciespartitionId[i]);

        SAFEFREE (speciespartitionId);
        numDefinedSpeciespartitions = 0;
        memAllocs[ALLOC_SPECIESPARTITIONS] = NO;
        memoryLetFree = YES;
        }

    if (memAllocs[ALLOC_TAXA] == YES)
        {
        if (taxaNames)
            {
            for (i = 0; i < taxonCount; i++)
                SAFEFREE (taxaNames[i]);
            }

        SAFEFREE (taxaNames);
        SAFEFREE (taxaInfo);
        SAFEFREE (tipCalibration);
        numTaxa = 0;
        memAllocs[ALLOC_TAXA] = NO;
        memoryLetFree = YES;
        }

    if (memAllocs[ALLOC_TMPSET] == YES)
        {
        SAFEFREE (tempSet);
        SAFEFREE (tempSetNeg);
        memAllocs[ALLOC_TMPSET] = NO;
        memoryLetFree = YES;
        }

    if (memAllocs[ALLOC_TAXASETS] == YES)
        {
        for (i = 0; i < numTaxaSets; i++)
            {
            SAFEFREE (taxaSetNames[i]);
            SAFEFREE (taxaSet[i]);
            }

        SAFEFREE (taxaSetNames);
        SAFEFREE (taxaSet);
        numTaxaSets = 0;
        memAllocs[ALLOC_TAXASETS] = NO;
        memoryLetFree = YES;
        }

    if (memAllocs[ALLOC_CONSTRAINTS] == YES)
        {
        for (i = 0; i < numDefinedConstraints; i++)
            {
            SAFEFREE (definedConstraint[i]);
            SAFEFREE (definedConstraintTwo[i]);
            SAFEFREE (definedConstraintPruned[i]);
            SAFEFREE (definedConstraintTwoPruned[i]);
            SAFEFREE (constraintNames[i]);
            }

        SAFEFREE (definedConstraint);
        SAFEFREE (definedConstraintTwo);
        SAFEFREE (definedConstraintsType);
        SAFEFREE (constraintNames);
        SAFEFREE (nodeCalibration);
        numDefinedConstraints = 0;
        SAFEFREE (tempActiveConstraints);
        memAllocs[ALLOC_CONSTRAINTS] = NO;
        memoryLetFree = YES;
        }

    if (numUserTrees > 0)
        {
        YvyraPrint ("%s   Deleting user trees\n", spacer);

        for (i = 0; i < numUserTrees; i++)
            {
            FreePolyTree (userTree[i]);
            userTree[i] = NULL;
            }

        numUserTrees = 0;
        }

    FreeCharacters();

    if (memoryLetFree == YES)
        YvyraPrint ("%s   Deleting previously defined taxa\n", spacer);

    /* reinitialize taxa variables */
    ResetTaxaFlags();

    return NO_ERROR;
}


int GetNumPartDivisions (int n)
{
    int         i, maxDiv, numDivs, *divFound;
    
    maxDiv = 0;
    for (i=0; i<numChar; i++)
        if (partitionId[i][n] > maxDiv)
            maxDiv = partitionId[i][n];

    divFound = (int *) SafeCalloc (maxDiv, sizeof(int));
    
    for (i=0; i<maxDiv; i++)
        divFound[i] = NO;
    
    for (i=0; i<numChar; i++)
        divFound[partitionId[i][n]] = YES;
        
    numDivs = 0;
    for (i=0; i<maxDiv; i++)
        if (divFound[i] == YES)
            numDivs++;
    
    free (divFound);

    return (numDivs + 1);
}


int GetToken (char *token, int *tokenType, char **sourceH)
{
    int             allNumbers, foundExp, foundExpSign;
    register char   *temp;
    char            *tempMax;
    
    (*tokenType) = 0;
    temp = token;
    tempMax = temp + CMD_STRING_LENGTH - 10;
    
    while (IsWhite(**sourceH) == 1 || IsWhite(**sourceH) == 2)
        {
        if (IsWhite(**sourceH) == 2)
            {
            *tokenType = RETURNSYMBOL;
            /* foundNewLine = YES;  Why is this commented out?? */
            /* YvyraPrint ("RETURN\n"); */
            }
        ++(*sourceH);
        }
    
    if (readWord == YES && **sourceH != '"')
        {
        if (**sourceH==';')
            {
            *temp++ = ';';
            *tokenType = SEMICOLON;
            }
        else
            {
            while (isgraph(**sourceH) && **sourceH!=';')
                {
                if (temp > tempMax)
                    {
                    *tokenType = NOTHING;
                    token[20]='\0';
                    YvyraPrint ("%s   Error while parsing a string. Token \"%s...[followed by at least %d  more characters]\" is too long.\n", spacer,token,tempMax-token-20);
                    YvyraPrint ("%s   Maximum allowed length of a token is %d\n", spacer,tempMax-token);
                    return (ERROR);
                    }
                *temp++ = *(*sourceH)++;
                }
            *tokenType = ALPHA;
            }
        *temp = '\0';
        readWord = NO;
        return (NO_ERROR);;
        }

    *tokenType = UNKNOWN_TOKEN_TYPE;
    if (IsIn(**sourceH,"="))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = EQUALSIGN;
        }
    else if (IsIn(**sourceH,";"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = SEMICOLON;
        }
    else if (IsIn(**sourceH,":"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = COLON;
        }
    else if (IsIn(**sourceH,","))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = COMMA;
        }
    else if (IsIn(**sourceH,"#"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = POUNDSIGN;
        }
    else if (IsIn(**sourceH,"("))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = LEFTPAR;
        }
    else if (IsIn(**sourceH,")"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = RIGHTPAR;
        }
    else if (IsIn(**sourceH,"{"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = LEFTCURL;
        }
    else if (IsIn(**sourceH,"}"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = RIGHTCURL;
        }
    else if (IsIn(**sourceH,"["))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = LEFTCOMMENT;
        }
    else if (IsIn(**sourceH,"]"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = RIGHTCOMMENT;
        }
    else if (IsIn(**sourceH,"?"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = QUESTIONMARK;
        }
    else if (IsIn(**sourceH,"-"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = DASH;
        }
    else if (IsIn(**sourceH,"$"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = DOLLAR;
        }
    else if (IsIn(**sourceH,"\"") && readWord == YES)
        {
        (*sourceH)++;
        while (**sourceH != '"' && **sourceH != '\0')
            {
            if (temp > tempMax)
                {
                *tokenType = NOTHING;
                token[20]='\0';
                YvyraPrint ("%s   Error while parsing a string. Token \"%s...[followed by at least %d  more characters]\" is too long.\n", spacer,token,tempMax-token-20);
                YvyraPrint ("%s   Maximum allowed length of a token is %d\n", spacer,tempMax-token);
                return (ERROR);
                }
            *temp++ = *((*sourceH)++);
            }
        *temp='\0';
        *tokenType = ALPHA;
        (*sourceH)++;
        readWord = NO;
        }
    else if (IsIn(**sourceH,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789."))
        {
        if (IsIn(**sourceH,"0123456789."))
            allNumbers = TRUE;
        else
            allNumbers = FALSE;
        foundExp = foundExpSign = FALSE;
        *temp++ = *(*sourceH)++;
        while (IsIn(**sourceH,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.-+"))
            {
            if (temp > tempMax)
                {
                *tokenType = NOTHING;
                token[20]='\0';
                YvyraPrint ("%s   Error while parsing a string. Token \"%s...[followed by at least %d  more characters]\" is too long.\n", spacer,token,tempMax-token-20);
                YvyraPrint ("%s   Maximum allowed length of a token is %d\n", spacer,tempMax-token);
                return (ERROR);
                }
            if (allNumbers == TRUE && !IsIn((*sourceH)[-1],"Ee") && **sourceH=='-')
                break;
            else if (allNumbers == TRUE && IsIn(**sourceH,"Ee") && foundExp == NO)
                foundExp = TRUE;
            else if (allNumbers == TRUE && IsIn(**sourceH,"+-") && IsIn((*sourceH)[-1],"Ee"))
                foundExpSign = TRUE; /* FIXME: Not used (from clang static analyzer) */
            else if (!IsIn(**sourceH,"0123456789."))
                allNumbers = FALSE;
            *temp++ = *(*sourceH)++;
            }
        if (allNumbers == TRUE)
            *tokenType = NUMBER;
        else
            *tokenType = ALPHA;
        }
    else if (IsIn(**sourceH,"*"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = ASTERISK;
        }
    else if (IsIn(**sourceH,"/"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = FORWARDSLASH;
        }
    else if (IsIn(**sourceH,"'\\'"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = BACKSLASH;
        }
    else if (IsIn(**sourceH,"!"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = EXCLAMATIONMARK;
        }
    else if (IsIn(**sourceH,"%"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = PERCENT;
        }
    else if (IsIn(**sourceH,"\""))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = QUOTATIONMARK;
        }
    else if (IsIn(**sourceH,"&"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = AMPERSAND;
        }
    else if (IsIn(**sourceH,"~+^@{}`><"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = WEIRD;
        }
    else if (IsIn(**sourceH,"|"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = VERTICALBAR;
        }

    *temp = '\0';
    return (NO_ERROR);
}


int GetUserHelp (char *helpTkn)
{
    int         i, j, k, tempInt;
    char        tempString[100];
    Model       *mp;
    
    if (!strcmp(helpTkn, "Begin"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Begin                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to format data or commands in the program. The correct   \n");
        YvyraPrint ("   usage is                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin <data or yvyra>;                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The valid uses of the \"begin\" command are                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("      begin yvyra;     (or 'begin mrbayes;' for backward compatibility)          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The \"data\" specifier is used to specify the beginning of a data block; your \n");
        YvyraPrint ("   character data should follow. For example, the following is an example of     \n");
        YvyraPrint ("   a data block for four taxa and ten DNA sites:                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        YvyraPrint ("         format datatype=dna;                                                    \n");
        YvyraPrint ("         matrix                                                                  \n");
        YvyraPrint ("         taxon_1  AACGATTCGT                                                     \n");
        YvyraPrint ("         taxon_2  AAGGATTCCA                                                     \n");
        YvyraPrint ("         taxon_3  AACGACTCCT                                                     \n");
        YvyraPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        YvyraPrint ("         ;                                                                       \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The other commands -- dimensions, format, and matrix -- are discussed         \n");
        YvyraPrint ("   in the appropriate help menu. The only thing to note here is that the         \n");
        YvyraPrint ("   block begins with a \"begin data\" command. The \"yvyra\" command block       \n");
        YvyraPrint ("   is used to enter commands specific to the yvyra program into the file.       \n");
        YvyraPrint ("   This allows you to automatically process commands on execution of the         \n");
        YvyraPrint ("   program. The following is a simple yvyra block:                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin yvyra;                                                               \n");
        YvyraPrint ("         charset first  = 1-10\\3;                                               \n");
        YvyraPrint ("         charset second = 2-10\\3;                                               \n");
        YvyraPrint ("         charset third  = 3-10\\3;                                               \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This yvyra block sets off the three \"charset\" commands, used to             \n");
        YvyraPrint ("   predefine some blocks of characters. The yvyra block can be very useful.      \n");
        YvyraPrint ("   For example, in this case, it would save you the time of typing the char-     \n");
        YvyraPrint ("   acter sets each time you executed the file. Also, note that every             \n");
        YvyraPrint ("   \"begin <data or yvyra>\" command ends with an \"end\". Finally, you can      \n");
        YvyraPrint ("   have so-called foreign blocks in the file. An example of a foreign block      \n");
        YvyraPrint ("   would be \"begin paup\". The program will simply skip this block. This is     \n");
        YvyraPrint ("   useful because it means that you can use the same file for yvyra, PAUP*     \n");
        YvyraPrint ("   or MacClade (although it isn't clear why you would want to use those other    \n");
        YvyraPrint ("   programs).                                                                    \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "End"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   End                                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to terminate a data or yvyra block. The correct        \n");
        YvyraPrint ("   usage is                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For more information on this, check the help for the \"begin\" command.       \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Endblock"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Endblock                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This is an older, deprecated version of \"End\", see that command.            \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Plot"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Plot                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command plots specified parameters in the .p file or one of the .p files \n");
        YvyraPrint ("   created during an MCMC analysis. An x-y graph of the parameter over the course\n");
        YvyraPrint ("   of the chain is created. The command can be useful for visually diagnosing    \n");
        YvyraPrint ("   convergence for many of the parameters of the phylogenetic model. The para-   \n");
        YvyraPrint ("   meter to be plotted is specified by the \"parameter\" option. Several para-   \n");
        YvyraPrint ("   meters can be plotted at once by using the \"match\" option, which has a      \n");
        YvyraPrint ("   default value of \"perfect\". For example, if you were to set \"parameter = pi\"\n");
        YvyraPrint ("   and \"match = consistentwith\", then all of the state frequency parameters    \n");
        YvyraPrint ("   would be plotted. You can also set \"match=all\", in which case all of the    \n");
        YvyraPrint ("   parameters are plotted.                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that the \"Sump\" command provides a different set of convergence diag-  \n");
        YvyraPrint ("   nostics tools that you may also want to explore. Unlike \"Plot\", \"Sump\" can\n");
        YvyraPrint ("   compare two or more parameter samples and will calculate convergence diagnos- \n");
        YvyraPrint ("   tics as well as parameter summaries for the pooled sample.                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Relburnin     -- If this option is set to 'Yes', then a proportion of the     \n");
        YvyraPrint ("                    samples will be discarded as burnin when creating the plot.  \n");
        YvyraPrint ("                    The proportion to be discarded is set with Burninfrac (see   \n");
        YvyraPrint ("                    Burninfrac below). When the Relburnin option is set to 'No', \n");
        YvyraPrint ("                    then a specific number of samples is discarded instead. This \n");
        YvyraPrint ("                    number is set by Burnin (see below). Note that the burnin    \n");
        YvyraPrint ("                    setting is shared across the 'comparetree', 'sump' and 'sumt'\n");
        YvyraPrint ("                    commands.                                                    \n");
        YvyraPrint ("   Burnin        -- Determines the number of samples (not generations) that will \n");
        YvyraPrint ("                    be discarded when summary statistics are calculated. The     \n");
        YvyraPrint ("                    value of this option is only relevant when Relburnin is set  \n");
        YvyraPrint ("                    to 'No'.                                                     \n");
        YvyraPrint ("   Burninfrac    -- Determines the fraction of samples that will be discarded    \n");
        YvyraPrint ("                    when creating a plot. The value of this parameter is only    \n");
        YvyraPrint ("                    relevant when Relburnin is set to 'Yes'. Example: A value of \n");
        YvyraPrint ("                    this option of 0.25 means that 25%% of the samples will be   \n");
        YvyraPrint ("                    discarded.                                                   \n");
        YvyraPrint ("   Filename      -- The name of the file to plot.                                \n");
        YvyraPrint ("   Parameter     -- Specification of parameters to be plotted. See above for     \n");
        YvyraPrint ("                    details.                                                     \n");
        YvyraPrint ("   Match         -- Specifies how to match parameter names to the Parameter      \n");
        YvyraPrint ("                    specification. See above for details.                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter       Options                      Current Setting                  \n");
        YvyraPrint ("   ------------------------------------------------------------                  \n");
        YvyraPrint ("   Relburnin       Yes/No                       %s                               \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        YvyraPrint ("   Burnin          <number>                     %d                               \n", chainParams.chainBurnIn);
        YvyraPrint ("   Burninfrac      <number>                     %1.2lf                           \n", chainParams.burninFraction);
        YvyraPrint ("   Filename        <name>                       %s                               \n", plotParams.plotFileName);
        YvyraPrint ("   Parameter       <name>                       %s                               \n", plotParams.parameter);
        YvyraPrint ("   Match           Perfect/Consistentwith/All   %s                               \n", plotParams.match);
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Dimensions"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Dimensions                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used in a data block to define the number of taxa and         \n");
        YvyraPrint ("   characters. The correct usage is                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      dimensions ntax=<number> nchar=<number>                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The dimensions must be the first command in a data block. The following       \n");
        YvyraPrint ("   provides an example of the proper use of this command:                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        YvyraPrint ("         format datatype=dna;                                                    \n");
        YvyraPrint ("         matrix                                                                  \n");
        YvyraPrint ("         taxon_1  AACGATTCGT                                                     \n");
        YvyraPrint ("         taxon_2  AAGGATTCCA                                                     \n");
        YvyraPrint ("         taxon_3  AACGACTCCT                                                     \n");
        YvyraPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        YvyraPrint ("         ;                                                                       \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Here, the dimensions command tells yvyra to expect a matrix with four       \n");
        YvyraPrint ("   taxa and 10 characters.                                                       \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Format"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Format                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used in a data block to define the format of the char-        \n");
        YvyraPrint ("   acter matrix. The correct usage is                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      format datatype=<name> ... <parameter>=<option>                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The format command must be the second command in a data block. The following  \n");
        YvyraPrint ("   provides an example of the proper use of this command:                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        YvyraPrint ("         format datatype=dna gap=-;                                              \n");
        YvyraPrint ("         matrix                                                                  \n");
        YvyraPrint ("         taxon_1  AACGATTCGT                                                     \n");
        YvyraPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        YvyraPrint ("         taxon_3  AACGACTCCT                                                     \n");
        YvyraPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        YvyraPrint ("         ;                                                                       \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Here, the format command tells yvyra to expect a matrix with DNA char-      \n");
        YvyraPrint ("   acters and with gaps coded as \"-\".                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The following are valid options for format:                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Datatype   -- This parameter MUST BE INCLUDED in the format command. More-    \n");
        YvyraPrint ("                 over, it must be the first parameter in the line. The           \n");
        YvyraPrint ("                 datatype command specifies what type of characters are          \n");
        YvyraPrint ("                 in the matrix. The following are valid options:                 \n");
        YvyraPrint ("                    Datatype = Dna: DNA states (A,C,G,T,R,Y,M,K,S,W,H,B,         \n");
        YvyraPrint ("                               V,D,N)                                            \n");
        YvyraPrint ("                    Datatype = Rna: DNA states (A,C,G,U,R,Y,M,K,S,W,H,B,         \n");
        YvyraPrint ("                               V,D,N)                                            \n");
        YvyraPrint ("                    Datatype = Protein: Amino acid states (A,R,N,D,C,Q,E,        \n");
        YvyraPrint ("                               G,H,I,L,K,M,F,P,S,T,W,Y,V)                        \n");
        YvyraPrint ("                    Datatype = Restriction: Restriction site (0,1) states        \n");
        YvyraPrint ("                    Datatype = Standard: Morphological (0,1) or (0,1,2)... states\n");
        YvyraPrint ("                    Datatype = Continuous: Real number valued states             \n");
        YvyraPrint ("                    Datatype = Mixed(<type>:<range>,...,<type>:<range>): A       \n");
        YvyraPrint ("                               mixture of the above datatypes. For example,      \n");
        YvyraPrint ("                               \"datatype=mixed(dna:1-100,protein:101-200)\"     \n");
        YvyraPrint ("                               would specify a mixture of DNA and amino acid     \n");
        YvyraPrint ("                               characters with the DNA characters occupying      \n");
        YvyraPrint ("                               the first 100 sites and the amino acid char-      \n");
        YvyraPrint ("                               acters occupying the last 100 sites.              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Interleave -- This parameter specifies whether the data matrix is in          \n");
        YvyraPrint ("                 interleave format. The valid options are \"Yes\" or \"No\",     \n");
        YvyraPrint ("                 with \"No\" as the default. An interleaved matrix looks like    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    format datatype=dna gap=- interleave=yes;                    \n");
        YvyraPrint ("                    matrix                                                       \n");
        YvyraPrint ("                    taxon_1  AACGATTCGT                                          \n");
        YvyraPrint ("                    taxon_2  AAGGAT--CA                                          \n");
        YvyraPrint ("                    taxon_3  AACGACTCCT                                          \n");
        YvyraPrint ("                    taxon_4  AAGGATTCCT                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    taxon_1  CCTGGTAC                                            \n");
        YvyraPrint ("                    taxon_2  CCTGGTAC                                            \n");
        YvyraPrint ("                    taxon_3  ---GGTAG                                            \n");
        YvyraPrint ("                    taxon_4  ---GGTAG                                            \n");
        YvyraPrint ("                    ;                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Gap        -- This parameter specifies the format for gaps. Note that         \n");
        YvyraPrint ("                 gap character can only be a single character and that it        \n");
        YvyraPrint ("                 cannot correspond to a standard state (e.g., A,C,G,T,R,Y,       \n");
        YvyraPrint ("                 M,K,S,W,H,B,V,D,N for nucleotide data).                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Missing    -- This parameter specifies the format for missing data. Note      \n");
        YvyraPrint ("                 that the missing character can only be a single character and   \n");
        YvyraPrint ("                 cannot correspond to a standard state (e.g., A,C,G,T,R,Y,       \n");
        YvyraPrint ("                 M,K,S,W,H,B,V,D,N for nucleotide data). This is often an        \n");
        YvyraPrint ("                 unnecessary parameter to set because many data types, such      \n");
        YvyraPrint ("                 as nucleotide or amino acid, already have a missing char-       \n");
        YvyraPrint ("                 acter specified. However, for morphological or restriction      \n");
        YvyraPrint ("                 site data, \"missing=?\" is often used to specify ambiguity     \n");
        YvyraPrint ("                 or unobserved data.                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Matchchar  -- This parameter specifies the matching character for the         \n");
        YvyraPrint ("                 matrix. For example,                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    format datatype=dna gap=- matchchar=.;                       \n");
        YvyraPrint ("                    matrix                                                       \n");
        YvyraPrint ("                    taxon_1  AACGATTCGT                                          \n");
        YvyraPrint ("                    taxon_2  ..G...--CA                                          \n");
        YvyraPrint ("                    taxon_3  .....C..C.                                          \n");
        YvyraPrint ("                    taxon_4  ..G.....C.                                          \n");
        YvyraPrint ("                    ;                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                 is equivalent to                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    format datatype=dna gap=-;                                   \n");
        YvyraPrint ("                    matrix                                                       \n");
        YvyraPrint ("                    taxon_1  AACGATTCGT                                          \n");
        YvyraPrint ("                    taxon_2  AAGGAT--CA                                          \n");
        YvyraPrint ("                    taxon_3  AACGACTCCT                                          \n");
        YvyraPrint ("                    taxon_4  AAGGATTCCT                                          \n");
        YvyraPrint ("                    ;                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The only non-standard NEXUS format option is the use of the \"mixed\",        \n");
        YvyraPrint ("   \"restriction\", and \"standard\" datatypes. Hence, if you use any of these   \n");
        YvyraPrint ("   datatype specifiers, a program like PAUP* or MacClade will report an error    \n");
        YvyraPrint ("   (as they should because yvyra is not strictly NEXUS compliant).             \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Matrix"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Matrix                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command specifies the actual data for the phylogenetic analysis.         \n");
        YvyraPrint ("   The character matrix should follow the dimensions and format commands         \n");
        YvyraPrint ("   in a data block. The matrix can have all of the characters for a taxon        \n");
        YvyraPrint ("   on a single line:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        YvyraPrint ("         format datatype=dna gap=-;                                              \n");
        YvyraPrint ("         matrix                                                                  \n");
        YvyraPrint ("         taxon_1  AACGATTCGT                                                     \n");
        YvyraPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        YvyraPrint ("         taxon_3  AACGACTCCT                                                     \n");
        YvyraPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        YvyraPrint ("         ;                                                                       \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or be in \"interleaved\" format:                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("         dimensions ntax=4 nchar=20;                                             \n");
        YvyraPrint ("         format datatype=dna gap=- interleave=yes;                               \n");
        YvyraPrint ("         matrix                                                                  \n");
        YvyraPrint ("         taxon_1  AACGATTCGT                                                     \n");
        YvyraPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        YvyraPrint ("         taxon_3  AACGACTCCT                                                     \n");
        YvyraPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("         taxon_1  TTTTCGAAGC                                                     \n");
        YvyraPrint ("         taxon_2  TTTTCGGAGC                                                     \n");
        YvyraPrint ("         taxon_3  TTTTTGATGC                                                     \n");
        YvyraPrint ("         taxon_4  TTTTCGGAGC                                                     \n");
        YvyraPrint ("         ;                                                                       \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that the taxon names must not have spaces. If you really want to         \n");
        YvyraPrint ("   indicate a space in a taxon name (perhaps between a genus and species         \n");
        YvyraPrint ("   name), then you might use an underline (\"_\"). There should be at            \n");
        YvyraPrint ("   least a single space after the taxon name, separating the name from           \n");
        YvyraPrint ("   the actual data on that line. There can be spaces between the char-           \n");
        YvyraPrint ("   acters.                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   If you have mixed data, then you specify all of the data in the same          \n");
        YvyraPrint ("   matrix. Here is an example that includes two different data types:            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      begin data;                                                                \n");
        YvyraPrint ("         dimensions ntax=4 nchar=20;                                             \n");
        YvyraPrint ("         format datatype=mixed(dna:1-10,standard:21-30) interleave=yes;          \n");
        YvyraPrint ("         matrix                                                                  \n");
        YvyraPrint ("         taxon_1  AACGATTCGT                                                     \n");
        YvyraPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        YvyraPrint ("         taxon_3  AACGACTCCT                                                     \n");
        YvyraPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("         taxon_1  0001111111                                                     \n");
        YvyraPrint ("         taxon_2  0111110000                                                     \n");
        YvyraPrint ("         taxon_3  1110000000                                                     \n");
        YvyraPrint ("         taxon_4  1000001111                                                     \n");
        YvyraPrint ("         ;                                                                       \n");
        YvyraPrint ("      end;                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The matrix command is terminated by a semicolon.                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Finally, just a note on data presentation. It is much easier for others       \n");
        YvyraPrint ("   to (1) understand your data and (2) repeat your analyses if you make          \n");
        YvyraPrint ("   your data clean, comment it liberally (using the square brackets), and        \n");
        YvyraPrint ("   embed the commands you used in a publication in the yvyra block.            \n");
        YvyraPrint ("   Remember that the data took a long time for you to collect. You might         \n");
        YvyraPrint ("   as well spend a little time making the data file look nice and clear to       \n");
        YvyraPrint ("   any that may later request the data for further analysis.                     \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Pairs"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Pairs                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to specify pairs of nucleotides. For example, your       \n");
        YvyraPrint ("   data may be RNA sequences with a known secondary structure of stems and       \n");
        YvyraPrint ("   loops. Substitutions in nucleotides involved in a Watson-Crick pairing        \n");
        YvyraPrint ("   in stems are not strictly independent; a change in one changes the prob-      \n");
        YvyraPrint ("   ability of a change in the partner. A solution to this problem is to          \n");
        YvyraPrint ("   expand the model around the pair of nucleotides in the stem. This             \n");
        YvyraPrint ("   command allows you to do this. The correct usage is:                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      pairs <NUC1>:<NUC2>, <NUC1>:<NUC2>,..., <NUC1>:<NUC2>;                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example,                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      pairs 30:56, 31:55, 32:54, 33:53, 34:52, 35:51, 36:50;                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   specifies pairings between nucleotides 30 and 56, 31 and 55, etc. Only        \n");
        YvyraPrint ("   nucleotide data (DNA or RNA) may be paired using this command. Note that      \n");
        YvyraPrint ("   in order for the program to actually implement a \"doublet\" model            \n");
        YvyraPrint ("   involving a 16 X 16 rate matrix, you must specify that the structure of       \n");
        YvyraPrint ("   the model is 16 X 16 using \"lset nucmodel=doublet\".                         \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Databreaks"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Databreaks                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to specify breaks in your input data matrix. Your        \n");
        YvyraPrint ("   data may be a mixture of genes or a mixture of different types of data.       \n");
        YvyraPrint ("   Some of the models implemented by yvyra account for nonindependence at      \n");
        YvyraPrint ("   adjacent characters. The autocorrelated gamma model, for example, allows      \n");
        YvyraPrint ("   rates at adjacent sites to be correlated. However, there is no way for        \n");
        YvyraPrint ("   such a model to tell whether two sites, adjacent in the matrix, are           \n");
        YvyraPrint ("   actually separated by many kilobases or megabases in the genome. The          \n");
        YvyraPrint ("   databreaks command allows you to specify such breaks. The correct             \n");
        YvyraPrint ("   usage is:                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      databreaks <break 1> <break 2> <break 3> ...                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example, say you have a data matrix of 3204 characters that include       \n");
        YvyraPrint ("   nucleotide data from three genes. The first gene covers characters 1 to       \n");
        YvyraPrint ("   970, the second gene covers characters 971 to 2567, and the third gene        \n");
        YvyraPrint ("   covers characters 2568 to 3204. Also, let's assume that the genes are         \n");
        YvyraPrint ("   not directly adjacent to one another in the genome, as might be likely        \n");
        YvyraPrint ("   if you have mitochondrial sequences. In this case, you can specify            \n");
        YvyraPrint ("   breaks between the genes using:                                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      databreaks 970 2567;                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The first break, between genes one and two, is after character 970 and        \n");
        YvyraPrint ("   the second break, between genes two and three, is after character 2567.       \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Acknowledgments"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Acknowledgments                                                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the authors' acknowledgments.                              \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "About"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   About                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command provides some general information about the program.             \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Version"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Version                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the release version of the program.                        \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Citations"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Citations                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows a thorough list of citations you may consider using        \n");
        YvyraPrint ("   when publishing the results of a yvyra analysis.                            \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmatrix"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Showmatrix                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the character matrix currently in memory.                  \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Speciespartition"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Speciespartition                                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Defines a partition of tips into species. The format for the speciespartition \n");
        YvyraPrint ("   command is                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      Speciespartition <name> = <species name>:<taxon list> ,...,<sp nm>:<tx lst>\n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The command enumerates comma separated list of pairs consisting of 'species   \n");
        YvyraPrint ("   name' and 'taxon list'. The 'taxon list' is a standard taxon list, as used by \n");
        YvyraPrint ("   the 'Taxset' command. This means that you can use either the index or the name\n");
        YvyraPrint ("   of a sequence ('taxon'). Ranges are specified using a dash, and a period can  \n");
        YvyraPrint ("   be used as a synonym of the last sequence in the matrix.                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example: speciespartition species = SpeciesA: 1, SpeciesB: 2-.            \n");
        YvyraPrint ("   Here, we name two species. SpeciesA is represented by a single sequence while \n");
        YvyraPrint ("   SpeciesB is represented by all remaining sequences in the matrix.             \n");
        YvyraPrint ("   Each sequence is specified by its row index in the data matrix.               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   As with ordinary partitioning you may define multiple species partitioning    \n");
        YvyraPrint ("   scheme. You have to use command 'set speciespartition' to enable use of one of\n");
        YvyraPrint ("   them.                                                                         \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Currently defined Speciespartitions:                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Number  Speciespartition name        Number of species                        \n");
        YvyraPrint ("   --------------------------------------------------------------------------    \n");
        for (i=0; i<numDefinedSpeciespartitions; i++)
            {
            tempInt=0;
            for (j=0; j<numTaxa; j++)
                {
                if (tempInt < speciespartitionId[j][i])
                    tempInt = speciespartitionId[j][i];
                }
            YvyraPrint ("   %4d    %-24.24s   %4d",i+1, speciespartitionNames[i], tempInt);
            YvyraPrint ("\n");
            }
        YvyraPrint ("                                                                                \n");
        YvyraPrint ("   --------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Constraint"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Constraint                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command defines a tree constraint. The format for the constraint         \n");
        YvyraPrint ("   command is                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      constraint <name> [hard|negative|partial] = <taxon list> [:<taxon list>]   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   There are three types of constraint implemented in yvyra. The type of the   \n");
        YvyraPrint ("   constraint is specified by using one of the three keywords 'hard', 'negative',\n");
        YvyraPrint ("   or 'partial' right after the name of the constraint. If no type is specified, \n");
        YvyraPrint ("   then the constraint is assumed to be 'hard'.                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   In a rooted tree, a 'hard' constraint forces the taxa in the list to form a   \n");
        YvyraPrint ("   monophyletic group. In an unrooted tree, the taxon split that separates the   \n");
        YvyraPrint ("   taxa in the list from other taxa is forced to be present. The interpretation  \n");
        YvyraPrint ("   of this depends on whether the tree is rooted on a taxon outside the list or  \n");
        YvyraPrint ("   a taxon in the list. If the outgroup is excluded , the taxa in the list are   \n");
        YvyraPrint ("   assumed to form a monophyletic group, but if the outgroup is included, the    \n");
        YvyraPrint ("   taxa that are not in the list are forced together.                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A 'negative' constraint bans all the trees that have the listed taxa in the   \n");
        YvyraPrint ("   same subtree. In other words, it is the opposite of a hard constraint.        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A 'partial' or backbone constraint is defined in terms of two sets of taxa    \n");
        YvyraPrint ("   separated by a colon character. The constraint forces all taxa in the first   \n");
        YvyraPrint ("   list to form a monophyletic group that does not include any taxon in the      \n");
        YvyraPrint ("   second list. Taxa that are not included in either list can be placed in any   \n");
        YvyraPrint ("   position on the tree, either inside or outside the constrained group. In an   \n");
        YvyraPrint ("   unrooted tree, the two taxon lists can be switched with each other with no    \n");
        YvyraPrint ("   effect. For a rooted tree, it is the taxa in the first list that have to be   \n");
        YvyraPrint ("   monophyletic, that is, these taxa must share a common ancestor not shared with\n");
        YvyraPrint ("   any taxon in the second list. The taxa in the second list may or may not fall \n");
        YvyraPrint ("   in a monophyletic group depending on the rooting of the tree.                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A list of taxa can be specified using a taxset, taxon names, taxon numbers, or\n");
        YvyraPrint ("   any combination of the above, separated by spaces. The constraint is treated  \n");
        YvyraPrint ("   as an absolute requirement of trees, that is, trees that are not compatible   \n");
        YvyraPrint ("   with the constraint have zero prior (and hence zero posterior) probability.   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   If you are interested in inferring ancestral states for a particular node,    \n");
        YvyraPrint ("   you need to 'hard' constrain that node first using the 'constraint' command.  \n");
        YvyraPrint ("   The same applies if you wish to calibrate an interior node in a dated         \n");
        YvyraPrint ("   analysis. For more information on how to infer ancestral states, see the help \n");
        YvyraPrint ("   for the 'report' command. For more on dating, see the 'calibrate' command.    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   It is important to note that simply defining a constraint using this          \n");
        YvyraPrint ("   command is not sufficient for the program to actually implement the           \n");
        YvyraPrint ("   constraint in an analysis. You must also enforce the constraints using        \n");
        YvyraPrint ("   'prset topologypr = constraints (<list of constraints>)'. For more infor-     \n");
        YvyraPrint ("   mation on this, see the help on the 'prset' command.                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Examples:                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      constraint myclade = Homo Pan Gorilla                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Defines a hard constraint forcing Homo, Pan, and Gorilla to form a mono-      \n");
        YvyraPrint ("   phyletic group or a split that does not include any other taxa.               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      constraint forbiddenclade negative = Homo Pan Gorilla                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Defines a negative constraint that associates all trees where Homo, Pan, and  \n");
        YvyraPrint ("   Gorilla form a monophyletic group with zero posterior probability. In other   \n");
        YvyraPrint ("   words, such trees will not be sampled during MCMC.                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      constraint backbone partial = Homo Gorilla : Mus                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Defines a partial constraint that keeps Mus outside of the clade defined by   \n");
        YvyraPrint ("   the most recent common ancestor of Homo and Gorilla. Other taxa are allowed to\n");
        YvyraPrint ("   sit anywhere in the tree. Note that this particular constraint is meaningless \n");
        YvyraPrint ("   in unrooted trees. yvyra does not assume anything about the position of the \n");
        YvyraPrint ("   outgroup unless it is explicitly included in the partial constraint. Therefore\n");
        YvyraPrint ("   a partial constraint must have at least two taxa on each side of the ':' to be\n");
        YvyraPrint ("   useful in analyses of unrooted trees. The case is different for rooted trees, \n");
        YvyraPrint ("   where it is sufficient for a partial constraint to have more than one taxon   \n");
        YvyraPrint ("   before the ':', as in the example given above, to constrain tree space.       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   To define a more complex constraint tree, simply combine constraints into a   \n");
        YvyraPrint ("   list when issuing the 'prset topologypr' command.                             \n");
        YvyraPrint ("                                                                                 \n");
        if (numDefinedConstraints > 0)
            {
            YvyraPrint ("   Currently defined constraints:                                                \n");
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   Number  Constraint name          type      Number of taxa in[:out]            \n");
            YvyraPrint ("   --------------------------------------------------------------------------    \n");       
            }
        for (i=0; i<numDefinedConstraints; i++)
            {
            strncpy (tempString, constraintNames[i], 22);
            YvyraPrint ("   %4d    %-22.22s   ",i+1, tempString);
            if (definedConstraintsType[i] == HARD)
                YvyraPrint ("hard      ");
            else if (definedConstraintsType[i] == PARTIAL)
                YvyraPrint ("partial   ");
            else
                {
                assert (definedConstraintsType[i] == NEGATIVE);
                YvyraPrint ("negative  ");
                }
            k = NumBits (definedConstraint[i], numTaxa/nBitsInALong + 1);
            YvyraPrint ("%d", k);
            if (definedConstraintsType[i] == PARTIAL)
                {
                k = NumBits (definedConstraintTwo[i], numTaxa/nBitsInALong + 1);
                YvyraPrint (":%d", k);
                }
            YvyraPrint ("\n");
            }
        YvyraPrint ("                                                                                \n");
        YvyraPrint ("   --------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Calibrate"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Calibrate                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command dates a terminal or interior node in the tree. The format is     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      calibrate <node_name> = <age_prior>                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   where <node_name> is the name of a defined interior constraint node or the    \n");
        YvyraPrint ("   name of a terminal node (tip) and <age_prior> is a prior probability distribu-\n");
        YvyraPrint ("   tion on the age of the node. The latter can either be a fixed date or a date  \n");
        YvyraPrint ("   drawn from one of the available prior probability distributions. In general,  \n");
        YvyraPrint ("   the available prior probability distributions are parameterized in terms of   \n");
        YvyraPrint ("   the expected mean age of the distribution to facilitate for users. Some dis-  \n");
        YvyraPrint ("   tributions put a positive probability on all ages above 0.0, while others in- \n");
        YvyraPrint ("   clude a minimum-age constraint and sometimes a maximum-age constraint. The    \n");
        YvyraPrint ("   available distributions and their parameters are:                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      calibrate <node_name> = fixed(<age>)                                       \n");
        YvyraPrint ("      calibrate <node_name> = uniform(<min_age>,<max_age>)                       \n");
        YvyraPrint ("      calibrate <node_name> = offsetexponential(<min_age>,<mean_age>)            \n");
        YvyraPrint ("      calibrate <node_name> = truncatednormal(<min_age>,<mean_age>,<stdev>)      \n");
        YvyraPrint ("      calibrate <node_name> = lognormal(<mean_age>,<stdev>)                      \n");
        YvyraPrint ("      calibrate <node_name> = offsetlognormal(<min_age>,<mean_age>,<stdev>)      \n");
        YvyraPrint ("      calibrate <node_name> = gamma(<mean_age>,<stdev>)                          \n");
        YvyraPrint ("      calibrate <node_name> = offsetgamma(<min_age>,<mean_age>,<stdev>)          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that mean_age is always the mean age and stdev the standard deviation of \n");
        YvyraPrint ("   the distribution measured in user-defined time units. This way of specifying  \n");
        YvyraPrint ("   the distribution parameters is often different from the parameterization used \n");
        YvyraPrint ("   elsewhere in the program. For instance, the standard parameters of the gamma  \n");
        YvyraPrint ("   distribution used by yvyra are shape (alpha) and rate (beta). If you want   \n");
        YvyraPrint ("   to use the standard parameterization, the conversions are as follows:         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exponential distribution: mean    = 1 / rate                               \n");
        YvyraPrint ("      gamma distribution:       mean    = alpha / beta                           \n");
        YvyraPrint ("                                st.dev. = square_root (alpha / beta^2)           \n");
        YvyraPrint ("      lognormal distribution:   mean    = exp (mean_log + st.dev._log^2/2)       \n");
        YvyraPrint ("                                st.dev. = square_root ((exp (st.dev._log^2) - 1) \n");
        YvyraPrint ("                                         * (exp (2*mean_log + st.dev._log^2))    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The truncated normal distribution is an exception in that the mean_age and    \n");
        YvyraPrint ("   stdev parameters are the mean and standard deviation of the underlying non-   \n");
        YvyraPrint ("   truncated normal distribution. The truncation will cause the modified distri- \n");
        YvyraPrint ("   bution to have a higher mean and lower standard deviation. The magnitude of   \n");
        YvyraPrint ("   that effect depends on how much of the tail of the distribution is removed.   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that previous to version 3.2.2, MrBayes used the standard rate parameter-\n");
        YvyraPrint ("   ization of the offset exponential. This should not cause a problem in most    \n");
        YvyraPrint ("   cases because the old parameterization will result in an error in more recent \n");
        YvyraPrint ("   versions of yvyra, and the likely source of the error is given in the error \n");
        YvyraPrint ("   message.                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For a practical example, assume that we had three fossil terminals named      \n");
        YvyraPrint ("   'FossilA', 'FossilB', and 'FossilC'. Assume further that we want to fix the   \n");
        YvyraPrint ("   age of FossilA to 100.0 million years, we think that FossilB is somewhere     \n");
        YvyraPrint ("   between 100.0 and 200.0 million years old, and that FossilC is at least 300.0 \n");
        YvyraPrint ("   million years old, possibly older but relatively unlikely to be more than     \n");
        YvyraPrint ("   400.0 million years old. Then we might use the commands:                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      calibrate FossilA = fixed(100) FossilB = uniform(100,200)                  \n");
        YvyraPrint ("      calibrate FossilC = offsetexponential(300,400)                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that it is possible to give more than one calibration for each           \n");
        YvyraPrint ("   'calibrate' statement. Thus, 'calibrate FossilA=<setting> FossilB=<setting>'  \n");
        YvyraPrint ("   would be a valid statement.                                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   To actually use the calibrations to obtain dated trees, you also need to set  \n");
        YvyraPrint ("   a clock model using relevant 'brlenspr' and 'nodeagepr' options of the 'prset'\n");
        YvyraPrint ("   command. You may also want to examine the 'clockvarpr' and 'clockratepr' op-  \n");
        YvyraPrint ("   tions. Furthermore, you need to activate the relevant constraint(s) using     \n");
        YvyraPrint ("   'topologypr', if you use any dated interior nodes in the tree.                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   You may wish to remove a calibration from an interior or terminal node, which \n");
        YvyraPrint ("   has previously been calibrated. You can do that using                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      calibrate <node_name> = unconstrained                                      \n");
        YvyraPrint ("                                                                                 \n");
        j = 0;
        for (i=0; i<numTaxa; i++)
            if (tipCalibration[i].prior != unconstrained)
                j++;
        for (i=0; i<numDefinedConstraints; i++)
            if (nodeCalibration[i].prior != unconstrained)
                j++;
        if (j > 0)
            {
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   Currently defined calibrations:                                               \n");
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   Node name                Type       Calibration                               \n");
            YvyraPrint ("   ------------------------------------------------------------------            \n");       
            for (i=0; i<numTaxa+numDefinedConstraints; i++)
                {
                if (i<numTaxa)
                    calibrationPtr = &tipCalibration[i];
                else
                    calibrationPtr = &nodeCalibration[i-numTaxa];
                if (calibrationPtr != NULL && calibrationPtr->prior != unconstrained)
                    {
                    if (i<numTaxa)
                        strncpy (tempString, taxaNames[i], 22);
                    else
                        strncpy (tempString, constraintNames[i-numTaxa], 22);
                    if (i<numTaxa)
                        YvyraPrint ("   %-22.22s   Terminal   %s\n", tempString, calibrationPtr->name);
                    else
                        YvyraPrint ("   %-22.22s   Interior   %s\n", tempString, calibrationPtr->name);
                    }
                }
            }
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmodel"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Showmodel                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the current model settings. The correct usage is           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      showmodel                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   After typing \"showmodel\", the modelling assumptions are shown on a          \n");
        YvyraPrint ("   partition-by-partition basis.                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Execute"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Execute                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command executes a file called <file name>. The correct usage is:        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      execute <file name>                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example,                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      execute replicase.nex                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   would execute the file named \"replicase.nex\". This file must be in the      \n");
        YvyraPrint ("   same directory as the executable.                                             \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Lset"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Lset                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command sets the parameters of the likelihood model. The likelihood      \n");
        YvyraPrint ("   function is the probability of observing the data conditional on the phylo-   \n");
        YvyraPrint ("   genetic model. In order to calculate the likelihood, you must assume a        \n");
        YvyraPrint ("   model of character change. This command lets you tailor the biological        \n");
        YvyraPrint ("   assumptions made in the phylogenetic model. The correct usage is              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      lset <parameter>=<option> ... <parameter>=<option>                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example, \"lset nst=6 rates=gamma\" would set the model to a general      \n");
        YvyraPrint ("   model of DNA substitution (the GTR) with gamma-distributed rate variation     \n");
        YvyraPrint ("   across sites.                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Applyto   -- This option allows you to apply the lset commands to specific    \n");
        YvyraPrint ("                partitions. This command should be the first in the list of      \n");
        YvyraPrint ("                commands specified in lset. Moreover, it only makes sense to     \n");
        YvyraPrint ("                be using this command if the data have been partitioned. A       \n");
        YvyraPrint ("                default partition is set on execution of a matrix. If the data   \n");
        YvyraPrint ("                are homogeneous (i.e., all of the same data type), then this     \n");
        YvyraPrint ("                partition will not subdivide the characters. Up to 30 other      \n");
        YvyraPrint ("                partitions can be defined, and you can switch among them using   \n");
        YvyraPrint ("                \"set partition=<partition name>\". Now, you may want to         \n");
        YvyraPrint ("                specify different models to different partitions of the data.    \n");
        YvyraPrint ("                Applyto allows you to do this. For example, say you have         \n");
        YvyraPrint ("                partitioned the data by codon position, and you want to apply    \n");
        YvyraPrint ("                a nst=2 model to the first two partitions and nst=6 to the       \n");
        YvyraPrint ("                last. This could be implemented in two uses of lset:             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                   lset applyto=(1,2) nst=2                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                   lset applyto=(3) nst=6                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                The first applies the parameters after \"applyto\" to the        \n");
        YvyraPrint ("                first and second partitions. The second lset applies nst=6       \n");
        YvyraPrint ("                to the third partition. You can also use applyto=(all), which    \n");
        YvyraPrint ("                attempts to apply the parameter settings to all of the data      \n");
        YvyraPrint ("                partitions. Importantly, if the option is not consistent with    \n");
        YvyraPrint ("                the data in the partition, the program will not apply the        \n");
        YvyraPrint ("                lset option to that partition.                                   \n");
        YvyraPrint ("   Nucmodel  -- This specifies the general form of the nucleotide substitution   \n");
        YvyraPrint ("                model. The options are \"4by4\" [the standard model of DNA       \n");
        YvyraPrint ("                substitution in which there are only four states (A,C,G,T/U)],   \n");
        YvyraPrint ("                \"doublet\" (a model appropriate for modelling the stem regions  \n");
        YvyraPrint ("                of ribosomal genes where the state space is the 16 doublets of   \n");
        YvyraPrint ("                nucleotides), \"codon\" (the substitution model is expanded      \n");
        YvyraPrint ("                around triplets of nucleotides--a codon), and \"Protein\"        \n");
        YvyraPrint ("                (triplets of nucleotides are translated to amino acids, which    \n");
        YvyraPrint ("                form the basis of the substitution model).                       \n");
        YvyraPrint ("   Nst       -- Sets the number of substitution types: \"1\" constrains all of   \n");
        YvyraPrint ("                the rates to be the same (e.g., a JC69 or F81 model); \"2\" all- \n");
        YvyraPrint ("                ows transitions and transversions to have potentially different  \n");
        YvyraPrint ("                rates (e.g., a K80 or HKY85 model); \"6\" allows all rates to    \n");
        YvyraPrint ("                be different, subject to the constraint of time-reversibility    \n");
        YvyraPrint ("                (e.g., a GTR model). Finally, 'nst' can be set to 'mixed', which \n");
        YvyraPrint ("                results in the Markov chain sampling over the space of all poss- \n");
        YvyraPrint ("                ible reversible substitution models, including the GTR model and \n");
        YvyraPrint ("                all models that can be derived from it model by grouping the six \n");
        YvyraPrint ("                rates in various combinations. This includes all the named models\n");
        YvyraPrint ("                above and a large number of others, with or without name.        \n");
        YvyraPrint ("   Code      -- Enforces the use of a particular genetic code. The default       \n");
        YvyraPrint ("                is the universal code. Other options include \"vertmt\" for      \n");
        YvyraPrint ("                vertebrate mitochondrial, \"invermt\", \"mycoplasma\", \"yeast\", \n");
        YvyraPrint ("                \"ciliate\", \"echinoderm\", \"euplotid\", and \"metmt\" (for    \n");
        YvyraPrint ("                metazoan mitochondrial except vertebrates).                      \n");
        YvyraPrint ("   Ploidy    -- Specifies the ploidy of the organism. Options are \"Haploid\",   \n");
        YvyraPrint ("                \"Diploid\" or \"Zlinked\". This option is used when a coalescent\n");
        YvyraPrint ("                prior is used on trees.                                          \n");
        YvyraPrint ("   Rates     -- Sets the model for among-site rate variation. In general, the    \n");
        YvyraPrint ("                rate at a site is considered to be an unknown random variable.   \n");
        YvyraPrint ("                The valid options are:                                           \n");
        YvyraPrint ("                * equal    -- No rate variation across sites.                    \n");
        YvyraPrint ("                * gamma    -- Gamma-distributed rates across sites. The rate     \n");
        YvyraPrint ("                              at a site is drawn from a gamma distribution.      \n");
        YvyraPrint ("                              The gamma distribution has a single parameter      \n");
        YvyraPrint ("                              that describes how much rates vary.                \n");
        YvyraPrint ("                * lnorm    -- Lognormal-distributed rates across sites. The      \n");
        YvyraPrint ("                              rate at a site is drawn from a lognormal           \n");
        YvyraPrint ("                              distribution, which has a single parameter, sigma  \n");
        YvyraPrint ("                              (SD in log scale) that describes how much rates    \n");
        YvyraPrint ("                              vary (mean in natural scale fixed to 1.0).         \n");
        YvyraPrint ("                * adgamma  -- Autocorrelated rates across sites. The marg-       \n");
        YvyraPrint ("                              inal rate distribution is gamma, but adjacent      \n");
        YvyraPrint ("                              sites have correlated rates.                       \n");
        YvyraPrint ("                * propinv  -- A proportion of the sites are invariable.          \n");
        YvyraPrint ("                * invgamma -- A proportion of the sites are invariable while     \n");
        YvyraPrint ("                              the rate for the remaining sites are drawn from    \n");
        YvyraPrint ("                              a gamma distribution.                              \n");
        YvyraPrint ("                * kmixture -- Site rates come from a mixture with k categories.  \n");
        YvyraPrint ("                              Category rates are drawn from an ordered flat      \n");
        YvyraPrint ("                              Dirichlet distribution with mean rather than sum   \n");
        YvyraPrint ("                              equal to 1.0.                                      \n");
        YvyraPrint ("                Note that MrBayes versions 2.0 and earlier supported options     \n");
        YvyraPrint ("                that allowed site specific rates (e.g., ssgamma). In versions    \n");
        YvyraPrint ("                3.0 and later, site specific rates are allowed, but set using    \n");
        YvyraPrint ("                the 'prset ratepr' command for each partition.                   \n");
        YvyraPrint ("   Ngammacat -- Sets the number of rate categories for the gamma distribution.   \n");
        YvyraPrint ("                The gamma distribution is continuous. However, it is virtually   \n");
        YvyraPrint ("                impossible to calculate likelihoods under the continuous gamma   \n");
        YvyraPrint ("                distribution. Hence, an approximation to the continuous gamma    \n");
        YvyraPrint ("                is used; the gamma distribution is broken into ncat categories   \n");
        YvyraPrint ("                of equal weight (1/ncat). The mean rate for each category rep-   \n");
        YvyraPrint ("                resents the rate for the entire category. This option allows     \n");
        YvyraPrint ("                you to specify how many rate categories to use when approx-      \n");
        YvyraPrint ("                imating the gamma. The approximation is better as ncat is inc-   \n");
        YvyraPrint ("                reased. In practice, \"ncat=4\" does a reasonable job of         \n");
        YvyraPrint ("                approximating the continuous gamma.                              \n");
        YvyraPrint ("   Nlnormcat -- Used to set the number of discrete categories used for the ap-   \n");
        YvyraPrint ("                proximation of the lognormal distribution, in the same way as    \n");
        YvyraPrint ("                the Ngammacat setting for the discrete gamma approximation.      \n");
        YvyraPrint ("                Default value is 4.                                              \n");
        YvyraPrint ("   Nmixtcat  -- Used to set the number of components in the k-mixture model of   \n");
        YvyraPrint ("                rate variation across sites. Default value is 4.                 \n");
        /* Temporarily disable this because of conflict with likelihood calculators.
           It should be renamed to samplerates when reintroduced.
        YvyraPrint ("   Usegibbs  -- Specifies whether site probabilities under the discrete gamma    \n");
        YvyraPrint ("                model of rate variation across sites will be summed across rate  \n");
        YvyraPrint ("                categories ('Usegibbs=No') or sampled using a Gibbs sampler      \n");
        YvyraPrint ("                ('Usegibbs=Yes'). The Gibbs sampling approach is much faster and \n");
        YvyraPrint ("                requires less memory but the likelihood of the sampled points    \n");
        YvyraPrint ("                will be considerably higher than with the standard approach of   \n");
        YvyraPrint ("                summing probabilities, so you need to be aware of this when com- \n");
        YvyraPrint ("                paring your results with those you obtain with other programs.   \n");
        YvyraPrint ("                Assume that you are using n rate categories in your discrete     \n");
        YvyraPrint ("                gamma distribution. Then the Gibbs approach is up to n times     \n");
        YvyraPrint ("                faster and requires 1/n as much memory as the standard method.   \n");
        YvyraPrint ("                Unfortunately, the state space also becomes larger so the chain  \n");
        YvyraPrint ("                may need more time to converge. The approach should work best    \n");
        YvyraPrint ("                for large trees, where the uncertainty concerning the best rate  \n");
        YvyraPrint ("                category for each site is negligible. Gibbs sampling cannot be   \n");
        YvyraPrint ("                used for the autocorrelated discrete gamma model, for standard   \n");
        YvyraPrint ("                data, or for restriction data. Also, yvyra will not use Gibbs  \n");
        YvyraPrint ("                sampling when you want to infer site rates.                      \n");
        YvyraPrint ("   Gibbsfreq -- Sets the frequency with which the rate categories of the discrete\n");
        YvyraPrint ("                gamma will be Gibbs sampled. In practice, we have found that a   \n");
        YvyraPrint ("                resampling frequency of every 100 MCMC generations works well for\n");
        YvyraPrint ("                reasonably long runs. The more frequent the Gibbs sampling, the  \n");
        YvyraPrint ("                slower the Gibbs sampling approach will be. If you have k rate   \n");
        YvyraPrint ("                categories and Gibbs sample them every n generations, then the   \n");
        YvyraPrint ("                time it takes to complete n generations will roughly be propor-  \n");
        YvyraPrint ("                tional to n+k. Compare this with the traditional approach of     \n");
        YvyraPrint ("                summing across the n rate categories in every generation, which  \n");
        YvyraPrint ("                requires time proportional to n*k. In practice, however, the     \n");
        YvyraPrint ("                speed difference is not quite as large as this.                  \n"); */
        YvyraPrint ("   Nbetacat  -- Sets the number of rate categories for the beta distribution.    \n");
        YvyraPrint ("                A symmetric beta distribution is used to model the stationary    \n");
        YvyraPrint ("                frequencies when morphological data are used. This option        \n");
        YvyraPrint ("                specifies how well the beta distribution will be approximated.   \n");
        YvyraPrint ("   Omegavar  -- Allows the nonsynonymous/synonymous rate ratio (omega) to vary   \n");
        YvyraPrint ("                across codons. Ny98 assumes that there are three classes, with   \n");
        YvyraPrint ("                potentially different omega values (omega1, omega2, omega3):     \n");
        YvyraPrint ("                omega2 = 1; 0 < omega1 < 1; and omega3 > 1. Like the Ny98 model, \n");
        YvyraPrint ("                the M3 model has three omega classes. However, their values are  \n");
        YvyraPrint ("                less constrained, with omega1 < omega2 < omega3. The default     \n");
        YvyraPrint ("                (omegavar = equal) has no variation on omega across sites.       \n");
        YvyraPrint ("   Covarion  -- This forces the use of a covarion-like model of substitution     \n");
        YvyraPrint ("                for nucleotide or amino acid data. The valid options are \"yes\" \n");
        YvyraPrint ("                and \"no\". The covarion model allows the rate at a site to      \n");
        YvyraPrint ("                change over its evolutionary history. Specifically, the site     \n");
        YvyraPrint ("                is either on or off. When it is off, no substitutions are poss-  \n");
        YvyraPrint ("                ible. When the process is on, substitutions occur according to   \n");
        YvyraPrint ("                a specified substitution model (specified using the other        \n");
        YvyraPrint ("                lset options).                                                   \n");
        YvyraPrint ("   Coding    -- This specifies how characters were sampled. If all site patterns \n");
        YvyraPrint ("                had the possibility of being sampled, then \"All\" should be     \n");
        YvyraPrint ("                specified (the default). Otherwise \"Variable\" (only variable   \n");
        YvyraPrint ("                characters had the possibility of being sampled), \"Informative\"\n");
        YvyraPrint ("                (only parsimony informative characters has the possibility of    \n");
        YvyraPrint ("                being sampled), \"Nosingletons\" (characters which are constant  \n");
        YvyraPrint ("                in all but one taxon were not sampled), \"Noabsencesites\" (char-\n");
        YvyraPrint ("                acters for which all taxa were coded as absent were not sampled),\n");
        YvyraPrint ("                \"Nopresencesites\" (characters for which all taxa were coded as \n");
        YvyraPrint ("                present were not sampled). \"All\" works for all data types.     \n");
        YvyraPrint ("                However, the others only work for morphological (All/Variable/   \n");
        YvyraPrint ("                Informative/Nosingletons) or restriction site (All/Variable/     \n");
        YvyraPrint ("                Informative/Nosingletons/Noabsencesites/Nopresencesites/         \n");
        YvyraPrint ("                Nosingletonpresence/Nosingletonabsence) data.                    \n");
        YvyraPrint ("  Statefrmod -- This option allows you to specify whether a \"stationary\"  \n");
        YvyraPrint ("                (= steady state) or a \"directional\" model of evolution should  \n");
        YvyraPrint ("                be used (the option \"mixed\" invokes a reversible jump over     \n");
        YvyraPrint ("                both alternatives). In the stationary (which is the standard)    \n");
        YvyraPrint ("                case, the state frequencies are assumed to be at equilibrium     \n");
        YvyraPrint ("                throughout the tree. If a directional model is chosen, then the  \n");
        YvyraPrint ("                state frequencies at the root are allowed to differ from the     \n");
        YvyraPrint ("                equilibrium frequencies. The directional and mixed models are    \n");
        YvyraPrint ("                currently only implemented for restriction data. Note that       \n");
        YvyraPrint ("                directional evolution means that the rooting of the tree matters.\n");
        YvyraPrint ("                Thus, although the tree is not a clock tree, it will have a root \n");
        YvyraPrint ("                under a directional model. When \"mixed\" is chosen, the chain   \n");
        YvyraPrint ("                samples the stationary state frequency model, with statefrmod=0  \n");
        YvyraPrint ("                indicating the stationary model and statefrmod=1 indicating the  \n");
        YvyraPrint ("                directional model.                                               \n");
        YvyraPrint ("   Parsmodel -- This forces calculation under the so-called parsimony model      \n");
        YvyraPrint ("                described by Tuffley and Steel (1998). The options are \"yes\"   \n");
        YvyraPrint ("                or \"no\". Note that the biological assumptions of this model    \n");
        YvyraPrint ("                are anything but parsimonious. In fact, this model assumes many  \n");
        YvyraPrint ("                more parameters than the next most complicated model implemented \n");
        YvyraPrint ("                in this program. If you really believe that the parsimony model  \n");
        YvyraPrint ("                makes the biological assumptions described by Tuffley and Steel, \n");
        YvyraPrint ("                then the parsimony method is miss-named.                         \n");
    /*  YvyraPrint ("   Augment   -- This allows the chain to consider the missing entries of         \n");
        YvyraPrint ("                the data matrix as random variables. A Gibbs sampler is          \n");
        YvyraPrint ("                used to sample states.                                           \n"); */
        YvyraPrint ("                                                                                 \n");
        if (numCurrentDivisions == 0)
            tempInt = 1;
        else
            tempInt = numCurrentDivisions;
        for (i=0; i<tempInt; i++)
            {
            if (numCurrentDivisions == 0)
                {
                YvyraPrint ("   Default model settings:                                                       \n");
                mp = &defaultModel;
                }
            else
                {
                YvyraPrint ("   Model settings for partition %d:                                              \n", i+1);
                mp = &modelParams[i];
                }
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   Parameter    Options                               Current Setting            \n");
            YvyraPrint ("   ------------------------------------------------------------------            \n");       
            YvyraPrint ("   Nucmodel     4by4/Doublet/Codon/Protein              %s                       \n", mp->nucModel);
            YvyraPrint ("   Nst          1/2/6/Mixed                             %s                       \n", mp->nst);
            YvyraPrint ("   Code         Universal/Vertmt/Invermt/Yeast/Mycoplasma/                       \n");
            YvyraPrint ("                Ciliate/Echinoderm/Euplotid/Metmt       %s                       \n", mp->geneticCode);
            YvyraPrint ("   Ploidy       Haploid/Diploid/Zlinked                 %s                       \n", mp->ploidy);
            YvyraPrint ("   Rates        Equal/Gamma/LNorm/Propinv/                                       \n");
            YvyraPrint ("                Invgamma/Adgamma/Kmixture               %s                       \n", mp->ratesModel);
            YvyraPrint ("   Ngammacat    <number>                                %d                       \n", mp->numGammaCats);
            YvyraPrint ("   Nlnormcat    <number>                                %d                       \n", mp->numLnormCats);
            YvyraPrint ("   Nmixtcat     <number>                                %d                       \n", mp->numMixtCats);
#if 0
/* Temporarily disable this because of conflict with likelihood calculators. It should be renamed to samplerates when reintroduced. */
            YvyraPrint ("   Usegibbs     Yes/No                                  %s                       \n", mp->useGibbs);
            YvyraPrint ("   Gibbsfreq    <number>                                %d                       \n", mp->gibbsFreq);
#endif
            YvyraPrint ("   Nbetacat     <number>                                %d                       \n", mp->numBetaCats);
            YvyraPrint ("   Omegavar     Equal/Ny98/M3                           %s                       \n", mp->omegaVar);
            YvyraPrint ("   Covarion     No/Yes                                  %s                       \n", mp->covarionModel);
            YvyraPrint ("   Coding       All/Variable/Informative/Nosingletons                            \n");
            YvyraPrint ("                Noabsencesites/Nopresencesites/                                  \n");
            YvyraPrint ("                Nosingletonabsence/Nosingletonpresence  %s                       \n", mp->codingString);
            YvyraPrint ("   Statefrmod   Stationary/Directional/Mixed            %s                       \n", mp->statefreqModel); //SK
            YvyraPrint ("   Parsmodel    No/Yes                                  %s                       \n", mp->parsModel);
        /*  YvyraPrint ("   Augment      No/Yes                                  %s                       \n", mp->augmentData); */
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   ------------------------------------------------------------------            \n");
            }
        }
    else if (!strcmp(helpTkn, "Prset"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Prset                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command sets the priors for the phylogenetic model. Remember that        \n");
        YvyraPrint ("   in a Bayesian analysis, you must specify a prior probability distribution     \n");
        YvyraPrint ("   for the parameters of the likelihood model. The prior distribution rep-       \n");
        YvyraPrint ("   resents your prior beliefs about the parameter before observation of the      \n");
        YvyraPrint ("   data. This command allows you to tailor your prior assumptions to a large     \n");
        YvyraPrint ("   extent.                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Applyto       -- This option allows you to apply the prset commands to        \n");
        YvyraPrint ("                    specific partitions. This command should be the first        \n");
        YvyraPrint ("                    in the list of commands specified in prset. Moreover, it     \n");
        YvyraPrint ("                    only makes sense to be using this command if the data        \n");
        YvyraPrint ("                    have been partitioned. A default partition is set on         \n");
        YvyraPrint ("                    execution of a matrix. If the data are homogeneous           \n");
        YvyraPrint ("                    (i.e., all of the same data type), then this partition       \n");
        YvyraPrint ("                    will not subdivide the characters. Up to 30 other part-      \n");
        YvyraPrint ("                    itions can be defined, and you can switch among them using   \n");
        YvyraPrint ("                    \"set partition=<partition name>\". Now, you may want to     \n");
        YvyraPrint ("                    specify different priors to different partitions of the      \n");
        YvyraPrint ("                    data. Applyto allows you to do this. For example, say        \n");
        YvyraPrint ("                    you have partitioned the data by codon position, and         \n");
        YvyraPrint ("                    you want to fix the statefreqs to equal for the first two    \n");
        YvyraPrint ("                    partitions but apply a flat Dirichlet prior to the state-    \n");
        YvyraPrint ("                    freqs of the last. This could be implemented in two uses of  \n");
        YvyraPrint ("                    prset:                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset applyto=(1,2) statefreqs=fixed(equal)               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset applyto=(3) statefreqs=dirichlet(1,1,1,1)           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The first applies the parameters after \"applyto\"           \n");
        YvyraPrint ("                    to the first and second partitions. The second prset         \n");
        YvyraPrint ("                    applies a flat Dirichlet to the third partition. You can     \n");
        YvyraPrint ("                    also use applyto=(all), which attempts to apply the para-    \n");
        YvyraPrint ("                    meter settings to all of the data partitions. Importantly,   \n");
        YvyraPrint ("                    if the option is not consistent with the data in the part-   \n");
        YvyraPrint ("                    ition, the program will not apply the prset option to        \n");
        YvyraPrint ("                    that partition.                                              \n");
        YvyraPrint ("   Tratiopr      -- This parameter sets the prior for the transition/trans-      \n");
        YvyraPrint ("                    version rate ratio (tratio). The options are:                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset tratiopr = beta(<number>, <number>)                 \n");
        YvyraPrint ("                       prset tratiopr = fixed(<number>)                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The program assumes that the transition and transversion     \n");
        YvyraPrint ("                    rates are independent gamma-distributed random variables     \n");
        YvyraPrint ("                    with the same scale parameter when beta is selected. If you  \n");
        YvyraPrint ("                    want a diffuse prior that puts equal emphasis on transition/ \n");
        YvyraPrint ("                    transversion rate ratios above 1.0 and below 1.0, then use a \n");
        YvyraPrint ("                    flat Beta, beta(1,1), which is the default. If you wish to   \n");
        YvyraPrint ("                    concentrate this distribution more in the equal-rates region,\n");
        YvyraPrint ("                    then use a prior of the type beta(x,x), where the magnitude  \n");
        YvyraPrint ("                    of x determines how much the prior is concentrated in the    \n");
        YvyraPrint ("                    equal rates region. For instance, a beta(20,20) puts more    \n");
        YvyraPrint ("                    probability on rate ratios close to 1.0 than a beta(1,1). If \n");
        YvyraPrint ("                    you think it is likely that the transition/transversion rate \n");
        YvyraPrint ("                    ratio is 2.0, you can use a prior of the type beta(2x,x),    \n");
        YvyraPrint ("                    where x determines how strongly the prior is concentrated on \n");
        YvyraPrint ("                    tratio values near 2.0. For instance, a beta(2,1) is much    \n");
        YvyraPrint ("                    more diffuse than a beta(80,40) but both have the expected   \n");
        YvyraPrint ("                    tratio 2.0 in the absence of data. The parameters of the     \n");
        YvyraPrint ("                    Beta can be interpreted as counts: if you have observed x    \n");
        YvyraPrint ("                    transitions and y transversions, then a beta(x+1,y+1) is a   \n");
        YvyraPrint ("                    good representation of this information. The fixed option    \n");
        YvyraPrint ("                    allows you to fix the tratio to a particular value.          \n");
        YvyraPrint ("   Revmatpr      -- This parameter sets the prior for the substitution rates     \n");
        YvyraPrint ("                    of the GTR model for nucleotide data. The options are:       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset revmatpr = dirichlet(<number>,<number>,...,<number>)\n");
        YvyraPrint ("                       prset revmatpr = fixed(<number>,<number>,...,<number>)    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The program assumes that the six substitution rates          \n");
        YvyraPrint ("                    are independent gamma-distributed random variables with the  \n");
        YvyraPrint ("                    same scale parameter when dirichlet is selected. The six     \n");
        YvyraPrint ("                    numbers in brackets each corresponds to a particular substi- \n");
        YvyraPrint ("                    tution type. Together, they determine the shape of the prior.\n");
        YvyraPrint ("                    The six rates are in the order A<->C, A<->G, A<->T, C<->G,   \n");
        YvyraPrint ("                    C<->T, and G<->T. If you want an uninformative prior you can \n");
        YvyraPrint ("                    use dirichlet(1,1,1,1,1,1), also referred to as a 'flat'     \n");
        YvyraPrint ("                    Dirichlet. This is the default setting. If you wish a prior  \n");
        YvyraPrint ("                    where the C<->T rate is 5 times and the A<->G rate 2 times   \n");
        YvyraPrint ("                    higher, on average, than the transversion rates, which are   \n");
        YvyraPrint ("                    all the same, then you should use a prior of the form        \n");
        YvyraPrint ("                    dirichlet(x,2x,x,x,5x,x), where x determines how much the    \n");
        YvyraPrint ("                    prior is focused on these particular rates. For more info,   \n");
        YvyraPrint ("                    see tratiopr. The fixed option allows you to fix the substi- \n");
        YvyraPrint ("                    tution rates to particular values.                           \n");
        YvyraPrint ("   Revratepr     -- This parameter sets the prior for each substitution rate of  \n");
        YvyraPrint ("                    the GTR model subspace when 'nst' is set to 'mixed' (see the \n");
        YvyraPrint ("                    'lset' command). The only option is                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset revratepr = symdir(<number>)                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    which will associate each independent rate in the rate matrix\n");
        YvyraPrint ("                    with a modified symmetric Dirichlet prior, where a singleton \n");
        YvyraPrint ("                    rate has the specified alpha parameter, while a rate that    \n");
        YvyraPrint ("                    applies to n pairwise substitution types has an alpha that is\n");
        YvyraPrint ("                    n times the specified number. The higher the specified num-  \n");
        YvyraPrint ("                    ber, the more focused the prior will be on equal rates. The  \n");
        YvyraPrint ("                    default value is 1, which gives an effect similar to a flat  \n");
        YvyraPrint ("                    Dirichlet.                                                   \n");
        YvyraPrint ("   Aamodelpr     -- This parameter sets the rate matrix for amino acid data.     \n");
        YvyraPrint ("                    You can either fix the model by specifying aamodelpr=fixed   \n");
        YvyraPrint ("                    (<model name>), where <model name> is 'poisson' (a glorified \n");
        YvyraPrint ("                    Jukes-Cantor model), 'jones', 'dayhoff', 'mtrev', 'mtmam',   \n");
        YvyraPrint ("                    'wag', 'rtrev', 'cprev', 'vt', 'blosum', 'lg', 'equalin'     \n");
        YvyraPrint ("                    (a glorified Felsenstein 1981 model), or 'gtr'. You can also \n");
        YvyraPrint ("                    average over the first ten models by specifying aamodelpr=   \n");
        YvyraPrint ("                    mixed. If you do so, the Markov chain will sample each model \n");
        YvyraPrint ("                    according to its probability. The sampled model is reported  \n");
        YvyraPrint ("                    as an index: poisson(0), jones(1), dayhoff(2), mtrev(3),     \n");
        YvyraPrint ("                    mtmam(4), wag(5), rtrev(6), cprev(7), vt(8), or blosum(9).   \n");
        YvyraPrint ("                    The 'Sump' command summarizes the MCMC samples and calculates\n");
        YvyraPrint ("                    the posterior probability estimate for each of these models. \n");
        YvyraPrint ("   Aarevmatpr    -- This parameter sets the prior for the substitution rates     \n");
        YvyraPrint ("                    of the GTR model for amino acid data. The options are:       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset aarevmatpr = dirichlet(<number>,<number>,...,<number>)\n");
        YvyraPrint ("                       prset aarevmatpr = fixed(<number>,<number>,...,<number>)  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The options are the same as those for 'Revmatpr' except that \n");
        YvyraPrint ("                    they are defined over the 190 rates of the time-reversible   \n");
        YvyraPrint ("                    GTR model for amino acids instead of over the 6 rates of the \n");
        YvyraPrint ("                    GTR model for nucleotides. The rates are in the order A<->R, \n");
        YvyraPrint ("                    A<->N, etc to Y<->V. In other words, amino acids are listed  \n");
        YvyraPrint ("                    in alphabetic order based on their full name. The first amino\n");
        YvyraPrint ("                    acid (Alanine) is then combined in turn with all amino acids \n");
        YvyraPrint ("                    following it in the list, starting with amino acid 2 (Argi-  \n");
        YvyraPrint ("                    nine) and finishing with amino acid 20 (Valine). The second  \n");
        YvyraPrint ("                    amino acid (Arginine) is then combined in turn with all amino\n");
        YvyraPrint ("                    acids following it, starting with amino acid 3 (Asparagine)  \n");
        YvyraPrint ("                    and finishing with amino acid 20 (Valine), and so on.        \n");
        YvyraPrint ("   Omegapr       -- This parameter specifies the prior on the nonsynonymous/     \n");
        YvyraPrint ("                    synonymous rate ratio. The options are:                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset omegapr = dirichlet(<number>,<number>)              \n");
        YvyraPrint ("                       prset omegapr = fixed(<number>)                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        YvyraPrint ("                    stitution model is set to codon using the lset command       \n");
        YvyraPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        YvyraPrint ("                    case when there is no variation in omega across sites (i.e., \n");
        YvyraPrint ("                    \"lset omegavar=equal\").                                    \n");
        YvyraPrint ("   Ny98omega1pr  -- This parameter specifies the prior on the nonsynonymous/     \n");
        YvyraPrint ("                    synonymous rate ratio for sites under purifying selection.   \n");
        YvyraPrint ("                    The options are:                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset Ny98omega1pr = beta(<number>,<number>)              \n");
        YvyraPrint ("                       prset Ny98omega1pr = fixed(<number>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        YvyraPrint ("                    stitution model is set to codon using the lset command       \n");
        YvyraPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        YvyraPrint ("                    case where omega varies across sites using the model of      \n");
        YvyraPrint ("                    Nielsen and Yang (1998) (i.e., \"lset omegavar=ny98\"). If   \n");
        YvyraPrint ("                    fixing the parameter, you must specify a number between      \n");
        YvyraPrint ("                    0 and 1.                                                     \n");
        YvyraPrint ("   Ny98omega3pr  -- This parameter specifies the prior on the nonsynonymous/     \n");
        YvyraPrint ("                    synonymous rate ratio for positively selected sites. The     \n");
        YvyraPrint ("                    options are:                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset Ny98omega3pr = uniform(<number>,<number>)           \n");
        YvyraPrint ("                       prset Ny98omega3pr = exponential(<number>)                \n");
        YvyraPrint ("                       prset Ny98omega3pr = fixed(<number>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        YvyraPrint ("                    stitution model is set to codon using the lset command       \n");
        YvyraPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        YvyraPrint ("                    case where omega varies across sites according to the        \n");
        YvyraPrint ("                    NY98 model. Note that if the NY98 model is specified         \n");
        YvyraPrint ("                    that this parameter must be greater than 1, so you should    \n");
        YvyraPrint ("                    not specify a uniform(0,10) prior, for example.              \n");
        YvyraPrint ("   M3omegapr     -- This parameter specifies the prior on the nonsynonymous/     \n");
        YvyraPrint ("                    synonymous rate ratios for all three classes of sites for    \n");
        YvyraPrint ("                    the M3 model. The options are:                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset M3omegapr = exponential                             \n");
        YvyraPrint ("                       prset M3omegapr = fixed(<number>,<number>,<number>)       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        YvyraPrint ("                    stitution model is set to codon using the lset command       \n");
        YvyraPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        YvyraPrint ("                    case where omega varies across sites using the M3 model of   \n");
        YvyraPrint ("                    Yang et al. (2000) (i.e., \"lset omegavar=M3\"). Under the   \n");
        YvyraPrint ("                    exponential prior, the four rates (dN1, dN2, dN3, and dS)    \n");
        YvyraPrint ("                    are all considered to be independent draws from the same     \n");
        YvyraPrint ("                    exponential distribution (the parameter of the exponential   \n");
        YvyraPrint ("                    does not matter, and so you don't need to specify it). The   \n");
        YvyraPrint ("                    rates dN1, dN2, and dN3 are taken to be the order statistics \n");
        YvyraPrint ("                    with dN1 < dN2 < dN3. These three rates are all scaled to    \n");
        YvyraPrint ("                    the same synonymous rate, dS. The other option is to simply  \n");
        YvyraPrint ("                    fix the three rate ratios to some values.                    \n");
        YvyraPrint ("   Codoncatfreqs -- This parameter specifies the prior on frequencies of sites   \n");
        YvyraPrint ("                    under purifying, neutral, and positive selection. The        \n");
        YvyraPrint ("                    options are:                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset codoncatfreqs = dirichlet(<num>,<num>,<num>)        \n");
        YvyraPrint ("                       prset codoncatfreqs = fixed(<number>,<number>,<number>)   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        YvyraPrint ("                    stitution model is set to codon using the lset command       \n");
        YvyraPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        YvyraPrint ("                    case where omega varies across sites using the models of     \n");
        YvyraPrint ("                    Nielsen and Yang (1998) (i.e., \"lset omegavar=ny98\")       \n");
        YvyraPrint ("                    or Yang et al. (2000) (i.e., \"lset omegavar=M3\")           \n");
        YvyraPrint ("                    Note that the sum of the three frequencies must be 1.        \n");
        YvyraPrint ("   Statefreqpr   -- This parameter specifies the prior on the state freq-        \n");
        YvyraPrint ("                    uencies. The options are:                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset statefreqpr = dirichlet(<number>)                   \n");
        YvyraPrint ("                       prset statefreqpr = dirichlet(<number>,...,<number>)      \n");
        YvyraPrint ("                       prset statefreqpr = fixed(equal)                          \n");
        YvyraPrint ("                       prset statefreqpr = fixed(empirical)                      \n");
        YvyraPrint ("                       prset statefreqpr = fixed(<number>,...,<number>)          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    For the dirichlet, you can specify either a single number    \n");
        YvyraPrint ("                    or as many numbers as there are states. If you specify a     \n");
        YvyraPrint ("                    single number, then the prior has all states equally         \n");
        YvyraPrint ("                    probable with a variance related to the single parameter     \n");
        YvyraPrint ("                    passed in.                                                   \n");
        YvyraPrint ("   Rootfreqpr    -- This prior is only available when the \"Directional\" model  \n");
        YvyraPrint ("                    was chosen as the Statefrmod in \"lset\". It specifies the   \n");
        YvyraPrint ("                    prior on the state freuencies at the root, in contrast to    \n");
        YvyraPrint ("                    the equilibrium state frequencies. The options are:          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset rootfreqpr = dirichlet(<number>)                    \n");
        YvyraPrint ("                       prset rootfreqpr = dirichlet(<number>,...,<number>)       \n");
        YvyraPrint ("                       prset rootfreqpr = fixed(<number>,...,<number>)           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    For the Dirichlet, you can specify either a single number    \n");
        YvyraPrint ("                    or as many numbers as there are states. If you specify a     \n");
        YvyraPrint ("                    single number, then the prior has all states equally         \n");
        YvyraPrint ("                    probable with a variance related to the single parameter     \n");
        YvyraPrint ("                    passed in.                                                   \n");
        YvyraPrint ("   Shapepr       -- This parameter specifies the prior for the gamma/lnorm shape \n");
        YvyraPrint ("                    parameter for among-site rate variation. The options are:    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset shapepr = uniform(<number>,<number>)                \n");
        YvyraPrint ("                       prset shapepr = exponential(<number>)                     \n");
        YvyraPrint ("                       prset shapepr = fixed(<number>)                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Pinvarpr      -- This parameter specifies the prior for the proportion of     \n");
        YvyraPrint ("                    invariable sites. The options are:                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset pinvarpr = uniform(<number>,<number>)               \n");
        YvyraPrint ("                       prset pinvarpr = fixed(<number>)                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    Note that the valid range for the parameter is between 0     \n");
        YvyraPrint ("                    and 1. Hence, \"prset pinvarpr=uniform(0,0.8)\" is valid     \n");
        YvyraPrint ("                    while \"prset pinvarpr=uniform(0,10)\" is not. The def-      \n");
        YvyraPrint ("                    ault setting is \"prset pinvarpr=uniform(0,1)\".             \n");
        YvyraPrint ("   Ratecorrpr    -- This parameter specifies the prior for the autocorrelation   \n");
        YvyraPrint ("                    parameter of the autocorrelated gamma distribution for       \n");
        YvyraPrint ("                    among-site rate variation. The options are:                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset ratecorrpr = uniform(<number>,<number>)             \n");
        YvyraPrint ("                       prset ratecorrpr = fixed(<number>)                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    Note that the valid range for the parameter is between -1    \n");
        YvyraPrint ("                    and 1. Hence, \"prset ratecorrpr=uniform(-1,1)\" is valid    \n");
        YvyraPrint ("                    while \"prset ratecorrpr=uniform(-11,10)\" is not. The       \n");
        YvyraPrint ("                    default setting is \"prset ratecorrpr=uniform(-1,1)\".       \n");
        YvyraPrint ("   Covswitchpr   -- This option sets the prior for the covarion switching        \n");
        YvyraPrint ("                    rates. The options are:                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset covswitchpr = uniform(<number>,<number>)            \n");
        YvyraPrint ("                       prset covswitchpr = exponential(<number>)                 \n");
        YvyraPrint ("                       prset covswitchpr = fixed(<number>,<number>)              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The covarion model has two rates: a rate from on to off      \n");
        YvyraPrint ("                    and a rate from off to on. The rates are assumed to have     \n");
        YvyraPrint ("                    independent priors that individually are either uniformly    \n");
        YvyraPrint ("                    or exponentially distributed. The other option is to         \n");
        YvyraPrint ("                    fix the switching rates, in which case you must specify      \n");
        YvyraPrint ("                    both rates. (The first number is off->on and the second      \n");
        YvyraPrint ("                    is on->off).                                                 \n");
        YvyraPrint ("   Symdirihyperpr - This option sets the prior for the stationary frequencies    \n");
        YvyraPrint ("                    of the states for morphological (standard) data. The         \n");
        YvyraPrint ("                    labelling of the states is somewhat arbitrary. For example,  \n");
        YvyraPrint ("                    the state \"1\" for different characters does not have the   \n");
        YvyraPrint ("                    same meaning. This is not true for DNA characters, for ex-   \n");
        YvyraPrint ("                    ample, where a \"G\" has the same meaning across characters. \n");
        YvyraPrint ("                    The fact that the labelling of morphological characters is   \n");
        YvyraPrint ("                    arbitrary makes it difficult to allow unequal character-     \n");
        YvyraPrint ("                    state frequencies. yvyra gets around this problem by       \n");
        YvyraPrint ("                    assuming that the states have a symmetric Dirichlet prior    \n");
        YvyraPrint ("                    (i.e. all Dirichlet parameters are equal). The variation in  \n");
        YvyraPrint ("                    the Dirichlet can be controlled by this parameter.           \n");
        YvyraPrint ("                    Symdirihyperpr specifies the distribution on the parameter   \n");
        YvyraPrint ("                    of the symmetric Dirichlet. The valid options are:           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset Symdirihyperpr = uniform(<number>,<number>)         \n");
        YvyraPrint ("                       prset Symdirihyperpr = exponential(<number>)              \n");
        YvyraPrint ("                       prset Symdirihyperpr = fixed(<number>)                    \n");
        YvyraPrint ("                       prset Symdirihyperpr = fixed(infinity)                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    If \"fixed(infinity)\" is chosen, the Dirichlet prior is     \n");
        YvyraPrint ("                    fixed such that all character states have equal frequency.   \n");
        YvyraPrint ("   Topologypr    -- This parameter specifies the prior probabilities of          \n");
        YvyraPrint ("                    phylogenies. The options are:                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset topologypr = uniform                                \n");
        YvyraPrint ("                       prset topologypr = speciestree                            \n");
        YvyraPrint ("                       prset topologypr = constraints(<list>)                    \n");
        YvyraPrint ("                       prset topologypr = fixed(<treename>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    If the prior is selected to be \"uniform\", the default,     \n");
        YvyraPrint ("                    then all possible trees are considered a priori equally      \n");
        YvyraPrint ("                    probable. The 'speciestree' option is used when the topology \n");
        YvyraPrint ("                    is constrained to fold inside a species tree together with   \n");
        YvyraPrint ("                    other (gene) trees. The constraints option allows you to     \n");
        YvyraPrint ("                    specify complicated prior probabilities on trees (constraints\n");
        YvyraPrint ("                    are discussed more fully in \"help constraint\"). Note that  \n");
        YvyraPrint ("                    you must specify a list of constraints that you wish to be   \n");
        YvyraPrint ("                    obeyed. The list can be either the constraints' name or      \n");
        YvyraPrint ("                    number. Finally, you can fix the topology to that of a user  \n");
        YvyraPrint ("                    tree defined in a trees block. Branch lengths will still be  \n");
        YvyraPrint ("                    sampled as usual on the fixed topology.                      \n");
        YvyraPrint ("   Brlenspr      -- This parameter specifies the prior probability dist-         \n");
        YvyraPrint ("                    ribution on branch lengths. The options are specified using: \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset brlenspr = <setting>                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    where <setting> is one of                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       unconstrained:uniform(<num>,<num>)                        \n");
        YvyraPrint ("                       unconstrained:exponential(<number>)                       \n");
        YvyraPrint ("                       unconstrained:twoexp(<num>,<num>)                         \n");
        YvyraPrint ("                       unconstrained:gammadir(<num>,<num>,<num>,<num>)           \n");
        YvyraPrint ("                       unconstrained:invgamdir(<num>,<num>,<num>,<num>)          \n");
        YvyraPrint ("                       clock:uniform                                             \n");
        YvyraPrint ("                       clock:birthdeath                                          \n");
        YvyraPrint ("                       clock:coalescence                                         \n");
        YvyraPrint ("                       clock:fossilization                                       \n");
        YvyraPrint ("                       clock:speciestree                                         \n");
        YvyraPrint ("                       fixed(<treename>)                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    Trees with unconstrained branch lengths are unrooted         \n");
        YvyraPrint ("                    whereas clock-constrained trees are rooted. The option       \n");
        YvyraPrint ("                    after the colon specifies the details of the probability     \n");
        YvyraPrint ("                    density of branch lengths. If you choose a birth-death       \n");
        YvyraPrint ("                    or coalescence prior, you may want to modify the details     \n");
        YvyraPrint ("                    of the parameters of those processes (speciation rate,       \n");
        YvyraPrint ("                    extinction rate and sample probability for the birth-death   \n");
        YvyraPrint ("                    prior; population size and clock rate parameter for the      \n");
        YvyraPrint ("                    coalescence prior). When gene trees are constrained to fold  \n");
        YvyraPrint ("                    inside species trees, the appropriate branch length prior is \n");
        YvyraPrint ("                    'clock:speciestree'. Under this model, it is possible to     \n");
        YvyraPrint ("                    control whether the population size is constant or variable  \n");
        YvyraPrint ("                    across the species tree using the 'popvarpr' setting.        \n");
        YvyraPrint ("                    Branch lengths can also be fixed but only if the topology is \n");
        YvyraPrint ("                    fixed.                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    For unconstrained branch lengths, yvyra offers five alter- \n");
        YvyraPrint ("                    native prior distributions. The first two are the simple     \n");
        YvyraPrint ("                    'uniform' and 'exponential' priors. The 'uniform' prior takes\n");
        YvyraPrint ("                    two parameters, the lower and upper bound of the uniform dis-\n");
        YvyraPrint ("                    tribution, respectively. The 'exponential' prior takes a sin-\n");
        YvyraPrint ("                    gle parameter, the rate of the exponential distribution. The \n");
        YvyraPrint ("                    mean of the exponential distribution is the inverse of the   \n");
        YvyraPrint ("                    rate. For instance, an 'exp(10)' distribution has an expected\n");
        YvyraPrint ("                    mean of 0.1.                                                 \n");
        YvyraPrint ("                    yvyra also offers three more complex prior distributions   \n");
        YvyraPrint ("                    on unconstrained branch lengths. The two-exponential prior   \n");
        YvyraPrint ("                    (Yang and Rannala 2005; Yang 2007) uses two different expo-  \n"); 
        YvyraPrint ("                    nential distributions, one for internal and one for external \n");
        YvyraPrint ("                    branch lengths. The two-exponential prior is invoked using   \n");
        YvyraPrint ("                    'twoexp(<r_I>,<r_E>)', where '<r_I>' is a number specifying  \n");
        YvyraPrint ("                    the rate of the exponential distribution on internal branch  \n");
        YvyraPrint ("                    lengths, while '<r_E>' is the rate for external branch       \n");
        YvyraPrint ("                    lengths. The prior mean for internal branch lengths is then  \n");
        YvyraPrint ("                    1/r_I, and for external ones is 1/r_E. For instance, to set  \n");
        YvyraPrint ("                    prior mean of internal branch lengths to 0.01, and external  \n");
        YvyraPrint ("                    ones to 0.1, use 'twoexp(100,10)'.                           \n");
        YvyraPrint ("                    The setting 'twoexp(10,10)' is equivalent to 'exp(10)'.      \n");
        YvyraPrint ("                    The compound Dirichlet priors 'gammadir(<a_T>,<b_T>,<a>,<c>)'\n");
        YvyraPrint ("                    and 'invgamdir(<a_T>,<b_T>,<a>,<c>)' specify a fairly diffuse\n");
        YvyraPrint ("                    prior on tree length 'T', and then partition the tree length \n");
        YvyraPrint ("                    into branch lengths according to a Dirichlet distribution    \n");
        YvyraPrint ("                    (Rannala et al. 2012). If 'T' is considered drawn from a     \n");
        YvyraPrint ("                    gamma distribution with parameters a_T and b_T, and with mean\n");
        YvyraPrint ("                    a_T/b_T, we recommend setting a_T = 1; if it is instead con- \n");
        YvyraPrint ("                    sidered drawn from an inverse gamma (invgamma) distribution  \n");
        YvyraPrint ("                    with parameters a_T and b_T, and with mean b_T/(a_T -1), then\n");
        YvyraPrint ("                    we recommend setting a_T = 3. In the latter case, b_T should \n");
        YvyraPrint ("                    be chosen so that the prior mean of T is reasonable for the  \n");
        YvyraPrint ("                    data. In the former case, setting b_T = 0.1 (corresponding to\n");
        YvyraPrint ("                    a mean tree length of 10) should be appropriate for a wide   \n");
        YvyraPrint ("                    range of tree lengths (at least in the interval 1 to 100).   \n");
        YvyraPrint ("                    The concentration parameter a of the Dirichlet distribution  \n");
        YvyraPrint ("                    is inversely related to the variance of the branch lengths,  \n");
        YvyraPrint ("                    while c is the ratio of the prior means for the internal and \n");
        YvyraPrint ("                    external branch lengths. The default setting, a = c = 1,     \n");
        YvyraPrint ("                    specifies a uniform Dirichlet distribution of branch lengths \n");
        YvyraPrint ("                    given the tree length. For instance, 'gammadir(1,0.1,1,1)'   \n");
        YvyraPrint ("                    specifies a compound Dirichlet prior on branch lengths, where\n");
        YvyraPrint ("                    tree length is associated with a gamma distribution with mean\n");
        YvyraPrint ("                    10, and branch length proportions are associated with a uni- \n");
        YvyraPrint ("                    form Dirichlet distribution (default).                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    For clock trees with calibrated external nodes (fossils),    \n");
        YvyraPrint ("                    yvyra also offers the fossilized birth-death prior:        \n");
        YvyraPrint ("                    'clock:fossilization'.                                       \n");
        YvyraPrint ("                    If 'SampleStrat' is set to 'fossiltip', it assumes that upon \n");
        YvyraPrint ("                    sampling the lineage is dead and won't produce descendants,  \n");
        YvyraPrint ("                    meaning each fossil sample is a tip. If 'SampleStrat' is set \n");
        YvyraPrint ("                    to 'random' (default), fossils are sampled serially along the\n");
        YvyraPrint ("                    birth-death tree (Stadler 2010), so they can be tips or an-  \n");
        YvyraPrint ("                    cestors. See 'Speciationpr', 'Extinctionpr', 'SampleStrat',  \n");
        YvyraPrint ("                    'Fossilizationpr' for more information.                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Treeagepr     -- This parameter specifies the prior probability distribution  \n");
        YvyraPrint ("                    on the tree age when a uniform or fossilization prior is used\n");
        YvyraPrint ("                    on the branch lengths of a clock tree.                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The options are:                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset treeagepr = <setting>                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    where <setting> is one of                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       fixed(<age>)                                              \n");
        YvyraPrint ("                       uniform(<min_age>,<max_age>)                              \n");
        YvyraPrint ("                       offsetexponential(<min_age>,<mean_age>)                   \n");
        YvyraPrint ("                       truncatednormal(<min_age>,<mean_age>,<st.dev.>)           \n");
        YvyraPrint ("                       lognormal(<mean_age>,<st.dev.>)                           \n");
        YvyraPrint ("                       offsetlognormal(<min_age>,<mean_age>,<st.dev.>)           \n");
        YvyraPrint ("                       gamma(<mean_age>,<st.dev.>)                               \n");
        YvyraPrint ("                       offsetgamma(<min_age>,<mean_age>,<st.dev.>)               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    These are the same options used for the 'Calibrate' command. \n");
        YvyraPrint ("                    Note that, unlike elsewhere in MrMayes, we always use the    \n");
        YvyraPrint ("                    mean and standard deviation of the resulting age distribution\n");
        YvyraPrint ("                    rather than the standard parameterization, if different. This\n");
        YvyraPrint ("                    is to facilitate for the users who want to focus on the in-  \n");
        YvyraPrint ("                    formation conveyed about the age. For those who wish to use  \n");
        YvyraPrint ("                    the standard parameterization, there are simple conversions  \n");
        YvyraPrint ("                    between the two. See the 'Calibrate' command for more infor- \n");
        YvyraPrint ("                    mation.                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The tree age is simply the age of the most recent common     \n");
        YvyraPrint ("                    ancestor of the tree. If the clock rate is fixed to 1.0,     \n");
        YvyraPrint ("                    which is the default, the tree age is equivalent to the      \n");
        YvyraPrint ("                    expected number of substitutions from the root to the tip of \n");
        YvyraPrint ("                    the tree, that is, tree height. The tree age prior ensures   \n");
        YvyraPrint ("                    that the joint probability for the uniform prior (or fossil- \n");
        YvyraPrint ("                    ization prior) model of branch lengths on a clock tree is    \n");
        YvyraPrint ("                    proper. The default setting is 'gamma(1,1)'. If the root node\n");
        YvyraPrint ("                    in the tree is calibrated, the root calibration replaces the \n");
        YvyraPrint ("                    tree age prior.                                              \n");
        YvyraPrint ("   Speciationpr  -- This parameter sets the prior on the net diversification rate\n");
        YvyraPrint ("                    i.e., (lambda - mu) in the birth-death model and the general \n");
        YvyraPrint ("                    case of fossilized birth-death (FBD) model; or (lambda - mu -\n");
        YvyraPrint ("                    psi) in the special case of the FBD model (fossiltip).       \n");
        YvyraPrint ("                    Values of this parameter are > 0. Prior options:             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset speciationpr = uniform(<number>,<number>)           \n");
        YvyraPrint ("                       prset speciationpr = exponential(<number>)                \n");
        YvyraPrint ("                       prset speciationpr = fixed(<number>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only relevant if the (fossil) birth-death  \n");
        YvyraPrint ("                    process is selected as the prior on branch lengths.          \n");
        YvyraPrint ("   Extinctionpr  -- This parameter sets the prior on the relative extinction rate\n");
        YvyraPrint ("                    (turnover), i.e., (mu / lambda) in the birth-death model and \n");
        YvyraPrint ("                    the general case of the FBD model; or (mu + psi) / lambda in \n");
        YvyraPrint ("                    the special case of the FBD model (fossiltip).               \n");
        YvyraPrint ("                    Values of this parameter are in range (0,1). Prior options:  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset extinctionpr = beta(<number>,<number>)              \n");
        YvyraPrint ("                       prset extinctionpr = fixed(<number>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only relevant if the (fossil) birth-death  \n");
        YvyraPrint ("                    process is selected as the prior on branch lengths.          \n");
        YvyraPrint (" Fossilizationpr -- This parameter sets the prior on the relative fossilization  \n");
        YvyraPrint ("                    rate, psi/(mu+psi), in the fossilized birth-death model.     \n");
        YvyraPrint ("                    Values of this parameter are in range (0,1). Prior options:  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset fossilizationpr = beta(<number>,<number>)           \n");
        YvyraPrint ("                       prset fossilizationpr = fixed(<number>)                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only relevant if the fossilized birth-death\n");
        YvyraPrint ("                    process is selected as the prior on branch lengths.          \n");
        YvyraPrint ("                    If SampleStrat is used to divide up time intervals, it sets  \n");
        YvyraPrint ("                    the i.i.d. prior for the parameter in each interval.         \n");
        YvyraPrint ("   Sampleprob    -- This parameter sets the fraction of extant taxa that are     \n");
        YvyraPrint ("                    sampled in the analysis. This is used with the birth-death   \n");
        YvyraPrint ("                    prior (Yang & Rannala 1997; Stadler 2009; Hohna et al. 2011),\n");
        YvyraPrint ("                    or the fossilized birth-death prior (Stadler 2010; Heath et  \n");
        YvyraPrint ("                    al. 2014; Zhang et al. 2016) on rooted timetrees.            \n");
        YvyraPrint ("                    Values are in range (0,1], with 1.0 for complete sampling.   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset sampleprob = <number>                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   SampleStrat   -- This parameter sets the strategy under which species were    \n");
        YvyraPrint ("                    sampled in the analysis. For the birth-death prior, 'birth-  \n");
        YvyraPrint ("                    death' (Hohna et al. 2011), three strategies: 'random',      \n");
        YvyraPrint ("                    'diversity' and 'cluster' sampling can be used for extant    \n");
        YvyraPrint ("                    taxa. No extinct sample (fossil) is allowed in this prior.   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    For data with extant and extinct samples, use 'prset brlenspr\n");
        YvyraPrint ("                    =clock:fossilization'. (Stadler 2010; Zhang et al. 2016)     \n");
        YvyraPrint ("                    For the fossilized birth-death prior, 'fossiltip' assumes    \n");
        YvyraPrint ("                    extant taxa are sampled randomly, extinct taxa (fossils) are \n");
        YvyraPrint ("                    sampled with constant rate, and upon sampling the lineage is \n");
        YvyraPrint ("                    dead and won't produce any descendant, so fossils are all at \n");
        YvyraPrint ("                    the tips. Except 'fossiltip', the following strategies allow \n");
        YvyraPrint ("                    fossils also being ancestors of other samples: 'random' (de- \n");
        YvyraPrint ("                    fault) assumes extant taxa are sampled uniformly at random   \n");
        YvyraPrint ("                    with probability rho; 'diversity' assumes extant taxa are    \n");
        YvyraPrint ("                    sampled with proportion rho to maximize diversity (Zhang et  \n");
        YvyraPrint ("                    al. 2016). Value of rho is set in Sampleprob (see above).    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    By default, the rates of speciation, extinction, and fossil  \n");
        YvyraPrint ("                    sampling are constant in the FBD model (Stadler 2010). To    \n");
        YvyraPrint ("                    allow rate variation through time (in a piecewise-constant   \n");
        YvyraPrint ("                    manner), append the number <s>: of shifts and corresponding  \n");
        YvyraPrint ("                    times <t_i> (in descending order) to the sampling strategy ( \n");
        YvyraPrint ("                    'random' or  'diversity'), in the order for relative fossil- \n");
        YvyraPrint ("                    sampling, net diversification and turnover rates. For example\n");
        YvyraPrint ("                      '4: 200 150 100 50, 4: 200 150 100 50, 4: 200 150 100 50'  \n");
        YvyraPrint ("                    specifies four shifts (5 epochs) for each rate, all shifting \n");
        YvyraPrint ("                    at the same time (note that the number and times can also be \n");
        YvyraPrint ("                    different for each rate, although not necessary).            \n");
        YvyraPrint ("                    Alternatively, simply providing only one series, e.g.,       \n");
        YvyraPrint ("                      'prset samplestrat=random 4: 200 150 100 50'               \n");
        YvyraPrint ("                    allows fossil-sampling rate to shift while keeping the       \n");
        YvyraPrint ("                    speciation and extinction rates constant.                    \n");
        YvyraPrint ("                    An implementation note for diversified sampling in the FBD   \n");
        YvyraPrint ("                    model: there is no fossil and no rate shift allowed between  \n");
        YvyraPrint ("                    x_cut (the cut-off time for sampling extant taxa) and the    \n");
        YvyraPrint ("                    present (Zhang et al. 2016).                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset samplestrat = random                                \n");
        YvyraPrint ("                       prset samplestrat = diversity                             \n");
        YvyraPrint ("                       prset samplestrat = cluster                               \n");
        YvyraPrint ("                       prset samplestrat = fossiltip                             \n");
        YvyraPrint ("                       prset samplestrat = random    <s>: ... <t_i> ...          \n");
        YvyraPrint ("                       prset samplestrat = diversity <s>: ... <t_i> ...          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Popsizepr     -- This parameter sets the prior on the population size compo-  \n");
        YvyraPrint ("                    nent of the coalescent parameter. The options are:           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset popsizepr = uniform(<number>,<number>)              \n");
        YvyraPrint ("                       prset popsizepr = lognormal(<number>,<number>)            \n");
        YvyraPrint ("                       prset popsizepr = normal(<number>,<number>)               \n");
        YvyraPrint ("                       prset popsizepr = gamma(<number>,<number>)                \n");
        YvyraPrint ("                       prset popsizepr = fixed(<number>)                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only relevant if the coalescence process is\n");
        YvyraPrint ("                    selected as the prior on branch lengths. Note that the set-  \n");
        YvyraPrint ("                    ting of 'ploidy' in 'lset' is important for how this para-   \n");
        YvyraPrint ("                    meter is interpreted.                                        \n");
        YvyraPrint ("   Popvarpr      -- In a gene tree - species tree model, this parameter deter-   \n");
        YvyraPrint ("                    mines whether the population size is the same for the entire \n");
        YvyraPrint ("                    species tree ('popvarpr = equal', the default), or varies    \n");
        YvyraPrint ("                    across branches of the species tree ('popvarpr=variable').   \n");
/*      YvyraPrint ("   Growthpr      -- This parameter sets the prior on the exponential growth      \n");
        YvyraPrint ("                    parameter of the coalescence process. The options are:       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset growthpr = uniform(<number>,<number>)               \n");
        YvyraPrint ("                       prset growthpr = exponential(<number>)                    \n");
        YvyraPrint ("                       prset growthpr = fixed(<number>)                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    This parameter is only relevant if the coalescence           \n");
        YvyraPrint ("                    process is selected as the prior on branch lengths.          \n"); */
        YvyraPrint ("   Nodeagepr     -- This parameter specifies the assumptions concerning the age  \n");
        YvyraPrint ("                    of the terminal and interior nodes in the tree. The default  \n");
        YvyraPrint ("                    model ('nodeagepr = unconstrained') assumes that all terminal\n");
        YvyraPrint ("                    nodes are of the same age while the age of interior nodes is \n");
        YvyraPrint ("                    unconstrained. The alternative ('nodeagepr = calibrated')    \n");
        YvyraPrint ("                    option derives a prior probability distribution on terminal  \n");
        YvyraPrint ("                    and interior node ages from the calibration settings (see    \n");
        YvyraPrint ("                    the 'calibrate' command). The 'nodeagepr' parameter is only  \n");
        YvyraPrint ("                    relevant for clock trees.                                    \n");
        YvyraPrint ("   Clockratepr   -- This parameter specifies the prior assumptions concerning the\n");
        YvyraPrint ("                    base substitution rate of the tree, measured in expected num-\n");
        YvyraPrint ("                    ber of substitutions per site per time unit. The default set-\n");
        YvyraPrint ("                    ting is 'Fixed(1.0)', which effectively means that the time  \n");
        YvyraPrint ("                    unit is the number of expected substitutions per site.       \n");
/*      YvyraPrint ("                    If you apply age constraints to the tree, the default setting\n");
        YvyraPrint ("                    changes automatically to 'Exponential(<x>)', where '<x>' (the\n");
        YvyraPrint ("                    rate of exponential) is ten times the age of the maximum age \n");
        YvyraPrint ("                    constraint. This will give you a very vague prior, which may \n");
        YvyraPrint ("                    or may not be adequate for your particular problem.          \n"); */
        YvyraPrint ("                    If you do not have any age calibrations in the tree, you can \n");
        YvyraPrint ("                    still calibrate the tree using 'Clockratepr'. For instance,  \n");
        YvyraPrint ("                    if you know that your sequence data evolve at a rate of 0.20 \n");
        YvyraPrint ("                    substitutions per million years, you might calibrate the tree\n");
        YvyraPrint ("                    by fixing the substitution rate to 0.20 using                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset clockratepr = fixed(0.20)                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    after which the tree will be calibrated using millions of    \n");
        YvyraPrint ("                    years as the unit.                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    You can also assign a prior probability distribution to the  \n");
        YvyraPrint ("                    substitution rate, accommodating the uncertainty of it.      \n");
        YvyraPrint ("                    When you calibrate the nodes, you should properly set this   \n");
        YvyraPrint ("                    prior to match the time unit of the calibrations.            \n");
        YvyraPrint ("                    You can choose among normal, lognormal, exponential and gamma\n");
        YvyraPrint ("                    distributions for this purpose. For instance, to assign a    \n");
        YvyraPrint ("                    normal distribution truncated at 0, so that only positive    \n");
        YvyraPrint ("                    values are allowed, and with mean 0.20 and standard deviation\n");
        YvyraPrint ("                    of 0.02, you would use                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset clockratepr = normal(0.20,0.02)                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    The lognormal distribution is parameterized in terms of the  \n");
        YvyraPrint ("                    mean and standard deviation on the log scale (natural logs). \n");
        YvyraPrint ("                    For instance,                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset clockratepr = lognormal(-1.61,0.10)                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    specifies a lognormal distribution with a mean of log values \n");
        YvyraPrint ("                    of -1.61 and a standard deviation of log values of 0.10. In  \n");
        YvyraPrint ("                    such a case, the mean value of the lognormal distribution is \n");
        YvyraPrint ("                    equal to e^(-1.61 + 0.10^2/2) = 0.20.                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    Note that the 'Clockratepr' parameter has no effect on non-  \n");
        YvyraPrint ("                    clock trees.                                                 \n");
        YvyraPrint ("   Clockvarpr    -- This parameter allows you to specify the type of clock you   \n");
        YvyraPrint ("                    are assuming. The default is 'strict', which corresponds to  \n");
        YvyraPrint ("                    the standard clock model where the evolutionary rate is      \n");
        YvyraPrint ("                    constant throughout the tree. For relaxed clock models, you  \n");
        YvyraPrint ("                    can use 'cpp', 'tk02', 'wn', 'igr', or 'iln'.                \n");
        YvyraPrint ("                    'cpp' invokes a relaxed clock model where the rate evolves   \n");
        YvyraPrint ("                    according to a Compound Poisson Process (CPP) (Huelsenbeck   \n");
        YvyraPrint ("                    et al., 2000).                                               \n");
        YvyraPrint ("                    'tk02' invokes the geometric Brownian Motion model described \n");
        YvyraPrint ("                    by Thorne and Kishino (2002). [autocorrelated lognormal]     \n");
        YvyraPrint ("                    'wn' invokes the white noise model (LePage et al., 2007)     \n");
        YvyraPrint ("                    where each branch has an independent rate drawn from a gamma \n");
        YvyraPrint ("                    distribution with variance proportional to the branch length.\n");
        YvyraPrint ("                    'igr' invokes the Independent Gamma Rate model. The differ-  \n");
        YvyraPrint ("                    ence from 'wn' is that these gamma distributions are i.i.d.  \n");
        YvyraPrint ("                    Thus, the variances do not depend on the branch lengths.     \n");
        YvyraPrint ("                    'iln' invokes the Independent Lognormal model. It is similar \n");
        YvyraPrint ("                    to 'igr' but the rates are i.i.d. lognormal distributions.   \n");
        YvyraPrint ("                    Each of the relaxed clock models has additional parameters   \n");
        YvyraPrint ("                    with priors. For the CPP model, it is 'cppratepr' and        \n");
        YvyraPrint ("                    'cppmultdevpr'; for the TK02 model, it is 'tk02varpr'; for   \n");
        YvyraPrint ("                    the WN model, it is 'wnvarpr'; for the IGR model, it is      \n");
        YvyraPrint ("                    'igrvarpr'; for the ILN model, it is 'ilnvarpr'.             \n");
        YvyraPrint ("                    The 'clockvarpr' parameter is only relevant for clock trees. \n");
        YvyraPrint ("   Cppratepr     -- This parameter allows you to specify a prior probability     \n");
        YvyraPrint ("                    distribution on the rate of the Poisson process generating   \n");
        YvyraPrint ("                    changes in the evolutionary rate in the CPP relaxed clock    \n");
        YvyraPrint ("                    model. You can either fix the rate or associate it with an   \n");
        YvyraPrint ("                    exponential prior using                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset cppratepr = fixed(<number>)                         \n");
        YvyraPrint ("                       prset cppratepr = exponential(<number>)                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    For instance, if you fix the rate to 2, then on a branch     \n");
        YvyraPrint ("                    with the length equal to one expressed in terms of average   \n");
        YvyraPrint ("                    expected number of substitution per site, you expect to see, \n"); 
        YvyraPrint ("                    on average, two rate-modifying events.                       \n");
        YvyraPrint ("                    If you put an exponential(0.1) on the rate, you will be      \n");
        YvyraPrint ("                    estimating the rate against a prior probability distribution \n");
        YvyraPrint ("                    where the expected rate is 10 (= 1/0.1).                     \n");
        YvyraPrint ("   Cppmultdevpr  -- This parameter allows you to specify the standard deviation  \n");
        YvyraPrint ("                    of the log-normal distribution from which the rate multi-    \n");
        YvyraPrint ("                    pliers of the CPP relaxed clock model are drawn. The standard\n");
        YvyraPrint ("                    deviation is given on the log scale. The default value of 1.0\n");
        YvyraPrint ("                    thus corresponds to rate multipliers varying from 0.37 (1/e) \n");
        YvyraPrint ("                    to 2.7 (e) when they are +/- one standard deviation from the \n");
        YvyraPrint ("                    expected mean. The expected mean of the logarithm of the mul-\n");
        YvyraPrint ("                    pliers is fixed to 0, ensuring that the expected mean rate is\n");
        YvyraPrint ("                    1.0. You can change the default value by using               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset cppmultdevpr = fixed(<number>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    where <number> is the standard deviation on the log scale.   \n");
        YvyraPrint ("   TK02varpr     -- This parameter allows you to specify the prior probability   \n");
        YvyraPrint ("                    distribution for the variance of the rate multiplier in the  \n");
        YvyraPrint ("                    Thorne-Kishino ('Brownian motion') relaxed clock model.      \n");
        YvyraPrint ("                    Specifically, the parameter specifies the rate at which the  \n");
        YvyraPrint ("                    variance increases with respect to the base rate of the      \n");
        YvyraPrint ("                    clock. If you have a branch of a length corresponding to 0.4 \n");
        YvyraPrint ("                    expected changes per site according to the base rate of the  \n");
        YvyraPrint ("                    clock, and the tk02var parameter has a value of 2.0, then the\n");
        YvyraPrint ("                    rate multiplier at the end of the branch will be drawn from a\n");
        YvyraPrint ("                    lognormal distribution with a variance of 0.4*2.0 (on the    \n");
        YvyraPrint ("                    linear, not the logarithm scale). The mean is the same as the\n");
        YvyraPrint ("                    rate multiplier at the start of the branch (again on the     \n");
        YvyraPrint ("                    linear scale).                                               \n");
        YvyraPrint ("                    You can set the parameter to a fixed value, or specify that  \n");
        YvyraPrint ("                    it is drawn from an exponential or uniform distribution:     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset tk02varpr = fixed(<number>)                         \n");
        YvyraPrint ("                       prset tk02varpr = exponential(<number>)                   \n");
        YvyraPrint ("                       prset tk02varpr = uniform(<number>,<number>)              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   WNvarpr       -- This parameter allows you to specify the prior distribution  \n");
        YvyraPrint ("                    for the variance of the independent branch rate (white noise)\n");
        YvyraPrint ("                    relaxed clock model. Specifically, the parameter specifies   \n");
        YvyraPrint ("                    the rate at which the variance increases with respect to the \n");
        YvyraPrint ("                    base rate of the clock. If you have a branch of a length     \n");
        YvyraPrint ("                    corresponding to 0.4 expected changes per site according to  \n");
        YvyraPrint ("                    the base rate of the clock, and the wnvar parameter has a    \n");
        YvyraPrint ("                    value of 2.0 , then the effective branch length will be drawn\n");
        YvyraPrint ("                    from a gamma distribution with a variance of 0.4*2.0.        \n");
        YvyraPrint ("                    You can set the parameter to a fixed value, or specify that  \n");
        YvyraPrint ("                    it is drawn from an exponential or uniform distribution:     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset wnvarpr = fixed(<number>)                           \n");
        YvyraPrint ("                       prset wnvarpr = exponential(<number>)                     \n");
        YvyraPrint ("                       prset wnvarpr = uniform(<number>,<number>)                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   IGRvarpr      -- This parameter allows you to specify the prior distribution  \n");
        YvyraPrint ("                    for the variance of the independent gamma rate (IGR) relaxed \n");
        YvyraPrint ("                    clock model. Specifically, the parameter specifies the rate  \n");
        YvyraPrint ("                    at which the variance increases with respect to the base rate\n");
        YvyraPrint ("                    of the clock. The difference from the white noise (WN) model \n");
        YvyraPrint ("                    is that the variance is the same for all branch rates and    \n");
        YvyraPrint ("                    does not depend on the branch length.  If igrvar parameter   \n");
        YvyraPrint ("                    has a value of 0.8 , then the rate multiplier will be drawn  \n");
        YvyraPrint ("                    from a gamma distribution with mean 1.0 and variance 0.8.    \n");
        YvyraPrint ("                    You can set the parameter to a fixed value, or specify that  \n");
        YvyraPrint ("                    it is drawn from an exponential or uniform distribution:     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset igrvarpr = fixed(<number>)                          \n");
        YvyraPrint ("                       prset igrvarpr = exponential(<number>)                    \n");
        YvyraPrint ("                       prset igrvarpr = uniform(<number>,<number>)               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ILNvarpr      -- This parameter allows you to specify the prior distribution  \n");
        YvyraPrint ("                    for the variance of the independent lognormal (ILN) relaxed  \n");
        YvyraPrint ("                    clock model. It is similar with IGR but each rate multiplier \n");
        YvyraPrint ("                    has a lognormal instead of a gamma distribution. If ilnvar   \n");
        YvyraPrint ("                    has a value of 0.8 , then the rate multiplier will be drawn  \n");
        YvyraPrint ("                    from a lognormal distribution with mean 1.0 and variance 0.8.\n");
        YvyraPrint ("                    You can set the parameter to a fixed value, or specify that  \n");
        YvyraPrint ("                    it is drawn from an exponential or uniform distribution:     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset ilnvarpr = fixed(<number>)                          \n");
        YvyraPrint ("                       prset ilnvarpr = exponential(<number>)                    \n");
        YvyraPrint ("                       prset ilnvarpr = uniform(<number>,<number>)               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Ratepr        -- This parameter allows you to specify the site specific rates \n");
        YvyraPrint ("                    model or any other model that allows different partitions to \n");
        YvyraPrint ("                    evolve at different rates. First, you must have defined a    \n");
        YvyraPrint ("                    partition of the characters. For example, you may define a   \n");
        YvyraPrint ("                    partition that divides the characters by codon position, if  \n");
        YvyraPrint ("                    you have DNA data. You can also divide your data using a     \n");
        YvyraPrint ("                    partition that separates different genes from each other.    \n");
        YvyraPrint ("                    The next step is to make the desired partition the active one\n");
        YvyraPrint ("                    using the set command. For example, if your partition is     \n");
        YvyraPrint ("                    called \"by_codon\", then you make that the active partition \n");
        YvyraPrint ("                    using \"set partition=by_codon\". Now that you have defined  \n");
        YvyraPrint ("                    and activated a partition, you can specify the rate multi-   \n");
        YvyraPrint ("                    pliers for the various partitions. The options are:          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                       prset ratepr = fixed                                      \n");
        YvyraPrint ("                       prset ratepr = variable                                   \n");
        YvyraPrint ("                       prset ratepr = dirichlet(<number>,<number>,...,<number>)  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                    If you specify \"fixed\", then the rate multiplier for       \n");
        YvyraPrint ("                    that partition is set to 1 (i.e., the rate is fixed to       \n");
        YvyraPrint ("                    the average rate across partitions). On the other hand,      \n");
        YvyraPrint ("                    if you specify \"variable\", then the rate is allowed to     \n");
        YvyraPrint ("                    vary across partitions subject to the constraint that the    \n");
        YvyraPrint ("                    average rate of substitution across the partitions is 1.     \n");
        YvyraPrint ("                    You must specify a variable rate prior for at least two      \n");
        YvyraPrint ("                    partitions, otherwise the option is not activated when       \n");
        YvyraPrint ("                    calculating likelihoods. The variable option automatically   \n");
        YvyraPrint ("                    associates the partition rates with a dirichlet(1,...,1)     \n");
        YvyraPrint ("                    prior. The dirichlet option is an alternative way of setting \n");
        YvyraPrint ("                    a partition rate to be variable, and also gives accurate     \n");
        YvyraPrint ("                    control of the shape of the prior. The parameters of the     \n");
        YvyraPrint ("                    Dirichlet are listed in the order of the partitions that the \n");
        YvyraPrint ("                    ratepr is applied to. For instance, \"prset applyto=(1,3,4)  \n");
        YvyraPrint ("                    ratepr = dirichlet(10,40,15)\" would set the Dirichlet para- \n");
        YvyraPrint ("                    meter 10 to partition 1, 40 to partition 3, and 15 to parti- \n");
        YvyraPrint ("                    tion 4. The Dirichlet distribution is applied to the weighted\n");
        YvyraPrint ("                    rates; that is, it weights the partition rates according to  \n");
        YvyraPrint ("                    the number of included characters in each partition.         \n");
        YvyraPrint ("   Generatepr    -- This parameter is similar to 'Ratepr' but applies to gene    \n");
        YvyraPrint ("                    trees in the multispecies coalescent, whereas 'Ratepr' app-  \n");
        YvyraPrint ("                    lies to partitions within genes.                             \n");
        YvyraPrint ("                                                                                 \n");
        if (numCurrentDivisions == 0)
            tempInt = 1;
        else
            tempInt = numCurrentDivisions;
        for (i=0; i<tempInt; i++)
            {
            if (numCurrentDivisions == 0)
                {
                YvyraPrint ("   Default model settings:                                                       \n");
                mp = &defaultModel;
                }
            else
                {
                YvyraPrint ("   Model settings for partition %d:                                              \n", i+1);
                mp = &modelParams[i];
                }
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   Parameter        Options                      Current Setting                 \n");
            YvyraPrint ("   ------------------------------------------------------------------            \n");

            YvyraPrint ("   Tratiopr         Beta/Fixed                   %s", mp->tRatioPr);
            if (!strcmp(mp->tRatioPr, "Beta"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->tRatioDir[0], mp->tRatioDir[1]);
            else
                YvyraPrint ("(%1.1lf)\n", mp->tRatioFix);

            YvyraPrint ("   Revmatpr         Dirichlet/Fixed              %s", mp->revMatPr);
            if (!strcmp(mp->revMatPr, "Dirichlet"))
                YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->revMatDir[0],
                mp->revMatDir[1], mp->revMatDir[2], mp->revMatDir[3],
                mp->revMatDir[4], mp->revMatDir[5]);
            else
                YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->revMatFix[0],
                mp->revMatFix[1], mp->revMatFix[2], mp->revMatFix[3],
                mp->revMatFix[4], mp->revMatFix[5]);

            YvyraPrint ("   Aamodelpr        Fixed/Mixed                  %s", mp->aaModelPr);
            if (!strcmp(mp->aaModelPr, "Fixed"))
                YvyraPrint ("(%s)\n", mp->aaModel);
            else
                YvyraPrint ("\n");

            YvyraPrint ("   Aarevmatpr       Dirichlet/Fixed              %s", mp->aaRevMatPr);
            if (!strcmp(mp->aaRevMatPr, "Dirichlet"))
                {
                for (j=1; j<190; j++)
                    if (AreDoublesEqual (mp->aaRevMatDir[0], mp->aaRevMatDir[j], 0.00001) == NO)
                        break;
                if (j==190)
                    YvyraPrint ("(%1.1lf,%1.1lf,...)\n", mp->aaRevMatDir[0], mp->aaRevMatDir[0]);
                else
                    YvyraPrint (" (use 'Showmodel' to see values set by user)\n");
                }
            else
                {
                for (j=1; j<190; j++)
                    if (AreDoublesEqual (mp->aaRevMatFix[0], mp->aaRevMatFix[j], 0.00001) == NO)
                        break;
                if (j==190)
                    YvyraPrint ("(%1.1lf,%1.1lf,...)\n", mp->aaRevMatFix[0], mp->aaRevMatFix[0]);
                else
                    YvyraPrint (" (use 'Showmodel' to see values set by user)\n");
                }

            YvyraPrint ("   Omegapr          Dirichlet/Fixed              %s", mp->omegaPr);
            if (!strcmp(mp->omegaPr, "Dirichlet"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->omegaDir[0], mp->omegaDir[1]);
            else
                YvyraPrint ("(%1.1lf)\n", mp->omegaFix);

            YvyraPrint ("   Ny98omega1pr     Beta/Fixed                   %s", mp->ny98omega1pr);
            if (!strcmp(mp->ny98omega1pr, "Beta"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->ny98omega1Beta[0], mp->ny98omega1Beta[1]);
            else if (!strcmp(mp->ny98omega1pr, "Fixed"))
                YvyraPrint ("(%1.1lf)\n", mp->ny98omega1Fixed);
                
            YvyraPrint ("   Ny98omega3pr     Uniform/Exponential/Fixed    %s", mp->ny98omega3pr);
            if (!strcmp(mp->ny98omega3pr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->ny98omega3Uni[0], mp->ny98omega3Uni[1]);
            else if (!strcmp(mp->ny98omega3pr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->ny98omega3Exp);
            else
                YvyraPrint ("(%1.1lf)\n", mp->ny98omega3Fixed);

            YvyraPrint ("   M3omegapr        Exponential/Fixed            %s", mp->m3omegapr);
            if (!strcmp(mp->m3omegapr, "Exponential"))
                YvyraPrint ("\n");
            else if (!strcmp(mp->m3omegapr, "Fixed"))
                YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", mp->m3omegaFixed[0], mp->m3omegaFixed[1], mp->m3omegaFixed[2]);
                
            YvyraPrint ("   Codoncatfreqs    Dirichlet/Fixed              %s", mp->codonCatFreqPr);
            if (!strcmp(mp->codonCatFreqPr, "Dirichlet"))
                YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", mp->codonCatDir[0], mp->codonCatDir[1], mp->codonCatDir[2]);
            else
                YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", mp->codonCatFreqFix[0], mp->codonCatFreqFix[1], mp->codonCatFreqFix[2]);

            YvyraPrint ("   Statefreqpr      Dirichlet/Fixed              %s", mp->stateFreqPr);
            if (!strcmp(mp->stateFreqPr, "Dirichlet"))
                {
                if (mp->dataType == DNA || mp->dataType == RNA)
                    {
                    if (!strcmp(mp->nucModel, "4by4"))
                        YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->stateFreqsDir[0], mp->stateFreqsDir[1],
                            mp->stateFreqsDir[2], mp->stateFreqsDir[3]);
                    else
                        YvyraPrint ("\n");
                    }
                else if (mp->dataType == RESTRICTION)
                    {
                    YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->stateFreqsDir[0], mp->stateFreqsDir[1]);
                    }
                else
                    YvyraPrint ("\n");
                }
            else if (!strcmp(mp->stateFreqPr, "Fixed"))
                {
                if (mp->dataType == DNA || mp->dataType == RNA)
                    {
                    if (!strcmp(mp->nucModel, "4by4"))
                        YvyraPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->stateFreqsFix[0], mp->stateFreqsFix[1],
                            mp->stateFreqsFix[2], mp->stateFreqsFix[3]);
                    else
                        YvyraPrint ("\n");
                    }
                else if (mp->dataType == RESTRICTION)
                    {
                    YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->stateFreqsFix[0], mp->stateFreqsFix[1]);
                    }
                else
                    YvyraPrint ("\n");
                }

            YvyraPrint ("   Shapepr          Uniform/Exponential/Fixed    %s", mp->shapePr);
            if (!strcmp(mp->shapePr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->shapeUni[0], mp->shapeUni[1]);
            else if (!strcmp(mp->shapePr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->shapeExp);
            else
                YvyraPrint ("(%1.1lf)\n", mp->shapeFix);

            YvyraPrint ("   Ratecorrpr       Uniform/Fixed                %s", mp->adGammaCorPr);
            if (!strcmp(mp->adGammaCorPr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->adgCorrUni[0], mp->adgCorrUni[1]);
            else
                YvyraPrint ("(%1.1lf)\n", mp->adgCorrFix);

            YvyraPrint ("   Pinvarpr         Uniform/Fixed                %s", mp->pInvarPr);
            if (!strcmp(mp->pInvarPr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->pInvarUni[0], mp->pInvarUni[1]);
            else
                YvyraPrint ("(%1.1lf)\n", mp->pInvarFix);

            YvyraPrint ("   Covswitchpr      Uniform/Exponential/Fixed    %s", mp->covSwitchPr);
            if (!strcmp(mp->covSwitchPr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->covswitchUni[0], mp->covswitchUni[1]);
            else if (!strcmp(mp->covSwitchPr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->covswitchExp);
            else
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->covswitchFix[0], mp->covswitchFix[1]);

            YvyraPrint ("   Symdirihyperpr   Uniform/Exponential/Fixed    %s", mp->symPiPr);
            if (!strcmp(mp->symPiPr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->symBetaUni[0], mp->symBetaUni[1]);
            else if (!strcmp(mp->covSwitchPr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->symBetaExp);
            else
                {
                if (mp->symBetaFix < 0)
                    YvyraPrint ("(Infinity)\n");
                else
                    YvyraPrint ("(%1.1lf)\n", mp->symBetaFix);
                }
            
            YvyraPrint ("   Topologypr       Uniform/Constraints/Fixed/   %s", mp->topologyPr);
            if (!strcmp(mp->topologyPr, "Constraints"))
                {
                YvyraPrint ("(");
                for (j=0; j<numDefinedConstraints; j++)
                    {
                    if (mp->activeConstraints[j] == YES)
                        {
                        YvyraPrint ("%d", j+1);
                        break;
                        }
                    }
               for (j++; j<numDefinedConstraints; j++)
                    {
                    if (mp->activeConstraints[j] == YES)
                        {
                        YvyraPrint (",%d", j+1);
                        }
                    }
                YvyraPrint (")\n");
                }
            else if (!strcmp(mp->topologyPr, "Fixed"))
                YvyraPrint ("(%s)\n", userTree[mp->topologyFix]->name);
            else
                YvyraPrint ("\n");
            YvyraPrint ("                    Speciestree                  \n");
            
            YvyraPrint ("   Brlenspr         Unconstrained/Clock/Fixed    %s", mp->brlensPr);
            if (!strcmp(mp->brlensPr, "Unconstrained"))
                {
                if (!strcmp(mp->unconstrainedPr, "Uniform"))
                    YvyraPrint (":Uni(%1.1lf,%1.1lf)\n", mp->brlensUni[0], mp->brlensUni[1]);
                else if (!strcmp(mp->unconstrainedPr, "GammaDir"))
                    YvyraPrint (":GammaDir(%1.1lf,%1.3lf,%1.1lf,%1.1lf)\n",
                                mp->brlensDir[0], mp->brlensDir[1], mp->brlensDir[2], mp->brlensDir[3]);
                else if (!strcmp(mp->unconstrainedPr, "invGamDir"))
                    YvyraPrint (":invGamDir(%1.1lf,%1.3lf,%1.1lf,%1.1lf)\n",
                                mp->brlensDir[0], mp->brlensDir[1], mp->brlensDir[2], mp->brlensDir[3]);
                else if (!strcmp(mp->unconstrainedPr, "twoExp"))
                    YvyraPrint (":twoExp(%1.1lf,%1.1lf)\n", mp->brlens2Exp[0], mp->brlens2Exp[1]);
                else
                    YvyraPrint (":Exp(%1.1lf)\n", mp->brlensExp);
                }
            else if (!strcmp(mp->brlensPr, "Clock"))
                {
                if (!strcmp(mp->clockPr,"Fixed"))
                    YvyraPrint (":%s(%s)\n", mp->clockPr, userTree[mp->brlensFix]->name);
                else
                    YvyraPrint (":%s\n", mp->clockPr);
                }
            else if (!strcmp(mp->brlensPr, "Fixed"))
                YvyraPrint ("(%s)\n", userTree[mp->brlensFix]->name);
            
            YvyraPrint ("   Treeagepr        Gamma/Uniform/Fixed/         %s\n", mp->treeAgePr.name);
            YvyraPrint ("                    Truncatednormal/Lognormal/   \n");
            YvyraPrint ("                    Offsetlognormal/Offsetgamma/ \n");
            YvyraPrint ("                    Offsetexponential            \n");
            
            YvyraPrint ("   Speciationpr     Uniform/Exponential/Fixed    %s", mp->speciationPr);
            if (!strcmp(mp->speciationPr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->speciationUni[0], mp->speciationUni[1]);
            else if (!strcmp(mp->speciationPr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->speciationExp);
            else
                YvyraPrint ("(%1.1lf)\n", mp->speciationFix);
            
            YvyraPrint ("   Extinctionpr     Beta/Fixed                   %s", mp->extinctionPr);
            if (!strcmp(mp->extinctionPr, "Beta"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->extinctionBeta[0], mp->extinctionBeta[1]);
            else if (!strcmp(mp->extinctionPr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->extinctionExp);
            else
                YvyraPrint ("(%1.1lf)\n", mp->extinctionFix);
            
            YvyraPrint ("   Fossilizationpr  Beta/Fixed                   %s", mp->fossilizationPr);
            if (!strcmp(mp->fossilizationPr, "Beta"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->fossilizationBeta[0], mp->fossilizationBeta[1]);
            else if (!strcmp(mp->fossilizationPr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->fossilizationExp);
            else
                YvyraPrint ("(%1.2lf)\n", mp->fossilizationFix);
            
            YvyraPrint ("   SampleStrat      Random/Diversity/Cluster/    %s\n", mp->sampleStrat);
            YvyraPrint ("                    FossilTip                    \n");
            
            YvyraPrint ("   Sampleprob       <number>                     %1.8lf\n", mp->sampleProb);
            
            YvyraPrint ("   Popsizepr        Lognormal/Gamma/Uniform/     %s", mp->popSizePr);
            if (!strcmp(mp->popSizePr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeUni[0], mp->popSizeUni[1]);
            else if (!strcmp(mp->popSizePr, "Lognormal"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeLognormal[0], mp->popSizeLognormal[1]);
            else if (!strcmp(mp->popSizePr, "Normal"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeNormal[0], mp->popSizeNormal[1]);
            else if (!strcmp(mp->popSizePr, "Gamma"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeGamma[0], mp->popSizeGamma[1]);
            else
                YvyraPrint ("(%1.1lf)\n", mp->popSizeFix);
            YvyraPrint ("                    Normal/Fixed                 \n");

            YvyraPrint ("   Popvarpr         Equal/Variable               %s\n", mp->popVarPr);

            /*
            YvyraPrint ("   Growthpr         Uniform/Exponential/         \n");
            YvyraPrint ("                    Fixed/Normal                 %s", mp->growthPr);
            if (!strcmp(mp->growthPr, "Uniform"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->growthUni[0], mp->growthUni[1]);
            else if (!strcmp(mp->growthPr, "Exponential"))
                YvyraPrint ("(%1.1lf)\n", mp->growthExp);
            else if (!strcmp(mp->growthPr, "Normal"))
                YvyraPrint ("(%1.1lf,%1.1lf)\n", mp->growthNorm[0], mp->growthNorm[1]);
            else
                YvyraPrint ("(%1.1lf)\n", mp->growthFix); 
            */

            YvyraPrint ("   Nodeagepr        Unconstrained/Calibrated     %s\n", mp->nodeAgePr);

            YvyraPrint ("   Clockratepr      Fixed/Normal/Lognormal/      %s", mp->clockRatePr);
            if (!strcmp(mp->clockRatePr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->clockRateFix);
            else if (!strcmp(mp->clockRatePr,"Exponential"))
                YvyraPrint ("(%1.2lf)\n", mp->clockRateExp);
            else if (!strcmp(mp->clockRatePr,"Normal"))
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->clockRateNormal[0], mp->clockRateNormal[1]);
            else if (!strcmp(mp->clockRatePr,"Lognormal"))
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->clockRateLognormal[0], mp->clockRateLognormal[1]);
            else
                {
                assert (!strcmp(mp->clockRatePr,"Gamma"));
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->clockRateGamma[0], mp->clockRateGamma[1]);
                }
            YvyraPrint ("                    Exponential/Gamma            \n");

            YvyraPrint ("   Clockvarpr       Strict/Cpp/TK02/WN/IGR/ILN   %s\n", mp->clockVarPr);

            YvyraPrint ("   Cppratepr        Fixed/Exponential            %s", mp->cppRatePr);
            if (!strcmp(mp->cppRatePr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->cppRateFix);
            else /* if (!strcmp(mp->cppRatePr,"Exponential")) */
                YvyraPrint ("(%1.2lf)\n", mp->cppRateExp);

            YvyraPrint ("   Cppmultdevpr     Fixed                        %s", mp->cppMultDevPr);
            YvyraPrint ("(%1.2lf)\n", mp->cppMultDevFix);

            YvyraPrint ("   TK02varpr        Fixed/Exponential/Uniform    %s", mp->tk02varPr);
            if (!strcmp(mp->tk02varPr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->tk02varFix);
            else if (!strcmp(mp->tk02varPr,"Exponential"))
                YvyraPrint ("(%1.2lf)\n", mp->tk02varExp);
            else
                {
                assert (!strcmp(mp->tk02varPr,"Uniform"));
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->tk02varUni[0], mp->tk02varUni[1]);
                }

            YvyraPrint ("   WNvarpr          Fixed/Exponential/Uniform    %s", mp->wnvarPr);
            if (!strcmp(mp->wnvarPr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->wnvarFix);
            else if (!strcmp(mp->wnvarPr,"Exponential"))
                YvyraPrint ("(%1.2lf)\n", mp->wnvarExp);
            else
                {
                assert (!strcmp(mp->wnvarPr,"Uniform"));
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->wnvarUni[0], mp->wnvarUni[1]);
                }

            YvyraPrint ("   IGRvarpr         Fixed/Exponential/Uniform    %s", mp->igrvarPr);
            if (!strcmp(mp->igrvarPr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->igrvarFix);
            else if (!strcmp(mp->igrvarPr,"Exponential"))
                YvyraPrint ("(%1.2lf)\n", mp->igrvarExp);
            else
                {
                assert (!strcmp(mp->igrvarPr,"Uniform"));
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->igrvarUni[0], mp->igrvarUni[1]);
                }
            
            YvyraPrint ("   ILNvarpr         Fixed/Exponential/Uniform    %s", mp->ilnvarPr);
            if (!strcmp(mp->ilnvarPr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->ilnvarFix);
            else if (!strcmp(mp->ilnvarPr,"Exponential"))
                YvyraPrint ("(%1.2lf)\n", mp->ilnvarExp);
            else
                {
                assert (!strcmp(mp->ilnvarPr,"Uniform"));
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->ilnvarUni[0], mp->ilnvarUni[1]);
                }
            
            YvyraPrint ("   Mixedvarpr       Fixed/Exponential/Uniform    %s", mp->mixedvarPr);
            if (!strcmp(mp->mixedvarPr, "Fixed"))
                YvyraPrint ("(%1.2lf)\n", mp->mixedvarFix);
            else if (!strcmp(mp->mixedvarPr,"Exponential"))
                YvyraPrint ("(%1.2lf)\n", mp->mixedvarExp);
            else
                {
                assert (!strcmp(mp->mixedvarPr,"Uniform"));
                YvyraPrint ("(%1.2lf,%1.2lf)\n", mp->mixedvarUni[0], mp->mixedvarUni[1]);
                }

            YvyraPrint ("   Ratepr           Fixed/Variable=Dirichlet     %s", mp->ratePr);
            if (!strcmp(mp->ratePr, "Dirichlet"))
                YvyraPrint ("(...,%1.1lf,...)\n", mp->ratePrDir);
            else
                YvyraPrint ("\n");

            YvyraPrint ("   Generatepr       Fixed/Variable=Dirichlet     %s", mp->generatePr);
            if (!strcmp(mp->generatePr, "Dirichlet"))
                YvyraPrint ("(...,%1.1lf,...)\n", mp->generatePrDir);
            else
                YvyraPrint ("\n");

            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   ------------------------------------------------------------------            \n");
            }
        }
    else if (!strcmp(helpTkn, "Ctype"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Ctype                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command sets the character ordering for standard-type data. The          \n");
        YvyraPrint ("   correct usage is:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      ctype <ordering>:<characters>                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The available options for the <ordering> specifier are:                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("     unordered    -- Movement directly from one state to another is              \n");
        YvyraPrint ("                     allowed in an instant of time.                              \n");
        YvyraPrint ("     ordered      -- Movement is only allowed between adjacent characters.       \n");
        YvyraPrint ("                     For example, perhaps only between 0 <-> 1 and 1 <-> 2       \n");
        YvyraPrint ("                     for a three state character ordered as 0 - 1 - 2.           \n");
        YvyraPrint ("     irreversible -- Rates of change for losses are 0.                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The characters to which the ordering is applied is specified in manner        \n");
        YvyraPrint ("   that is identical to commands such as \"include\" or \"exclude\". For         \n");
        YvyraPrint ("   example,                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      ctype ordered: 10 23 45                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   defines characters 10, 23, and 45 to be of type ordered. Similarly,           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      ctype irreversible: 54-67  71-92                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   defines characters 54 to 67 and characters 71 to 92 to be of type             \n");
        YvyraPrint ("   irreversible. You can use the \".\" to denote the last character, and         \n");
        YvyraPrint ("   \"all\" to denote all of the characters. Finally, you can use the             \n");
        YvyraPrint ("   specifier \"\\\" to apply the ordering to every n-th character or             \n");
        YvyraPrint ("   you can use predefined charsets to specify the character.                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Only one ordering can be used on any specific application of ctype.           \n");
        YvyraPrint ("   If you want to apply different orderings to different characters, then        \n");
        YvyraPrint ("   you need to use ctype multiple times. For example,                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      ctype ordered: 1-50                                                        \n");
        YvyraPrint ("      ctype irreversible: 51-100                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   sets characters 1 to 50 to be ordered and characters 51 to 100 to be          \n");
        YvyraPrint ("   irreversible.                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The ctype command is only sensible with morphological (here called            \n");
        YvyraPrint ("   \"standard\") characters. The program ignores attempts to apply char-         \n");
        YvyraPrint ("   acter orderings to other types of characters, such as DNA characters.         \n");

        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Propset"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Propset                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command allows the user to change the details of the MCMC samplers       \n");
        YvyraPrint ("   (moves) that update the state of the chain. The usage is:                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      propset  <move_name>$<tuning-parameter>=<value>                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Assume we have a branch-length parameter called 'V{all}', which is sampled by \n");
        YvyraPrint ("   the move 'Nodeslider(V{all})' (note that the parameter name is included in    \n");
        YvyraPrint ("   the move name). This move has two tuning parameters: (1) 'prob', the relative \n");
        YvyraPrint ("   proposal probability (a weight defining its probability relative to other     \n");
        YvyraPrint ("   moves); and (2) 'lambda', the tuning  parameter of the branch length          \n");
        YvyraPrint ("   multiplier. A list of the tuning parameters is available by using 'Showmoves' \n");
        YvyraPrint ("   (see below). To change the relative proposal probability to 20 and the tuning \n");
        YvyraPrint ("   parameter to 0.7, use:                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      propset Nodeslider(V{all})$prob=20 Nodeslider(V{all})$lambda=0.7           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This change would apply to all chains in all runs. It is also possible to set \n");
        YvyraPrint ("   the tuning parameters of individual runs and chains using the format:         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      propset  <move_name>$<tuning-parameter>(<run>,<chain>)=<value>             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   where <run> and <chain> are the index numbers of the run and chain for which  \n");
        YvyraPrint ("   you want to change the value. If you leave out the index of the run, the      \n");
        YvyraPrint ("   change will apply to all runs; if you leave out the index of the chain, the   \n");
        YvyraPrint ("   change will similarly apply to all chains. To switch off the                  \n");
        YvyraPrint ("   Nodeslider(V{all}) move in chain 2 of all runs, use:                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      propset  Nodeslider(V{all})$prob(,2)=0                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   It is important to note that all moves are not available until the model has  \n");
        YvyraPrint ("   been completely defined. Any change to the model will cause all proposal      \n");
        YvyraPrint ("   tuning parameters to return to their default values. To see a list of all the \n");
        YvyraPrint ("   moves that are currently switched on for the model, use 'showmoves'. You can  \n");
        YvyraPrint ("   also see other available moves by using 'showmoves allavailable=yes'. A list  \n");
        YvyraPrint ("   of the moves for each parameter in the model is available by using the command\n");
        YvyraPrint ("   'Showparams'. If you change proposal probabilities, make sure that all        \n");
        YvyraPrint ("   parameters that are not fixed in your model have at least one move switched   \n");
        YvyraPrint ("   on.                                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   One word of warning: You should be extremely careful when modifying any       \n");
        YvyraPrint ("   of the chain parameters using 'propset'. It is quite possible to completely   \n");
        YvyraPrint ("   wreck any hope of achieving convergence by inappropriately setting the        \n");
        YvyraPrint ("   tuning parameters. In general, you want to set move tuning parameters such    \n");
        YvyraPrint ("   that the acceptance rate of the move is intermediate (we suggest targeting    \n");
        YvyraPrint ("   the range 10%% to 70%% acceptance, if possible). If the acceptance rate is    \n");
        YvyraPrint ("   outside of this range, the MCMC chain will probably not sample that parameter \n");
        YvyraPrint ("   very efficiently. The acceptance rates for all moves in the cold chain(s) are \n");
        YvyraPrint ("   summarized at the end of each run in the screen output. The acceptance rates  \n");
        YvyraPrint ("   (potentially for all chains, cold and heated) are also printed to the .mcmc   \n");
        YvyraPrint ("   file if MCMC convergence diagnostics are turned on (using 'Mcmc' or 'Mcmcp'). \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Log"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Log                                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command allows output to the screen to also be output to a file.         \n");
        YvyraPrint ("   The usage is:                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      log start/stop filename=<name> append/replace                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The options are:                                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Start/Stop     -- Starts or stops logging of output to file.                  \n");
        YvyraPrint ("   Append/Replace -- Either append to or replace existing file.                  \n");
        YvyraPrint ("   Filename       -- Name of log file (currently, the name of the log            \n");
        YvyraPrint ("                     file is \"%s\").\n", logFileName);
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Translate"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Translate                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used by yvyra to specify the mapping between taxon names    \n");
        YvyraPrint ("   and taxon numbers in a Nexus tree file. For instance,                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      translate                                                                  \n");
        YvyraPrint ("         1 Homo,                                                                 \n");
        YvyraPrint ("         2 Pan,                                                                  \n");
        YvyraPrint ("         3 Gorilla,                                                              \n");
        YvyraPrint ("         4 Hylobates;                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   establishes that the taxon labeled 1 in the trees that follow is Homo, the    \n");
        YvyraPrint ("   taxon labeled 2 is Pan, etc.                                                  \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Usertree"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Usertree                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command allows you to specify a user tree. The user tree can then be     \n");
        YvyraPrint ("   used as a starting tree for a MCMC analysis. The format for the command is    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      usertree = <tree in Newick format>                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example,                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      usertree = (A,B,(C,D))                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   specifies an unrooted tree of four species. Note that the program re-         \n");
        YvyraPrint ("   quires that trees are binary (i.e., strictly bifurcating). Hence, there       \n");
        YvyraPrint ("   can be only one three-way split, as shown in the example. If the tree         \n");
        YvyraPrint ("   is not binary, the program will return an error.                              \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Mcmc"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Mcmc                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command starts the Markov chain Monte Carlo (MCMC) analysis. The         \n");
        YvyraPrint ("   posterior probability of phylogenetic trees (and other parameters of the      \n");
        YvyraPrint ("   substitution model) cannot be determined analytically. Instead, MCMC is       \n");
        YvyraPrint ("   used to approximate the posterior probabilities of trees by drawing           \n");
        YvyraPrint ("   (dependent) samples from the posterior distribution. This program can         \n");
        YvyraPrint ("   implement a variant of MCMC called \"Metropolis-coupled Markov chain Monte    \n");
        YvyraPrint ("   Carlo\", or MCMCMC for short. Basically, \"Nchains\" are run, with            \n");
        YvyraPrint ("   Nchains - 1 of them heated. The chains are labelled 1, 2, ..., Nchains.       \n");
        YvyraPrint ("   The heat that is applied to the i-th chain is B = 1 / (1 + temp X i). B       \n");
        YvyraPrint ("   is the power to which the posterior probability is raised. When B = 0, all    \n");
        YvyraPrint ("   trees have equal probability and the chain freely visits trees. B = 1 is      \n");
        YvyraPrint ("   the \"cold\" chain (or the distribution of interest). MCMCMC can mix          \n");
        YvyraPrint ("   better than ordinary MCMC; after all of the chains have gone through          \n");
        YvyraPrint ("   one cycle, two chains are chosen at random and an attempt is made to          \n");
        YvyraPrint ("   swap the states (with the probability of a swap being determined by the       \n");
        YvyraPrint ("   Metropolis et al. equation). This allows the chain to potentially jump        \n");
        YvyraPrint ("   a valley in a single bound. The correct usage is                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      mcmc <parameter> = <value> ... <parameter> = <value>                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example,                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      mcmc ngen=100000 nchains=4 temp=0.5                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   performs a MCMCMC analysis with four chains with the temperature set to       \n");
        YvyraPrint ("   0.5. The chains would be run for 100,000 cycles.                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Ngen         -- This option sets the number of cycles for the MCMC alg-       \n");
        YvyraPrint ("                   orithm. This should be a big number as you want the chain     \n");
        YvyraPrint ("                   to first reach stationarity, and then remain there for        \n");
        YvyraPrint ("                   enough time to take lots of samples.                          \n");
        YvyraPrint ("   Nruns        -- How many independent analyses are started simultaneously.     \n");
        YvyraPrint ("   Nchains      -- How many chains are run for each analysis for the MCMCMC      \n");
        YvyraPrint ("                   variant. The default is 4: 1 cold chain and 3 heated chains.  \n");
        YvyraPrint ("                   If Nchains is set to 1, yvyra will use regular MCMC sam-    \n");
        YvyraPrint ("                   pling, without heating.                                       \n");
        YvyraPrint ("   Temp         -- The temperature parameter for heating the chains. The higher  \n");
        YvyraPrint ("                   the temperature, the more likely the heated chains are to     \n");
        YvyraPrint ("                   move between isolated peaks in the posterior distribution.    \n");
        YvyraPrint ("                   However, excessive heating may lead to very low acceptance    \n");
        YvyraPrint ("                   rates for swaps between different chains. Before changing the \n");
        YvyraPrint ("                   default setting, however, note that the acceptance rates of   \n");
        YvyraPrint ("                   swaps tend to fluctuate during the burn-in phase of the run.  \n");
        YvyraPrint ("   Reweight     -- Here, you specify three numbers, that respectively represent  \n");
        YvyraPrint ("                   the percentage of characters to decrease in weight, the       \n");
        YvyraPrint ("                   percentage of characters to increase in weight, and the       \n");
        YvyraPrint ("                   increment. An increase/decrease in weight is achieved by      \n");
        YvyraPrint ("                   replicating/removing a character in the matrix. This is       \n");
        YvyraPrint ("                   only done to non-cold chains. The format for this parameter   \n");
        YvyraPrint ("                   is \"reweight=(<number>,<number>)\" or \"reweight=(<number>,  \n");
        YvyraPrint ("                   <number>,<number>)\".                                         \n");
        YvyraPrint ("   Swapfreq     -- This specifies how often swaps of states between chains are   \n");
        YvyraPrint ("                   attempted. You must be running at least two chains for this   \n");
        YvyraPrint ("                   option to be relevant. The default is Swapfreq=1, resulting   \n");
        YvyraPrint ("                   in Nswaps (see below) swaps being tried each generation of    \n");
        YvyraPrint ("                   the run. If Swapfreq is set to 10, then Nswaps swaps will be  \n");
        YvyraPrint ("                   tried every tenth generation of the run.                      \n");
        YvyraPrint ("   Nswaps       -- The number of swaps tried for each swapping generation of the \n");
        YvyraPrint ("                   chain (see also Swapfreq).                                    \n");
        YvyraPrint ("   Samplefreq   -- This specifies how often the Markov chain is sampled. You     \n");
        YvyraPrint ("                   can sample the chain every cycle, but this results in very    \n");
        YvyraPrint ("                   large output files. Thinning the chain is a way of making     \n");
        YvyraPrint ("                   these files smaller and making the samples more independent.  \n");
        YvyraPrint ("   Printfreq    -- This specifies how often information about the chain is       \n");
        YvyraPrint ("                   printed to the screen.                                        \n");
        YvyraPrint ("   Printall     -- If set to NO, only cold chains in a MCMC analysis are printed \n");
        YvyraPrint ("                   to screen. If set to YES, both cold and heated chains will be \n");
        YvyraPrint ("                   output. This setting only affects the printing to screen, it  \n");
        YvyraPrint ("                   does not change the way values are written to file.           \n");
        YvyraPrint ("   Printmax     -- The maximum number of chains to print to screen.              \n");
        YvyraPrint ("   Mcmcdiagn    -- Determines whether acceptance ratios of moves and swaps will  \n");
        YvyraPrint ("                   be printed to file. The file will be named similarly to the   \n");
        YvyraPrint ("                   '.p' and '.t' files, but will have the ending '.mcmc'. If     \n");
        YvyraPrint ("                   more than one independent analysis is run simultaneously (see \n");
        YvyraPrint ("                   Nruns below), convergence diagnostics for tree topology will  \n");
        YvyraPrint ("                   also be printed to this file. The convergence diagnostic used \n");
        YvyraPrint ("                   is the average standard deviation in partition frequency      \n");
        YvyraPrint ("                   values across independent analyses. The Burnin setting (see   \n");
        YvyraPrint ("                   below) determines how many samples will be discarded as burnin\n");
        YvyraPrint ("                   before calculating the partition frequencies. The Minpartfreq \n");
        YvyraPrint ("                   setting (see below) determines the minimum partition frequency\n");
        YvyraPrint ("                   required for a partition to be included in the calculation. As\n");
        YvyraPrint ("                   the independent analyses approach stationarity (converge), the\n");
        YvyraPrint ("                   value of the diagnostic is expected to approach zero.         \n");
        YvyraPrint ("   Diagnfreq    -- The number of generations between the calculation of MCMC     \n");
        YvyraPrint ("                   diagnostics (see Mcmcdiagn above).                            \n");
        YvyraPrint ("   Diagnstat    -- The statistic to use for run-time convergence diagnostics.    \n");
        YvyraPrint ("                   Choices are 'Avgstddev' for average standard deviation of     \n");
        YvyraPrint ("                   split frequencies and 'Maxstddev' for maximum standard devia- \n");
        YvyraPrint ("                   tion of split frequencies.                                    \n");
        YvyraPrint ("   Savetrees    -- If you are using a relative burnin for run-time convergence   \n");
        YvyraPrint ("                   diagnostics, tree samples need to be deleted from split       \n");
        YvyraPrint ("                   frequency counters as the cut-off point for the burnin moves  \n");
        YvyraPrint ("                   during the run. If 'Savetrees' is set to 'No', tree samples   \n");
        YvyraPrint ("                   to be discarded are read back in from file. If 'Savetrees' is \n");
        YvyraPrint ("                   set to 'Yes', the tree samples to be removed will be stored   \n");
        YvyraPrint ("                   in the internal memory instead. This can use up a lot of      \n");
        YvyraPrint ("                   memory in large analyses.                                     \n");
        YvyraPrint ("   Minpartfreq  -- The minimum frequency required for a partition to be included \n");
        YvyraPrint ("                   in the calculation of the topology convergence diagnostic. The\n");
        YvyraPrint ("                   partition is included if the minimum frequency is reached in  \n");
        YvyraPrint ("                   at least one of the independent tree samples that are com-    \n");
        YvyraPrint ("                   pared.                                                        \n");
        YvyraPrint ("   Allchains    -- If this option is set to YES, acceptance ratios for moves are \n");
        YvyraPrint ("                   recorded for all chains, cold or heated. By default, only the \n");
        YvyraPrint ("                   acceptance ratios for the cold chain are recorded.            \n");
        YvyraPrint ("   Allcomps     -- If this option is set to YES, topological convergence diag-   \n");
        YvyraPrint ("                   nostics are calculated over all pairwise comparisons of runs. \n");
        YvyraPrint ("                   If it is set to NO, only the overall value is reported.       \n");
        YvyraPrint ("   Relburnin    -- If this option is set to YES, then a proportion of the sampled\n");
        YvyraPrint ("                   values will be discarded as burnin when calculating the con-  \n");
        YvyraPrint ("                   vergence diagnostic. The proportion to be discarded is set    \n");
        YvyraPrint ("                   with Burninfrac (see below). When the Relburnin option is set \n");
        YvyraPrint ("                   to NO, then a specific number of samples will be discarded    \n");
        YvyraPrint ("                   instead. This number is set by Burnin (see below).            \n");
        YvyraPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
        YvyraPrint ("                   be discarded when convergence diagnostics are calculated.     \n");
        YvyraPrint ("                   The value of this option is only relevant when Relburnin is   \n");
        YvyraPrint ("                   set to NO.                                                    \n");
        YvyraPrint ("   BurninFrac   -- Determines the fraction of samples that will be discarded     \n");
        YvyraPrint ("                   when convergence diagnostics are calculated. The value of     \n");
        YvyraPrint ("                   this option is only relevant when Relburnin is set to YES.    \n");
        YvyraPrint ("                   Example: A value for this option of 0.25 means that 25%% of   \n");
        YvyraPrint ("                   the samples will be discarded.                                \n");
        YvyraPrint ("   Stoprule     -- If this option is set to NO, then the chain is run the number \n");
        YvyraPrint ("                   of generations determined by Ngen. If it is set to YES, and   \n");
        YvyraPrint ("                   topological convergence diagnostics are calculated (Mcmcdiagn \n");
        YvyraPrint ("                   is set to YES), then the chain will be stopped before the pre-\n");
        YvyraPrint ("                   determined number of generations if the convergence diagnostic\n");
        YvyraPrint ("                   falls below the stop value.                                   \n");
        YvyraPrint ("   Stopval      -- The critical value for the topological convergence diagnostic.\n");
        YvyraPrint ("                   Only used when Stoprule and Mcmcdiagn are set to yes, and     \n");
        YvyraPrint ("                   more than one analysis is run simultaneously (Nruns > 1).     \n");
        YvyraPrint ("   Checkpoint   -- If this parameter is set to 'Yes', all the current parameter  \n");
        YvyraPrint ("                   values of all chains will be printed to a check-pointing file \n");
        YvyraPrint ("                   every 'Checkfreq' generation of the analysis. The file will be\n");
        YvyraPrint ("                   named <Filename>.ckp and allows you to restart the analysis   \n");
        YvyraPrint ("                   from the last check point. This can be handy if you are       \n");
        YvyraPrint ("                   running a long analysis and want to extend it, or if there is \n");
        YvyraPrint ("                   a risk that a long analysis will be inadvertently interrupted \n");
        YvyraPrint ("                   by hardware failure or other factors that are out of your     \n");
        YvyraPrint ("                   control.                                                      \n");
        YvyraPrint ("   Checkfreq    -- The number of generations between check-pointing. See the     \n");
        YvyraPrint ("                   'Checkpoint' parameter above for more information.            \n");
        YvyraPrint ("   Filename     -- The name of the files that will be generated. Two files       \n");
        YvyraPrint ("                   are generated: \"<Filename>.t\" and \"<Filename>.p\".         \n");
        YvyraPrint ("                   The .t file contains the trees whereas the .p file con-       \n");
        YvyraPrint ("                   tains the sampled values of the parameters.                   \n");
        YvyraPrint ("   Startparams  -- The starting values for the model parameters are set to       \n");
        YvyraPrint ("                   arbitrary or random values when the parameters are created.   \n");
        YvyraPrint ("                   These starting values can be altered using the 'Startvals'    \n");
        YvyraPrint ("                   command. The 'Startparams=reset' option allows you to reset   \n");
        YvyraPrint ("                   the starting values to the default at the start of the ana-   \n");
        YvyraPrint ("                   lysis, overriding any previous user-defined starting values.  \n");
        YvyraPrint ("                   Under the default option, 'current', the chains will use the  \n");
        YvyraPrint ("                   current starting values.                                      \n");
        YvyraPrint ("   Starttree    -- The starting tree(s) for the chain can either be randomly     \n");
        YvyraPrint ("                   selected or user-defined. It might be a good idea to          \n");
        YvyraPrint ("                   start from randomly chosen trees; convergence seems           \n");
        YvyraPrint ("                   likely if independently run chains, each of which             \n");
        YvyraPrint ("                   started from different random trees, converge to the same     \n");
        YvyraPrint ("                   answer. If you want the chain to start from user-defined      \n");
        YvyraPrint ("                   trees instead, you first need to read in your tree(s) from a  \n");
        YvyraPrint ("                   Nexus file with a 'trees' block, and then you need to set the \n");
        YvyraPrint ("                   starting tree(s) using the 'Startvals' command. Finally, you  \n");
        YvyraPrint ("                   need to make sure that 'Starttree' is set to 'current'. If    \n");
        YvyraPrint ("                   you do not set the starting tree(s), the chains will start    \n");
        YvyraPrint ("                   with random trees. Setting 'Starttree' to 'random' causes     \n");
        YvyraPrint ("                   new starting trees to be drawn randomly at the start of the   \n");
        YvyraPrint ("                   run, overwriting any previous user-defined starting trees.    \n");
        YvyraPrint ("   Nperts       -- This is the number of random perturbations to apply to the    \n");
        YvyraPrint ("                   user starting tree. This allows you to have something         \n");
        YvyraPrint ("                   between completely random and user-defined trees start        \n");
        YvyraPrint ("                   the chain.                                                    \n");
        YvyraPrint ("   Data         -- When Data is set to NO, the chain is run without data. This   \n");
        YvyraPrint ("                   should be used only for examining induced priors. DO NOT SET  \n");
        YvyraPrint ("                   'DATA' TO 'NO' UNLESS YOU KNOW WHAT YOU ARE DOING!            \n");
        YvyraPrint ("   Ordertaxa    -- Determines whether taxa should be ordered before trees are    \n");
        YvyraPrint ("                   printed to file. If set to 'Yes', terminals in the sampled    \n");
        YvyraPrint ("                   trees will be reordered to match the order of the taxa in the \n");
        YvyraPrint ("                   data matrix as closely as possible. By default, trees will be \n");
        YvyraPrint ("                   printed without reordering of taxa.                           \n");
        YvyraPrint ("   Append       -- Set this to 'Yes' to append the results of the current run to \n");
        YvyraPrint ("                   a previous run. yvyra will first read in the results of the \n");
        YvyraPrint ("                   previous run (number of generations and sampled splits) and   \n");
        YvyraPrint ("                   will then continue that run where you left it off. Make sure  \n");
        YvyraPrint ("                   that the output file names used in the previous run are the   \n");
        YvyraPrint ("                   same as those in the current run.                             \n");
        YvyraPrint ("   Autotune     -- Set this to 'Yes' to autotune the proposals that change       \n");
        YvyraPrint ("                   substitution model parameters. When set to 'No', the tuning   \n");
        YvyraPrint ("                   parameters are fixed to their starting values. Note that the  \n");
        YvyraPrint ("                   autotuning occurs independently for each chain. The target    \n");
        YvyraPrint ("                   acceptance rate for each move can be changed using the        \n");
        YvyraPrint ("                   'Propset' command.                                            \n");
        YvyraPrint ("   Tunefreq     -- When a proposal has been tried 'Tunefreq' times, its tuning   \n");
        YvyraPrint ("                   parameter is adjusted to reach the target acceptance rate     \n");
        YvyraPrint ("                   if 'Autotune' is set to 'Yes'.                                \n");
        YvyraPrint ("                                                                                 \n");
        PrintSettings ("Mcmc");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Mcmcp"))
        {
        // PrintYesNo (chainParams.saveBrlens, yesNoStr);
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Mcmcp                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command sets the parameters of the Markov chain Monte Carlo (MCMC)       \n");
        YvyraPrint ("   analysis without actually starting the chain. This command is identical       \n");
        YvyraPrint ("   in all respects to Mcmc, except that the analysis will not start after        \n");
        YvyraPrint ("   this command is issued. For more details on the options, check the help       \n");
        YvyraPrint ("   menu for Mcmc.\n");
        YvyraPrint ("                                                                                 \n");
        PrintSettings ("Mcmc");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Ss"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Ss                                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to start stepping-stone sampling, which is an efficient  \n");
        YvyraPrint ("   and accurate method for estimating the marginal likelihood of the currently   \n");
        YvyraPrint ("   specified model. It is considerably more accurate than the harmonic mean of   \n");
        YvyraPrint ("   the likelihoods from a standard MCMC run on the model (calculated by the      \n");
        YvyraPrint ("   'Sump' command) but it requires a separate MCMC-like run. To be more specific,\n");
        YvyraPrint ("   stepping-stone sampling uses importance sampling to estimate each ratio in a  \n");
        YvyraPrint ("   series of discrete steps bridging the posterior and prior distributions.      \n");
        YvyraPrint ("   The importance distributions that are used are called power posterior distri- \n");
        YvyraPrint ("   butions, and are defined as prior*(likelihood^beta). By varying beta from 1 to\n");
        YvyraPrint ("   0, we get a series of distributions that connect the posterior (beta = 1) to  \n");
        YvyraPrint ("   the prior (beta = 0).                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The power posterior distributions are sampled using MCMC. First, we start a   \n");
        YvyraPrint ("   standard MCMC chain on the posterior distribution, and let it run until we    \n");
        YvyraPrint ("   have reached the criterion specified by the 'Burninss' option. After this, we \n");
        YvyraPrint ("   step through the power posterior distributions until we reach the prior dis-  \n");
        YvyraPrint ("   tribution. In each of the 'Nsteps' steps, we sample from a new power poster-  \n");
        YvyraPrint ("   ior distribution with a distinct beta value. The beta values correspond to    \n");
        YvyraPrint ("   'Nsteps' evenly spaced quantiles in a Beta distribution with the parameters   \n");
        YvyraPrint ("   'Alpha' and 1.0. For the first sampling step, the beta value is equal to the  \n");
        YvyraPrint ("   last quantile, i.e., it is close to 1.0. For each successive step, the beta   \n");
        YvyraPrint ("   value takes on the value of the next quantile, in decreasing order, until it  \n");
        YvyraPrint ("   reaches the value of 0.0. If you change value of 'FromPrior' from default 'No'\n");
        YvyraPrint ("   to 'Yes' then the direction of power posterior change during SS analysis is   \n");
        YvyraPrint ("   opposite to the one described above, i.e. we start from sampling prior and    \n");
        YvyraPrint ("   finish close to posterior.                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The 'Ss' procedure uses the same machinery as the standard 'Mcmc' algorithm,  \n");
        YvyraPrint ("   and shares most of its parameters with the 'Mcmc' and 'Mcmcp' commands. All   \n");
        YvyraPrint ("   'Mcmc' parameters, except those related to burnin, have the same meaning and  \n");
        YvyraPrint ("   usage in the 'Ss' command as they have in the 'Mcmc' command. The 'Mcmc'      \n");
        YvyraPrint ("   burnin parameters are used to set up burnin within each step. The 'Ss' command\n");
        YvyraPrint ("   also uses its own burnin parameter, 'Burninss' (see below for details). The   \n");
        YvyraPrint ("   'Ss' command also has its own parameters for specifying the number of steps   \n");
        YvyraPrint ("   and the shape of the Beta distribution from which the beta values are computed\n");
        YvyraPrint ("   (see below).                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that the 'Ngen' parameter of 'Mcmc' is used to set the maximum number of \n");
        YvyraPrint ("   generations processed, including both the burnin and the following steps in   \n");
        YvyraPrint ("   the stepping-stone sampling phase. For instance, assume that 'Burninss' is set\n");
        YvyraPrint ("   to '-1', 'Nsteps' to '49', 'Ngen' to '1000000' and 'Samplefreq' to '1000'.    \n");
        YvyraPrint ("   We will then get 1,000 samples in total (1,000,000 / 1,000). These will fall  \n");
        YvyraPrint ("   into 50 bins, one of which represents the burnin and is discarded. Each step  \n");
        YvyraPrint ("   in the algorithm will thus be represented by 20 samples.                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   More information on 'Mcmc' parameters is available in the help for the 'Mcmc' \n");
        YvyraPrint ("   and 'Mcmcp' commands. Only the exclusive 'Ss' parameters are listed below.    \n");
        YvyraPrint ("   These can only be set up using the 'Ss' command, while the parameters shared  \n");
        YvyraPrint ("   with 'Mcmc' and 'Mcmcp' can also be set up using those commands.              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The correct usage is                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      ss <parameter>=<value> ... <parameter>=<value>                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that a command:                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      ss <setting parameters shared with mcmc> <setting exclusive ss parameters> \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   would be equivalent to executing two commands:                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("     mcmcp <setting parameters shared with mcmc>;                                \n");
        YvyraPrint ("     ss <setting exclusive ss parameters>;                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For more information on the stepping-stone algorithm, see:                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Xie, W., P. O. Lewis, Y. Fan, L. Kuo, and M.-H. Chen. 2011. Improving marginal\n");
        YvyraPrint ("      likelihood estimation for Bayesian phylogenetic model selection. Systematic\n");
        YvyraPrint ("      Biology 60:150-160.                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Available options:                                                            \n");
        YvyraPrint ("   (NB: Only exclusive ss parameters listed here. For additional parameters, see \n");
        YvyraPrint ("        help on 'Mcmc' or 'Mcmcp'.                                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Alpha        -- The beta values used in the stepping-stone sampling procedure \n");
        YvyraPrint ("                   correspond to evenly spaced quantiles from a Beta('Alpha',1.0)\n");
        YvyraPrint ("                   distribution. The parameter 'Alpha' determines the skewness of\n");
        YvyraPrint ("                   the beta values. If 'Alpha' is set to '1.0', the beta values  \n");
        YvyraPrint ("                   would be spaced uniformly on the interval (0.0,1.0). However, \n");
        YvyraPrint ("                   better results are obtained if the beta values are skewed.    \n");
        YvyraPrint ("                   Empirically, it was observed that 'Alpha' values in the range \n");
        YvyraPrint ("                   of 0.3 to 0.5 produce the most accurate results.              \n");
        YvyraPrint ("   Burninss     -- Fixed number of samples discarded before sampling of the first\n");
        YvyraPrint ("                   step starts. 'Burninss' can be specified using either a pos-  \n");
        YvyraPrint ("                   itive or a negative number. If the number is positive, it is  \n");
        YvyraPrint ("                   interpreted as the number of samples to discard as burnin. If \n");
        YvyraPrint ("                   the number is negative, its absolute value is interpreted as  \n");
        YvyraPrint ("                   the length of the burnin in terms of the length of each of the\n");
        YvyraPrint ("                   following steps in the stepping-stone algorithm. For instance,\n");
        YvyraPrint ("                   a value of '-1' means that the length of the burnin is the    \n");
        YvyraPrint ("                   same as the length of each of the subsequent steps.           \n");
        YvyraPrint ("   Nsteps       -- Number of steps in the stepping-stone algorithm. Typically, a \n");
        YvyraPrint ("                   number above 30 is sufficient for accurate results.           \n");
        YvyraPrint ("   FromPrior    -- If it is set to 'Yes', it indicates that in the first step we \n"); 
        YvyraPrint ("                   sample from the prior, with each consecutive step we sample   \n");
        YvyraPrint ("                   closer to the posterior. 'No' indicates the opposite direction\n");
        YvyraPrint ("                   of power posterior change, i.e. in the first step we sample   \n");
        YvyraPrint ("                   close to the posterior, and with each consecutive step we     \n");
        YvyraPrint ("                   sample closer to the prior.                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter          Options               Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Alpha              <number>              %1.2lf\n", chainParams.alphaSS);
        YvyraPrint ("   BurninSS           <number>              %d\n", chainParams.burninSS);
        YvyraPrint ("   Nsteps             <number>              %d\n", chainParams.numStepsSS);
        YvyraPrint ("   FromPrior           Yes/No               %s                                   \n", chainParams.startFromPriorSS == YES ? "Yes" : "No");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
else if (!strcmp(helpTkn, "Ssp"))
        {
        // PrintYesNo (chainParams.saveBrlens, yesNoStr);
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Ssp                                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command sets the parameters of the stepping-stone sampling               \n");
        YvyraPrint ("   analysis without actually starting the chain. This command is identical       \n");
        YvyraPrint ("   in all respects to Ss, except that the analysis will not start after          \n");
        YvyraPrint ("   this command is issued. For more details on the options, check the help       \n");
        YvyraPrint ("   menu for Ss.\n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter          Options               Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Alpha              <number>              %1.2lf\n", chainParams.alphaSS);
        YvyraPrint ("   BurninSS           <number>              %d\n", chainParams.burninSS);
        YvyraPrint ("   Nsteps             <number>              %d\n", chainParams.numStepsSS);
        YvyraPrint ("   FromPrior           Yes/No               %s                                   \n", chainParams.startFromPriorSS == YES ? "Yes" : "No");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
else if (!strcmp(helpTkn, "Set"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Set                                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to set some general features of the model or program     \n");
        YvyraPrint ("   behavior. The correct usage is                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      set <parameter>=<value> ... <parameter>=<value>                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Available options:                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Seed         -- Sets the seed number for the random number generator. The     \n");
        YvyraPrint ("                   random number seed is initialized haphazardly at the beg-     \n");
        YvyraPrint ("                   inning of each yvyra session. This option allows you to     \n");
        YvyraPrint ("                   set the seed to some specific value, thereby allowing you     \n");
        YvyraPrint ("                   to exactly repeat an analysis. If the analysis uses swapping  \n");
        YvyraPrint ("                   between cold and heated chains, you must also set the swap    \n");
        YvyraPrint ("                   seed (see below) to exactly repeat the analysis.              \n");
        YvyraPrint ("   Swapseed     -- Sets the seed used for generating the swapping sequence       \n");
        YvyraPrint ("                   when Metropolis-coupled heated chains are used. This seed     \n");
        YvyraPrint ("                   is initialized haphazardly at the beginning of each yvyra   \n");
        YvyraPrint ("                   session. This option allows you to set the seed to some       \n");
        YvyraPrint ("                   specific value, thereby allowing you to exactly repeat a      \n");
        YvyraPrint ("                   swap sequence. See also the 'Seed' option.                    \n");
        YvyraPrint ("   Dir          -- The working directory. Specifies the absolute or relative path\n");
        YvyraPrint ("                   to the working directory. If left empty, the working directory\n");
        YvyraPrint ("                   is the current directory.                                     \n");
        YvyraPrint ("   Partition    -- Set this option to a valid partition id, either the number or \n");
        YvyraPrint ("                   name of a defined partition, to enforce a specific partition- \n");
        YvyraPrint ("                   ing of the data. When a data matrix is read in, a partition   \n");
        YvyraPrint ("                   called \"Default\" is automatically created. It divides the   \n");
        YvyraPrint ("                   data into one part for each data type. If you only have one   \n");
        YvyraPrint ("                   data type, DNA for instance, the default partition will not   \n");
        YvyraPrint ("                   divide up the data at all. The default partition is always    \n");
        YvyraPrint ("                   the first partition, so 'set partition=1' is the same as      \n");
        YvyraPrint ("                   'set partition=default'.                                      \n");
        YvyraPrint ("   Speciespartition -- Set this option to a valid speciespartition id, either the\n");
        YvyraPrint ("                   number or name of a defined speciespartition, to enforce a    \n");
        YvyraPrint ("                   specific partitioning of taxa to species. When a data matrix  \n");
        YvyraPrint ("                   is read in, a speciespartition called \"Default\" is auto-    \n");
        YvyraPrint ("                   matically created. It assigns one taxon for each species. The \n"); 
        YvyraPrint ("                   default speciespartition is always the first speciespartition,\n");
        YvyraPrint ("                   so 'set speciespartition=1' is the same as                    \n");
        YvyraPrint ("                   'set speciespartition=default'.                               \n");
        YvyraPrint ("   Autoclose    -- If autoclose is set to 'yes', then the program will not prompt\n");
        YvyraPrint ("                   you during the course of executing a file. This is particular-\n");
        YvyraPrint ("                   ly useful when you run yvyra in batch mode.                 \n");
        YvyraPrint ("   Nowarnings   -- If nowarnings is set to yes, then the program will not prompt \n");
        YvyraPrint ("                   you when overwriting or appending an output file that is al-  \n");
        YvyraPrint ("                   ready present. If 'nowarnings=no' (the default setting), then \n");
        YvyraPrint ("                   the program prompts the user before overwriting output files. \n");
        YvyraPrint ("   Autoreplace  -- When nowarnings is set to yes, then yvyra will by default   \n");
        YvyraPrint ("                   overwrite output files that already exists. This may cause    \n");
        YvyraPrint ("                   irrecoverable loss of previous results if you have not removed\n");
        YvyraPrint ("                   or renamed the files from previous runs. To override this be- \n");
        YvyraPrint ("                   havior, set autoreplace to no, in which case new output will  \n");
        YvyraPrint ("                   be appended to existing files instead.                        \n");
        YvyraPrint ("   Quitonerror  -- If quitonerror is set to yes, then the program will quit when \n");
        YvyraPrint ("                   an error is encountered, after printing an error message. If  \n");
        YvyraPrint ("                   quitonerror is set to no (the default setting), then the      \n");
        YvyraPrint ("                   program will wait for additional commands from the command    \n");
        YvyraPrint ("                   line after the error message is printed.                      \n");
        YvyraPrint ("   Scientific   -- Set this option to 'Yes' to write sampled values to file in   \n");
        YvyraPrint ("                   scientific format and to 'No' to write them in fixed format.  \n");
        YvyraPrint ("                   Fixed format is easier for humans to read but you risk losing \n");
        YvyraPrint ("                   precision for small numbers. For instance, sampled values that\n");
        YvyraPrint ("                   are less than 1E-6 will print to file as '0.000000' if fixed  \n");
        YvyraPrint ("                   format is used and 'precision' is set to 6.                   \n");
        YvyraPrint ("   Precision    -- Precision allows you to set the number of decimals to be prin-\n");
        YvyraPrint ("                   ted when sampled values are written to file. Precision must be\n");
        YvyraPrint ("                   in the range 3 to 15.                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter          Options               Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Seed               <number>              %ld                                  \n", globalSeed);
        YvyraPrint ("   Swapseed           <number>              %ld                                  \n", swapSeed);
        YvyraPrint ("   Dir                <name>                \"%s\"\n", workingDir);
        if (defMatrix == YES)
            YvyraPrint ("   Partition          <name>                %s\n", partitionNames[partitionNum]);
        else
            YvyraPrint ("   Partition          <name>                \"\"\n");
        if (defTaxa == YES)
            YvyraPrint ("   Speciespartition   <name>                %s\n", speciespartitionNames[speciespartitionNum]);
        else
            YvyraPrint ("   Speciespartition   <name>                \"\"\n");
        YvyraPrint ("   Autoclose          Yes/No                %s                                   \n", autoClose == YES ? "Yes" : "No");
        YvyraPrint ("   Nowarnings         Yes/No                %s                                   \n", noWarn == YES ? "Yes" : "No");
        YvyraPrint ("   Autoreplace        Yes/No                %s                                   \n", autoOverwrite == YES ? "Yes" : "No");
        YvyraPrint ("   Quitonerror        Yes/No                %s                                   \n", quitOnError == YES ? "Yes" : "No");
        YvyraPrint ("   Scientific         Yes/No                %s                                   \n", scientific == YES ? "Yes" : "No");
        YvyraPrint ("   Precision          <number>              %d                                   \n", precision);
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Charset"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Charset                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command defines a character set. The format for the charset command      \n"); 
        YvyraPrint ("   is                                                                            \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      charset <name> = <character numbers>                                       \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example, \"charset first_pos = 1-720\\3\" defines a character set         \n");
        YvyraPrint ("   called \"first_pos\" that includes every third site from 1 to 720.            \n");
        YvyraPrint ("   The character set name cannot have any spaces in it. The slash (\\)           \n");
        YvyraPrint ("   is a nifty way of telling the program to assign every third (or               \n");
        YvyraPrint ("   second, or fifth, or whatever) character to the character set.                \n");
        YvyraPrint ("   This option is best used not from the command line, but rather as a           \n");
        YvyraPrint ("   line in the yvyra block of a file. Note that you can use \".\" to           \n");
        YvyraPrint ("   stand in for the last character (e.g., charset 1-.\\3).                       \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Outgroup"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Outgroup                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command assigns a taxon to the outgroup. The correct usage is:           \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      outgroup <number>/<taxon name>                                             \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example, \"outgroup 3\" assigns the third taxon in the matrix to be       \n");
        YvyraPrint ("   the outgroup. Similarly, \"outgroup Homo_sapiens\" assigns the taxon          \n");
        YvyraPrint ("   \"Homo_sapiens\" to be the outgroup (assuming that there is a taxon named     \n");
        YvyraPrint ("   \"Homo_sapiens\" in the matrix). Only a single taxon can be assigned to       \n");
        YvyraPrint ("   be the outgroup.                                                              \n");
        YvyraPrint ("                                                                                 \n");
        if (defTaxa == YES)
            YvyraPrint ("   Current outgroup: %s (taxon no. %d)\n", taxaNames[outGroupNum], outGroupNum+1);
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showusertrees"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Showusertrees                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the currently defined user trees. The correct usage        \n");
        YvyraPrint ("   is \"showusertrees\".                                                         \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmcmctrees"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Showmcmctrees                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the current trees used by the Markov chains.               \n");
        YvyraPrint ("   is \"showmcmctrees\".                                                         \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Deroot"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Deroot                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command deroots the user tree. If the tree is already unrooted, a        \n");
        YvyraPrint ("   warning is issued. The correct usage is \"deroot\".                           \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Root"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Root                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command roots the tree. If the tree is already rooted, a warning         \n");
        YvyraPrint ("   is issued. The tree is rooted at the midpoint between the outgroup species    \n");
        YvyraPrint ("   and the ingroup species. The correct usage is \"root\".                       \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Taxset"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Taxset                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command defines a taxon set. The format for the taxset command           \n"); 
        YvyraPrint ("   is                                                                            \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      taxset <name> = <taxon names or numbers>                                   \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example, \"taxset apes = Homo Pan Gorilla Orang gibbon\" defines a        \n");
        YvyraPrint ("   taxon set called \"apes\" that includes five taxa (namely, apes).             \n");
        YvyraPrint ("   You can assign up to 30 taxon sets. This option is best used                  \n");
        YvyraPrint ("   not from the command line but rather as a line in the yvyra block           \n");
        YvyraPrint ("   of a file.                                                                    \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Taxlabels"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Taxlabels                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command defines taxon labels. It could be used within taxa block.        \n"); 
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Charstat"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Charstat                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the status of all the characters. The correct usage        \n");
        YvyraPrint ("   is                                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      charstat                                                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   After typing \"charstat\", the character number, whether it is excluded       \n");
        YvyraPrint ("   or included, and the partition identity are shown. The output is paused       \n");
        YvyraPrint ("   every 100 characters. This pause can be turned off by setting autoclose       \n");
        YvyraPrint ("   to \"yes\" (set autoclose=yes).                                               \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Taxastat"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Taxastat                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the status of all the taxa. The correct usage is           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      taxastat                                                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   After typing \"taxastat\", the taxon number, name, and whether it is          \n");
        YvyraPrint ("   excluded or included are shown.                                               \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Partition"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Partition                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command allows you to specify a character partition. The format for      \n"); 
        YvyraPrint ("   this command is                                                               \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      partition <name> = <num parts>:<chars in first>, ...,<chars in last>       \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example, \"partition by_codon = 3:1st_pos,2nd_pos,3rd_pos\" specifies     \n"); 
        YvyraPrint ("   a partition called \"by_codon\" which consists of three parts (first,         \n"); 
        YvyraPrint ("   second, and third codon positions). Here, we are assuming that the sites      \n"); 
        YvyraPrint ("   in each partition were defined using the charset command. You can specify     \n"); 
        YvyraPrint ("   a partition without using charset as follows:                                 \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      partition by_codon = 3:1 4 6 9 12,2 5 7 10 13,3 6 8 11 14                  \n"); 
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   However, we recommend that you use the charsets to define a set of char-      \n"); 
        YvyraPrint ("   acters and then use these predefined sets when defining the partition.        \n"); 
        YvyraPrint ("   Also, it makes more sense to define a partition as a line in the yvyra      \n"); 
        YvyraPrint ("   block than to issue the command from the command line (then again, you        \n"); 
        YvyraPrint ("   may be a masochist, and want to do extra work).                               \n"); 
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Exclude"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Exclude                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command excludes characters from the analysis. The correct usage is      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exclude <number> <number> <number>                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or                                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exclude <number> - <number>                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or                                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exclude <charset>                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or some combination thereof. Moreover, you can use the specifier \"\\\" to    \n");
        YvyraPrint ("   exclude every nth character. For example, the following                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exclude 1-100\\3                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   would exclude every third character. As a specific example,                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exclude 2 3 10-14 22                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   excludes sites 2, 3, 10, 11, 12, 13, 14, and 22 from the analysis. Also,      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      exclude all                                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   excludes all of the characters from the analysis. Excluding all characters    \n");
        YvyraPrint ("   does not leave you much information for inferring phylogeny.                  \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Include"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Include                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command includes characters that were previously excluded from the       \n");
        YvyraPrint ("   analysis. The correct usage is                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      include <number> <number> <number>                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or                                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      include <number> - <number>                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or                                                                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      include <charset>                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   or some combination thereof. Moreover, you can use the specifier \"\\\" to    \n");
        YvyraPrint ("   include every nth character. For example, the following                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      include 1-100\\3                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   would include every third character. As a specific example,                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      include 2 3 10-14 22                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   includes sites 2, 3, 10, 11, 12, 13, 14, and 22 from the analysis. Also,      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      include all                                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   includes all of the characters in the analysis. Including all of the          \n");
        YvyraPrint ("   characters (even if many of them are bad) is a very total-evidence-like       \n");
        YvyraPrint ("   thing to do. Doing this will make a certain group of people very happy.       \n");
        YvyraPrint ("   On the other hand, simply using this program would make those same people     \n");
        YvyraPrint ("   unhappy.                                                                      \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Delete"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Delete                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command deletes taxa from the analysis. The correct usage is:            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      delete <name and/or number and/or taxset> ...                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A list of the taxon names or taxon numbers (labelled 1 to Ntax in the order   \n");
        YvyraPrint ("   in the matrix) or taxset(s) can be used.  For example, the following:         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      delete 1 2 Homo_sapiens                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   deletes taxa 1, 2, and the taxon labelled Homo_sapiens from the analysis.     \n");
        YvyraPrint ("   You can also use \"all\" to delete all of the taxa. For example,              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      delete all                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   deletes all of the taxa from the analysis. Of course, a phylogenetic anal-    \n");
        YvyraPrint ("   ysis that does not include any taxa is fairly uninteresting.                  \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Restore"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Restore                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command restores taxa to the analysis. The correct usage is:             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      restore <name and/or number and/or taxset> ...                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A list of the taxon names or taxon numbers (labelled 1 to Ntax in the order   \n");
        YvyraPrint ("   in the matrix) or taxset(s) can be used.  For example, the following:         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      restore 1 2 Homo_sapiens                                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   restores taxa 1, 2, and the taxon labelled Homo_sapiens to the analysis.      \n");
        YvyraPrint ("   You can also use \"all\" to restore all of the taxa. For example,             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      restore all                                                                \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   restores all of the taxa to the analysis.                                     \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Quit"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Quit                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command quits the program. The correct usage is:                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      quit                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   It is a very easy command to use properly.                                    \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Disclaimer"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Disclaimer                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the disclaimer for the program. In short, the disclaimer   \n");
        YvyraPrint ("   states that the authors are not responsible for any silly things you may do   \n");
        YvyraPrint ("   to your computer or any unforeseen but possibly nasty things the computer     \n");
        YvyraPrint ("   program may inadvertently do to you.                                          \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Unlink"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Unlink                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command unlinks model parameters across partitions of the data. The      \n");
        YvyraPrint ("   correct usage is:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      unlink <parameter name> = (<all> or <partition list>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A little background is necessary to understand this command. Upon exe-        \n");
        YvyraPrint ("   cution of a file, a default partition is set up. This partition refer-        \n");
        YvyraPrint ("   enced either by its name (\"default\") or number (0). If your data are        \n");
        YvyraPrint ("   all of one type, then this default partition does not actually divide up      \n");
        YvyraPrint ("   your characters. However, if your datatype is mixed, then the default         \n");
        YvyraPrint ("   partition contains as many divisions as there are datatypes in your           \n");
        YvyraPrint ("   character matrix. Of course, you can also define other partitions, and        \n");
        YvyraPrint ("   switch among them using the set command (\"set partition=<name/number>\").    \n");
        YvyraPrint ("   Importantly, you can also assign model parameters to individual part-         \n");
        YvyraPrint ("   itions or to groups of them using the \"applyto\" option in lset and          \n");
        YvyraPrint ("   prset. When the program attempts to perform an analysis, the model is         \n");
        YvyraPrint ("   set for individual partitions. If the same parameter applies to differ-       \n");
        YvyraPrint ("   partitions and if that parameter has the same prior, then the program         \n");
        YvyraPrint ("   will link the parameters: that is, it will use a single value for the         \n");
        YvyraPrint ("   parameter. The program's default, then, is to strive for parsimony.           \n");
        YvyraPrint ("   However, there are lots of cases where you may want unlink a parameter        \n");
        YvyraPrint ("   across partitions. For example, you may want a different transition/          \n");
        YvyraPrint ("   transversion rate ratio to apply to different partitions. This command        \n");
        YvyraPrint ("   allows you to unlink the parameters, or to make them different across         \n");
        YvyraPrint ("   partitions. The converse of this command is \"link\", which links to-         \n");
        YvyraPrint ("   gether parameters that were previously told to be different. The list         \n");
        YvyraPrint ("   of parameters that can be unlinked includes:                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      Tratio          -- Transition/transversion rate ratio                      \n");
        YvyraPrint ("      Revmat          -- Substitution rates of GTR model                         \n");
        YvyraPrint ("      Omega           -- Nonsynonymous/synonymous rate ratio                     \n");
        YvyraPrint ("      Statefreq       -- Character state frequencies                             \n");
        YvyraPrint ("      Shape           -- Gamma/LNorm shape parameter                             \n");
        YvyraPrint ("      Pinvar          -- Proportion of invariable sites                          \n");
        YvyraPrint ("      Correlation     -- Correlation parameter of autodiscrete gamma             \n");
        YvyraPrint ("      Ratemultiplier  -- Rate multiplier for partitions                          \n");
        YvyraPrint ("      Switchrates     -- Switching rates for covarion model                      \n");
        YvyraPrint ("      Topology        -- Topology of tree                                        \n");
        YvyraPrint ("      Brlens          -- Branch lengths of tree                                  \n");
        YvyraPrint ("      Popsize         -- Population size for coalescence process                 \n");
        YvyraPrint ("      Growthrate      -- Growth rate of coalescence process                      \n"); 
        YvyraPrint ("      Aamodel         -- Aminoacid rate matrix                                   \n"); 
        YvyraPrint ("      Cpprate         -- Rate of Compound Poisson Process (CPP)                  \n"); 
        YvyraPrint ("      Cppmultdev      -- Standard dev. of CPP rate multipliers (log scale)       \n"); 
        YvyraPrint ("      Cppevents       -- CPP events                                              \n"); 
        YvyraPrint ("      TK02var         -- Variance increase in TK02 relaxed clock model           \n"); 
        YvyraPrint ("      WNvar           -- Variance increase in WN relaxed clock model             \n");
        YvyraPrint ("      IGRvar          -- Variance in IGR relaxed clock model                     \n");
        YvyraPrint ("      ILNvar          -- Variance in ILN relaxed clock model                     \n");
    //  YvyraPrint ("      TK02branchrates -- Branch rates of TK02  relaxed clock model               \n");
    //  YvyraPrint ("      WNbranchrates   -- Branch rates of WN    relaxed clock model               \n");
    //  YvyraPrint ("      IGRbranchrates  -- Branch rates of IGR   relaxed clock model               \n");
    //  YvyraPrint ("      ILNbranchrates  -- Branch rates of ILN   relaxed clock model               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example,                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      unlink shape=(all)                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   unlinks the gamma/lnorm shape parameter across all partitions of the data.    \n");
        YvyraPrint ("   You can use \"showmodel\" to see the current linking status of the            \n");
        YvyraPrint ("   characters.                                                                   \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Link"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Link                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command links model parameters across partitions of the data. The        \n");
        YvyraPrint ("   correct usage is:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      link <parameter name> = (<all> or <partition list>)                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The list of parameters that can be linked includes:                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      Tratio          -- Transition/transversion rate ratio                      \n");
        YvyraPrint ("      Revmat          -- Substitution rates of GTR model                         \n");
        YvyraPrint ("      Omega           -- Nonsynonymous/synonymous rate ratio                     \n");
        YvyraPrint ("      Statefreq       -- Character state frequencies                             \n");
        YvyraPrint ("      Shape           -- Gamma/LNorm shape parameter                             \n");
        YvyraPrint ("      Pinvar          -- Proportion of invariable sites                          \n");
        YvyraPrint ("      Correlation     -- Correlation parameter of autodiscrete gamma             \n");
        YvyraPrint ("      Ratemultiplier  -- Rate multiplier for partitions                          \n");
        YvyraPrint ("      Switchrates     -- Switching rates for covarion model                      \n");
        YvyraPrint ("      Topology        -- Topology of tree                                        \n");
        YvyraPrint ("      Brlens          -- Branch lengths of tree                                  \n");
        YvyraPrint ("      Popsize         -- Population size for coalescence process                 \n");
        YvyraPrint ("      Growthrate      -- Growth rate of coalescence process                      \n");
        YvyraPrint ("      Aamodel         -- Aminoacid rate matrix                                   \n");
        YvyraPrint ("      Cpprate         -- Rate of Compound Poisson Process (CPP)                  \n");
        YvyraPrint ("      Cppmultdev      -- Standard dev. of CPP rate multipliers (log scale)       \n");
        YvyraPrint ("      Cppevents       -- CPP events                                              \n");
        YvyraPrint ("      TK02var         -- Variance increase in TK02 relaxed clock model           \n");
        YvyraPrint ("      WNvar           -- Variance increase in WN relaxed clock model             \n");
        YvyraPrint ("      IGRvar          -- Variance in IGR relaxed clock model                     \n");
        YvyraPrint ("      ILNvar          -- Variance in ILN relaxed clock model                     \n");
    //  YvyraPrint ("      TK02branchrates -- Branch rates of TK02  relaxed clock model               \n");
    //  YvyraPrint ("      WNbranchrates   -- Branch rates of WN    relaxed clock model               \n");
    //  YvyraPrint ("      IGRbranchrates  -- Branch rates of IGR   relaxed clock model               \n");
    //  YvyraPrint ("      ILNbranchrates  -- Branch rates of ILN   relaxed clock model               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   For example,                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      link shape=(all)                                                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   links the gamma/lnorm shape parameter across all partitions of the data.      \n");
        YvyraPrint ("   You can use \"showmodel\" to see the current linking status of the            \n");
        YvyraPrint ("   characters. For more information on this command, see the help menu           \n");
        YvyraPrint ("   for link's converse, unlink (\"help unlink\");                                \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Help"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Help                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command provides useful information on the use of this program. The      \n");
        YvyraPrint ("   correct usage is                                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      help                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   which gives a list of all available commands with a brief description of      \n");
        YvyraPrint ("   each or                                                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      help <command>                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   which gives detailed information on the use of <command>.                     \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Sump"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Sump                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   During an MCMC analysis, yvyra prints the sampled parameter values to one or\n");
        YvyraPrint ("   more tab-delimited text files, one for each independent run in your analysis. \n");
        YvyraPrint ("   The command 'Sump' summarizes the information in this parameter file or these \n");
        YvyraPrint ("   parameter files. By default, the root of the parameter file name(s) is assumed\n");
        YvyraPrint ("   to be the name of the last matrix-containing nexus file. yvyra also remem-  \n");
        YvyraPrint ("   bers the number of independent runs in the last analysis that you set up, re- \n");
        YvyraPrint ("   gardless of whether you actually ran it. For instance, if there were two in-  \n");
        YvyraPrint ("   dependent runs, which is the initial setting when you read in a new matrix,   \n");
        YvyraPrint ("   yvyra will assume that there are two parameter files with the endings       \n");
        YvyraPrint ("   '.run1.p' and '.run2.p'. You can change the root of the file names and the    \n");
        YvyraPrint ("   number of runs using the 'Filename' and 'Nruns' settings.                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   When you invoke the 'Sump' command, three items are output: (1) a generation  \n");
        YvyraPrint ("   plot of the likelihood values; (2) estimates of the marginal likelihood of    \n");
        YvyraPrint ("   the model; and (3) a table with the mean, variance, and 95 percent credible   \n");
        YvyraPrint ("   interval for the sampled parameters. All three items are output to screen.    \n");
        YvyraPrint ("   The table of marginal likelihoods is also printed to a file with the ending   \n");
        YvyraPrint ("   '.lstat' and the parameter table to a file with the ending '.pstat'. For some \n");
        YvyraPrint ("   model parameters, there may also be a '.mstat' file.                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   When running 'Sump' you typically want to discard a specified number or       \n");
        YvyraPrint ("   fraction of samples from the beginning of the chain as the burn in. This is   \n");
        YvyraPrint ("   done using the same mechanism used by the 'Mcmc' command. That is, if you     \n");
        YvyraPrint ("   run an MCMC analysis with a relative burn in of 25 %% of samples for con-     \n");
        YvyraPrint ("   vergence diagnostics, then the same burn in will be used for a subsequent     \n");
        YvyraPrint ("   sump command, unless a different burn in is specified. That is, issuing       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   sump                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   immediately after 'Mcmc', will result in using the same burn in settings as   \n");
        YvyraPrint ("   for the 'Mcmc' command. All burnin settings are reset to default values every \n");
        YvyraPrint ("   time a new matrix is read in, namely relative burnin ('relburnin=yes') with   \n");
        YvyraPrint ("   25 %% of samples discarded ('burninfrac = 0.25').                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Relburnin    -- If this option is set to 'Yes', then a proportion of the      \n");
        YvyraPrint ("                   samples will be discarded as burnin when calculating summary  \n");
        YvyraPrint ("                   statistics. The proportion to be discarded is set with        \n");
        YvyraPrint ("                   'Burninfrac' (see below). When the 'Relburnin' option is set  \n");
        YvyraPrint ("                   to 'No', then a specific number of samples is discarded       \n");
        YvyraPrint ("                   instead. This number is set by 'Burnin' (see below). Note that\n");
        YvyraPrint ("                   the burnin setting is shared across the 'Sumt', 'Sump', and   \n");
        YvyraPrint ("                   'Mcmc' commands.                                              \n");
        YvyraPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
        YvyraPrint ("                   be discarded when summary statistics are calculated. The      \n");
        YvyraPrint ("                   value of this option is only applicable when 'Relburnin' is   \n");
        YvyraPrint ("                   set to 'No'.                                                  \n");
        YvyraPrint ("   Burninfrac   -- Determines the fraction of samples that will be discarded when\n");
        YvyraPrint ("                   summary statistics are calculated. The setting only takes     \n");
        YvyraPrint ("                   effect if 'Relburnin' is set to 'Yes'.                        \n");
        YvyraPrint ("   Nruns        -- Determines how many '.p' files from independent analyses that \n");
        YvyraPrint ("                   will be summarized. If Nruns > 1 then the names of the files  \n");
        YvyraPrint ("                   are derived from 'Filename' by adding '.run1.p', '.run2.p',   \n");
        YvyraPrint ("                   etc. If Nruns=1, then the single file name is obtained by     \n");
        YvyraPrint ("                   adding '.p' to 'Filename'.                                    \n");
        YvyraPrint ("   Filename     -- The name of the file to be summarized. This is the base of the\n");
        YvyraPrint ("                   file name to which endings are added according to the current \n");
        YvyraPrint ("                   setting of the 'Nruns' parameter. If 'Nruns' is 1, then only  \n");
        YvyraPrint ("                   '.p' is added to the file name. Otherwise, the endings will   \n");
        YvyraPrint ("                   be '.run1.p', '.run2.p', etc.                                 \n");
        YvyraPrint ("   Outputname   -- Base name of the file(s) to which 'Sump' results will be      \n");
        YvyraPrint ("                   printed.                                                      \n");
        YvyraPrint ("   Hpd          -- Determines whether credibility intervals will be given as the \n");
        YvyraPrint ("                   region of Highest Posterior Density ('Yes') or as the interval\n");
        YvyraPrint ("                   containing the median 95 %% of sampled values ('No').         \n");
        YvyraPrint ("   Minprob      -- Determines the minimum probability of submodels to be included\n");
        YvyraPrint ("                   in summary statistics. Only applicable to models that explore \n");
        YvyraPrint ("                   submodel spaces, like 'nst=mixed' and 'aamodelpr=mixed'.      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter       Options                  Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        YvyraPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        YvyraPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        YvyraPrint ("   Nruns           <number>                 %d                                   \n", sumpParams.numRuns);
        if (sumpParams.numRuns == 1)
            YvyraPrint ("   Filename        <name>                   %s<.p>\n", sumpParams.sumpFileName);
        else
            YvyraPrint ("   Filename        <name>                   %s<.run<i>.p>\n", sumpParams.sumpFileName);
        YvyraPrint ("   Outputname      <name>                   %s<.pstat etc>\n", sumpParams.sumpOutfile);
        YvyraPrint ("   Hpd             Yes/No                   %s                                   \n", sumpParams.HPD == YES ? "Yes" : "No");
        YvyraPrint ("   Minprob         <number>                 %1.3lf                               \n", sumpParams.minProb);
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Sumss"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Sumss                                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command summarizes results of stepping stone analyses. It is a tool to   \n");
        YvyraPrint ("   investigate the obtained results, and to help find the proper step burn-in.   \n");
        YvyraPrint ("   To get more help information on stepping-stone analyses, use 'help ss'.       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   During stepping-stone analysis, yvyra collects the sampled likelihoods in   \n");
        YvyraPrint ("   order to estimate the marginal likelihood at the end. It also prints the sam- \n");
        YvyraPrint ("   pled parameter values to one or more tab-delimited text files, one for each   \n");
        YvyraPrint ("   independent run in your analysis. The command 'Sumss' summarizes likelihood   \n");
        YvyraPrint ("   values stored in these parameter files and calculates marginal likelihood es- \n");
        YvyraPrint ("   timates. The names of the files that are summarized are exactly the same as   \n");
        YvyraPrint ("   the names of the files used for the 'sump' command. In fact, the 'filename'   \n");
        YvyraPrint ("   setting is a shared setting for the 'sump' and 'sumss' commands. That is, if  \n");
        YvyraPrint ("   you change the setting in one of the commands, it would change the setting in \n");
        YvyraPrint ("   the other command as well.                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   When you invoke the 'Sumss' command, three items are output: (1) 'Step contri-\n");
        YvyraPrint ("   bution table' - summarizes the contribution of each step to the overall esti- \n");
        YvyraPrint ("   mate; (2) 'Step plot' - plot of the likelihood values for the initial burn-in \n");
        YvyraPrint ("   phase or a chosen step in the stepping-stone algorithm; (3) 'Joined plot' -   \n");
        YvyraPrint ("   summarizes sampling across all steps in the algorithm.                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Step contribution table                                                       \n");
        YvyraPrint ("   The printed table is similar to the one output to the .ss file. The main pur- \n");
        YvyraPrint ("   pose of the table is to summarize marginal likelihood for different values of \n");
        YvyraPrint ("   the step burn-in after the stepping stone  analysis has finished. The burn-in \n");
        YvyraPrint ("   is controlled by the 'Relburnin', 'Burnin' and 'Burninfrac' settings.         \n");
        YvyraPrint ("   Note that during stepping-stone analyses, step contributions to marginal      \n");
        YvyraPrint ("   likelihood are calculated based on all generations excluding burn-in. 'Sumss' \n");
        YvyraPrint ("   on the other hand makes estimates based only on the sampled generations. This \n");
        YvyraPrint ("   may lead to slight difference in results compared to the one printed to the   \n");
        YvyraPrint ("   .ss file.                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Step plot                                                                     \n");
        YvyraPrint ("   The main objective of the plot is to provide a close look at a given step     \n");
        YvyraPrint ("   in the analysis. Which step is printed here is defined by the 'Steptoplot'    \n");
        YvyraPrint ("   setting.  The plot could be used to inspect if the chosen step burn-in is     \n");
        YvyraPrint ("   appropriate for the given step. It could also be used to check if the initial \n");
        YvyraPrint ("   burn-in phase has converged. Note that the amount of discarded samples is     \n");
        YvyraPrint ("   controlled by the 'Discardfrac' setting, and not by the ordinary burn-in      \n");
        YvyraPrint ("   settings.                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Joined plot                                                                   \n");
        YvyraPrint ("   Different steps sample from different power posterior distributions. When we  \n");
        YvyraPrint ("   switch from one distribution to another, it takes some number of generations  \n");
        YvyraPrint ("   before the chain settles at the correct stationary distribution. This lag is  \n");
        YvyraPrint ("   called a 'temperature lag' and if the corresponding samples are not removed,  \n");
        YvyraPrint ("   it will result in a biased estimate. It is difficult to determine the lag be- \n");
        YvyraPrint ("   forehand, but yvyra allows you to explore different step burn-in settings   \n");
        YvyraPrint ("   after you have finished the stepping-stone algorithm, without having to rerun \n");
        YvyraPrint ("   the whole analysis. The 'Joined plot' helps to facilitate the choice of the   \n");
        YvyraPrint ("   right step burn-in. The plot summarizes samples across all steps and gives you\n");
        YvyraPrint ("   a quick overview of the whole analysis.                                       \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Specifically, the following procedure is used to obtain the joined plot. Each \n");
        YvyraPrint ("   step has the same number N of samples taken. We number each sample 1 to N     \n");
        YvyraPrint ("   within steps according to the order in which the samples are taken. The first \n"); 
        YvyraPrint ("   sample in each step is numbered 1, and the last sample is N. For each number i\n");
        YvyraPrint ("   in [1,..., N], we sum up log likelihoods for all samples numbered i across all\n");
        YvyraPrint ("   steps. The joined plot is a graph of the step number versus the normalized    \n");
        YvyraPrint ("   sums we get in the procedure describe above. This directly visualizes the tem-\n");
        YvyraPrint ("   perature lag and allows you to select the appropriate step burn-in.           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Ideally, after you discard the appropriate step burn-in, the graph should     \n");
        YvyraPrint ("   appear as white noise around the estimated value. If you see an increasing or \n");
        YvyraPrint ("   decreasing tendency in the beginning of the graph, you should increase the    \n");
        YvyraPrint ("   step burn-in. If you see an increasing or decreasing tendency across the whole\n");
        YvyraPrint ("   graph, then the initial burn-in phase was not long enough. In this case, you  \n");
        YvyraPrint ("   need to rerun the analysis with a longer initial burn-in.                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   To make it easier to observe tendencies in the plotted graph you can choose   \n");
        YvyraPrint ("   different levels of curve smoothing. If 'Smoothing' is set to k, it means that\n");
        YvyraPrint ("   for each step i we take an average over step i and k neighboring samples in   \n");
        YvyraPrint ("   both directions, i.e., the k-smoothed estimate for step i is an average over  \n");
        YvyraPrint ("   values for steps [i-k,...,i+k].                                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Allruns      -- If set to 'Yes', it forces all runs to be printed on the same \n");
        YvyraPrint ("                   graph when drawing joined and step plots. If set to 'No', each\n");
        YvyraPrint ("                   run is printed on a separate plot.                            \n");
        YvyraPrint ("   Askmore      -- Long analyses may produce huge .p files. Reading in them may  \n");
        YvyraPrint ("                   take several minutes. If you want to investigate different    \n");
        YvyraPrint ("                   aspects of your analyses, it could be very inconvenient to    \n");
        YvyraPrint ("                   wait for several minutes each time you want to get a new sum- \n");
        YvyraPrint ("                   mary for different settings. If you set 'Askmore' to 'YES',   \n");
        YvyraPrint ("                   sumss will read .p files only once. After responding to the   \n");
        YvyraPrint ("                   original query, it will interactively ask you if you wish to  \n");
        YvyraPrint ("                   produce more tables and plots for different settings of       \n");
        YvyraPrint ("                   'Burnin' or 'Smoothing' (see below).                          \n");
        YvyraPrint ("   Relburnin    -- If this option is set to 'Yes', then a proportion of the      \n");
        YvyraPrint ("                   samples from each step will be discarded as burnin when calcu-\n");
        YvyraPrint ("                   lsting summary statistics. The proportion to be discarded is  \n");
        YvyraPrint ("                   set with 'Burninfrac' (see below). When the 'Relburnin' option\n");
        YvyraPrint ("                   is set to 'No', then a specific number of samples is discarded\n");
        YvyraPrint ("                   instead. This number is set by 'Burnin'. Note that the burnin \n");
        YvyraPrint ("                   settings --- 'Relburnin', 'Burnin', and 'Burninfrac' --- are  \n");
        YvyraPrint ("                   shared across the 'Sumt', 'Sump', 'Sumss' and 'Mcmc' commands.\n");
        YvyraPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
        YvyraPrint ("                   be discarded from each step when summary statistics are calcu-\n");
        YvyraPrint ("                   lated. The value of this option is only applicable when       \n");
        YvyraPrint ("                   'Relburnin' is set to 'No'.                                   \n");
        YvyraPrint ("   Burninfrac   -- Determines the fraction of samples that will be discarded from\n");
        YvyraPrint ("                   each step when summary statistics are calculated. The setting \n");
        YvyraPrint ("                   only takes effect if 'Relburnin' is set to 'Yes'.             \n");
        YvyraPrint ("   Discardfrac  -- Determines the fraction of samples that will be discarded when\n");
        YvyraPrint ("                   a step plot is printed. It is similar to the 'Burninfrac' set-\n");
        YvyraPrint ("                   ting, but unlike 'Burninfrac' it is used only for better vis- \n");
        YvyraPrint ("                   ualization of the step plot. It has no effect on the number of\n");
        YvyraPrint ("                   samples discarded during marginal likelihood computation.     \n");
        YvyraPrint ("   Filename     -- The name of the file to be summarized. This is the base of the\n");
        YvyraPrint ("                   file name to which endings are added according to the current \n");
        YvyraPrint ("                   setting of the 'Nruns' parameter. If 'Nruns' is 1, then only  \n");
        YvyraPrint ("                   '.p' is added to the file name. Otherwise, the endings will   \n");
        YvyraPrint ("                   be '.run1.p', '.run2.p', etc. Note that the 'Filename' setting\n");
        YvyraPrint ("                   is shared with 'sump' command.                                \n");
        YvyraPrint ("   Nruns        -- Determines how many '.p' files from independent analyses that \n");
        YvyraPrint ("                   will be summarized. If Nruns > 1 then the names of the files  \n");
        YvyraPrint ("                   are derived from 'Filename' by adding '.run1.p', '.run2.p',   \n");
        YvyraPrint ("                   etc. If Nruns=1, then the single file name is obtained by     \n");
        YvyraPrint ("                   adding '.p' to 'Filename'.                                    \n");
        YvyraPrint ("   Steptoplot   -- Defines which step will be printed in the step plot.If the    \n");
        YvyraPrint ("                   value is set to 0, then the initial sample from the posterior \n");
        YvyraPrint ("                   will be used.                                                 \n");
        YvyraPrint ("   Smoothing    -- Determines smoothing of the joined plot (see above). A value  \n");
        YvyraPrint ("                   equal to 0 results in no smoothing.                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter       Options                  Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Allruns         Yes/No                   %s                                   \n", sumssParams.allRuns == YES ? "Yes" : "No");
        YvyraPrint ("   Askmore         Yes/No                   %s                                   \n", sumssParams.askForMorePlots == YES ? "Yes" : "No");
        YvyraPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        YvyraPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        YvyraPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        YvyraPrint ("   Discardfrac     <number>                 %1.2lf                               \n", sumssParams.discardFraction);
        if (sumpParams.numRuns == 1)
            YvyraPrint ("   Filename        <name>                   %s<.p>\n", sumpParams.sumpFileName);
        else
            YvyraPrint ("   Filename        <name>                   %s<.run<i>.p>\n", sumpParams.sumpFileName);        
        YvyraPrint ("   Nruns           <number>                 %d                                   \n", sumpParams.numRuns);
        YvyraPrint ("   Steptoplot      <number>                 %d                                   \n", sumssParams.stepToPlot);
        YvyraPrint ("   Smoothing       <number>                 %d                                   \n", sumssParams.smoothing);
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Comparetree"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Comparetree                                                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command compares the trees in two files, called \"filename1\" and        \n");
        YvyraPrint ("   \"filename2\". It will output a bivariate plot of the split frequencies       \n");
        YvyraPrint ("   as well as plots of the tree distance as a function of the generation. The    \n");
        YvyraPrint ("   plots can be used to get a quick indication of whether two runs have con-     \n");
        YvyraPrint ("   verged onto the same set of trees. The \"Comparetree\" command will also      \n");
        YvyraPrint ("   produce a \".pairs\" file and a \".dists\" file (these file endings are added \n");
        YvyraPrint ("   to the end of the \"Outputname\"). The \".pairs\" file contains the paired    \n");
        YvyraPrint ("   split frequencies from the two tree samples; the \".dists\" file contains the \n");
        YvyraPrint ("   tree distance values.                                                         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Note that the \"Sumt\" command provides a different set of convergence diag-  \n");
        YvyraPrint ("   nostics tools that you may also want to explore. Unlike \"Comparetree\",      \n");
        YvyraPrint ("   \"Sumt\" can compare more than two tree samples and will calculate consensus  \n");
        YvyraPrint ("   trees and split frequencies from the pooled samples.                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Relburnin     -- If this option is set to 'Yes', then a proportion of the     \n");
        YvyraPrint ("                    samples will be discarded as burnin when calculating summary \n");
        YvyraPrint ("                    statistics. The proportion to be discarded is set with       \n");
        YvyraPrint ("                    Burninfrac (see below). When the Relburnin option is set to  \n");
        YvyraPrint ("                    'No', then a specific number of samples is discarded instead.\n");
        YvyraPrint ("                    This number is set by Burnin (see below). Note that the      \n");
        YvyraPrint ("                    burnin setting is shared with the 'Mcmc', 'Sumt', 'Sump' and \n");
        YvyraPrint ("                    'Plot' commands.                                             \n");
        YvyraPrint ("   Burnin        -- Determines the number of samples (not generations) that will \n");
        YvyraPrint ("                    be discarded when summary statistics are calculated. The     \n");
        YvyraPrint ("                    value of this option is only relevant when Relburnin is set  \n");
        YvyraPrint ("                    to 'No'.                                                     \n");
        YvyraPrint ("   BurninFrac    -- Determines the fraction of samples that will be discarded    \n");
        YvyraPrint ("                    when summary statistics are calculated. The value of this    \n");
        YvyraPrint ("                    option is only relevant when Relburnin is set to 'Yes'.      \n");
        YvyraPrint ("                    Example: A value for this option of 0.25 means that 25%% of  \n");
        YvyraPrint ("                    the samples will be discarded.                               \n");
        YvyraPrint ("   Minpartfreq   -- The minimum probability of partitions to include in summary  \n");
        YvyraPrint ("                    statistics.                                                  \n");
        YvyraPrint ("   Filename1     -- The name of the first tree file to compare.                  \n");
        YvyraPrint ("   Filename2     -- The name of the second tree file to compare.                 \n");
        YvyraPrint ("   Outputname    -- Name of the file to which 'comparetree' results will be      \n");
        YvyraPrint ("                    printed.                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter       Options                  Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        YvyraPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        YvyraPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        YvyraPrint ("   Minpartfreq     <number>                 %1.2lf                               \n", comptreeParams.minPartFreq);
        YvyraPrint ("   Filename1       <name>                   %s                                   \n", comptreeParams.comptFileName1);
        YvyraPrint ("   Filename2       <name>                   %s                                   \n", comptreeParams.comptFileName2);
        YvyraPrint ("   Outputname      <name>                   %s                                   \n", comptreeParams.comptOutfile);
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Sumt"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Sumt                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used to produce summary statistics for trees sampled during   \n");
        YvyraPrint ("   a Bayesian MCMC analysis. You can either summarize trees from one individual  \n");
        YvyraPrint ("   analysis, or trees coming from several independent analyses. In either case,  \n");
        YvyraPrint ("   all the sampled trees are read in and the proportion of the time any single   \n");
        YvyraPrint ("   taxon bipartition (split) is found is counted. The proportion of the time that\n");
        YvyraPrint ("   the bipartition is found is an approximation of the posterior probability of  \n");
        YvyraPrint ("   the bipartition. (Remember that a taxon bipartition is defined by removing a  \n");
        YvyraPrint ("   branch on the tree, dividing the tree into those taxa to the left and right   \n");
        YvyraPrint ("   of the removed branch. This set is called a taxon bipartition.) The branch    \n");
        YvyraPrint ("   length of the bipartition is also recorded, if branch lengths have been saved \n");
        YvyraPrint ("   to file. The result is a list of the taxon bipartitions found, the frequency  \n");
        YvyraPrint ("   with which they were found, the posterior probability of the bipartition      \n");
        YvyraPrint ("   and, the mean and variance of the branch lengths or node depths, and various  \n");
        YvyraPrint ("   other statistics.                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The key to the partitions is output to a file with the suffix '.parts'. The   \n");
        YvyraPrint ("   summary statistics pertaining to bipartition probabilities are output to a    \n");
        YvyraPrint ("   file with the suffix '.tstat', and the statistics pertaining to branch or node\n");
        YvyraPrint ("   parameters are output to a file with the suffix '.vstat'.                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   A consensus tree is also printed to a file with the suffix '.con.tre' and     \n");
        YvyraPrint ("   printed to the screen as a cladogram, and as a phylogram if branch lengths    \n");
        YvyraPrint ("   have been saved. The consensus tree is either a 50 percent majority rule tree \n");
        YvyraPrint ("   or a majority rule tree showing all compatible partitions. If branch lengths  \n");
        YvyraPrint ("   have been recorded during the run, the '.con.tre' file will contain a consen- \n");
        YvyraPrint ("   sus tree with branch lengths and interior nodes labelled with support values. \n");
        YvyraPrint ("   By default, the consensus tree will also contain other summary information in \n");
        YvyraPrint ("   a format understood by the program 'FigTree'. To use a simpler format under-  \n");
        YvyraPrint ("   stood by other tree-drawing programs, such as 'TreeView', set 'Conformat' to  \n");
        YvyraPrint ("   'Simple'.                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   yvyra also produces a file with the ending \".trprobs\" that contains a list\n");
        YvyraPrint ("   of all the trees that were found during the MCMC analysis, sorted by their    \n");
        YvyraPrint ("   probabilities. This list of trees can be used to construct a credible set of  \n");
        YvyraPrint ("   trees. For example, if you want to construct a 95 percent credible set of     \n");
        YvyraPrint ("   trees, you include all of those trees whose cumulative probability is less    \n");
        YvyraPrint ("   than or equal to 0.95. You have the option of displaying the trees to the     \n");
        YvyraPrint ("   screen using the \"Showtreeprobs\" option. The default is to not display the  \n");
        YvyraPrint ("   trees to the screen; the number of different trees sampled by the chain can   \n");
        YvyraPrint ("   be quite large. If you are analyzing a large set of taxa, you may actually    \n");
        YvyraPrint ("   want to skip the calculation of tree probabilities entirely by setting        \n");
        YvyraPrint ("   'Calctreeprobs' to 'No'.                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   When calculating summary statistics you probably want to skip those trees that\n");
        YvyraPrint ("   were sampled in the initial part of the run, the so-called burn-in period. The\n");
        YvyraPrint ("   number of skipped samples is controlled by the 'Relburnin', 'Burnin', and     \n");
        YvyraPrint ("   'Burninfrac' settings, just as for the 'Mcmc' command. Since version 3.2.0,   \n");
        YvyraPrint ("   the burn-in settings are shared across the 'Sumt', 'Sump' and 'Mcmc' commands.\n");
        YvyraPrint ("   That is, changing the burn-in setting for one command will change the settings\n");
        YvyraPrint ("   for subsequent calls to any of the other commands.                            \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   If you are summarizing the trees sampled in several independent analyses,     \n");
        YvyraPrint ("   such as those resulting from setting the 'Nruns' option of the 'Mcmc' command \n");
        YvyraPrint ("   to a value larger than 1, yvyra will also calculate convergence diagnostics \n");
        YvyraPrint ("   for the sampled topologies and branch lengths. These values can help you      \n");
        YvyraPrint ("   determine whether it is likely that your chains have converged.               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The 'Sumt' command expands the 'Filename' according to the current values of  \n");
        YvyraPrint ("   the 'Nruns' and 'Ntrees' options. For instance, if both 'Nruns' and 'Ntrees'  \n");
        YvyraPrint ("   are set to 1, 'Sumt' will try to open a file named '<Filename>.t'. If 'Nruns' \n");
        YvyraPrint ("   is set to 2 and 'Ntrees' to 1, then 'Sumt' will open two files, the first     \n");
        YvyraPrint ("   named '<Filename>.run1.t' and the second '<Filename>.run2.t', etc. By default,\n");
        YvyraPrint ("   the 'Filename' option is set such that 'Sumt' automatically summarizes all the\n");
        YvyraPrint ("   results from your immediately preceding 'Mcmc' command. You can also use the  \n");
        YvyraPrint ("   'Sumt' command to summarize tree samples in older analyses. If you want to do \n");
        YvyraPrint ("   that, remember to first read in a matrix so that yvyra knows what taxon     \n");
        YvyraPrint ("   names to expect in the trees. Then set the 'Nruns', 'Ntrees' and 'Filename'   \n");
        YvyraPrint ("   options appropriately if they differ from the yvyra defaults.               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Relburnin     -- If this option is set to YES, then a proportion of the       \n");
        YvyraPrint ("                    samples will be discarded as burnin when calculating summary \n");
        YvyraPrint ("                    statistics. The proportion to be discarded is set with       \n");
        YvyraPrint ("                    Burninfrac (see below). When the Relburnin option is set to  \n");
        YvyraPrint ("                    NO, then a specific number of samples is discarded instead.  \n");
        YvyraPrint ("                    This number is set by Burnin (see below). Note that the      \n");
        YvyraPrint ("                    burnin setting is shared across the 'Sumt', 'Sump', and      \n");
        YvyraPrint ("                    'Mcmc' commands.                                             \n");
        YvyraPrint ("   Burnin        -- Determines the number of samples (not generations) that will \n");
        YvyraPrint ("                    be discarded when summary statistics are calculated. The     \n");
        YvyraPrint ("                    value of this option is only relevant when Relburnin is set  \n");
        YvyraPrint ("                    to NO.                                                       \n");
        YvyraPrint ("   BurninFrac    -- Determines the fraction of samples that will be discarded    \n");
        YvyraPrint ("                    when summary statistics are calculated. The value of this    \n");
        YvyraPrint ("                    option is only relevant when Relburnin is set to YES.        \n");
        YvyraPrint ("                    Example: A value for this option of 0.25 means that 25%% of  \n");
        YvyraPrint ("                    the samples will be discarded.                               \n");
        YvyraPrint ("   Nruns         -- Determines how many '.t' files from independent analyses that\n");
        YvyraPrint ("                    will be summarized. If Nruns > 1 then the names of the files \n");
        YvyraPrint ("                    are derived from 'Filename' by adding '.run1.t', '.run2.t',  \n");
        YvyraPrint ("                    etc. If Nruns=1 and Ntrees=1 (see below), then only '.t' is  \n");
        YvyraPrint ("                    added to 'Filename'.                                         \n");
        YvyraPrint ("   Ntrees        -- Determines how many trees there are in the sampled model. If \n");
        YvyraPrint ("                    'Ntrees' > 1 then the names of the files are derived from    \n");
        YvyraPrint ("                    'Filename' by adding '.tree1.t', '.tree2.t', etc. If there   \n");
        YvyraPrint ("                    are both multiple trees and multiple runs, the filenames will\n");
        YvyraPrint ("                    be '<Filename>.tree1.run1.t', '<Filename>.tree1.run2.t', etc.\n");
        YvyraPrint ("   Filename      -- The name of the file(s) to be summarized. This is the base of\n");
        YvyraPrint ("                    the file name, to which endings are added according to the   \n");
        YvyraPrint ("                    current settings of the 'Nruns' and 'Ntrees' options.        \n");
        YvyraPrint ("   Minpartfreq   -- The minimum probability of partitions to include in summary  \n");
        YvyraPrint ("                    statistics.                                                  \n");
        YvyraPrint ("   Contype       -- Type of consensus tree. 'Halfcompat' results in a 50%% major-\n");
        YvyraPrint ("                    ity rule tree, 'Allcompat' adds all compatible groups to such\n");
        YvyraPrint ("                    a tree.                                                      \n");
        YvyraPrint ("   Conformat     -- Format of consensus tree. The 'Figtree' setting results in a \n");
        YvyraPrint ("                    consensus tree formatted for the program FigTree, with rich  \n");
        YvyraPrint ("                    summary statistics. The 'Simple' setting results in a simple \n");
        YvyraPrint ("                    consensus tree written in a format read by a variety of pro- \n");
        YvyraPrint ("                    grams.                                                       \n");
        YvyraPrint ("   Outputname    -- Base name of the file(s) to which 'sumt' results will be     \n");
        YvyraPrint ("                    printed. The default is the same as 'Filename'.              \n");
        YvyraPrint ("   Calctreeprobs -- Determines whether tree probabilities should be calculated.  \n");
        YvyraPrint ("   Showtreeprobs -- Determines whether tree probabilities should be displayed on \n");
        YvyraPrint ("                    screen.                                                      \n");
        YvyraPrint ("   Hpd           -- Determines whether credibility intervals will be given as the\n");
        YvyraPrint ("                    region of Highest Posterior Density ('Yes') or as the inter- \n");
        YvyraPrint ("                    val containing the median 95 %% of sampled values ('No').    \n");
        YvyraPrint ("   Savebrparams  -- Set this option to 'yes' to save all sampled branch and node \n");
        YvyraPrint ("                    parameter values to a separate file with the filename ending \n");
        YvyraPrint ("                    in '.brparams'. All partitions with a posterior probability  \n");
        YvyraPrint ("                    larger than Minbrparamfreq will be included.                 \n");
        YvyraPrint ("   Minbrparamfreq -- The minimum probability of partitions for which to save     \n");
        YvyraPrint ("                     parameter values to file if 'Savebrparams' is set to 'yes'. \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Current settings:                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter       Options                  Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        YvyraPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        YvyraPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        YvyraPrint ("   Nruns           <number>                 %d                                   \n", sumtParams.numRuns);
        YvyraPrint ("   Ntrees          <number>                 %d                                   \n", sumtParams.numTrees);
        if (sumtParams.numRuns == 1 && sumtParams.numTrees == 1)
            YvyraPrint ("   Filename        <name>                   %s<.t>\n", sumtParams.sumtFileName);
        else if (sumtParams.numRuns == 1 && sumtParams.numTrees > 1)
            YvyraPrint ("   Filename        <name>                   %s<.tree<i>.t>\n", sumtParams.sumtFileName);
        else if (sumtParams.numRuns > 1 && sumtParams.numTrees == 1)
            YvyraPrint ("   Filename        <name>                   %s<.run<i>.t>\n", sumtParams.sumtFileName);
        else if (sumtParams.numRuns > 1 && sumtParams.numTrees > 1)
            YvyraPrint ("   Filename        <name>                   %s<.tree<i>.run<i>.t>\n", sumtParams.sumtFileName);
        YvyraPrint ("   Minpartfreq     <number>                 %1.2lf                               \n", sumtParams.minPartFreq);
        YvyraPrint ("   Contype         Halfcompat/Allcompat     %s\n", sumtParams.sumtConType);
        YvyraPrint ("   Conformat       Figtree/Simple           %s                                   \n", sumtParams.consensusFormat == SIMPLE ? "Simple" : "Figtree");
        YvyraPrint ("   Outputname      <name>                   %s<.parts etc>\n", sumtParams.sumtOutfile);
        YvyraPrint ("   Calctreeprobs   Yes/No                   %s                                   \n", sumtParams.calcTreeprobs == YES ? "Yes" : "No");
        YvyraPrint ("   Showtreeprobs   Yes/No                   %s                                   \n", sumtParams.showSumtTrees == YES ? "Yes" : "No");
        YvyraPrint ("   Hpd             Yes/No                   %s                                   \n", sumtParams.HPD == YES ? "Yes" : "No");
        YvyraPrint ("   Savebrparams    Yes/No                   %s                                   \n", sumtParams.saveBrParams == YES ? "Yes" : "No");
        YvyraPrint ("   Minbrparamfreq  <number>                 %1.2lf                               \n", sumtParams.minBrParamFreq);
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Tree"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Tree                                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command is used by yvyra to write trees to a nexus tree file. Trees    \n");
        YvyraPrint ("   are written in the Newick format. For instance,                               \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      tree ((1,2),3,4);                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   describes an unrooted tree with taxa 1 and 2 being more closely related to    \n");
        YvyraPrint ("   each other than to taxa 3 and 4. If branch lengths are saved to file, they    \n");
        YvyraPrint ("   are given after a colon sign immediately following the terminal taxon or the  \n");
        YvyraPrint ("   interior node they refer to. An example of an unrooted tree with branch       \n");
        YvyraPrint ("   lengths is:                                                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      tree ((1:0.064573,2:0.029042):0.041239,3:0.203988,4:0.187654);             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Trees that are rooted (clock trees) are written with a basal dichotomy        \n");
        YvyraPrint ("   instead of a basal trichotomy. If the tree described above had been rooted    \n");
        YvyraPrint ("   on the branch leading to taxon 4, it would have been represented as:          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      tree (((1,2),3),4);                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Report"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Report                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command allows you to control how the posterior distribution is          \n");
        YvyraPrint ("   reported. For rate parameters, it allows you to choose among several popular  \n");
        YvyraPrint ("   parameterizations. The report command also allows you to request printing of  \n");
        YvyraPrint ("   some model aspects that are usually not reported. For instance, if a node is  \n");
        YvyraPrint ("   constrained in the analysis, yvyra can print the probabilities of the       \n");
        YvyraPrint ("   ancestral states at that node. Similarly, if there is rate variation in the   \n");
        YvyraPrint ("   model, yvyra can print the inferred site rates, and if there is omega varia-\n");
        YvyraPrint ("   tion, yvyra can print the inferred omega (positive selection) values for    \n");
        YvyraPrint ("   each codon. In a complex model with several partitions, each partition is     \n");
        YvyraPrint ("   controlled separately using the same 'Applyto' mechanism as in the 'Lset' and \n");
        YvyraPrint ("   'Prset' commands.                                                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Options:                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Applyto   -- This option allows you to apply the report commands to specific  \n");
        YvyraPrint ("                partitions. This command should be the first in the list of      \n");
        YvyraPrint ("                commands specified in 'report'.                                  \n");
        YvyraPrint ("                For example,                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                   report applyto=(1,2) tratio=ratio                             \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                   report applyto=(3) tratio=dirichlet                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("                would result in the transition and transversion rates of the     \n");
        YvyraPrint ("                first and second partitions in the model being reported as a     \n");
        YvyraPrint ("                ratio and the transition and transversion rates of the third     \n");
        YvyraPrint ("                partition being reported as proportions of the rate sum (the     \n");
        YvyraPrint ("                Dirichlet parameterization).                                     \n");
        YvyraPrint ("   Tratio    -- This specifies the report format for the transition and trans-   \n");
        YvyraPrint ("                version rates of a nucleotide substitution model with nst=2.     \n");
        YvyraPrint ("                If 'ratio' is selected, the rates will be reported as a ratio    \n");
        YvyraPrint ("                (transition rate/transversion rate). If 'dirichlet' is selected, \n");
        YvyraPrint ("                the transition and transversion rates will instead be reported   \n");
        YvyraPrint ("                as proportions of the rate sum. For example, if the transition   \n");
        YvyraPrint ("                rate is three times the transversion rate and 'ratio' is selec-  \n");
        YvyraPrint ("                ted, this will reported as a single value, '3.0'. If 'dirichlet' \n");
        YvyraPrint ("                is selected instead, the same rates will be reported using two   \n");
        YvyraPrint ("                values, '0.75 0.25'. The sum of the Dirichlet values is always 1.\n");
        YvyraPrint ("                Although the Dirichlet format may be unfamiliar to some users,   \n");
        YvyraPrint ("                it is more convenient for specifying priors than the ratio       \n");
        YvyraPrint ("                format.                                                          \n");
        YvyraPrint ("   Revmat    -- This specifies the report format for the substitution rates of   \n");
        YvyraPrint ("                a GTR substitution model for nucleotide or amino acid data. If   \n");
        YvyraPrint ("                'ratio' is selected, the rates will be reported scaled to the    \n");
        YvyraPrint ("                G-T rate (for nucleotides) or the Y-V rate (for amino acids). If \n");
        YvyraPrint ("                'dirichlet' is specified instead, the rates are reported as pro- \n");
        YvyraPrint ("                portions of the rate sum. For instance, assume that the C-T rate \n");
        YvyraPrint ("                is twice the A-G rate and four times the transversion rates,     \n");
        YvyraPrint ("                which are equal. If the report format is set to 'ratio', this    \n");
        YvyraPrint ("                would be reported as '1.0 2.0 1.0 1.0 4.0 1.0' since the rates   \n");
        YvyraPrint ("                are reported in the order rAC, rAG, rAT, rCG, rCT, rGT and scaled\n");
        YvyraPrint ("                relative to the last rate, the G-T rate. If 'dirichlet' is selec-\n");
        YvyraPrint ("                ted instead, the same rates would have been reported as '0.1 0.2 \n");
        YvyraPrint ("                0.1 0.1 0.4 0.1' since the rates are now scaled so that they sum \n");
        YvyraPrint ("                to 1.0. The Dirichlet format is the parameterization used for    \n");
        YvyraPrint ("                formulating priors on the rates.                                 \n");
        YvyraPrint ("   Ratemult  -- This specifies the report format used for the rate multiplier of \n");
        YvyraPrint ("                different model partitions. Three formats are available. If      \n");
        YvyraPrint ("                'scaled' is selected, then rates are scaled such that the mean   \n");
        YvyraPrint ("                rate per site across partitions is 1.0. If 'ratio' is chosen,    \n");
        YvyraPrint ("                the rates are scaled relative to the rate of the first parti-    \n");
        YvyraPrint ("                tion. Finally, if 'dirichlet' is chosen, the rates are given as  \n");
        YvyraPrint ("                proportions of the rate sum. The latter is the format used       \n");
        YvyraPrint ("                when formulating priors on the rate multiplier.                  \n");
        YvyraPrint ("   Tree      -- This specifies the report format used for the tree(s). Two op-   \n");
        YvyraPrint ("                tions are available. 'Topology' results in only the topology     \n");
        YvyraPrint ("                being printed to file, whereas 'brlens' causes branch lengths to \n");
        YvyraPrint ("                to be printed as well.                                           \n");
        YvyraPrint ("   Ancstates -- If this option is set to 'yes', yvyra will print the pro-      \n");
        YvyraPrint ("                bability of the ancestral states at all constrained nodes. Typ-  \n");
        YvyraPrint ("                ically, you are interested in the ancestral states of only a few \n");
        YvyraPrint ("                characters and only at one node in the tree. To perform such     \n");
        YvyraPrint ("                an analysis, first define and enforce a topology constraint      \n");
        YvyraPrint ("                using 'constraint' and 'prset topologypr = constraints (...)'.   \n");
        YvyraPrint ("                Then put the character(s) of interest in a separate partition and\n");
        YvyraPrint ("                set yvyra to report the ancestral states for that partition.   \n");
        YvyraPrint ("                For instance, if the characters of interest are in partition 2,  \n");
        YvyraPrint ("                use 'report applyto=(2) ancstates=yes' to force yvyra to print \n");
        YvyraPrint ("                the probability of the ancestral states of those characters at   \n");
        YvyraPrint ("                the constrained node to the '.p' file.                           \n");
        YvyraPrint ("   Siterates -- If this option is set to 'yes' and the relevant model has rate   \n");
        YvyraPrint ("                variation across sites, then the site rates, weighted over rate  \n");
        YvyraPrint ("                categories, will be reported to the '.p' file.                   \n");
        YvyraPrint ("   Possel    -- If this option is set to 'yes' and the relevant model has omega  \n");
        YvyraPrint ("                variation across sites, the probability that each model site     \n");
        YvyraPrint ("                (codon in this case) is positively selected will be written to   \n");
        YvyraPrint ("                file.                                                            \n");
        YvyraPrint ("   Siteomega -- If this option is set to 'yes' and the relevant model has omega  \n");
        YvyraPrint ("                variation across sites, the weighted omega value (over omega     \n");
        YvyraPrint ("                categories) for each model site will be reported to file.        \n");
        YvyraPrint ("                                                                                 \n");
        if (numCurrentDivisions == 0)
            tempInt = 1;
        else
            tempInt = numCurrentDivisions;
        for (i=0; i<tempInt; i++)
            {
            if (numCurrentDivisions == 0)
                {
                YvyraPrint ("   Default report settings:                                                       \n");
                mp = &defaultModel;
                }
            else
                {
                YvyraPrint ("   Current report settings for partition %d:                                              \n", i+1);
                mp = &modelParams[i];
                }
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   Parameter       Options                  Current Setting                      \n");
            YvyraPrint ("   --------------------------------------------------------                      \n");
            YvyraPrint ("   Tratio          Ratio/Dirichlet          %s                                   \n", mp->tratioFormat);
            YvyraPrint ("   Revmat          Ratio/Dirichlet          %s                                   \n", mp->revmatFormat);
            YvyraPrint ("   Ratemult        Scaled/Ratio/Dirichlet   %s                                   \n", mp->ratemultFormat);
            YvyraPrint ("   Tree            Brlens/Topology          %s                                   \n", mp->treeFormat);
            YvyraPrint ("   Ancstates       Yes/No                   %s                                   \n", mp->inferAncStates);
            YvyraPrint ("   Siterates       Yes/No                   %s                                   \n", mp->inferSiteRates);
            YvyraPrint ("   Possel          Yes/No                   %s                                   \n", mp->inferPosSel);
            YvyraPrint ("   Siteomega       Yes/No                   %s                                   \n", mp->inferSiteOmegas);
            YvyraPrint ("                                                                                 \n");
            YvyraPrint ("   ------------------------------------------------------------------            \n");       
            }
        }
    else if (!strcmp(helpTkn, "Manual"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Manual                                                                        \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command allows you to generate a text file containing help information   \n");
        YvyraPrint ("   on all the available commands. This text file can be used as an up-to-date    \n");
        YvyraPrint ("   command reference. You can set the name of the text file using the            \n");
        YvyraPrint ("   \"filename\" option; the default is \"commref_mb<version>.txt\".              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Parameter       Options                  Current Setting                      \n");
        YvyraPrint ("   --------------------------------------------------------                      \n");
        YvyraPrint ("   Filename        <name>                   %s                                   \n", manFileName);
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmoves"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Showmoves                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows the MCMC samplers (moves) that are switched on for the     \n");
        YvyraPrint ("   parameters in the current model. The basic usage is                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      showmoves                                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   If you want to see all available moves, use                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      showmoves allavailable=yes                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   If you want to change any of the tuning parameters for the moves, use the     \n");
        YvyraPrint ("   'propset' command.                                                            \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showparams"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Showparams                                                                    \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This command shows all of the parameters in the current model. The basic      \n");
        YvyraPrint ("   usage is                                                                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      showparams                                                                 \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The parameters are listed together with their priors, the available moves,    \n");
        YvyraPrint ("   and the current value(s), which will be used as the starting values in the    \n");
        YvyraPrint ("   next MCMC analysis.                                                           \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Startvals"))
        {
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        YvyraPrint ("   Startvals                                                                     \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   Use this command to change the current values for parameters in your model.   \n");
        YvyraPrint ("   These values will be used as the starting values in the next MCMC analysis.   \n");
        YvyraPrint ("   The basic format is:                                                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      startvals <param>=(<value_1>,<value_2>,...,<value_n>)                      \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   for all substitution model parameters. The format is slightly different for   \n");
        YvyraPrint ("   parameters that are written to a tree file:                                   \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      startvals <param>=<tree_name>                                              \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   This version of the command will look for a tree with the specified name      \n");
        YvyraPrint ("   among the trees read in previously when parsing a tree block. The information \n");
        YvyraPrint ("   stored in that tree will be used to set the starting value of the parameter.  \n");
        YvyraPrint ("   The parameters that are set using this mechanism include topology and branch  \n");
        YvyraPrint ("   length parameters, as well as relaxed clock branch rates, CPP events and      \n");
        YvyraPrint ("   CPP branch rate multipliers.                                                  \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   The above versions of the command will set the value for all runs and chains. \n");
        YvyraPrint ("   You can also set the value for an individual run and chain by using the format\n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      startvals <param>(<run>,<chain>)=(<value_1>,...)                           \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   where <run> is the index of the run and <chain> the index of the chain. If    \n");
        YvyraPrint ("   the run index is omitted, the values will be changed for all runs. Similarly, \n");
        YvyraPrint ("   if the chain index is omitted, all chains will be set to the specified value. \n");
        YvyraPrint ("   For example, if we wanted to set the values of the stationary frequency       \n");
        YvyraPrint ("   parameter pi{1} to (0.1,0.1,0.4,0.4) for all chains in run 1, and to          \n");
        YvyraPrint ("   (0.3,0.3,0.2,0.2) for chain 3 of run 2, we would use                          \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("      startvals pi{1}(1,)=(0.1,0.1,0.4,0.4) pi{1}(2,3)=(0.3,0.3,0.2,0.2)         \n");
        YvyraPrint ("                                                                                 \n");
        YvyraPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else
        {
        return (ERROR);
        }
        
    return (NO_ERROR);
}


/* IsAmbig: This function returns YES if character is set as ambiguous
   either by using parenthetic notation or by ambiguity codes. It returns
   NO if character is unambiguous, missing or gapped */ 
int IsAmbig (int charCode, int dType)
{
    if (dType == DNA || dType == RNA || dType == STANDARD || dType == RESTRICTION || dType == PROTEIN)
        {
        if (charCode != MISSING && charCode != GAP)
            if (NBits(charCode) > 1)
                return (YES);
        }
    else if (dType == CONTINUOUS)
        {
        /* do nothing, these cannot be partly ambiguous */
        }
    else
        {
        YvyraPrint ("Unknown datatype in \"IsAmbig\"\n", spacer);
        }

    return (NO);
}


int IsArgValid (char *tk, char *validArg)
{
    int         i, j, k, tkLen, targetLen, numDiff, numStrMatches;
    char        tempStr[100];
    ParmInfoPtr p;

    p = paramPtr;
    tkLen = (int) strlen(tk);

    numStrMatches = i = j = 0;
    do
        {
        if (p->valueList[i] == '|' || p->valueList[i] == '\0')
            {
            tempStr[j++] = '\0';
            targetLen = (int) strlen(tempStr);
            if (tkLen <= targetLen)
                {
                numDiff = 0;
                for (k=0; k<tkLen; k++)
                    if (ChangeCase(tk[k]) != ChangeCase(tempStr[k]))
                        numDiff++;
                if (numDiff == 0)
                    {
                    numStrMatches++;
                    strcpy (validArg, tempStr);
                    }
                }
            j = 0;
            }
        else
            tempStr[j++] = p->valueList[i];
        i++;
        }
    while (p->valueList[i] != '\0');
        
    if (numStrMatches == 0)
        {
        YvyraPrint ("%s   No valid match for argument \"%s\"\n", spacer, tk);
        return (ERROR);
        }
    else if (numStrMatches == 1)
        {
        return (NO_ERROR);
        }
    else
        {
        YvyraPrint ("%s   Argument \"%s\" is ambiguous\n", spacer, tk);
        return (ERROR);
        }
}


int IsIn (char ch, char *s)
{
    while (*s)
        {
        if (*s++ == ch)
            return 1;
        }
    return 0;
}


int IsMissing (int charCode, int dType)
{
    if (dType == DNA || dType == RNA)
        {
        if (charCode == 15 || charCode == 16)
            return (YES);
        }
    else if (dType == STANDARD || dType == PROTEIN)
        {
        if (charCode == MISSING || charCode == GAP)
            return (YES);
        }
    else if (dType == RESTRICTION)
        {
        if (charCode == 3 || charCode == 4)
            return (YES);
        }
    else if (dType == CONTINUOUS)
        {

        }
    else
        {
        YvyraPrint ("Unknown datatype in \"IsMissing\"\n", spacer);
        }
    return (NO);
}


int IsSame (char *s1, char *s2)
{
    int         i, nDiff, isIdentical, len;
    
    isIdentical = YES;
    if (strlen(s1) != strlen(s2))
        isIdentical = NO; /* strings cannot be identical because they are different lengths */
    
    /* now, we go through both strings, one character at a time, to see if
       any are different */
    if (strlen(s1) > strlen(s2))
        len = (int) strlen(s2);
    else
        len = (int) strlen(s1);
    i = nDiff = 0;
    while (i < len)
        {
        if (tolower(s1[i]) != tolower(s2[i]))
            nDiff++;
        i++;
        }
    if (nDiff == 0 && isIdentical == YES)
        return (SAME);
    else if (nDiff == 0 && isIdentical == NO)
        return (CONSISTENT_WITH);
    else
        return (DIFFERENT);
}


int IsWhite (char c)
{
    if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
        {
        if (c == '\n' || c == '\r')
            return 2;
        return 1;
        }
    return 0;
}


int NucID (char nuc)
{
    char        n;
    
    if (nuc == 'U' || nuc == 'u')
        n = 'T';
    else
        n = nuc;

    if (n == 'A' || n == 'a')
        {
        return 1;
        }
    else if (n == 'C' || n == 'c')
        {
        return 2;
        }
    else if (n == 'G' || n == 'g')
        {
        return 4;
        }
    else if (n == 'T' || n == 't')
        {
        return 8;
        }
    else if (n == 'R' || n == 'r')
        {
        return 5;
        }
    else if (n == 'Y' || n == 'y')
        {
        return 10;
        }
    else if (n == 'M' || n == 'm')
        {
        return 3;
        }
    else if (n == 'K' || n == 'k')
        {
        return 12;
        }
    else if (n == 'S' || n == 's')
        {
        return 6;
        }
    else if (n == 'W' || n == 'w')
        {
        return 9;
        }
    else if (n == 'H' || n == 'h')
        {
        return 11;
        }
    else if (n == 'B' || n == 'b')
        {
        return 14;
        }
    else if (n == 'V' || n == 'v')
        {
        return 7;
        }
    else if (n == 'D' || n == 'd')
        {
        return 13;
        }
    else if (n == 'N' || n == 'n')
        {
        return 15;
        }
    else if (n == gapId)
        {
        return GAP;
        }
    else if (n == missingId)
        {
        return MISSING;
        }
    else
        return -1;
}


/*-------| ParseCommand |------------------------------------------------
|
|   This function is used to parse a file. The expected format is:
|   
|      command parameter=value parameter=value ... ;
|
|   For example, the following is a valid line for this parser:
|
|      lset nst=2;
|
|   In some cases, however, the format is:
|
|      command stuff more stuff ... ;
|
|   For example, when reading a data file, the matrix command might be:
|
|      matrix
|         taxon_1 data
|         taxon_2 data
|         taxon_3 data
|         ;
|
|   Like before, the command and all of the stuff for that command are
|   terminated by a semicolon.
|
*/
int ParseCommand (char *s)
{
    int             rc, tokenType, inError, numMatches, skipCmd;
    char            errStr[100];

    numMatches = 0;     /* Avoid gcc warnings (actually set in call to FindValidCommand) */
    cmdStr = s;
    tokenP = &s[0];

#   if defined (ECHO_PROCESSED_COMMANDS)
        YvyraPrint ("Currently processing command: %s\n", s);
#   endif
    
    inError = skipCmd = NO;
    do
        {
        /* Get the next token. A token is a valid word in a line. Token type is defined in "bayes.h". */
        if (GetToken (token, &tokenType, &tokenP))
            {
            inError = YES; 
            break;
            }
        if (strlen(token) > 0 || tokenType == ALPHA)
            {
#           if defined (SHOW_TOKENS)
            YvyraPrint ("%s\n", token);
#           endif
            if (tokenType == LEFTCOMMENT)
                {
                /* If the token is a left comment "[", then we don't want to
                   actually process commands until we find a right comment.  */
                /* The exception is if readComment is set to YES, in which case
                   we will leave it to the parser functions to decide on whether
                   they want to read the comment or not */
                if (readComment == NO || inComment == YES)
                    {
                    inComment = YES;
                    numComments++;
                    }
                }
            if (inComment == NO && inForeignBlock == NO)
                {
                if (tokenType != SEMICOLON)
                    {
                    /* If the token is not a semicolon, then we will be processing 
                       either a command or a parameter. */
                    if (expecting == Expecting(COMMAND))
                        {
                        /* We are expecting to find a command (defined above in "commands[]"). Find the 
                           correct command and set a pointer to that command. */
                        commandPtr = NULL;
                        if (FindValidCommand (token, &numMatches) == ERROR)
                            {
                            /* We couldn't find the command or the user did not specify enough letters
                               to unambiguously determine the command. The command pointer (commandPtr)
                               is NULL. */
                            if (numMatches == 0)    
                                YvyraPrint ("%s   Could not find command \"%s\"\n", spacer, token);
                            else 
                                YvyraPrint ("%s   Ambiguous command \"%s\"\n", spacer, token);
                            inError = YES;
                            }
                        else
                            {
                            /* We did find a valid command. Set what we are expecting to see next. */
                            expecting = commandPtr->expect;
                            
                            /* Check to see if we have one of the so-called special cases in which a 
                               command is not necessarily followed by a parameter (e.g., matrix). If we
                               do have a special case, then we want to set the parameter pointer (paramPtr)
                               appropriately. In this case, simply go to the first parameter in the parmList. */
                            if (commandPtr->specialCmd == YES)
                                {
                                isFirstMatrixRead = YES;
                                foundFirst = NO;
                                paramPtr = paramTable + commandPtr->parmList[0];
                                }
                            if (strcmp(commandPtr->string, "Execute")==0)
                                {
                                /* set the tokenizer to recognize quoted strings */
                                readWord = YES;
                                }
                            }
                        }
                    else 
                        {
                        /* We are expecting to find a parameter or a value for the parameter, not a command. */
                        if ((expecting & Expecting(PARAMETER)) == Expecting(PARAMETER) && 
                            (expecting & Expecting(tokenType)) != Expecting(tokenType))
                            {
                            /* Specifically, if we are here, we need to go through the parameter list,
                               checking to see if the token is a valid parameter. */
                            expecting = (expecting & Expecting(PARAMETER));
                            if (FindValidParam (token, &numMatches) == ERROR)
                                {
                                /* The token is not a valid parameter. */
                                if (numMatches == 0)
                                    YvyraPrint ("%s   Could not find parameter \"%s\"\n", spacer, token);
                                else 
                                    YvyraPrint ("%s   Ambiguous parameter \"%s\"\n", spacer, token);
                                inError = YES;
                                }
                            else
                                {
                                /* The token is a valid parameter. Call the appropriate function ("DoXxxxParm"). */
                                if ((paramPtr->fp)(paramPtr->string, token) == ERROR)
                                    {
                                    if (strcmp("Xxxxxxxxxx", paramPtr->string))
                                        YvyraPrint ("%s   Error when setting parameter \"%s\" (1)\n", spacer, paramPtr->string);
                                    inError = YES;
                                    }
                                }
                            }
                        else
                            {
                            /* Otherwise, we are expecting a value for the parameter. Call the appropriate function ("DoXxxxParm"). */
                            if ((expecting & Expecting(tokenType)) != 0)
                                expecting = (expecting & Expecting(tokenType));
                            if ((Expecting(tokenType) & expecting) == Expecting(tokenType))
                                {
                                if ((paramPtr->fp)(paramPtr->string, token) == ERROR)
                                    {
                                    if (strcmp("Xxxxxxxxxx", paramPtr->string))
                                        YvyraPrint ("%s   Error when setting parameter \"%s\" (2)\n", spacer, paramPtr->string);
                                    inError = YES;
                                    }
                                }
                            else
                                {
                                inError = YES;
                                WhatVariableExp (expecting, errStr);
                                YvyraPrint ("%s   Expecting '%s'\n", spacer, errStr+1);  /* there will be an initial space in errStr so print from pos 1 */
                                if (numOpenExeFiles > 0)
                                    YvyraPrint ("%s   Instead found '%s' in command '%s'\n",
                                        spacer, token, commandPtr->string);
                                else
                                    YvyraPrint ("%s   Instead found '%s' in command '%s' at position %d\n",
                                        spacer, token, commandPtr->string, tokenP - cmdStr - strlen(token)+1);
                                }
                            }
                        }
                    }
                else
                    {
                    /* The token is a semicolon. This means that we are at the end of processing one command. We
                       need to clean things up. We do this by calling the finishing function ("DoXxxx"). */
                    if ((Expecting(SEMICOLON) & expecting) == Expecting(SEMICOLON))
                        {
                        if (commandPtr->cmdFxnPtr != NULL)
                            {
                            /* Finish up the command here. */
                            rc = (commandPtr->cmdFxnPtr) ();
                            if (rc  == ERROR || rc == ABORT)
                                {
                                if (rc == ABORT)
                                    {
                                    YvyraPrint ("   Mcmc run aborted\n");
                                    }
                                else if (rc == SKIP_COMMAND)
                                    {
                                    YvyraPrint ("   Cancelled execution of command\n");
                                    skipCmd = YES;
                                    }
                                else
                                    {
                                    YvyraPrint ("%s   Error in command \"%s\"\n", spacer, commandPtr->string);
                                    inError = YES;
                                    }
                                }
                            }
                        /* if the user typed "quit", then we want to bail out of this loop, with a NO_ERROR_QUIT */
                        if (!strcmp(commandPtr->string, "Quit"))
                            return (NO_ERROR_QUIT);
                        expecting = Expecting(COMMAND);
                        }
                    else
                        {
                        inError = YES;
                        WhatVariableExp (expecting, errStr);
                        YvyraPrint ("%s   Expecting %s\n", spacer, errStr);
                        }
                    }
                }
            /* Check to see if a comment is terminated. A comment can either be a right comment "]" or, if we were in a foreign nexus block
               (e.g., a "paup" block) the terminating comment will be "end". */
            if (tokenType == RIGHTCOMMENT)
                {
                if (inComment == NO && readComment == NO)
                    {
                    YvyraPrint ("%s   Found \"]\", without having previously found \"[\"\n", spacer);
                    inError = YES; 
                    }
                else if (inComment == NO && readComment == YES)
                    {
                    /* This is OK, we just pass through and rely on the command to handle the RIGHTCOMMENT */
                    }
                else
                    {
                    numComments--;
                    if (numComments == 0)
                        inComment = NO;
                    }
                }
            if ((IsSame(token, "end") == SAME || IsSame(token, "endblock") == SAME) && inForeignBlock == YES)
                {
                strcpy (spacer, "");
                inForeignBlock = NO;
                }
            }
        
        } while ((*token || tokenType == ALPHA) && inError == NO && skipCmd == NO);
        
    if (inError == YES)
        {
        readComment = NO;   /* reset this in case it is set to YES in command and we get an error exit */
        return (ERROR);
        }
    else
        return (NO_ERROR);
}


void PrintSettings (char *command)
{
    char yesNoStr[20];

    if (!strcmp(command,"Mcmc"))
        {
        YvyraPrint ("   Parameter       Options               Current Setting                         \n");
        YvyraPrint ("   -----------------------------------------------------                         \n");
        YvyraPrint ("   Ngen            <number>              %d                                      \n", chainParams.numGen);
        YvyraPrint ("   Nruns           <number>              %d                                      \n", chainParams.numRuns);
        YvyraPrint ("   Nchains         <number>              %d                                      \n", chainParams.numChains);
        YvyraPrint ("   Temp            <number>              %lf                                     \n", chainParams.chainTemp);
        YvyraPrint ("   Reweight        <number>,<number>     %1.2lf,%1.2lf                           \n", chainParams.weightScheme[0], chainParams.weightScheme[1]);
        YvyraPrint ("   Swapfreq        <number>              %d                                      \n", chainParams.swapFreq);
        YvyraPrint ("   Nswaps          <number>              %d                                      \n", chainParams.numSwaps);
        YvyraPrint ("   Samplefreq      <number>              %d                                      \n", chainParams.sampleFreq);
        YvyraPrint ("   Printfreq       <number>              %d                                      \n", chainParams.printFreq);
        PrintYesNo (chainParams.printAll, yesNoStr);
        YvyraPrint ("   Printall        Yes/No                %s                                      \n", yesNoStr);
        YvyraPrint ("   Printmax        <number>              %d                                      \n", chainParams.printMax);
        PrintYesNo (chainParams.mcmcDiagn, yesNoStr);
        YvyraPrint ("   Mcmcdiagn       Yes/No                %s                                      \n", yesNoStr);
        YvyraPrint ("   Diagnfreq       <number>              %d                                      \n", chainParams.diagnFreq);
        if (chainParams.diagnStat == AVGSTDDEV)
            strcpy (yesNoStr, "Avgstddev");
        else
            strcpy (yesNoStr, "Maxstddev");
        YvyraPrint ("   Diagnstat       Avgstddev/Maxstddev   %s                                     \n", yesNoStr);
        YvyraPrint ("   Minpartfreq     <number>              %1.2lf                                 \n", chainParams.minPartFreq);
        PrintYesNo (chainParams.allChains, yesNoStr);
        YvyraPrint ("   Allchains       Yes/No                %s                                     \n", yesNoStr);
        PrintYesNo (chainParams.allComps, yesNoStr);
        YvyraPrint ("   Allcomps        Yes/No                %s                                     \n", yesNoStr);
        PrintYesNo (chainParams.relativeBurnin, yesNoStr);
        YvyraPrint ("   Relburnin       Yes/No                %s                                     \n", yesNoStr);
        YvyraPrint ("   Burnin          <number>              %d                                     \n", chainParams.chainBurnIn);
        YvyraPrint ("   Burninfrac      <number>              %1.2lf                                 \n", chainParams.burninFraction);
        PrintYesNo (chainParams.stopRule, yesNoStr);
        YvyraPrint ("   Stoprule        Yes/No                %s                                     \n", yesNoStr);
        YvyraPrint ("   Stopval         <number>              %1.2lf                                 \n", chainParams.stopVal);
        PrintYesNo (chainParams.saveTrees, yesNoStr);
        YvyraPrint ("   Savetrees       Yes/No                %s                                     \n", yesNoStr);
        PrintYesNo (chainParams.checkPoint, yesNoStr);
        YvyraPrint ("   Checkpoint      Yes/No                %s                                     \n", yesNoStr);
        YvyraPrint ("   Checkfreq       <number>              %d                                     \n", chainParams.checkFreq);
        YvyraPrint ("   Filename        <name>                %s.<p/t>\n", chainParams.chainFileName);
        YvyraPrint ("   Startparams     Current/Reset         %s                                     \n", chainParams.startParams);
        YvyraPrint ("   Starttree       Current/Random/       %s                                     \n", chainParams.startTree);
        YvyraPrint ("                   Parsimony                                                    \n");
        YvyraPrint ("   Nperts          <number>              %d                                     \n", chainParams.numStartPerts);
        PrintYesNo (chainParams.runWithData, yesNoStr);
        YvyraPrint ("   Data            Yes/No                %s                                     \n", yesNoStr);
        YvyraPrint ("   Ordertaxa       Yes/No                %s                                     \n", chainParams.orderTaxa == YES? "Yes" : "No");
        YvyraPrint ("   Append          Yes/No                %s                                     \n", chainParams.append == YES? "Yes" : "No");
        YvyraPrint ("   Autotune        Yes/No                %s                                     \n", chainParams.autotune == YES? "Yes" : "No");
        YvyraPrint ("   Tunefreq        <number>              %d                                     \n", chainParams.tuneFreq);
        YvyraPrint ("                                                                                \n");
        }
}


void PrintYesNo (int yn, char s[4])
{
    if (yn == YES)
        strcpy (s, "Yes");
    else
        strcpy (s, "No");
}


int ProtID (char aa)
{
    if (aa == 'A' || aa == 'a')      /* Ala */
        {
        return 1;
        }
    else if (aa == 'R' || aa == 'r') /* Arg */
        {
        return 2;
        }
    else if (aa == 'N' || aa == 'n') /* Asn */
        {
        return 4;
        }
    else if (aa == 'D' || aa == 'd') /* Asp */
        {
        return 8;
        }
    else if (aa == 'C' || aa == 'c') /* Cys */
        {
        return 16;
        }
    else if (aa == 'Q' || aa == 'q') /* Gln */
        {
        return 32;
        }
    else if (aa == 'E' || aa == 'e') /* Glu */
        {
        return 64;
        }
    else if (aa == 'G' || aa == 'g') /* Gly */
        {
        return 128;
        }
    else if (aa == 'H' || aa == 'h') /* His */
        {
        return 256;
        }
    else if (aa == 'I' || aa == 'i') /* Ile */
        {
        return 512;
        }
    else if (aa == 'L' || aa == 'l') /* Leu */
        {
        return 1024;
        }
    else if (aa == 'K' || aa == 'k') /* Lys */
        {
        return 2048;
        }
    else if (aa == 'M' || aa == 'm') /* Met */
        {
        return 4096;
        }
    else if (aa == 'F' || aa == 'f') /* Phe */
        {
        return 8192;
        }
    else if (aa == 'P' || aa == 'p') /* Pro */
        {
        return 16384;
        }
    else if (aa == 'S' || aa == 's') /* Ser */
        {
        return 32768;
        }
    else if (aa == 'T' || aa == 't') /* Thr */
        {
        return 65536;
        }
    else if (aa == 'W' || aa == 'w') /* Trp */
        {
        return 131072;
        }
    else if (aa == 'Y' || aa == 'y') /* Tyr */
        {
        return 262144;
        }
    else if (aa == 'V' || aa == 'v') /* Val */
        {
        return 524288;
        }
    else if (aa == 'X' || aa == 'x') /* Nonidentified */
        {
        return MISSING;
        }
    else if (aa == gapId)
        {
        return GAP;
        }
    else if (aa == missingId)
        {
        return MISSING;
        }
    else
        return -1;
}


int RemoveLastFromString (char *s1)
{
    int     i, j, numPrev, numRemoved;
    
    /* We remove the last name from the string simply by deleting the last "|". */
       
    i = numPrev = 0;
    while (s1[i] != '\0')
        {
        if (s1[i] == '|')
            numPrev++;
        i++;
        }
        
    i = j = numRemoved = 0;
    while (s1[i] != '\0')
        {
        if (s1[i] == '|')
            j++;
        if (numPrev == j)
            {
            s1[i] = ' ';
            numRemoved++;
            break;
            }
        i++;
        }

    if (numRemoved != 1)
        {
        YvyraPrint ("%s   Could not find name to remove\n", spacer);
        return (ERROR);
        }

    return (NO_ERROR);
}


int MBResID (char nuc)
{
    char        n;
    
    n = nuc;

    if (n == '0' || n == 'a' || n == 'A')
        {
        return 1;
        }
    else if (n == '1' || n == 'b' || n == 'B')
        {
        return 2;
        }
    else if (n == gapId)
        {
        return GAP;
        }
    else if (n == missingId)
        {
        return MISSING;
        }
    else
        return -1;
}


/* Reset character flags */
void ResetCharacterFlags (void)
{
    /* reset all characters flags */
    numChar              = 0;                        /* number of defined characters                  */
    defChars             = NO;                       /* flag for whether number of characters is known*/
    defMatrix            = NO;                       /* flag for whether matrix is successful read    */
    matrixHasPoly        = NO;                       /* flag for whether matrix has polymorphisms     */
    isInAmbig            = NO;                       /* flag for whether the parser is within ()      */
    isInPoly             = NO;                       /* flag for whether the parser is within {}      */
    defPartition         = NO;                       /* flag for whether character partition is read  */
    defPairs             = NO;                       /* flag indicating whether pairs have been defined */
    numDefinedPartitions = 0;                        /* number of defined partitions                  */
    partitionNum         = 0;                        /* partition number currently enforced           */
    numCurrentDivisions  = 0;                        /* number of partitions of data                  */
    numCharSets          = 0;                        /* holds number of character sets                */
    numDivisions         = 1;                        /* holds number of partitions                    */
    isMixed              = NO;                       /* are data mixed ?                              */
    dataType             = NONE;                     /* holds datatype                                */
    matchId              = '\0';                     /* no default for match character                */
    gapId                = '\0';                     /* no default for gap character                  */
    missingId            = '\0';                     /* no default for missing characters             */
}


/* Reset taxa flags */
void ResetTaxaFlags (void)
{
    numTaxa                 = 0;                         /* number of taxa                                */
    numNamedTaxa            = 0;                         /* number of named taxa                          */
    defTaxa                 = NO;                        /* flag for whether number of taxa is known      */
    isTaxsetDef             = NO;                        /* is a taxlabels set defined                    */
    numDefinedConstraints   = 0;                         /* holds number of defined constraints           */
    definedConstraint       = NULL;
    definedConstraintTwo    = NULL;
    definedConstraintPruned       = NULL;
    definedConstraintTwoPruned    = NULL;
    constraintNames         = NULL;
    nodeCalibration         = NULL;
    tempActiveConstraints   = NULL;                      /* holds temp info on active constraints         */
    outGroupNum             = 0;                         /* default outgroup                              */
    numTaxaSets             = 0;                         /* holds number of taxa sets                     */
}


/* SetPartition: Set model partition */
int SetPartition (int part)
{
    int     i, j;
    
    /* Free space for modelParams and modelSettings */
    if (memAllocs[ALLOC_MODEL] == YES)
        {
        for (i=0; i<numCurrentDivisions; i++)
          free (modelParams[i].activeConstraints);
        free (modelParams);
        free (modelSettings);
        modelParams = NULL;
        modelSettings = NULL;
        memAllocs[ALLOC_MODEL] = NO;
        }

    /* Set model partition */
    partitionNum = part;
    numCurrentDivisions = 0;

    /* Set numCurrentDivisions to maximum division a character belongs to in partition part */
    for (i=0; i<numChar; i++)
        {
        j = partitionId[i][part];
        if (j > numCurrentDivisions)
            numCurrentDivisions = j;
        }

    /* Allocate space for partition models */
    modelParams = (Model *) SafeCalloc (numCurrentDivisions, sizeof (Model));
    modelSettings = (ModelInfo *) SafeCalloc (numCurrentDivisions, sizeof (ModelInfo));
    if (!modelParams || !modelSettings)
        {
        YvyraPrint ("%s   Could not allocate modelParams or modelSettings\n", spacer);
        if (modelParams)
            free (modelParams);
        if (modelSettings)
            free (modelSettings);
        return (ERROR);
        }
    memAllocs[ALLOC_MODEL] = YES;

    numVars = (int *) SafeRealloc ((void *) numVars, 3 * (size_t)numCurrentDivisions * sizeof(int));
    tempLinkUnlinkVec = numVars + numCurrentDivisions;
    activeParts       = numVars + 2*numCurrentDivisions;

    tempNum = (YFlt *) SafeRealloc ((void *) tempNum, 6 * sizeof(YFlt));

    activeParams[0] = (int *) SafeRealloc ((void *) (activeParams[0]), (size_t)NUM_LINKED * (size_t)numCurrentDivisions * sizeof(int));
    for (i=1; i<NUM_LINKED; i++)
        activeParams[i] = activeParams[0] + i*numCurrentDivisions;
 
    linkTable[0] = (int *) SafeRealloc ((void *) (linkTable[0]), 3 * (size_t)NUM_LINKED * (size_t)numCurrentDivisions * sizeof(int));
    tempLinkUnlink[0] = linkTable[0] + NUM_LINKED*numCurrentDivisions;
    for (i=1; i<NUM_LINKED; i++)
        {
        linkTable[i]      = linkTable[0] + i*numCurrentDivisions;
        tempLinkUnlink[i] = tempLinkUnlink[0] + i*numCurrentDivisions;
        }

    return (NO_ERROR);
}


/* SetSpeciespartition: Set speciespartition */
int SetSpeciespartition (int part)
{
    int     i, j;
    
    /* Set model partition */
    speciespartitionNum = part;
    numSpecies = 0;

    /* Set numSpecies to maximum species a taxon belongs to in partition part */
    for (i=0; i<numTaxa; i++)
        {
        j = speciespartitionId[i][part];
        if (j > numSpecies)
            numSpecies = j;
        }

    return (NO_ERROR);
}


int SetTaxaFromTranslateTable (void)
{
    int     i;

    if (numTaxa != 0)
        return ERROR;

    for (i=0; i<numTranslates; i++)
        {
        if (strlen(transFrom[i])>99)
            {
            YvyraPrint ("%s   Taxon name %s is too long. Maximum 99 characters is allowed.\n", spacer, transFrom[i]);
            return (ERROR);
            }
        AddString(&taxaNames, numTaxa, transFrom[i]);
        numTaxa++;
        }
    
    return NO_ERROR;
}


void SetUpParms (void)
{
    ParmInfoPtr p = paramTable;

    PARAM   (0, "NEXUS",          DoNexusParm,       "NEXUS|\0");
    PARAM   (1, "Data",           DoBeginParm,       "\0");
    PARAM   (2, "Mrbayes",        DoBeginParm,       "\0");
    PARAM   (3, "Trees",          DoBeginParm,       "\0");
    PARAM   (4, "Ntax",           DoDimensionsParm,  "\0");
    PARAM   (5, "Nchar",          DoDimensionsParm,  "\0");
    PARAM   (6, "Interleave",     DoFormatParm,      "Yes|No|\0");
    PARAM   (7, "Datatype",       DoFormatParm,      "Dna|Rna|Protein|Restriction|Standard|Continuous|Mixed|\0");
    PARAM   (8, "Gap",            DoFormatParm,      "\0");
    PARAM   (9, "Missing",        DoFormatParm,      "\0");
    PARAM  (10, "Matchchar",      DoFormatParm,      "\0");
    PARAM  (11, "MatrixInfo",     DoMatrixParm,      "\0");
    PARAM  (12, "Filename",       DoExecuteParm,     "\0");
    PARAM  (13, "Autoclose",      DoSetParm,         "Yes|No|\0");
    PARAM  (14, "Partition",      DoSetParm,         "\0");
    PARAM  (15, "Xxxxxxxxxx",     DoCharsetParm,     "\0");
    PARAM  (16, "Xxxxxxxxxx",     DoPartitionParm,   "\0");
    PARAM  (17, "Seed",           DoMcmcParm,        "\0");
    PARAM  (18, "Ngen",           DoMcmcParm,        "\0");
    PARAM  (19, "Samplefreq",     DoMcmcParm,        "\0");
    PARAM  (20, "Printfreq",      DoMcmcParm,        "\0");
    PARAM  (21, "Nchains",        DoMcmcParm,        "\0");
    PARAM  (22, "Temp",           DoMcmcParm,        "\0");
    PARAM  (23, "Filename",       DoMcmcParm,        "\0");
    PARAM  (24, "Burnin",         DoMcmcParm,        "\0");
    PARAM  (25, "Starttree",      DoMcmcParm,        "Random|Current|User|Parsimony|NJ|\0");
    PARAM  (26, "Nperts",         DoMcmcParm,        "\0");
    PARAM  (27, "Savebrlens",     DoMcmcParm,        "Yes|No|\0");
    PARAM  (28, "Nucmodel",       DoLsetParm,        "4by4|Doublet|Codon|Protein|\0");
    PARAM  (29, "Nst",            DoLsetParm,        "1|2|6|Mixed|\0");
    PARAM  (30, "Aamodel",        DoLsetParm,        "Poisson|Equalin|Jones|Dayhoff|Mtrev|Mtmam|Wag|Rtrev|Cprev|Vt|Blosum|Blossum|LG|\0");
    PARAM  (31, "Parsmodel",      DoLsetParm,        "Yes|No|\0");
    PARAM  (32, "Omegavar",       DoLsetParm,        "Equal|Ny98|M3|M10|\0");
    PARAM  (33, "Code",           DoLsetParm,        "Universal|Vertmt|Invermt|Mycoplasma|Yeast|Ciliate|Echinoderm|Euplotid|Metmt|\0");
    PARAM  (34, "Coding",         DoLsetParm,        "All|Variable|Informative|Nosingletons|Noabsencesites|Nopresencesites|Nosingletonpresence|Nosingletonabsence|\0");
    PARAM  (35, "Seqerror",       DoPrsetParm,       "\0");
    PARAM  (36, "Tratiopr",       DoPrsetParm,       "Beta|Fixed|\0");
    PARAM  (37, "Revmatpr",       DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM  (38, "Omegapr",        DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM  (39, "Statefreqpr",    DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM  (40, "Ngammacat",      DoLsetParm,        "\0");
    PARAM  (41, "Shapepr",        DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (42, "Ratecorrpr",     DoPrsetParm,       "Uniform|Fixed|\0");
    PARAM  (43, "Pinvarpr",       DoPrsetParm,       "Uniform|Fixed|\0");
    PARAM  (44, "Covswitchpr",    DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (45, "Xxxxxxxxxx",     DoExcludeParm,     "\0");
    PARAM  (46, "Xxxxxxxxxx",     DoIncludeParm,     "\0");
    PARAM  (47, "Xxxxxxxxxx",     DoDeleteParm,      "\0");
    PARAM  (48, "Xxxxxxxxxx",     DoRestoreParm,     "\0");
    PARAM  (49, "Xxxxxxxxxx",     DoTaxasetParm,     "\0");
    PARAM  (50, "Xxxxxxxxxx",     DoHelpParm,        "\0");
    PARAM  (51, "Applyto",        DoLsetParm,        "\0");
    PARAM  (52, "Rates",          DoLsetParm,        "Equal|Gamma|LNorm|Propinv|Invgamma|Adgamma|Kmixture|\0");
    PARAM  (53, "Covarion",       DoLsetParm,        "Yes|No|\0");
    PARAM  (54, "Applyto",        DoPrsetParm,       "\0");
    PARAM  (55, "Tratio",         DoLinkParm,        "\0");
    PARAM  (56, "Revmat",         DoLinkParm,        "\0");
    PARAM  (57, "Omega",          DoLinkParm,        "\0");
    PARAM  (58, "Statefreq",      DoLinkParm,        "\0");
    PARAM  (59, "Shape",          DoLinkParm,        "\0");
    PARAM  (60, "Pinvar",         DoLinkParm,        "\0");
    PARAM  (61, "Correlation",    DoLinkParm,        "\0");
    PARAM  (62, "Ratemultiplier", DoLinkParm,        "\0");
    PARAM  (63, "Switchrates",    DoLinkParm,        "\0");
    PARAM  (64, "Symdirihyperpr", DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (65, "Xxxxxxxxxx",     DoCtypeParm,       "\0");
    PARAM  (66, "Xxxxxxxxxx",     DoConstraintParm,  "\0");
    PARAM  (67, "Topologypr",     DoPrsetParm,       "Uniform|Constraints|Fixed|Speciestree|\0");
    PARAM  (68, "Brlenspr",       DoPrsetParm,       "Unconstrained|Clock|Relaxedclock|Fixed|\0");
    PARAM  (69, "Speciationpr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (70, "Extinctionpr",   DoPrsetParm,       "Beta|Exponential|Fixed|\0");
    PARAM  (71, "Popsizepr",      DoPrsetParm,       "Lognormal|Uniform|Gamma|Normal|Fixed|\0");
    PARAM  (72, "Topology",       DoLinkParm,        "\0");
    PARAM  (73, "Brlens",         DoLinkParm,        "\0");
    PARAM  (74, "Speciationrate", DoLinkParm,        "\0");
    PARAM  (75, "Extinctionrate", DoLinkParm,        "\0");
    PARAM  (76, "Popsize",        DoLinkParm,        "\0");
    PARAM  (77, "Ratepr",         DoPrsetParm,       "Variable|Dirichlet|Fixed|\0");
    PARAM  (78, "Xxxxxxxxxx",     DoOutgroupParm,    "\0");
    PARAM  (79, "Xxxxxxxxxx",     DoTreeParm,        "\0");
    PARAM  (80, "Filename",       DoSumtParm,        "\0");
    PARAM  (81, "Burnin",         DoSumtParm,        "\0");
    PARAM  (82, "Contype",        DoSumtParm,        "Halfcompat|Allcompat|\0");
    PARAM  (83, "Xxxxxxxxxx",     DoTranslateParm,   "\0");
    PARAM  (84, "Swapfreq",       DoMcmcParm,        "\0");
    PARAM  (85, "Start",          DoLogParm,         "\0");
    PARAM  (86, "Stop",           DoLogParm,         "\0");
    PARAM  (87, "Filename",       DoLogParm,         "\0");
    PARAM  (88, "Append",         DoLogParm,         "\0");
    PARAM  (89, "Replace",        DoLogParm,         "\0");
    PARAM  (90, "Nbetacat",       DoLsetParm,        "\0");
    PARAM  (91, "Augment",        DoLsetParm,        "Yes|No|\0");
    PARAM  (92, "Xxxxxxxxxx",     DoPairsParm,       "\0");
    PARAM  (93, "Xxxxxxxxxx",     DoBreaksParm,      "\0");
    PARAM  (94, "Nowarnings",     DoSetParm,         "Yes|No|\0");
    PARAM  (95, "Showtreeprobs",  DoSumtParm,        "Yes|No|\0");
    PARAM  (96, "Filename",       DoSumpParm,        "\0");
    PARAM  (97, "Burnin",         DoSumpParm,        "\0");
    PARAM  (98, "Reweight",       DoMcmcParm,        "\0");
    PARAM  (99, "Noop",           DoMcmcParm,        "\0");
    PARAM (100, "Ny98omega1pr",   DoPrsetParm,       "Beta|Fixed|\0");
    PARAM (101, "Ny98omega3pr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM (102, "Codoncatfreqs",  DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM (103, "Sampleprob",     DoPrsetParm,       "\0");
    PARAM (104, "Aamodelpr",      DoPrsetParm,       "Fixed|Mixed|\0");
    PARAM (105, "Aamodel",        DoLinkParm,        "\0");
    PARAM (106, "Filename",       DoPlotParm,        "\0");
    PARAM (107, "Parameter",      DoPlotParm,        "\0");
    PARAM (108, "Match",          DoPlotParm,        "Perfect|Consistentwith|All|\0");
    PARAM (109, "Burnin",         DoPlotParm,        "\0");
    PARAM (110, "Brownscalepr",   DoPrsetParm,       "Uniform|Gamma|Fixed|\0");
    PARAM (111, "Browncorrpr",    DoPrsetParm,       "Uniform|Fixed|\0");
    PARAM (112, "Pbf",            DoMcmcParm,        "Yes|No|\0");
    PARAM (113, "Pbfinitburnin",  DoMcmcParm,        "\0");
    PARAM (114, "Pbfsamplefreq",  DoMcmcParm,        "\0");
    PARAM (115, "Pbfsampletime",  DoMcmcParm,        "\0");
    PARAM (116, "Pbfsampleburnin",DoMcmcParm,        "\0");
    PARAM (117, "Growthpr",       DoPrsetParm,       "Uniform|Exponential|Fixed|Normal|\0");
    PARAM (118, "Growthrate",     DoLinkParm,        "\0");
    PARAM (119, "Xxxxxxxxxx",     DoCalibrateParm,   "Unconstrained|Fixed|Uniform|Offsetexponential|Truncatednormal|Lognormal|Offsetlognormal|Gamma|Offsetgamma|\0");
    PARAM (120, "Calwaitpr",      DoPrsetParm,       "Exponential|Fixed|\0");     /* not used but leave it in to not destroy mapping to commands */
    PARAM (121, "M3omegapr",      DoPrsetParm,       "Exponential|Fixed|\0");
    PARAM (122, "Applyto",        DoReportParm,      "\0");
    PARAM (123, "Tratio",         DoReportParm,      "Dirichlet|Ratio|\0");
    PARAM (124, "Revmat",         DoReportParm,      "Dirichlet|Ratio|\0");
    PARAM (125, "Ratemult",       DoReportParm,      "Dirichlet|Scaled|Ratio|\0");
    PARAM (126, "Filename",       DoManualParm,      "\0");
    PARAM (127, "Filename1",      DoCompareTreeParm, "\0");
    PARAM (128, "Filename2",      DoCompareTreeParm, "\0");
    PARAM (129, "Outputname",     DoCompareTreeParm, "\0");
    PARAM (130, "Burnin",         DoCompareTreeParm, "\0");
    PARAM (131, "Ploidy",         DoLsetParm,        "Haploid|Diploid|Zlinked|\0");
    PARAM (132, "Swapadjacent",   DoMcmcParm,        "Yes|No|\0");
    PARAM (133, "Treeagepr",      DoPrsetParm,       "Fixed|Uniform|Offsetexponential|Truncatednormal|Lognormal|Offsetlognormal|Gamma|Offsetgamma|\0");
    PARAM (134, "Ancstates",      DoReportParm,      "Yes|No|\0");
    PARAM (135, "Siterates",      DoReportParm,      "Yes|No|\0");
    PARAM (136, "Possel",         DoReportParm,      "Yes|No|\0");
    PARAM (137, "Plot",           DoSumpParm,        "Yes|No|\0");
    PARAM (138, "Table",          DoSumpParm,        "Yes|No|\0");
    PARAM (139, "Minprob",        DoSumpParm,        "\0");
    PARAM (140, "Printtofile",    DoSumpParm,        "Yes|No|\0");
    PARAM (141, "Outputname",     DoSumpParm,        "\0");
    PARAM (142, "Redirect",       DoMcmcParm,        "Yes|No|\0");
    PARAM (143, "Swapseed",       DoMcmcParm,        "\0");
    PARAM (144, "Runidseed",      DoMcmcParm,        "\0");
    PARAM (145, "Quitonerror",    DoSetParm,         "Yes|No|\0");
    PARAM (146, "Savebrparams",   DoSumtParm,        "Yes|No|\0");
    PARAM (147, "Minbrparamfreq", DoSumtParm,        "\0");
    PARAM (148, "Minpartfreq",    DoMcmcParm,        "\0");
    PARAM (149, "Allchains",      DoMcmcParm,        "Yes|No|\0");
    PARAM (150, "Mcmcdiagn",      DoMcmcParm,        "Yes|No|\0");
    PARAM (151, "Diagnfreq",      DoMcmcParm,        "\0");
    PARAM (152, "Nruns",          DoMcmcParm,        "\0");
    PARAM (153, "Stoprule",       DoMcmcParm,        "Yes|No|\0");
    PARAM (154, "Stopval",        DoMcmcParm,        "\0");
    PARAM (155, "Relburnin",      DoMcmcParm,        "Yes|No|\0");
    PARAM (156, "Burninfrac",     DoMcmcParm,        "\0");
    PARAM (157, "Allcomps",       DoMcmcParm,        "Yes|No|\0");
    PARAM (158, "Printall",       DoMcmcParm,        "Yes|No|\0");
    PARAM (159, "Printmax",       DoMcmcParm,        "\0");
    PARAM (160, "Data",           DoMcmcParm,        "Yes|No|\0");
    PARAM (161, "Nruns",          DoSumpParm,        "\0");
    PARAM (162, "Allruns",        DoSumpParm,        "Yes|No|\0");
    PARAM (163, "Nruns",          DoSumtParm,        "\0");
    PARAM (164, "Ntrees",         DoSumtParm,        "\0");
    PARAM (165, "Calctreeprobs",  DoSumtParm,        "Yes|No|\0");
    PARAM (166, "Ordertaxa",      DoMcmcParm,        "Yes|No|\0");
    PARAM (167, "Ordertaxa",      DoSumtParm,        "Yes|No|\0");
    PARAM (168, "Aarevmatpr",     DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM (169, "Nswaps",         DoMcmcParm,        "\0");
    PARAM (170, "Autoreplace",    DoSetParm,         "Yes|No|\0");
    PARAM (171, "Npthreads",      DoSetParm,         "\0");
    PARAM (172, "Cppratepr",      DoPrsetParm,       "Fixed|Exponential|\0");
    PARAM (173, "Cppmultdevpr",   DoPrsetParm,       "Fixed|\0");
    PARAM (174, "TK02varpr",      DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (175, "Tfile",          DoSumtParm,        "\0");
    PARAM (176, "Pfile",          DoSumpParm,        "\0");
    PARAM (177, "Autocomplete",   DoSumtParm,        "Yes|No|\0");
    PARAM (178, "Autocomplete",   DoSumpParm,        "Yes|No|\0");
    PARAM (179, "Userlevel",      DoSetParm,         "Standard|Developer|\0");
    PARAM (180, "Allavailable",   DoShowmovesParm,   "Yes|No|\0");
    PARAM (181, "Seed",           DoSetParm,         "\0");
    PARAM (182, "Swapseed",       DoSetParm,         "\0");
    PARAM (183, "Clockratepr",    DoPrsetParm,       "Fixed|Normal|Lognormal|Exponential|Gamma|\0");
    PARAM (184, "Nodeagepr",      DoPrsetParm,       "Unconstrained|Calibrated|\0");
    PARAM (185, "Clockvarpr",     DoPrsetParm,       "Strict|Cpp|TK02|WN|IGR|ILN|Mixed|\0");
    PARAM (186, "Xxxxxxxxxx",     DoPropsetParm,     "\0");
    PARAM (187, "Xxxxxxxxxx",     DoStartvalsParm,   "\0");
    PARAM (188, "Usegibbs",       DoLsetParm,        "Yes|No|\0");
    PARAM (189, "Gibbsfreq",      DoLsetParm,        "\0");
    PARAM (190, "Checkpoint",     DoMcmcParm,        "Yes|No|\0");
    PARAM (191, "Checkfreq",      DoMcmcParm,        "\0");
    PARAM (192, "Tree",           DoReportParm,      "Topology|Brlens|\0");
    PARAM (193, "Cpprate",        DoLinkParm,        "\0");
    PARAM (194, "Cppmultdev",     DoLinkParm,        "\0");
    PARAM (195, "Cppevents",      DoLinkParm,        "\0");
    PARAM (196, "TK02var",        DoLinkParm,        "\0");
    PARAM (197, "TK02branchrates",DoLinkParm,        "\0");
    PARAM (198, "Savetrees",      DoMcmcParm,        "Yes|No|\0");
    PARAM (199, "Diagnstat",      DoMcmcParm,        "Avgstddev|Maxstddev|\0");
    PARAM (200, "Startparams",    DoMcmcParm,        "Reset|Current|\0");
    PARAM (201, "Characters",     DoBeginParm,       "\0");
    PARAM (202, "Startingtrees",  DoMcmcParm,        "\0");
    PARAM (203, "Xxxxxxxxxx",     DoUserTreeParm,    "\0");
    PARAM (204, "Outputname",     DoSumtParm,        "\0");
    PARAM (205, "Table",          DoSumtParm,        "Yes|No|\0");
    PARAM (206, "Summary",        DoSumtParm,        "Yes|No|\0");
    PARAM (207, "Consensus",      DoSumtParm,        "Yes|No|\0");
    PARAM (208, "Minpartfreq",    DoSumtParm,        "\0");
    PARAM (209, "Relburnin",      DoSumtParm,        "Yes|No|\0");
    PARAM (210, "Burninfrac",     DoSumtParm,        "\0");
    PARAM (211, "Relburnin",      DoSumpParm,        "Yes|No|\0");
    PARAM (212, "Burninfrac",     DoSumpParm,        "\0");
    PARAM (213, "Append",         DoMcmcParm,        "Yes|No|\0");
    PARAM (214, "Autotune",       DoMcmcParm,        "Yes|No|\0");
    PARAM (215, "Tunefreq",       DoMcmcParm,        "\0");
    PARAM (216, "Scientific",     DoSetParm,         "Yes|No|\0");
    PARAM (217, "Siteomega",      DoReportParm,      "Yes|No|\0");
    PARAM (218, "IGRvarpr",       DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (219, "Symbols",        DoFormatParm,      "\0");
    PARAM (220, "Equate",         DoFormatParm,      "\0");
    PARAM (221, "Relburnin",      DoCompareTreeParm, "Yes|No|\0");
    PARAM (222, "Burninfrac",     DoCompareTreeParm, "\0");
    PARAM (223, "Minpartfreq",    DoCompareTreeParm, "\0");
    PARAM (224, "Relburnin",      DoPlotParm,        "Yes|No|\0");
    PARAM (225, "Burninfrac",     DoPlotParm,        "\0");
    PARAM (226, "Taxa",           DoBeginParm,       "\0");
    PARAM (227, "Xxxxxxxxxx",     DoBeginParm,       "\0");
    PARAM (228, "Xxxxxxxxxx",     DoTaxlabelsParm,   "\0");
    PARAM (229, "Dir",            DoSetParm,         "\0");
    PARAM (230, "Conformat",      DoSumtParm,        "Figtree|Simple|\0");
    PARAM (231, "Hpd",            DoSumpParm,        "Yes|No|\0");
    PARAM (232, "Hpd",            DoSumtParm,        "Yes|No|\0");
    PARAM (233, "Usebeagle",      DoSetParm,         "Yes|No|\0");
    PARAM (234, "Beagledevice",   DoSetParm,         "Cpu|Gpu|\0");
    PARAM (235, "Beagleprecision",DoSetParm,         "Single|Double|\0");
    PARAM (236, "Beaglesse",      DoSetParm,         "Yes|No|\0");
    PARAM (237, "Beagleopenmp",   DoSetParm,         "Yes|No|\0"); /* not in use */
    PARAM (238, "Beaglethreads",  DoSetParm,         "Yes|No|\0");
    PARAM (239, "Beaglescaling",  DoSetParm,         "Always|Dynamic|\0");
    PARAM (240, "Beaglefreq",     DoSetParm,         "\0");
    PARAM (241, "Popvarpr",       DoPrsetParm,       "Equal|Variable|\0");
    PARAM (242, "IGRvar",         DoLinkParm,        "\0");
    PARAM (243, "IGRbranchrates", DoLinkParm,        "\0");
    PARAM (244, "Xxxxxxxxxx",     DoSpeciespartitionParm,   "\0");
    PARAM (245, "Speciespartition",  DoSetParm,      "\0");
    PARAM (246, "Revratepr",      DoPrsetParm,       "Symdir|\0");
    PARAM (247, "Samplestrat",    DoPrsetParm,       "Random|Diversity|Cluster|FossilTip|\0");
    PARAM (248, "Burninss",       DoSsParm,          "\0");
    PARAM (249, "Nsteps",         DoSsParm,          "\0");
    PARAM (250, "Alpha",          DoSsParm,          "\0");
    PARAM (251, "WNvarpr",        DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (252, "WNvar",          DoLinkParm,        "\0");
    PARAM (253, "WNbranchrates",  DoLinkParm,        "\0");
    PARAM (254, "ILNvarpr",       DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (255, "ILNvar",         DoLinkParm,        "\0");
    PARAM (256, "ILNbranchlens",  DoLinkParm,        "\0");
    PARAM (257, "FromPrior",      DoSsParm,          "Yes|No|\0");
    PARAM (258, "Filename",       DoSumSsParm,       "\0");
    PARAM (259, "Burnin",         DoSumSsParm,       "\0");
    PARAM (260, "Nruns",          DoSumSsParm,       "\0");
    PARAM (261, "Allruns",        DoSumSsParm,       "Yes|No|\0");
    PARAM (262, "Askmore",        DoSumSsParm,       "Yes|No|\0");
    PARAM (263, "Relburnin",      DoSumSsParm,       "Yes|No|\0");
    PARAM (264, "Burninfrac",     DoSumSsParm,       "\0");
    PARAM (265, "Discardfrac",    DoSumSsParm,       "\0");
    PARAM (266, "Smoothing",      DoSumSsParm,       "\0");
    PARAM (267, "Steptoplot",     DoSumSsParm,       "\0");
    PARAM (268, "Precision",      DoSetParm,         "\0");
    PARAM (269, "Fossilizationpr",   DoPrsetParm,    "Beta|Exponential|Fixed|\0");
    PARAM (270, "Browncorr",      DoLinkParm,        "\0");
    PARAM (271, "Generatepr",     DoPrsetParm,       "Variable|Fixed|\0");
    PARAM (272, "Mixedvarpr",     DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (273, "Mixedvar",       DoLinkParm,        "\0");
    PARAM (274, "Mixedbrchrates", DoLinkParm,        "\0");
    PARAM (275, "Beagleresource", DoSetParm,         "\0");
    PARAM (276, "Nlnormcat",      DoLsetParm,        "\0");
    PARAM (277, "Nmixtcat",       DoLsetParm,        "\0");
    PARAM (278, "Beaglethreadcount",  DoSetParm,     "\0");
    PARAM (279, "Beaglefloattips",DoSetParm,         "Yes|No|\0");
    PARAM (280, "Statefreqmodel", DoLsetParm,        "Stationary|Directional|Mixed|\0"); //SK
    PARAM (281, "Rootfreqpr",     DoPrsetParm,       "Dirichlet|Fixed|\0"); //SK
    PARAM (282, "Statefrmod",     DoLsetParm,        "Stationary|Directional|Mixed|\0"); //SK
    PARAM (283, "Sitelikes",     DoReportParm,      "Yes|No|\0");
    PARAM (284, "Xxxxxxxxxx",   DoUsertypeParm,    "\0");
    PARAM (285, "Xxxxxxxxxx",   DoUsertypeParm,    "\0");
    PARAM (286, "Xxxxxxxxxx",   DoWtsetParm,       "\0");
    PARAM (287, "Xxxxxxxxxx",   DoTipweightsParm,  "\0");
    PARAM (288, "Xxxxxxxxxx",   DoChardiagParm,    "\0");
    PARAM (289, "Xxxxxxxxxx",   DoSensitivityParm, "\0");
    PARAM (290, "Xxxxxxxxxx",   DoCharlabelsParm,  "\0");

    /* NOTE: If a change is made to the parameter table, make certain you change
            NUMPARAMS (now 291; one more than last index) at the top of this file. */
    /* CmdType commands[] */
}


void ShowNodes (TreeNode *p, int indent, int isThisTreeRooted)
{
    if (p != NULL)
        {
        printf ("   ");
        if (p->left == NULL && p->right == NULL && p->anc != NULL)
            {
            printf ("%*cN %d (l=%d r=%d a=%d) %1.15lf (%s) isDated=%d ",
            indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->length, p->label, p->isDated);
            }
        else if (p->left != NULL && p->right == NULL && p->anc == NULL)
            {
            if (isThisTreeRooted == NO)
                {
                if (p->label[0] == '\0' || p->label[0] == '\n' || p->label[0] == ' ')
                    printf ("%*cN %d (l=%d r=%d a=%d) (---) ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
                else
                    printf ("%*cN %d (l=%d r=%d a=%d) (%s) ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->label);
                }
            else
                {
                printf ("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ",
                indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
                }
            }
        else
            {
            if (p->anc != NULL)
                {
                if (p->anc->anc == NULL && isThisTreeRooted == YES)
                    printf ("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
                else    
                    printf ("%*cN %d (l=%d r=%d a=%d) %1.15lf ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->length);
                }
            }
        if (isThisTreeRooted == YES)
            printf ("depth=%1.15lf\n", p->nodeDepth);
        else
            printf ("\n");
        ShowNodes (p->left,  indent + 2, isThisTreeRooted);
        ShowNodes (p->right, indent + 2, isThisTreeRooted);
        }
}


int StandID (char nuc)
{
    /* Note that if you change how many states are recognized, you need
       to look at IsMissing */
    char n = nuc;

    if (n == '0')
        {
        return 1;
        }
    else if (n == '1')
        {
        return 2;
        }
    else if (n == '2')
        {
        return 4;
        }
    else if (n == '3')
        {
        return 8;
        }
    else if (n == '4')
        {
        return 16;
        }
    else if (n == '5')
        {
        return 32;
        }
    else if (n == '6')
        {
        return 64;
        }
    else if (n == '7')
        {
        return 128;
        }
    else if (n == '8')
        {
        return 256;
        }
    else if (n == '9')
        {
        return 512;
        }
    else if (n == 'A' || n == 'a')
        {
        return 1024;
        }
    else if (n == 'B' || n == 'b')
        {
        return 2048;
        }
    else if (n == 'C' || n == 'c')
        {
        return 4096;
        }
    else if (n == 'D' || n == 'd')
        {
        return 8192;
        }
    else if (n == 'E' || n == 'e')
        {
        return 16384;
        }
    else if (n == 'F' || n == 'f')
        {
        return 32768;
        }
    else if (n == 'G' || n == 'g')
        {
        return 65536;
        }
    else if (n == 'H' || n == 'h')
        {
        return 131072;
        }
    else if (n == 'I' || n == 'i')
        {
        return 262144;
        }
    else if (n == 'J' || n == 'j')
        {
        return 524288;
        }
    else if (n == 'K' || n == 'k')
        {
        return 1048576;
        }
    else if (n == 'L' || n == 'l')
        {
        return 2097152;
        }
    else if (n == 'M' || n == 'm')
        {
        return 4194304;
        }
    else if (n == 'N' || n == 'n')
        {
        return 8388608;
        }
    else if (n == missingId)
        {
        return MISSING;
        }
    else if (n == gapId)
        {
        return GAP;
        }
    else
        return -1;
}
int StateCode_Std (int n)
{
    /* max 24 states: 0-9 A-N */
    if (n <= 9 && n >= 0)
        return '0' + n;
    else if (n == 10)
        return 'A';
    else if (n == 11)
        return 'B';
    else if (n == 12)
        return 'C';
    else if (n == 13)
        return 'D';
    else if (n == 14)
        return 'E';
    else if (n == 15)
        return 'F';
    else if (n == 16)
        return 'G';
    else if (n == 17)
        return 'H';
    else if (n == 18)
        return 'I';
    else if (n == 19)
        return 'J';
    else if (n == 20)
        return 'K';
    else if (n == 21)
        return 'L';
    else if (n == 22)
        return 'M';
    else if (n == 23)
        return 'N';
    else return '?';
}


void WhatVariableExp (BitsLong exp, char *st)
{
    int         n;
    
    strcpy (st, "");
    n = 0;
    if (exp == 0)
        strcat(st, " nothing");
    else
        {
        if ((exp & Expecting(COMMAND)) == Expecting(COMMAND))
            {
            strcat(st, " command");
            n++;
            }
        if ((exp & Expecting(PARAMETER)) == Expecting(PARAMETER))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " parameter");
            n++;
            }
        if ((exp & Expecting(EQUALSIGN)) == Expecting(EQUALSIGN))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " =");
            n++;
            }
        if ((exp & Expecting(COLON)) == Expecting(COLON))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " :");
            n++;
            }
        if ((exp & Expecting(SEMICOLON)) == Expecting(SEMICOLON))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ;");
            n++;
            }
        if ((exp & Expecting(COMMA)) == Expecting(COMMA))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ,");
            n++;
            }
        if ((exp & Expecting(POUNDSIGN)) == Expecting(POUNDSIGN))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " #");
            n++;
            }
        if ((exp & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ?");
            n++;
            }
        if ((exp & Expecting(DASH)) == Expecting(DASH))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " -");
            n++;
            }
        if ((exp & Expecting(LEFTPAR)) == Expecting(LEFTPAR))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " (");
            n++;
            }
        if ((exp & Expecting(RIGHTPAR)) == Expecting(RIGHTPAR))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " )");
            n++;
            }
        if ((exp & Expecting(LEFTCOMMENT)) == Expecting(LEFTCOMMENT))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " [");
            n++;
            }
        if ((exp & Expecting(RIGHTCOMMENT)) == Expecting(RIGHTCOMMENT))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ]");
            n++;
            }
        if ((exp & Expecting(ALPHA)) == Expecting(ALPHA))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " <name>");
            n++;
            }
        if ((exp & Expecting(NUMBER)) == Expecting(NUMBER))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " <number>");
            n++;
            }
        if ((exp & Expecting(RETURNSYMBOL)) == Expecting(RETURNSYMBOL))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " return");
            n++;
            }
        if ((exp & Expecting(ASTERISK)) == Expecting(ASTERISK))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " *");
            n++;
            }
        if ((exp & Expecting(BACKSLASH)) == Expecting(BACKSLASH))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " /");
            n++;
            }
        if ((exp & Expecting(BACKSLASH)) == Expecting(BACKSLASH))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " \\");
            n++;
            }
        if ((exp & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " !");
            n++;
            }
        if ((exp & Expecting(PERCENT)) == Expecting(PERCENT))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " %");
            n++;
            }
        if ((exp & Expecting(LEFTCURL)) == Expecting(LEFTCURL))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " {");
            n++;
            }
        if ((exp & Expecting(RIGHTCURL)) == Expecting(RIGHTCURL))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " }");
            n++;
            }
        if ((exp & Expecting(WEIRD)) == Expecting(WEIRD))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " <whatever>");
            n++;
            }
        if ((exp & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " |");
            n++;
            }
        if ((exp & Expecting(UNKNOWN_TOKEN_TYPE)) == Expecting(UNKNOWN_TOKEN_TYPE))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " no clue");
            n++;
            }
        }
}
char WhichStand (int x)
{
    if (x == 1)
        return ('0');
    else if (x == 2)
        return ('1');
    else if (x == 4)
        return ('2');
    else if (x == 8)
        return ('3');
    else if (x == 16)
        return ('4');
    else if (x == 32)
        return ('5');
    else if (x == 64)
        return ('6');
    else if (x == 128)
        return ('7');
    else if (x == 256)
        return ('8');
    else if (x == 512)
        return ('9');
    else if (x == 1024)
        return ('a');
    else if (x == 2048)
        return ('b');
    else if (x == 4096)
        return ('c');
    else if (x == 8192)
        return ('d');
    else if (x == 16384)
        return ('e');
    else if (x == 32768)
        return ('f');
    else if (x == 65536)
        return ('g');
    else if (x == 131072)
        return ('h');
    else if (x == 262144)
        return ('i');
    else if (x == 524288)
        return ('j');
    else if (x == 1048576)
        return ('k');
    else if (x == 2097152)
        return ('l');
    else if (x == 4194304)
        return ('m');
    else if (x == 8388608)
        return ('n');
    else if (x > 0 && x < 8388608)
        return ('*');
    else if (x == MISSING)
        return ('?');
    else if (x == GAP)
        return ('-');
    else 
        return (' ');
}


YFlt WhichCont (int x)
{
    return ((YFlt)(x / 1000.0));
}

#ifndef COMMAND_H_
#define COMMAND_H_

int         AddString (char ***list, int len, char *token);
BitsLong    Expecting (int y);
int         CheckString (char **list, int len, char *token, int *matchIndex);
int         CheckStringValidity (char *s);
int         DoExecute (void);
int         FreeMatrix (void);
int         GetToken (char *token, int *tokenType, char **sourceH);
int         FindValidCommand (char *tk, int *numMatches);
int         IsArgValid (char *s, char *validArg);
int         IsIn (char ch, char *s);
int         IsSame (char *s1, char *s2);
int         IsWhite (char c);
int         ParseCommand (char *s);
void        ResetCharacterFlags (void);
void        ResetTaxaFlags (void);
int         RootUserTree (TreeNode *p);
void        SetUpParms (void);
void        ShowNodes (TreeNode *p, int indent, int isThisTreeRooted);
int         ShowTree (Tree *t);
int         StateCode_Std (int n);
char        WhichStand (int x);

#endif  /* COMMAND_H_ */

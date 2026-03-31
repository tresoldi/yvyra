#ifndef CHARWEIGHT_H_
#define CHARWEIGHT_H_

#include "bayes.h"

/*
 *  charweight.h: Estimated character weight sampling.
 *
 *  Characters with weight: estimated get their weights sampled
 *  during MCMC via log-normal multiplier proposals with Gamma priors.
 */

#define MAX_ESTIMATED_WEIGHTS  200

/* Initialize weight estimation system. Called before MCMC. */
int  InitCharWeights (void);

/* Free weight estimation resources. Called after MCMC. */
void FreeCharWeights (void);

/* Propose a weight change for one random estimated-weight character.
   Returns the log acceptance ratio (lnPriorRatio + lnProposalRatio).
   Updates the weight in place; caller must accept/reject. */
YFlt ProposeCharWeight (int chain, RandLong *seed);

/* Accept or reject the last weight proposal. */
void AcceptCharWeight (void);
void RejectCharWeight (void);

/* Get the current weight multiplier for compressed character c.
   Returns 1.0 for fixed-weight characters. */
YFlt GetCharWeightMultiplier (int compCharIdx);

/* Print weight parameter header for .p file */
void PrintCharWeightHeader (FILE *fp);

/* Print weight parameter values for .p file */
void PrintCharWeightValues (FILE *fp);

/* Are there any estimated weights? */
extern int hasEstimatedWeights;

#endif /* CHARWEIGHT_H_ */

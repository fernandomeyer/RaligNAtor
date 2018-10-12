/*Copyright (C) 2013  Fernando Meyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <pthread.h>

#include "align.h"
#include "redblack/redblack.h"
#include "alignThread.h"

// Global variables for chaining
Chainparam *chainParamGlob;
const void *rbvalCR;
struct rbtree *rbCR; // Tree for storing the output from chain2dim
int matchSeqId;
unsigned int *uiSeqEndPosGlob;
bool bForwardSeqGlob;

typedef struct {
	unsigned long chaincounter;
} Counter;

static void chaining2gooutput(void *data, const GtChain2Dimmatchtable *matchtable, const GtChain2Dim *chain) {
	unsigned long idx, chainlength;
	int /*offset,*/ chainscore;
	Counter *counter = (Counter *) data;
	Chain *chainingResult = NULL;
	GtChain2Dimmatchvalues value;

	chainlength = gt_chain_chainlength(chain);
	chainscore  = gt_chain_chainscore(chain);

	if (counter->chaincounter >= chainParamGlob->topkscoring) {
		return;
	}

	if (chainlength < chainParamGlob->minchainlength || chainscore < chainParamGlob->minchainscore) {
		return;
	}

	if ((chainingResult = (Chain *) malloc(sizeof(Chain))) == NULL) {
		fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	chainingResult->seqId        = matchSeqId;
	chainingResult->chainLength  = gt_chain_chainlength(chain);
	chainingResult->chainScore   = chainscore;
	chainingResult->forwardSeq   = bForwardSeqGlob;
	chainingResult->fragmentInfo = NULL;
	chainingResult->bStoredInReverseOrder = gt_chain_storedinreverseorder(chain);

	if ((chainingResult->fragmentInfo = calloc(chainingResult->chainLength, sizeof(ChainFragment))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"fragmentInfo\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	for (idx = 0; idx < chainlength; idx++) {
		gt_chain_extractchainelem(&value, matchtable, chain, idx);
		chainingResult->fragmentInfo[idx].startpos[0] = value.startpos[0];
		chainingResult->fragmentInfo[idx].endpos[0]   = value.endpos[0];
		chainingResult->fragmentInfo[idx].weight      = value.weight;
		chainingResult->fragmentInfo[idx].match       = (Match *) value.whatever;

		chainingResult->fragmentInfo[idx].startpos[1] = value.startpos[1];
		chainingResult->fragmentInfo[idx].endpos[1]   = value.endpos[1];
	}

	// Insert chaining result into rb tree
	rbvalCR = rbsearch((void *) chainingResult, rbCR);
	// if(rbvalCR == NULL) identical chain is already in the tree

	counter->chaincounter++;
}

int** compute2DPositions(MultiPattern *multiPattern, bool bForwardStrand) {
	Pattern *pattern = multiPattern->pattern;
	int numPatterns = multiPattern->iNumPatterns;
	int **patternArtificialInt, ii;

	patternArtificialInt = (int **) calloc(numPatterns, sizeof(int*));

	if (multiPattern->bUsesStartPosParam) {

		for (ii = 0; ii < numPatterns; ii++) {
			patternArtificialInt[ii] = (int *) calloc(2, sizeof(int));
			patternArtificialInt[ii][0] = pattern[ii].startpos;
			patternArtificialInt[ii][1] = patternArtificialInt[ii][0] + pattern[ii].iLength - 1;
		}

	} else if (bForwardStrand) {

		for (ii = 0; ii < numPatterns; ii++) {
			patternArtificialInt[ii] = (int *) calloc(2, sizeof(int));
			if(ii == 0) {
				patternArtificialInt[ii][0] = 0;
				patternArtificialInt[ii][1] = pattern[ii].iLength - 1;
			} else {
				patternArtificialInt[ii][0] = patternArtificialInt[ii - 1][1] + 1;
				patternArtificialInt[ii][1] = patternArtificialInt[ii][0] + pattern[ii].iLength - 1;
			}
		}

	} else {

		for (ii = numPatterns - 1; ii >= 0; ii--) {
			patternArtificialInt[numPatterns - 1 - ii] = (int *) calloc(2, sizeof(int));
			if((numPatterns - 1 - ii) == 0) {
				patternArtificialInt[numPatterns - 1 - ii][0] = 0;
				patternArtificialInt[numPatterns - 1 - ii][1] = pattern[ii].iLength - 1;
			} else {
				patternArtificialInt[numPatterns - 1 - ii][0] = patternArtificialInt[numPatterns - 2 - ii][1] + 1;
				patternArtificialInt[numPatterns - 1 - ii][1] = patternArtificialInt[numPatterns - 1 - ii][0] + pattern[ii].iLength - 1;
			}
		}
	}

	return patternArtificialInt;
}

void chainAll(void *vrb, void *vrbCR, unsigned int *uiNumMatchesSeq, SearchParam *searchParam,
		MultiPattern *multiPattern, AffixArray *affixArray, bool bFreeMatches, bool bForwardSeq) {
	Pattern *patternP = multiPattern->pattern;
	int iNumPatterns = multiPattern->iNumPatterns;

	chainParamGlob  = searchParam->chainParam;
	bForwardSeqGlob = bForwardSeq;
	uiSeqEndPosGlob = affixArray->multiSeq->seqEndPos;

	RBLIST *rblist;
	struct rbtree *rb = vrb;
	Match *match;
	int **patternArtificialInt, i;
	unsigned int uiOffset;

	SearchParam searchParamLocal;
	rbtree *rbLocal = NULL;

	rbCR = vrbCR;

	bool hasErr;
	GtChain2Dim *chain = NULL;
	GtChain2Dimmode *chainmode = NULL;

	const unsigned int presortdim = 1U;

	if ((rblist = rbopenlist(rb)) == NULL) {
		fprintf(stderr, "Insufficient memory from rbopenlist().\n");
		exit(1);
	}

	patternArtificialInt = compute2DPositions(multiPattern, bForwardSeq);

	GtChain2Dimmatchtable *matchtable = NULL;
	unsigned int uiCount = 0;
	unsigned int uiPrevSeqId = UINT_MAX; // purpose: error checking


	//Matches must be sorted by sequence for chaining. If this is not the case, create a local rbtree
	if (searchParam->cPrintMatchesByScore > 0 && !searchParam->bPrintMatchesBySeq) {
		searchParamLocal.bPrintMatchesBySeq   = true;
		searchParamLocal.cPrintMatchesByScore = 0;

		if ((rbLocal = rbinit(compareRB, &searchParamLocal)) == NULL) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		while((match = (Match *) rbreadlist(rblist))) {
			rbsearch((void *) match, rbLocal);
		}
		rbcloselist(rblist);

		if ((rblist = rbopenlist(rbLocal)) == NULL) {
			fprintf(stderr, "Insufficient memory from rbopenlist().\n");
			exit(1);
		}
	}


	while((match = (Match *) rbreadlist(rblist))) {

		if (uiCount == 0) {
			if (uiPrevSeqId == match->iSeqId) {
				fprintf(stderr, "The number of matches stored in uiNumMatchesSeq are incorrect.\n");
				exit(1);
			}
			matchtable = gt_chain_matchtable_new(uiNumMatchesSeq[match->iSeqId]);
		}

		GtChain2Dimmatchvalues fragment;

		if (bForwardSeq) {
			fragment.startpos[0] = patternArtificialInt[match->patternStrandDirection->pattern->iId][0];
			fragment.endpos[0]   = patternArtificialInt[match->patternStrandDirection->pattern->iId][1];
		} else {
			fragment.startpos[0] = patternArtificialInt[iNumPatterns - 1 - match->patternStrandDirection->pattern->iId][0];
			fragment.endpos[0]   = patternArtificialInt[iNumPatterns - 1 - match->patternStrandDirection->pattern->iId][1];
		}

		fragment.weight      = patternP[match->patternStrandDirection->pattern->iId].weight - match->iCost;
		fragment.whatever    = match;

		if(match->iSeqId > 0)
			uiOffset = uiSeqEndPosGlob[match->iSeqId - 1] + 1;
		else
			uiOffset = 0;
		fragment.startpos[1] = match->uiPos - uiOffset + 1;
		fragment.endpos[1]   = match->uiEndPos - uiOffset + 1;

		gt_chain_matchtable_add(matchtable, &fragment);

		uiCount++;
		if (uiNumMatchesSeq[match->iSeqId] == uiCount) {
			hasErr = false;
			gt_chain_fillthegapvalues(matchtable);
			gt_chain_applyweight(chainParamGlob->weightfactor, matchtable);
			gt_chain_possiblysortmatches(NULL, matchtable, presortdim);
			chain = gt_chain_chain_new();
			chainmode = gt_chain_chainmode_new(chainParamGlob->maxgap,
					!chainParamGlob->islocal,
					chainParamGlob->globalparm,
					chainParamGlob->islocal,
					chainParamGlob->localparm,
					NULL);

			if (chainmode == NULL) {
				hasErr = true;
				break;
			}
			if (!hasErr) {
				Counter counter;
				counter.chaincounter = 0;

				matchSeqId = match->iSeqId;

				gt_chain_fastchaining(chainmode,
						chain,
						matchtable,
						true,
						presortdim,
						true,
						chaining2gooutput,
						&counter,
						NULL);
			}

			gt_chain_chain_delete(chain);
			gt_chain_chainmode_delete(chainmode);
			gt_chain_matchtable_delete(matchtable);

			uiCount = 0;
		}
		uiPrevSeqId = match->iSeqId;
		if (bFreeMatches) {
			rbdelete((void *)&match, rb);
			free(match);
		}
	}

	rbcloselist(rblist);

	for (i = 0; i < iNumPatterns; i++) {
		free(patternArtificialInt[i]);
	}
	free(patternArtificialInt);

	if (rbLocal != NULL)
		rbdestroy(rbLocal);
}

void reportChainingResults(void *vrb, Chainparam *chainparam, AffixArray *affixArray) {
	struct rbtree *rb = vrb;
	RBLIST *rblist;
	Chain *chainingResult;
	int ii, ij, ij2, iSeqPos, iDescLegth, iNumChains = 0;
	char *cv;
	bool bShowChain = chainparam->show;
	Match *match;
	bool bForwardSeq;
	bool **iupacTable = affixArray->alphabet->iupacTable;
	unsigned int uiOffset;
	int iNumSeqs = affixArray->multiSeq->numSeqs;

	if((rblist=rbopenlist(rb)) == NULL) {
		fprintf(stderr, "Insufficient memory from rbopenlist().\n");
		exit(1);
	}

	bool *bSeqHasChain = (bool *) calloc(iNumSeqs, sizeof(bool));

	if (!chainparam->show2) fprintf(stderr, "%c[sequence]                                	[chain score]	[chain length]	[strand]\n", LINESYMBOL);
	while((chainingResult=(Chain *)rbreadlist(rblist))) {
		bSeqHasChain[chainingResult->seqId] = true;

		if (chainparam->show2) continue;

		bForwardSeq = chainingResult->forwardSeq;
		iNumChains++;
		iDescLegth = strlen(affixArray->multiSeq->seqDescription[chainingResult->seqId]);
		printf(">");
		for(ii = 0; ii < 40; ii++) {
			if(ii < iDescLegth)
				printf("%c", affixArray->multiSeq->seqDescription[chainingResult->seqId][ii]);
			else
				printf(" ");
		}
		if(ii < iDescLegth)
			printf("...");
		else
			printf("   ");
		printf("	");
		printf("%d		", chainingResult->chainScore);
		printf("%d 		", chainingResult->chainLength);
		printf("%c\n", bForwardSeq ? '+' : '-');

		if (bShowChain && chainingResult->fragmentInfo != NULL) {
			if (bForwardSeq)
				for (ii = 0; ii < chainingResult->chainLength; ii++) {
					printf("%d %d ", chainingResult->fragmentInfo[ii].startpos[0], chainingResult->fragmentInfo[ii].endpos[0]);
					printf("%d %d ", chainingResult->fragmentInfo[ii].startpos[1], chainingResult->fragmentInfo[ii].endpos[1]);
					printf("%d\n", chainingResult->fragmentInfo[ii].weight);
				}
			else
				for (ii = chainingResult->chainLength; ii > 0; ii--) {
					printf("%d %d ", chainingResult->fragmentInfo[ii - 1].startpos[0], chainingResult->fragmentInfo[ii - 1].endpos[0]);
					printf("%d %d ", chainingResult->fragmentInfo[ii - 1].startpos[1], chainingResult->fragmentInfo[ii - 1].endpos[1]);
					printf("%d\n", chainingResult->fragmentInfo[ii - 1].weight);
				}

			if (bForwardSeq)
				for (ii = 0; ii < chainingResult->chainLength; ii++) {
					match = chainingResult->fragmentInfo[ii].match;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpInsertion)
							printf("%c", '-');
						else if (match->cSeqOperations[ij] == OpDeletion || match->cSeqOperations[ij] == OpReplacement)
							printf("%c", match->patternStrandDirection->structure[ij2++]);
						else
							printf("ERROR");
					}
					printf(" ");
				}
			else
				for (ii = chainingResult->chainLength; ii > 0; ii--) {
					match = chainingResult->fragmentInfo[ii - 1].match;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpInsertion)
							printf("%c", '-');
						else if (match->cSeqOperations[ij] == OpDeletion || match->cSeqOperations[ij] == OpReplacement)
							printf("%c", match->patternStrandDirection->structure[ij2++]);
						else
							printf("ERROR");
					}
					printf(" ");
				}
			printf("\n");

			if (bForwardSeq)
				for (ii = 0; ii < chainingResult->chainLength; ii++) {
					match = chainingResult->fragmentInfo[ii].match;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpInsertion)
							printf("%c", '-');
						else if (match->cSeqOperations[ij] == OpDeletion || match->cSeqOperations[ij] == OpReplacement)
							printf("%c", c(match->patternStrandDirection->seq[ij2++]));
						else
							printf("ERROR");
					}
					printf(" ");
				}
			else
				for (ii = chainingResult->chainLength; ii > 0; ii--) {
					match = chainingResult->fragmentInfo[ii - 1].match;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpInsertion)
							printf("%c", '-');
						else if (match->cSeqOperations[ij] == OpDeletion || match->cSeqOperations[ij] == OpReplacement)
							printf("%c", c(match->patternStrandDirection->seq[ij2++]));
						else
							printf("ERROR");
					}
					printf(" ");
				}
			printf("\n");

			// Print additional line
			if (bForwardSeq)
				for (ii = 0; ii < chainingResult->chainLength; ii++) {
					match = chainingResult->fragmentInfo[ii].match;
					cv = affixArray->multiSeq->convSequences + match->uiPos;
					iSeqPos = 0;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpInsertion || match->cSeqOperations[ij] == OpDeletion)
							printf("%c", ' ');
						else if (match->cSeqOperations[ij] == OpReplacement)
							printf("%c", match->patternStrandDirection->seq[ij2] == cv[iSeqPos] ? '|' : (iupacTable[(unsigned int) match->patternStrandDirection->seq[ij2]][(unsigned int) cv[iSeqPos]] ? '+' : ' '));
						else
							printf("ERROR");
						if (match->cSeqOperations[ij] == OpDeletion || match->cSeqOperations[ij] == OpReplacement)
							ij2++;
						if (match->cSeqOperations[ij] == OpInsertion || match->cSeqOperations[ij] == OpReplacement)
							iSeqPos++;
					}
					printf(" ");
				}
			else
				for (ii = chainingResult->chainLength; ii > 0; ii--) {
					match = chainingResult->fragmentInfo[ii - 1].match;
					cv = affixArray->multiSeq->convSequences + match->uiPos;
					iSeqPos = 0;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpInsertion || match->cSeqOperations[ij] == OpDeletion)
							printf("%c", ' ');
						else if (match->cSeqOperations[ij] == OpReplacement)
							printf("%c", match->patternStrandDirection->seq[ij2] == cv[iSeqPos] ? '|' : (iupacTable[(unsigned int) match->patternStrandDirection->seq[ij2]][(unsigned int) cv[iSeqPos]] ? '+' : ' '));
						else
							printf("ERROR");
						if (match->cSeqOperations[ij] == OpDeletion || match->cSeqOperations[ij] == OpReplacement)
							ij2++;
						if (match->cSeqOperations[ij] == OpInsertion || match->cSeqOperations[ij] == OpReplacement)
							iSeqPos++;
					}
					printf(" ");
				}
			printf("\n");

			if (bForwardSeq)
				for (ii = 0; ii < chainingResult->chainLength; ii++) {
					match = chainingResult->fragmentInfo[ii].match;
					cv = affixArray->multiSeq->convSequences + match->uiPos;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpDeletion)
							printf("%c", '-');
						else if (match->cSeqOperations[ij] == OpInsertion || match->cSeqOperations[ij] == OpReplacement)
							printf("%c", c(cv[ij2++]));
						else
							printf("ERROR");
					}
					printf(" ");
				}
			else
				for (ii = chainingResult->chainLength; ii > 0; ii--) {
					match = chainingResult->fragmentInfo[ii - 1].match;
					cv = affixArray->multiSeq->convSequences + match->uiPos;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						if (match->cSeqOperations[ij] == OpDeletion)
							printf("%c", '-');
						else if (match->cSeqOperations[ij] == OpInsertion || match->cSeqOperations[ij] == OpReplacement)
							printf("%c", c(cv[ij2++]));
						else
							printf("ERROR");
					}
					printf(" ");
				}
			printf("\n");

			if (bForwardSeq)
				for (ii = 0; ii < chainingResult->chainLength; ii++) {
					match = chainingResult->fragmentInfo[ii].match;
					cv = affixArray->multiSeq->convSequences + match->uiPos;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						printf("%c", match->cArcOperations[ij]);
					}
					printf(" ");
				}
			else
				for (ii = chainingResult->chainLength; ii > 0; ii--) {
					match = chainingResult->fragmentInfo[ii - 1].match;
					cv = affixArray->multiSeq->convSequences + match->uiPos;
					for (ij = 0, ij2 = 0; ij < match->uiNumOperations; ij++) {
						printf("%c", match->cArcOperations[ij]);
					}
					printf(" ");
				}
			printf("\n");

			free(chainingResult->fragmentInfo);
		}
		free(chainingResult);
	}

	if (!chainparam->show2) fprintf(stderr, "\nTotal number of chains: %d\n\n", iNumChains);

	//Print sequences with chains
	if (chainparam->show2) {
		for (ii = 0; ii < iNumSeqs; ii++) {
			if (bSeqHasChain[ii]) {

				if(ii > 0)
					uiOffset = affixArray->multiSeq->seqEndPos[ii - 1] + 1;
				else
					uiOffset = 0;

				cv = affixArray->multiSeq->convSequences + uiOffset;

				iSeqPos = 0;
				printf(">%s\n", affixArray->multiSeq->seqDescription[ii]);
				while (cv[iSeqPos] != $) {
					printf("%c", c(cv[iSeqPos++]));
				}
				printf("\n");
			}
		}
		printf("\n");
	}

	free(bSeqHasChain);

	rbcloselist(rblist);
}

void printMatchInOneLine(Match *match, AffixArray *affixArray) {
	unsigned int uiMatchingPos, uiOffset, uitemp, uitemp2, uik = 0;
	unsigned int uiNumOperations = match->uiNumOperations;
	char *cSeqOperations = match->cSeqOperations;

	if(match->iSeqId > 0)
		uiOffset = affixArray->multiSeq->seqEndPos[match->iSeqId - 1] + 1;
	else
		uiOffset = 0;
	uiMatchingPos = match->uiPos - uiOffset + 1;

	printf("%d	%d	%d	", match->patternStrandDirection->pattern->iId + 1, match->iSeqId + 1, uiMatchingPos);

	char *cu;
	char *cu_str;
	char *cv = affixArray->multiSeq->convSequences + match->uiPos;

	cu     = match->patternStrandDirection->seq;
	cu_str = match->patternStrandDirection->structure;

	// Print pattern
	uitemp2 = 0;
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperations[uitemp] == OpInsertion)
			printf("%c", '-');
		else if (cSeqOperations[uitemp] == OpDeletion || cSeqOperations[uitemp] == OpReplacement)
			printf("%c", c(cu[uitemp2++]));
		else
			printf("ERROR");
	}
	printf("	");

	// Print structure
	uitemp2 = 0;
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperations[uitemp] == OpInsertion)
			printf("%c", '-');
		else if (cSeqOperations[uitemp] == OpDeletion || cSeqOperations[uitemp] == OpReplacement)
			printf("%c", cu_str[uitemp2++]);
		else
			printf("ERROR");
	}
	printf("	");

	// Print database sequence
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperations[uitemp] == OpDeletion)
			printf("%c", '-');
		else if (cSeqOperations[uitemp] == OpInsertion || cSeqOperations[uitemp] == OpReplacement)
			printf("%c", c(cv[uik++]));
		else
			printf("ERROR");
	}
	printf("	");

	// Print database structure
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		printf("%c", match->cArcOperations[uitemp]);
	}
	printf("	");


	printf("%d	%d	%s", match->iCost, match->iScore, match->bForwardStrand ? "+" : "-");
}

void printMatch(Match *match, AffixArray *affixArray) {
	unsigned int uiMatchingPos, uiOffset;

	if(match->iSeqId > 0)
		uiOffset = affixArray->multiSeq->seqEndPos[match->iSeqId - 1] + 1;
	else
		uiOffset = 0;
	uiMatchingPos = match->uiPos - uiOffset + 1;

	printf("\n Pattern desc.: %s, Matching pos. = %d, Strand = %s, Score = %d, Edist = %d\n\n", match->patternStrandDirection->pattern->desc, uiMatchingPos, match->patternStrandDirection->bForwardStrand ? "forward" : "reverse", match->iScore, match->iCost);
	printAlignment(match, affixArray->multiSeq->convSequences + match->uiPos, 0, affixArray);
}

int compareRB(const void *pa, const void *pb, const void *config) {
	Match *ma = (Match *)pa;
	Match *mb = (Match *)pb;
	int iCmpResult;

	char cPrintMatchesByScore = ((SearchParam *) config)->cPrintMatchesByScore;
	bool bPrintMatchesBySeq   = ((SearchParam *) config)->bPrintMatchesBySeq;

	if (cPrintMatchesByScore > 0) {
		if (bPrintMatchesBySeq) {
			if (ma->iSeqId < mb->iSeqId) return -1;
			if (ma->iSeqId > mb->iSeqId) return 1;
		}

		if (cPrintMatchesByScore == 1) {
			if (ma->patternStrandDirection->pattern->iId < mb->patternStrandDirection->pattern->iId) return -1;
			if (ma->patternStrandDirection->pattern->iId > mb->patternStrandDirection->pattern->iId) return 1;

			if (ma->iScore > mb->iScore) return -1;
			if (ma->iScore < mb->iScore) return 1;

		} else if (cPrintMatchesByScore == 2) {
			if (ma->patternStrandDirection->pattern->iId < mb->patternStrandDirection->pattern->iId) return -1;
			if (ma->patternStrandDirection->pattern->iId > mb->patternStrandDirection->pattern->iId) return 1;

			if (ma->iScore < mb->iScore) return -1;
			if (ma->iScore > mb->iScore) return 1;
		}
	}

	if (ma->iSeqId < mb->iSeqId) return -1;
	if (ma->iSeqId > mb->iSeqId) return 1;

	if (ma->uiPos < mb->uiPos) return -1;
	if (ma->uiPos > mb->uiPos) return 1;

	if (ma->uiEndPos < mb->uiEndPos) return -1;
	if (ma->uiEndPos > mb->uiEndPos) return 1;

	if (ma->iCost < mb->iCost) return -1;
	if (ma->iCost > mb->iCost) return 1;

	if (ma->uiNumOperations < mb->uiNumOperations) return -1;
	if (ma->uiNumOperations > mb->uiNumOperations) return 1;

	if (ma->patternStrandDirection->pattern->iId < mb->patternStrandDirection->pattern->iId) return -1;
	if (ma->patternStrandDirection->pattern->iId > mb->patternStrandDirection->pattern->iId) return 1;

	if (cPrintMatchesByScore == 0) {
		if (ma->iScore > mb->iScore) return -1;
		if (ma->iScore < mb->iScore) return 1;
	}

	if ((iCmpResult = strcmp(ma->cSeqOperations, mb->cSeqOperations)) < 0)
		return -1;
	else if (iCmpResult > 0)
		return 1;

	return 0;
}

int compareRBCR(const void *pa, const void *pb, const void *config) {
	Chain *ma = (Chain *)pa;
	Chain *mb = (Chain *)pb;
	int i=0, longestChain;

	if(ma->chainScore > mb->chainScore) return -1;
	if(ma->chainScore < mb->chainScore) return 1;

	if(ma->seqId < mb->seqId) return -1;
	if(ma->seqId > mb->seqId) return 1;

	if(ma->forwardSeq > mb->forwardSeq) return -1;
	if(ma->forwardSeq < mb->forwardSeq) return 1;

	longestChain = ma->chainLength > mb->chainLength ? ma->chainLength : mb->chainLength;
	while (i < longestChain) {
		if (ma->fragmentInfo[i].startpos[1] > mb->fragmentInfo[i].startpos[1]) return -1;
		if (ma->fragmentInfo[i].startpos[1] < mb->fragmentInfo[i].startpos[1]) return 1;
		i++;
	}

	return 0;
}

void printSortedMatches(void *vrbF, void *vrbR, AffixArray *affixArray, SearchParam *searchParam) {
	int i;
	Match *match;
	RBLIST *rblist = NULL;
	rbtree *rb;

	if (vrbF == NULL && vrbR == NULL)
		return;

	for (i = 0; i < 2; i++) {
		if (i == 0 && vrbF != NULL)
			rb = vrbF;
		else if (i == 1 && vrbR != NULL)
			rb = vrbR;
		else
			continue;

		if ((rblist = rbopenlist(rb))==NULL) {
			fprintf(stderr, "Insufficient memory from rbopenlist() - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}

		while((match = (Match *)rbreadlist(rblist))) {
			pthread_mutex_lock(&mutexPrintSTDOUT);
			if (!searchParam->bPrintTable) {
				printf(">%s\n", affixArray->multiSeq->seqDescription[match->iSeqId]);
			}

			if (searchParam->bPrintTable) {
				printMatchInOneLine(match, affixArray);
			} else {
				printMatch(match, affixArray);
			}
			printf("\n");
			pthread_mutex_unlock(&mutexPrintSTDOUT);
		}
	}

	if (rblist != NULL)
		rbcloselist(rblist);
}

#define swap(x, a, b) { tmpmatch = match[a]; \
		match[a] = match[b]; match[b] = tmpmatch; }

void sortMatchesByPos(Match **match, int left, int right) {
	int i, last;
	Match *tmpmatch;

	if(left >= right)
		return;

	i = (left + right) >> 1; //i = (left + right)/2;
	swap(match, left, i);
	last = left;
	for(i = left+1; i <= right; i++){
		if(match[i]->uiPos < match[left]->uiPos) {
			++last;
			swap(array, last, i);
		}
	}
	swap(match, left, last);
	sortMatchesByPos(match, left, last-1);
	sortMatchesByPos(match, last+1, right);
}

unsigned int removeDuplicateMatches(Match **match, const int iMatches) {
	unsigned int i, j;

	if (iMatches == 0)
		return 0;
	if (iMatches == 1)
		return 1;

	j = 0;
	for (i = 1; i < iMatches; i++) {
		if (match[i]->uiPos <= match[j]->uiEndPos && match[i]->uiEndPos >= match[j]->uiPos && match[i]->iSeqId == match[j]->iSeqId) {
			if (match[i]->iScore > match[j]->iScore) {
				free(match[j]->cSeqOperations);
				free(match[j]->cArcOperations);
				free(match[j]);
				match[j] = match[i];
				match[i] = NULL;
			} else {
				free(match[i]->cSeqOperations);
				free(match[i]->cArcOperations);
				free(match[i]);
				match[i] = NULL;
			}
		} else {
			j++;
			if (i > j) {
				match[j] = match[i];
				match[i] = NULL;
			}
		}
	}

	return j + 1; // Number of remaining matches
}

void processAlignment_MT(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, const int iAbsoluteMatchPos,
		char *cSeqOperations, char *cArcOperations, const int iMinLength, const int iMaxLength, AffixArray *affixArray, ThreadVars *threadVars) {

	Pattern *pattern = patternStrandDirection->pattern;
	int ij, iIndels, iDeletions;
	int im = pattern->iLength;
	int iMatrixLine = patternStrandDirection->patternStructures->iLineIndex[im];
	unsigned int uiNumOperations;
	Match *match;

	SearchParam *searchParam   = threadVars->searchParam;
	unsigned int *uiNumMatches = &threadVars->uiNumMatchesPthread;
	Match **matches            = threadVars->matches;

	for (ij = iMinLength; ij <= iMaxLength; ij++) {
		if (P(edist, 0, iMatrixLine, ij) <= pattern->uiThreshold) {
			uiNumOperations = 0;
			iDeletions = iIndels = 0;

			getEditOperations(cSeqOperations, cArcOperations, operation, trace, &iIndels, &iDeletions, pattern->uiIndels, patternStrandDirection->patternStructures, &uiNumOperations, 0, im, ij, true);

			if (iIndels > pattern->uiIndels)
				continue;

			if (searchParam->bPrintOutMatches || searchParam->chainParam->isactive || searchParam->bFilterOverlaps) {
				cSeqOperations[uiNumOperations] = 0;

				match = (Match *) malloc(sizeof(Match));
				match->patternStrandDirection = patternStrandDirection;
				match->iSeqId     = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, iAbsoluteMatchPos);
				match->uiPos      = iAbsoluteMatchPos;
				match->uiEndPos   = match->uiPos + uiNumOperations - iDeletions - 1;

				match->cSeqOperations = (char *) calloc(uiNumOperations + 1, sizeof(char));
				strncpy(match->cSeqOperations, cSeqOperations, uiNumOperations);

				match->cArcOperations = (char *) calloc(uiNumOperations + 1, sizeof(char));
				strncpy(match->cArcOperations, cArcOperations, uiNumOperations);

				match->iCost           = P(edist, 0, iMatrixLine, ij);
				match->iScore          = pattern->weight - match->iCost;
				match->bForwardStrand  = patternStrandDirection->bForwardStrand;
				match->uiNumOperations = uiNumOperations;

				if (!searchParam->chainParam->show2 && !searchParam->bPrintMatchesBySeq && searchParam->cPrintMatchesByScore == 0 && !searchParam->bFilterOverlaps) {
					pthread_mutex_lock(&mutexPrintSTDOUT);
					if (searchParam->bPrintTable) {
						printMatchInOneLine(match, affixArray);
						printf("\n");
					} else if (searchParam->bPrintOutMatches) {
						printf(">%s\n", affixArray->multiSeq->seqDescription[match->iSeqId]);
						printMatch(match, affixArray);
						printf("\n");
					}
					pthread_mutex_unlock(&mutexPrintSTDOUT);
				}

				if (searchParam->bFilterOverlaps || searchParam->chainParam->isactive || searchParam->bPrintMatchesBySeq || searchParam->cPrintMatchesByScore > 0) {
					if((*uiNumMatches) > 0 && (*uiNumMatches) % BUFFER1 == 0) {
						if((threadVars->matches = matches = (Match **) realloc(matches, ((*uiNumMatches) + BUFFER1) * sizeof(Match *))) == NULL) {
							fprintf(stderr,"Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
							exit(1);
						}
					}
					matches[*uiNumMatches] = match;
				} else {
					free(match->cSeqOperations);
					free(match->cArcOperations);
					free(match);
				}
			}
			(*uiNumMatches)++;
		}
	}

}

void processAlignment_S_MT(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, const int iAbsoluteMatchPos,
		char *cSeqOperations, char *cArcOperations, const int iMinLength, const int iMaxLength, AffixArray *affixArray, ThreadVars *threadVars) {

	Pattern *pattern = patternStrandDirection->pattern;
	int ij, iIndels, iDeletions;
	int im = pattern->iLength;
	int iMatrixLine = patternStrandDirection->patternStructures->iLineIndex[im];
	unsigned int uiNumOperations;
	Match *match;

	SearchParam *searchParam   = threadVars->searchParam;
	unsigned int *uiNumMatches = &threadVars->uiNumMatchesPthread;
	Match **matches            = threadVars->matches;

	for (ij = iMinLength; ij <= iMaxLength; ij++) {
		if (S(edist, 0, iMatrixLine, ij) <= pattern->uiThreshold) {
			uiNumOperations = 0;
			iDeletions = iIndels = 0;

			getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, &iIndels, &iDeletions, pattern->uiIndels, patternStrandDirection->patternStructures, &uiNumOperations, 0, im, ij, true);

			// If searching online and match belongs to the previous partition...
			if (threadVars->iStartPos > 0 &&
					(iAbsoluteMatchPos + uiNumOperations - iDeletions - 1) <= threadVars->threadsManager->threadVars[threadVars->iThreadID - 1].iEndPos)
				continue;

			if (iIndels > pattern->uiIndels)
				continue;

			if (searchParam->bPrintOutMatches || searchParam->chainParam->isactive || searchParam->bFilterOverlaps) {
				cSeqOperations[uiNumOperations] = 0;

				match = (Match *) malloc(sizeof(Match));
				match->patternStrandDirection = patternStrandDirection;
				match->iSeqId     = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, iAbsoluteMatchPos);
				match->uiPos      = iAbsoluteMatchPos;
				match->uiEndPos   = match->uiPos + uiNumOperations - iDeletions - 1;

				match->cSeqOperations = (char *) calloc(uiNumOperations + 1, sizeof(char));
				strncpy(match->cSeqOperations, cSeqOperations, uiNumOperations);

				match->cArcOperations = (char *) calloc(uiNumOperations + 1, sizeof(char));
				strncpy(match->cArcOperations, cArcOperations, uiNumOperations);

				match->iCost           = S(edist, 0, iMatrixLine, ij);
				match->iScore          = pattern->weight - match->iCost;
				match->bForwardStrand  = patternStrandDirection->bForwardStrand;
				match->uiNumOperations = uiNumOperations;

				if (!searchParam->chainParam->show2 && !searchParam->bPrintMatchesBySeq && searchParam->cPrintMatchesByScore == 0 && !searchParam->bFilterOverlaps) {
					pthread_mutex_lock(&mutexPrintSTDOUT);
					if (searchParam->bPrintTable) {
						printMatchInOneLine(match, affixArray);
						printf("\n");
					} else if (searchParam->bPrintOutMatches) {
						printf(">%s\n", affixArray->multiSeq->seqDescription[match->iSeqId]);
						printMatch(match, affixArray);
						printf("\n");
					}
					pthread_mutex_unlock(&mutexPrintSTDOUT);
				}

				if (searchParam->bFilterOverlaps || searchParam->chainParam->isactive || searchParam->bPrintMatchesBySeq || searchParam->cPrintMatchesByScore > 0) {
					if((*uiNumMatches) > 0 && (*uiNumMatches) % BUFFER1 == 0) {
						if((threadVars->matches = matches = (Match **) realloc(matches, ((*uiNumMatches) + BUFFER1) * sizeof(Match *))) == NULL) {
							fprintf(stderr,"Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
							exit(1);
						}
					}
					matches[*uiNumMatches] = match;
				} else {
					free(match->cSeqOperations);
					free(match->cArcOperations);
					free(match);
				}
			}
			(*uiNumMatches)++;
		}
	}

}


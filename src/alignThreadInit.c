/*Copyright (C) 2013  Fernando Meyer, Ole Eigenbrod, Clemens Hoebsch

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "align.h"
#include "alignThread.h"
#include "redblack/redblack.h"

int alignGlobalInit(MultiPattern *multiPattern, AffixArray *affixArray, SearchParam *searchParam) {
	int ii;
	bool bSeqSearched, bRevSeqSearched;
	Pattern *pattern;
	PatternStrandDirection *patternStrandDirection;

	bSeqSearched = bRevSeqSearched = false;

	for (ii = 0; ii < multiPattern->iNumPatterns;) {
		pattern = &multiPattern->pattern[ii];

		if (pattern->forwardStrand != NULL && !bSeqSearched) {
			patternStrandDirection = pattern->forwardStrand;
			bSeqSearched = true;

		} else if (pattern->reverseStrand != NULL && !bRevSeqSearched) {
			patternStrandDirection = pattern->reverseStrand;
			bRevSeqSearched = true;

		} else {
			ii++;
			bSeqSearched = bRevSeqSearched = false;
			continue;
		}
		patternStrandDirection->matches = NULL;

		if (!searchParam->chainParam->show2) {
			fprintf(stderr, "Costs: Replacement     = %d\n", pattern->cost->iCostReplacement);
			fprintf(stderr, "       Deletion        = %d\n", pattern->cost->iCostDeletion);
			fprintf(stderr, "       Arc-breaking    = %d\n", pattern->cost->iCostArcBreak);
			fprintf(stderr, "       Arc-altering    = %d\n", pattern->cost->iCostArcAltering);
			fprintf(stderr, "       Arc-removing    = %d\n\n", pattern->cost->iCostArcRemoving);
			fflush(stderr);
		}

		alignGlobal(patternStrandDirection, affixArray, searchParam);
	}

	return EXIT_SUCCESS;
}

void threadMain(MultiPattern *multiPattern, AffixArray *affixArray,
		bool **compRules, SearchParam *searchParam, int iNumThreads, int iNumPartitions) {

	if (searchParam->cVariant == V_GLOBAL) {
		alignGlobalInit(multiPattern, affixArray, searchParam);
		return;
	}

	pthread_t * threads;
	int *matchArray;
	int nPerPart = 1;
	//int *patternSetup;
	int ii, ii2;
	int iNumPatterns = multiPattern->iNumPatterns;
	unsigned int uiTotalNumMatches = 0;
	ThreadsManager threadsManager;
	ThreadVars *threadVars;
	RBPointers *rbPointers;
	struct timeval start, end;
	Match *match;

	pthread_mutex_init(&mutexThreadCounter, NULL);
	pthread_mutex_init(&mutexPrintSTDOUT, NULL);
	pthread_mutex_init(&mutexPrintSTDERROR, NULL);
	pthread_cond_init(&condThreadCount, NULL);

	threadsManager.iNumThreads               = iNumThreads;
	threadsManager.iThreadCounter            = 0;
	threadsManager.uiNumMatchesPerPatternFor = 0;
	threadsManager.uiNumMatchesPerPatternRev = 0;
	threadsManager.matchesPerPatternFor      = NULL;
	threadsManager.matchesPerPatternRev      = NULL;

	if (searchParam->bSearchForwardString && searchParam->bSearchReverseString)
		threadsManager.iNumJobs = multiPattern->iNumPatterns * iNumPartitions * 2;
	else
		threadsManager.iNumJobs = multiPattern->iNumPatterns * iNumPartitions;

	if((threads = malloc(sizeof(pthread_t) * (threadsManager.iNumJobs))) == NULL){
		fprintf(stderr, "Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	if((matchArray = calloc(iNumPatterns, sizeof(int)))== NULL){
		fprintf(stderr, "Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	rbPointers = initRBPointers(searchParam, affixArray);
	threadsManager.threadVars = threadVars = initThreads(&threadsManager, iNumThreads, iNumPartitions, multiPattern,
			searchParam, affixArray, rbPointers, nPerPart);

	gettimeofday(&start, NULL);

	for (ii = 0; ii < threadsManager.iNumJobs; ii++) {
		pthread_mutex_lock(&mutexThreadCounter);
		while (threadsManager.iThreadCounter >= iNumThreads) {
			pthread_cond_wait(&condThreadCount, &mutexThreadCounter);
		}
		threadsManager.iThreadCounter++;
		pthread_mutex_unlock(&mutexThreadCounter);

		pthread_create(&threads[ii], NULL, alignThreaded, (void *) &threadVars[ii]);
	}

	for (ii = 0; ii < threadsManager.iNumJobs; ii++) { //TODO: join only threads of the same pattern
		pthread_join(threads[ii], NULL);
	}

	if (searchParam->bFilterOverlaps || searchParam->chainParam->isactive || searchParam->bPrintMatchesBySeq || searchParam->cPrintMatchesByScore > 0) {
		for (ii = 0; ii < iNumPatterns; ii++) {
			if (searchParam->bFilterOverlaps) {
				if (searchParam->bSearchForwardString) {
					sortMatchesByPos(threadsManager.matchesPerPatternFor[ii], 0, threadsManager.uiNumMatchesPerPatternFor[ii] - 1);
					threadsManager.uiNumMatchesPerPatternFor[ii] = removeDuplicateMatches(threadsManager.matchesPerPatternFor[ii], threadsManager.uiNumMatchesPerPatternFor[ii]);
				}
				if (searchParam->bSearchReverseString) {
					sortMatchesByPos(threadsManager.matchesPerPatternRev[ii], 0, threadsManager.uiNumMatchesPerPatternRev[ii] - 1);
					threadsManager.uiNumMatchesPerPatternRev[ii] = removeDuplicateMatches(threadsManager.matchesPerPatternRev[ii], threadsManager.uiNumMatchesPerPatternRev[ii]);
				}
			}

			if (searchParam->bSearchForwardString) {
				for (ii2 = 0; ii2 < threadsManager.uiNumMatchesPerPatternFor[ii]; ii2++) {
					match = threadsManager.matchesPerPatternFor[ii][ii2];
					if (rbsearch((void *) threadsManager.matchesPerPatternFor[ii][ii2], rbPointers->rbF) != NULL) {
						searchParam->uiNumMatchesSeq[match->iSeqId]++;
					}
				}
			}
			if (searchParam->bSearchReverseString) {
				for (ii2 = 0; ii2 < threadsManager.uiNumMatchesPerPatternRev[ii]; ii2++) {
					match = threadsManager.matchesPerPatternRev[ii][ii2];
					if (rbsearch((void *) threadsManager.matchesPerPatternRev[ii][ii2], rbPointers->rbR) != NULL) {
						searchParam->uiNumMatchesReverseSeq[match->iSeqId]++;
					}
				}
			}
		}
	}

	uiTotalNumMatches = 0;
	for (ii = 0; ii < iNumPatterns; ii++) {
		if (searchParam->bSearchForwardString)
			uiTotalNumMatches += threadsManager.uiNumMatchesPerPatternFor[ii];
		if (searchParam->bSearchReverseString)
			uiTotalNumMatches += threadsManager.uiNumMatchesPerPatternRev[ii];

		if (!searchParam->chainParam->show2 && searchParam->cVariant != V_GLOBAL) {
			if (searchParam->bSearchForwardString)
				fprintf(stderr, "Number of matches for pattern %s in the forward strand(s): %d\n", multiPattern->pattern[ii].desc, threadsManager.uiNumMatchesPerPatternFor[ii]);
			if (searchParam->bSearchReverseString)
				fprintf(stderr, "Number of matches for pattern %s in the reverse complement strand(s): %d\n", multiPattern->pattern[ii].desc, threadsManager.uiNumMatchesPerPatternRev[ii]);
		}


		free(threadVars[ii * iNumPartitions].patternStrandDirection->matches); //todo: check for both strands
	}

	if (!searchParam->chainParam->show2 && searchParam->cVariant != V_GLOBAL)
		fprintf(stderr, "Total number of matches: %d\n", uiTotalNumMatches);

	gettimeofday(&end, NULL);
	if (!searchParam->chainParam->show2) fprintf(stderr, "Time: %.4f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );

	if (!searchParam->chainParam->show2 && searchParam->bPrintOutMatches && (searchParam->bPrintMatchesBySeq || searchParam->cPrintMatchesByScore != 0 || searchParam->bFilterOverlaps)) {
		if (searchParam->bPrintTable) fprintf(stderr, "[pattern id]	[target seq. id]	[matching pos.]	[pattern string]	[pattern structure]	[target string]	[target structure]	[edist]	[score]	[strand]\n");
		printSortedMatches(rbPointers->rbF, rbPointers->rbR, affixArray, searchParam);
	}
	chain(searchParam, rbPointers, multiPattern, affixArray);

	free(threads);
	free(matchArray);
	free(rbPointers); //todo: free related pointers

	pthread_mutex_destroy(&mutexThreadCounter);
	pthread_mutex_destroy(&mutexPrintSTDOUT);
	pthread_mutex_destroy(&mutexPrintSTDERROR);

	pthread_cond_destroy(&condThreadCount);

	freeThreads(threadVars, iNumPatterns, iNumPartitions);
}

void *alignThreaded(void * tmanageIn) { //called by each thread
	ThreadVars * threadVars;
	threadVars = (ThreadVars *) tmanageIn;

	SearchParam *searchParam = threadVars->searchParam;
	PatternStrandDirection *patternStrandDirection = threadVars->patternStrandDirection;
	Pattern *pattern = patternStrandDirection->pattern;

	if (threadVars->patternStrandDirection->bForwardStrand) {
		if (threadVars->startIdx == 0 && threadVars->iStartPos == 0) {
			if (!searchParam->chainParam->show2 && searchParam->cVariant != V_GLOBAL) {
				pthread_mutex_lock(&mutexPrintSTDERROR);
				fprintf(stderr, "\n%cSearching for pattern %s in the forward sequence(s)... \n", LINESYMBOL, pattern->desc);
				fflush(stdout);
				pthread_mutex_unlock(&mutexPrintSTDERROR);
			}
		}
	} else {
		if (threadVars->startIdx == 0 && threadVars->iStartPos == 0) {
			if (!searchParam->chainParam->show2 && searchParam->cVariant != V_GLOBAL) {
				pthread_mutex_lock(&mutexPrintSTDERROR);
				fprintf(stderr, "\n%cSearching for pattern %s in the reverse complement sequence(s)... \n", LINESYMBOL, pattern->desc);
				fflush(stdout);
				pthread_mutex_unlock(&mutexPrintSTDERROR);
			}
		}
	}

	if(threadVars->startIdx == 0 && threadVars->iStartPos == 0){
		if (!searchParam->chainParam->show2) {
			pthread_mutex_lock(&mutexPrintSTDERROR);
			if (searchParam->cVariant != V_GLOBAL) {
				fprintf(stderr, "Cost threshold (edist) = %d\n", pattern->uiThreshold);
				fprintf(stderr, "Max. allowed indels    = %d\n", pattern->uiIndels);
				fprintf(stderr, "Min./Max. match length = %d / %d\n", pattern->iLength > pattern->uiIndels ? pattern->iLength - pattern->uiIndels : 1, pattern->iLength + pattern->uiIndels);
				fprintf(stderr, "Max. match score       = %d\n", pattern->weight);
			}
			fprintf(stderr, "Costs: Replacement     = %d\n", pattern->cost->iCostReplacement);
			fprintf(stderr, "       Deletion        = %d\n", pattern->cost->iCostDeletion);
			fprintf(stderr, "       Arc-breaking    = %d\n", pattern->cost->iCostArcBreak);
			fprintf(stderr, "       Arc-altering    = %d\n", pattern->cost->iCostArcAltering);
			fprintf(stderr, "       Arc-removing    = %d\n\n", pattern->cost->iCostArcRemoving);
			fflush(stdout);
			pthread_mutex_unlock(&mutexPrintSTDERROR);
		}

		if (!searchParam->chainParam->show2 && !searchParam->bFilterOverlaps && searchParam->bPrintOutMatches && searchParam->bPrintTable && !searchParam->bPrintMatchesBySeq && searchParam->cPrintMatchesByScore == 0) {
			pthread_mutex_lock(&mutexPrintSTDERROR);
			fprintf(stderr, "[pattern id]	[target seq. id]	[matching pos.]	[pattern string]	[pattern structure]	[target string]	[target structure]	[edist]	[score]	[strand]\n");
			pthread_mutex_unlock(&mutexPrintSTDERROR);
		}
	}

	switch (searchParam->cVariant) {
	case V_LESA:
		alignESABasedThreaded(threadVars);
		break;
	case V_LGSLINK:
		generalSLinkAlignNEWThreaded(threadVars);
		break;
	case V_SCAN:
		alignOnlineScanMainThreaded(threadVars);
		break;
	case V_SCANLA:
		alignOnlineScanMainLA_NEWThreaded(threadVars);
		break;
	}

	pthread_mutex_lock(&mutexThreadCounter);
	threadVars->threadsManager->iThreadCounter--;
	pthread_mutex_unlock(&mutexThreadCounter);
	pthread_cond_broadcast(&condThreadCount); //submit new counter value to waiting patterns

	pthread_exit(NULL);
}

int *initCheckpoints(AffixArray *affixArray){
	int iRemoveOverlapsInterval;
	int ii;
	static int iRemoveControl[FILTERFREQ];

	iRemoveOverlapsInterval = (affixArray->length / FILTERFREQ) + 1;

	iRemoveControl[0] = iRemoveOverlapsInterval;
	for (ii = 1; ii < FILTERFREQ; ii++) {
		iRemoveControl[ii] = iRemoveControl[ii - 1] + iRemoveOverlapsInterval;
	}
	return iRemoveControl;
}

ThreadVars *initThreads(ThreadsManager *threadsManager, int numThreads, int iNumPartitions,
		MultiPattern *multiPattern, SearchParam *searchParam,
		AffixArray *affixArray, RBPointers *rbPointers, int nPerPart) {

	ThreadVars *threadVars;
	int ii, j, k = 0, k2, *iFilterOverlapsControl;
	bool *bVisitedSuffixTab;
	int iEndOverDollarSymbol = 0;
	Mutex *mutex;

	iFilterOverlapsControl = initCheckpoints(affixArray);

	char *cSeqs;
	int imax, iAlignmentLength, icount;

	if (searchParam->bSearchForwardString && searchParam->bSearchReverseString) {
		if((threadVars = malloc(sizeof(ThreadVars) * (multiPattern->iNumPatterns * iNumPartitions * 2))) == NULL){
			fprintf(stderr, "Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	} else {
		if((threadVars = malloc(sizeof(ThreadVars) * (multiPattern->iNumPatterns * iNumPartitions))) == NULL){
			fprintf(stderr, "Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	}

	if (searchParam->cVariant == V_LESA || searchParam->cVariant == V_LGSLINK){
		nPerPart = ((affixArray->length - affixArray->multiSeq->numSeqs + 1)/ iNumPartitions);
	} else if(searchParam->cVariant == V_SCAN || searchParam->cVariant == V_SCANLA){
		nPerPart = ((affixArray->length) /  iNumPartitions);
	}

	if (searchParam->bSearchForwardString) {
		threadsManager->matchesPerPatternFor = (Match ***) calloc(multiPattern->iNumPatterns, sizeof(Match **));
		threadsManager->uiNumMatchesPerPatternFor = (unsigned int *) calloc(multiPattern->iNumPatterns, sizeof(unsigned int));
	}
	if (searchParam->bSearchReverseString) {
		threadsManager->matchesPerPatternRev = (Match ***) calloc(multiPattern->iNumPatterns, sizeof(Match **));
		threadsManager->uiNumMatchesPerPatternRev = (unsigned int *) calloc(multiPattern->iNumPatterns, sizeof(unsigned int));
	}

	for (ii = 0; ii < multiPattern->iNumPatterns; ii++) {
		if(searchParam->cVariant == V_LGSLINK) {
			bVisitedSuffixTab = (bool *) calloc(affixArray->length, sizeof(bool));
		} else {
			bVisitedSuffixTab = NULL;
			iEndOverDollarSymbol = 0;
		}

		if (searchParam->bSearchForwardString)
			threadsManager->matchesPerPatternFor[ii] = (Match **) malloc(sizeof(Match *));
		else
			threadsManager->matchesPerPatternRev[ii] = (Match **) malloc(sizeof(Match *));

		mutex = (Mutex *) malloc(sizeof(Mutex));

		pthread_mutex_init(&mutex->mutexMatchCount, NULL);
		pthread_mutex_init(&mutex->mutexTree, NULL);

		imax = multiPattern->pattern[ii].iLength + multiPattern->pattern[ii].uiIndels;
		iAlignmentLength = imax > affixArray->length ? affixArray->length : imax;

		for (j = 0; j < iNumPartitions; j++, k++) {
			threadVars[k].iThreadID = k;
			threadVars[k].threadsManager = threadsManager;

			if (searchParam->bSearchForwardString) {
				threadVars[k].patternStrandDirection    = multiPattern->pattern[ii].forwardStrand;
			} else {
				threadVars[k].patternStrandDirection = multiPattern->pattern[ii].reverseStrand;
			}

			threadVars[k].patternStrandDirection->matches = NULL;

			threadVars[k].affixArray          = affixArray;
			threadVars[k].searchParam         = searchParam;
			threadVars[k].rbPointers          = rbPointers;
			threadVars[k].iFilterOverlapsControl = iFilterOverlapsControl;
			threadVars[k].bVisitedSuffix      = bVisitedSuffixTab;

			threadVars[k].uiNumMatchesPthread = 0;
			threadVars[k].mutex               = mutex;

			threadVars[k].matches            = (Match **) calloc(BUFFER1, sizeof(Match *));

			if (searchParam->cVariant == V_LESA || searchParam ->cVariant == V_LGSLINK){
				threadVars[k].iStartPos = 0;
				threadVars[k].startIdx = j * nPerPart;
				if (j == (iNumPartitions - 1)) {
					threadVars[k].endIdx = affixArray->length - affixArray->multiSeq->numSeqs;
				} else {
					threadVars[k].endIdx = ((j + 1) * nPerPart) - 1;
				}

			} else if (searchParam->cVariant == V_SCAN || searchParam->cVariant == V_SCANLA){
				threadVars[k].startIdx = 0;
				threadVars[k].iStartPos = j * nPerPart - iEndOverDollarSymbol;

				cSeqs = affixArray->multiSeq->convSequences + threadVars[k].iStartPos;
				if (*(cSeqs) == $) threadVars[k].iStartPos++;

				if (j == (iNumPartitions - 1)) {
					threadVars[k].iEndPos = affixArray->length - 2;
				} else {
					threadVars[k].iEndPos = ((j + 1) * nPerPart) - 1;

					cSeqs = affixArray->multiSeq->convSequences + threadVars[k].iEndPos;
					icount = 0;

					while ((cSeqs - icount) != affixArray->multiSeq->convSequences &&
							*(cSeqs - icount++) != $ && icount < iAlignmentLength);

					if (icount == iAlignmentLength || (cSeqs - icount) == affixArray->multiSeq->convSequences) {
						icount = 0;
						while (*(cSeqs + ++icount) != $ && icount < iAlignmentLength);

						threadVars[k].iEndPos += icount - 1;

						iEndOverDollarSymbol = 0;
					} else {
						threadVars[k].iEndPos  -= icount;
						iEndOverDollarSymbol = icount - 1;
					}

				}

			}
		}
	}

	k2 = k;
	k  = 0;
	if (searchParam->bSearchForwardString && searchParam->bSearchReverseString) {
		for (ii = 0; ii < multiPattern->iNumPatterns; ii++) {
			if(searchParam->cVariant == V_LGSLINK) {
				bVisitedSuffixTab = (bool *) calloc(affixArray->length, sizeof(bool));
			} else {
				bVisitedSuffixTab = NULL;
				iEndOverDollarSymbol = 0;
			}

			threadsManager->matchesPerPatternRev[ii] = (Match **) malloc(sizeof(Match *));

			mutex = (Mutex *) malloc(sizeof(Mutex));

			pthread_mutex_init(&mutex->mutexMatchCount, NULL);
			pthread_mutex_init(&mutex->mutexTree, NULL);

			imax = multiPattern->pattern[ii].iLength + multiPattern->pattern[ii].uiIndels;
			iAlignmentLength = imax > affixArray->length ? affixArray->length : imax;

			for (j = 0; j < iNumPartitions; j++, k2++, k++) {
				threadVars[k2].iThreadID = k2;
				threadVars[k2].threadsManager = threadsManager;

				threadVars[k2].patternStrandDirection    = multiPattern->pattern[ii].reverseStrand;

				threadVars[k2].patternStrandDirection->matches = NULL;

				threadVars[k2].affixArray          = affixArray;
				threadVars[k2].searchParam         = searchParam;
				threadVars[k2].rbPointers          = rbPointers;
				threadVars[k2].iFilterOverlapsControl = iFilterOverlapsControl;
				threadVars[k2].bVisitedSuffix      = bVisitedSuffixTab;

				threadVars[k2].uiNumMatchesPthread = 0;
				threadVars[k2].mutex               = mutex;

				threadVars[k2].matches = (Match **) calloc(BUFFER1, sizeof(Match *));

				threadVars[k2].iStartPos = threadVars[k].iStartPos;
				threadVars[k2].iEndPos   = threadVars[k].iEndPos;
				threadVars[k2].startIdx  = threadVars[k].startIdx;
				threadVars[k2].endIdx    = threadVars[k].endIdx;
			}
		}
	}

	//fprintf(stderr, "{%d}\n", threadVars->matchesPerPatternRev == NULL);
	//exit(0);

	return threadVars;
}
/*free tmanage struct*/
void freeThreads(ThreadVars * threadVars, int iNumPatterns, int iNumPartitions) {
	int ii, iLastThread;
	SearchParam *searchParam = threadVars->searchParam;
	ThreadsManager *threadsManager = threadVars->threadsManager;

	for (ii = 0; ii < iNumPatterns; ii++) {
		free(threadVars[ii * iNumPartitions].bVisitedSuffix);

		pthread_mutex_destroy(&threadVars[ii * iNumPartitions].mutex->mutexMatchCount);
		pthread_mutex_destroy(&threadVars[ii * iNumPartitions].mutex->mutexTree);
		free(threadVars[ii * iNumPartitions].mutex);

		if (searchParam->bSearchForwardString) {
			free(threadsManager->matchesPerPatternFor[ii]);
		} else if (searchParam->bSearchReverseString) {
			free(threadsManager->matchesPerPatternRev[ii]);
		}
	}
	if (searchParam->bSearchForwardString && searchParam->bSearchReverseString) {
		iLastThread = iNumPatterns * iNumPartitions;
		for (ii = 0; ii < iNumPatterns; ii++) {
			free(threadVars[iLastThread + (ii * iNumPartitions)].bVisitedSuffix);

			pthread_mutex_destroy(&threadVars[iLastThread + (ii * iNumPartitions)].mutex->mutexMatchCount);
			pthread_mutex_destroy(&threadVars[iLastThread + (ii * iNumPartitions)].mutex->mutexTree);
			free(threadVars[iLastThread + (ii * iNumPartitions)].mutex);
		}
	}

	if (searchParam->bSearchForwardString) {
		free(threadsManager->matchesPerPatternFor);
		free(threadsManager->uiNumMatchesPerPatternFor);
	}
	if (searchParam->bSearchReverseString) {
		free(threadsManager->matchesPerPatternRev);
		free(threadsManager->uiNumMatchesPerPatternRev);
	}

	free(threadVars);
}

RBPointers *initRBPointers(SearchParam *searchParam, AffixArray *affixArray) {
	rbtree *rbF, *rbR, *rbCR = NULL;
	RBPointers *rbPointers = (RBPointers *) malloc(sizeof(RBPointers));

	if (searchParam->chainParam->isactive) {
		if ((rbCR = rbinit(compareRBCR, NULL )) == NULL ) { // Tree for chaining report
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	if (searchParam->bSearchForwardString
			&& (searchParam->chainParam->isactive
					|| searchParam->bPrintMatchesBySeq
					|| searchParam->cPrintMatchesByScore != 0
					|| searchParam->bFilterOverlaps)) {
		if ((searchParam->uiNumMatchesSeq = (unsigned int *) calloc(affixArray->multiSeq->numSeqs, sizeof(unsigned int))) == NULL ) {
			fprintf(stderr, "Memory allocation failed for \"uiNumMatchesSeq\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((rbF = rbinit(compareRB, searchParam)) == NULL ) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	} else {
		searchParam->uiNumMatchesSeq = NULL;
		rbF = NULL;
	}

	if (searchParam->bSearchReverseString
			&& (searchParam->chainParam->isactive
					|| searchParam->bPrintMatchesBySeq
					|| searchParam->cPrintMatchesByScore != 0
					|| searchParam->bFilterOverlaps)) {
		if ((searchParam->uiNumMatchesReverseSeq = (unsigned int *) calloc(affixArray->multiSeq->numSeqs, sizeof(unsigned int))) == NULL ) {
			fprintf(stderr, "Memory allocation failed for \"uiNumMatchesReverseSeq\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((rbR = rbinit(compareRB, searchParam)) == NULL ) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	} else {
		searchParam->uiNumMatchesReverseSeq = NULL;
		rbR = NULL;
	}

	rbPointers->rbCR = rbCR;
	rbPointers->rbF = rbF;
	rbPointers->rbR = rbR;

	return rbPointers;
}

void chain(SearchParam *searchParam, RBPointers *rbPointers, MultiPattern *multiPattern, AffixArray *affixArray) {
	struct timeval start_chain, end_chain;

	if (searchParam->chainParam->isactive) {
		if (!searchParam->chainParam->show2) fprintf(stderr, "\n%cChaining matches... ", LINESYMBOL);
		fflush(stdout);

		gettimeofday(&start_chain, NULL);

		if (searchParam->bSearchForwardString)
			chainAll(rbPointers->rbF, rbPointers->rbCR, searchParam->uiNumMatchesSeq,
					searchParam, multiPattern, affixArray, false,
					true);
		if (searchParam->bSearchReverseString)
			chainAll(rbPointers->rbR, rbPointers->rbCR,
					searchParam->uiNumMatchesReverseSeq,
					searchParam, multiPattern, affixArray, false,
					false);

		if (!searchParam->chainParam->show2) fprintf(stderr, "done\n");

		gettimeofday(&end_chain, NULL);

		if (!searchParam->chainParam->show2)
			fprintf(stderr, "Time:     %.4f ms\n\n", (double) ((end_chain.tv_sec - start_chain.tv_sec) * 1000 + (end_chain.tv_usec - start_chain.tv_usec) / 1000.0));

		reportChainingResults(rbPointers->rbCR, searchParam->chainParam, affixArray);
	}
}

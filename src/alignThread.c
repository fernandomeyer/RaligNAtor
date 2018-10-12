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

void deliverMatches(ThreadVars *threadVars) {
	if (threadVars->uiNumMatchesPthread == 0) return;

	Match **matches;
	unsigned int uiNumMatches;

	int iPatternID                  = threadVars->patternStrandDirection->pattern->iId;
	Match **matchesThread           = threadVars->matches;
	unsigned int uiNumMatchesThread = threadVars->uiNumMatchesPthread;
	SearchParam *searchParam        = threadVars->searchParam;

	if (threadVars->patternStrandDirection->bForwardStrand) {
		uiNumMatches = threadVars->threadsManager->uiNumMatchesPerPatternFor[iPatternID];
		threadVars->threadsManager->uiNumMatchesPerPatternFor[iPatternID] += uiNumMatchesThread;
		matches = threadVars->threadsManager->matchesPerPatternFor[iPatternID];
	} else {
		uiNumMatches = threadVars->threadsManager->uiNumMatchesPerPatternRev[iPatternID];
		threadVars->threadsManager->uiNumMatchesPerPatternRev[iPatternID] += uiNumMatchesThread;
		matches = threadVars->threadsManager->matchesPerPatternRev[iPatternID];
	}

	threadVars->uiNumMatchesPthread = 0;

	if (searchParam->bFilterOverlaps || searchParam->chainParam->isactive || searchParam->bPrintMatchesBySeq || searchParam->cPrintMatchesByScore > 0) {
		if((matches = (Match **) realloc(matches, (uiNumMatches + uiNumMatchesThread) * sizeof(Match *))) == NULL) {
			fprintf(stderr,"Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
			exit(1);
		}

		if (threadVars->patternStrandDirection->bForwardStrand) {
			threadVars->threadsManager->matchesPerPatternFor[iPatternID] = matches;
			memcpy(threadVars->threadsManager->matchesPerPatternFor[iPatternID] + uiNumMatches, matchesThread, uiNumMatchesThread * sizeof(Match *));
		} else {
			threadVars->threadsManager->matchesPerPatternRev[iPatternID] = matches;
			memcpy(threadVars->threadsManager->matchesPerPatternRev[iPatternID] + uiNumMatches, matchesThread, uiNumMatchesThread * sizeof(Match *));
		}

		free(threadVars->matches);
		if((threadVars->matches = (Match **) calloc(BUFFER1, sizeof(Match *))) == NULL) {
			fprintf(stderr,"Memory allocation failed. File: %s Line: %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	}

}

void *generalSLinkAlignNEWThreaded(void * tmanageIn) {
	ThreadVars * threadVars;
	threadVars = (ThreadVars *) tmanageIn;

	int startIdx = threadVars->startIdx;
	int endIdx   = threadVars->endIdx;

	AffixArray *affixArray = threadVars->affixArray;

	PatternStrandDirection *patternStrandDirection = threadVars->patternStrandDirection;
	Pattern *pattern = patternStrandDirection->pattern;

	SearchParam *searchParam = threadVars->searchParam;
	bool *bVisitedSuffix     = threadVars->bVisitedSuffix;

	int ii, ij, ik, ikStart=0, ikEnd=0, iRegion, iLastRegion=0, iLastRegion2=0, iLastRegionCheck=0, iLastRegionCheck2=0;
	int iLcpCheck, iDepthEnded, iDepthEnded2 = 0, iRoundCount = 0;
	int iSuffix, iSuffix2, iSuffixToMark, iLcp = 0, iLcp2 = 0, iReadingDepth = 0;
	int iUnalignedPrefixLength, iLink, iLinkDown, iLinkUp, iNumRemainingColumnsToMove, iR, ijGoal;
	bool bSuccess, bFollowedSLink;
	LR LR;
	char *seqs, *cSeqOperations, *cArcOperations;
	PatternRegion *patternRegion;
	Mutex *mutex = threadVars->mutex;

	int *iSufArray = affixArray->esa->xarray;
	int iLength = affixArray->length;
	int iNumNonEmptySuffixes = iLength - affixArray->multiSeq->numSeqs + 1;

	int im = pattern->iLength;

	int imax = im + pattern->uiIndels;
	int imin = im > pattern->uiIndels ? im - pattern->uiIndels : 1;

	cSeqOperations = (char *) calloc(imax + 1, sizeof(char));
	cArcOperations = (char *) calloc(imax + 1, sizeof(char));

	unsigned int ***edist, ***edist2, ***trace, ***trace2;
	char ***operation, ***operation2;

	unsigned int **edistTmp, **traceTmp;
	char **operationTmp;

	edistTmp     = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	traceTmp     = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	operationTmp = (char **) calloc(imax + 1, sizeof(char *));

	edist     = uiNew3DMatrixX(im, imax);
	operation = cNew3DMatrixX(im, imax);
	trace     = uiNew3DMatrixX(im, imax);

	/*pthread_mutex_lock(&mutexPatternSetup);
	while (threadVars->patternSetup[threadVars->pattern->iId] != 1) {
		pthread_cond_wait(&condPatternSetup, &mutexPatternSetup);
	}
	pthread_mutex_unlock(&mutexPatternSetup);*/

	initializeRowsAndColumnsRegions_S(edist, operation, trace, pattern->cost, patternStrandDirection->structure, im, imax, patternStrandDirection->patternStructures);

	edist2     = uiNew3DMatrixX(im, imax);
	operation2 = cNew3DMatrixX(im, imax);
	trace2     = uiNew3DMatrixX(im, imax);
	initializeRowsAndColumnsRegions_S(edist2, operation2, trace2, pattern->cost, patternStrandDirection->structure, im, imax, patternStrandDirection->patternStructures);

	ComputedEntries *compEntries  = newComputedEntriesMatrix(patternStrandDirection, true);
	ComputedEntries *compEntries2 = newComputedEntriesMatrix(patternStrandDirection, true);

	PatternRegion *patternRegionZero = patternStrandDirection->patternStructures->patternRegionInsideOut;

	int *iRoundComputedPerRow     = (int *) calloc(imax + 1, sizeof(int));
	int *iRoundComputedPerRowInit = (int *) calloc(imax + 1, sizeof(int));
	int *iEmptyArray = (int *) calloc(imax + 1, sizeof(int));

	int iRemoveControlIndex = 0;
	/*Setting iRemoveControl(checkpoints) for no-overlaps option*/
	for (ii = 0; ii < FILTERFREQ; ii++){
		if(threadVars->iFilterOverlapsControl[ii] > startIdx){
			iRemoveControlIndex = ii;
			break;
		}
	}

	if (searchParam->bUseSequenceBasedFilter)
		applySequenceBasedFilter(affixArray, patternStrandDirection, bVisitedSuffix, startIdx, endIdx);

	iSuffix2 = startIdx;
	while (iSuffix2 < endIdx) {
		bFollowedSLink = false;
		iRoundCount = 0;

		while (bVisitedSuffix[iSufArray[iSuffix2]] && ++iSuffix2 < iNumNonEmptySuffixes) {
			if ((iLcp = lcpvalue(iSuffix2)) < iLcp2)
				iLcp2 = iLcp;
		}

		iSuffix = iSuffix2;

		/*while (!bVisitedSuffix[iSufArray[iSuffix]]) {
			pthread_mutex_lock(&mutexVisitedSuffix);
			bVisitedSuffix[iSufArray[iSuffix]] = true;
			pthread_mutex_unlock(&mutexVisitedSuffix);*/
			//bVisitedSuffix[iSufArray[iSuffix]] = true;
		while (true) {
			/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
			if (!bVisitedSuffix[iSufArray[iSuffix]]) {
				bVisitedSuffix[iSufArray[iSuffix]] = true;
			} else {
				pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
				break;
			}
			pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/
			if (!__sync_bool_compare_and_swap(&bVisitedSuffix[iSufArray[iSuffix]], false, true))
				break;

			seqs = affixArray->multiSeq->convSequences + iSufArray[iSuffix];

			bSuccess = false;

			if (bFollowedSLink) {

				iReadingDepth--;
				while (*(seqs + iReadingDepth) != $ && ++iReadingDepth < imax);
				if (iReadingDepth < imin) break;
				//pthread_mutex_lock(&mutexPattern);
				iLastRegion = alignESA_LA_NEW_LastColumnS1(edist2,
						operation2,
						trace2,
						patternStrandDirection,
						compEntries2,
						iEmptyArray,
						seqs,
						iReadingDepth,
						imax,
						affixArray,
						&LR,
						iRoundCount,
						iRoundComputedPerRow,
						&bSuccess,
						true,
						iLastRegion,
						&iLastRegionCheck2);
				iDepthEnded = LR.iR;
				//pthread_mutex_unlock(&mutexPattern);

				ikEnd--;

			} else {


				if (iLcp2 < imax) {
					iReadingDepth = iLcp2;
					while (*(seqs + iReadingDepth) != $ && ++iReadingDepth < imax);
					if (iReadingDepth < imin) break;
				} else {
					iReadingDepth = imax;
				}

				memcpy(iRoundComputedPerRow, iRoundComputedPerRowInit, (imax + 1) * sizeof(int));
				//pthread_mutex_lock(&mutexPattern);
				iLastRegionCheck2 = iLastRegionCheck = iLastRegion = iLastRegion2 = alignESA_LA_NEW_S(edist,
						operation,
						trace,
						patternStrandDirection,
						compEntries,
						iEmptyArray,
						seqs,
						iReadingDepth,
						imax,
						affixArray,
						iLcp2,
						iDepthEnded2,
						iLastRegion2,
						iRoundComputedPerRow,
						&bSuccess,
						&LR,
						true);
				iDepthEnded = iDepthEnded2 = LR.iR;
				//pthread_mutex_unlock(&mutexPattern);
				ikEnd = patternRegionZero[iLastRegionCheck].ikBigEnd;

			}


			iRoundCount++;

			if (!bSuccess) {
				//bVisitedSuffix[iSufArray[iSuffix]] = true;

				iSuffixToMark = iSuffix + 1;

				iLcpCheck = iDepthEnded + pattern->uiIndels;
				iLcp = lcpvalue(iSuffixToMark);

				if (bFollowedSLink) {
					while (iLcp >= iLcpCheck && (bVisitedSuffix[iSufArray[iSuffixToMark]] = true) && ++iSuffixToMark < iLength && !bVisitedSuffix[iSufArray[iSuffixToMark]]) {
						iLcp = lcpvalue(iSuffixToMark);
					}
					/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
					while (iLcp >= iLcpCheck) {
						bVisitedSuffix[iSufArray[iSuffixToMark]] = true;
						if (++iSuffixToMark >= iLength || bVisitedSuffix[iSufArray[iSuffixToMark]]) break;
						iLcp = lcpvalue(iSuffixToMark);
					}
					pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/

					if ((iSuffix - 1) >= 0) {
						iSuffixToMark = iSuffix - 1;
						iLcp = lcpvalue(iSuffixToMark + 1);

						while (iLcp >= iLcpCheck && (bVisitedSuffix[iSufArray[iSuffixToMark]] = true) && --iSuffixToMark >= 0 && !bVisitedSuffix[iSufArray[iSuffixToMark]]) {
							iLcp = lcpvalue(iSuffixToMark + 1);
						}
						/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
						while (iLcp >= iLcpCheck) {
							bVisitedSuffix[iSufArray[iSuffixToMark]] = true;
							if (--iSuffixToMark < 0 || bVisitedSuffix[iSufArray[iSuffixToMark]]) break;
							iLcp = lcpvalue(iSuffixToMark + 1);
						}
						pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/
					}
				} else {
					iLcp2 = iLcp;
					while (iLcp >= iLcpCheck && (bVisitedSuffix[iSufArray[iSuffixToMark]] = true) && ++iSuffixToMark < iLength) {
						if ((iLcp = lcpvalue(iSuffixToMark)) < iLcp2)
							iLcp2 = iLcp;
					}
					/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
					while (iLcp >= iLcpCheck) {
						bVisitedSuffix[iSufArray[iSuffixToMark]] = true;
						if (++iSuffixToMark >= iLength) break;
						if ((iLcp = lcpvalue(iSuffixToMark)) < iLcp2) iLcp2 = iLcp;
					}
					pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/
					iSuffix2 = iSuffixToMark;
				}

				if ((iUnalignedPrefixLength = LR.iL - pattern->uiIndels - 1) > 0) {

					iLink = affixArray->esa->xsufinv[iSufArray[iSuffix] + iUnalignedPrefixLength];
					iLcpCheck = LR.iR + pattern->uiIndels - iUnalignedPrefixLength;

					iLinkDown = iLink + 1;
					while (iLinkDown < iLength && lcpvalue(iLinkDown) >= iLcpCheck) {
						if (iSufArray[iLinkDown] - iUnalignedPrefixLength >= 0) {
							//pthread_mutex_lock(&mutex->mutexVisitedSuffix);
							bVisitedSuffix[iSufArray[iLinkDown] - iUnalignedPrefixLength] = true;
							//pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
						}
						iLinkDown++;
					}

					iLinkUp = iLink - 1;
					while (iLinkUp >= 0 && lcpvalue(iLinkUp + 1) >= iLcpCheck) {
						if (iSufArray[iLinkUp] - iUnalignedPrefixLength >= 0) {
							//pthread_mutex_lock(&mutex->mutexVisitedSuffix);
							bVisitedSuffix[iSufArray[iLinkUp] - iUnalignedPrefixLength] = true;
							//pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
						}
						iLinkUp--;
					}
				}

			} else {

				/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
				if (bVisitedSuffix[iSufArray[iSuffix]]) {
					pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
					break;
				}
				bVisitedSuffix[iSufArray[iSuffix]] = true;
				pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/

				iLcpCheck = (iDepthEnded + pattern->uiIndels) > imax ? imax : iDepthEnded + pattern->uiIndels;

				if (bFollowedSLink) {
					iSuffixToMark = iSuffix;

					processAlignment_S_MT(edist2,
							operation2,
							trace2,
							patternStrandDirection,
							iSufArray[iSuffixToMark],
							cSeqOperations,
							cArcOperations,
							imin,
							iReadingDepth,
							affixArray,
							//searchParam,
							//&uiNumMatches, //&(tmanage->uiNumMatchesPattern[tmanage->pattern->iId]),
							threadVars);

					if (++iSuffixToMark < iLength) {
						iLcp = lcpvalue(iSuffixToMark);
						while (iLcp >= iLcpCheck) {
							/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
							if (bVisitedSuffix[iSufArray[iSuffixToMark]]) {
								pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
								break;
							}
							bVisitedSuffix[iSufArray[iSuffixToMark]] = true;
							pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/

							//bool __sync_bool_compare_and_swap (type *ptr, type oldval, type newval, ...)
							if (!__sync_bool_compare_and_swap(&bVisitedSuffix[iSufArray[iSuffixToMark]], false, true))
								break;

							processAlignment_S_MT(edist2,
									operation2,
									trace2,
									patternStrandDirection,
									iSufArray[iSuffixToMark],
									cSeqOperations,
									cArcOperations,
									imin,
									iReadingDepth,
									affixArray,
									//searchParam,
									//&uiNumMatches, //&(tmanage->uiNumMatchesPattern[tmanage->pattern->iId]),
									threadVars);

							if (++iSuffixToMark >= iLength) {
								break;
							}
							iLcp = lcpvalue(iSuffixToMark);
						}
					}

					if ((iSuffixToMark = iSuffix - 1) >= 0) {
						iLcp = lcpvalue(iSuffixToMark + 1);

						while (iLcp >= iLcpCheck) {
							/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
							if (bVisitedSuffix[iSufArray[iSuffixToMark]]) {
								pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
								break;
							}
							bVisitedSuffix[iSufArray[iSuffixToMark]] = true;
							pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/
							if (!__sync_bool_compare_and_swap(&bVisitedSuffix[iSufArray[iSuffixToMark]], false, true))
								break;

							processAlignment_S_MT(edist2,
									operation2,
									trace2,
									patternStrandDirection,
									iSufArray[iSuffixToMark],
									cSeqOperations,
									cArcOperations,
									imin,
									iReadingDepth,
									affixArray,
									//searchParam,
									//&uiNumMatches, //&(tmanage->uiNumMatchesPattern[tmanage->pattern->iId]),
									threadVars);

							if (--iSuffixToMark < 0) {
								break;
							}
							iLcp = lcpvalue(iSuffixToMark + 1);
						}
					}
					//pthread_mutex_unlock(&mutex->mutexVisitedSuffix);

					//pthread_mutex_unlock(&mutexPattern);

				} else {

					//pthread_mutex_lock(&mutexPattern);

					processAlignment_S_MT(edist,
							operation,
							trace,
							patternStrandDirection,
							iSufArray[iSuffix],
							cSeqOperations,
							cArcOperations,
							imin,
							iReadingDepth,
							affixArray,
							threadVars);

					//if (++iSuffix < iLength) {
					if (++iSuffix <= endIdx) {
						iLcp = lcpvalue(iSuffix);

						while (iLcp >= iLcpCheck) {
							/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
							if (bVisitedSuffix[iSufArray[iSuffix]]) {
								pthread_mutex_unlock(&mutex->mutexVisitedSuffix);
								break;
							}
							bVisitedSuffix[iSufArray[iSuffix]] = true;
							pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/

							if (!__sync_bool_compare_and_swap(&bVisitedSuffix[iSufArray[iSuffix]], false, true))
								break;

							processAlignment_S_MT(edist,
									operation,
									trace,
									patternStrandDirection,
									iSufArray[iSuffix],
									cSeqOperations,
									cArcOperations,
									imin,
									iReadingDepth,
									affixArray,
									threadVars);

							//if (++iSuffix >= iLength) {
							if (++iSuffix > endIdx) {
								break;
							}
							iLcp = lcpvalue(iSuffix);
						}
					}
					iSuffix2 = --iSuffix;

					if (searchParam->bFilterOverlaps && iSuffix > threadVars->iFilterOverlapsControl[iRemoveControlIndex]) {
						iRemoveControlIndex++;

						pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
						deliverMatches(threadVars);
						if (threadVars->patternStrandDirection->bForwardStrand) {
							sortMatchesByPos(threadVars->threadsManager->matchesPerPatternFor[pattern->iId], 0,
									threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] - 1);
							threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternFor[pattern->iId],
									threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId]);
						} else {
							sortMatchesByPos(threadVars->threadsManager->matchesPerPatternRev[pattern->iId], 0,
									threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] - 1);
							threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternRev[pattern->iId],
									threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId]);
						}
						pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);
					}
				}

			}

			iSuffix = affixArray->esa->xsufinv[iSufArray[iSuffix] + 1];

			/*pthread_mutex_lock(&mutex->mutexVisitedSuffix);
			if (bVisitedSuffix[iSufArray[iSuffix]])
				bVisitedSuffix1 = true;
			else
				bVisitedSuffix1 = false;
			pthread_mutex_unlock(&mutex->mutexVisitedSuffix);*/

			if (((iSuffix < startIdx) || (iSuffix > iSuffix2)) && !bVisitedSuffix[iSufArray[iSuffix]] && *(seqs + imin) != $) {

				if (bFollowedSLink) {
					ikStart = 0;
					if (patternRegionZero[iLastRegionCheck2].ikBigEnd > ikEnd)
						ikEnd = patternRegionZero[iLastRegionCheck2].ikBigEnd;

					if (ikStart == 0) ikStart = 1;
					ik = ikStart - 1;
					memcpy(edistTmp, edist2[ik], sizeof(unsigned int **) * (imax - ik)); // Backup first matrices
					memcpy(traceTmp, trace2[ik], sizeof(unsigned int **) * (imax - ik));
					memcpy(operationTmp, operation2[ik], sizeof(char **) * (imax - ik));

					ij = imax - ikStart;
					for (ik = ikStart; ik <= ikEnd; ik++) {
						memcpy(edist2[ik - 1] + 1, edist2[ik] + 1, sizeof(unsigned int **) * (imax - ik)); // Copy all columns
						memcpy(trace2[ik - 1] + 1, trace2[ik] + 1, sizeof(unsigned int **) * (imax - ik));
						memcpy(operation2[ik - 1] + 1, operation2[ik] + 1, sizeof(char **) * (imax - ik));

						edist2[ik][imax - ik] = edistTmp[ij]; // Copy a column of backed-up matrix to the last column of another matrix
						trace2[ik][imax - ik] = traceTmp[ij];
						operation2[ik][imax - ik] = operationTmp[ij];
						ij--;

						S(edist2, ik, 0, imax - ik) = S(edist2, ik, 0, imax - ik - 1) + pattern->cost->iCostDeletion;
					}

					iNumRemainingColumnsToMove = imax - ikEnd - 1;
					if (iNumRemainingColumnsToMove > 0) {
						memcpy(edist2[ikEnd] + 1, edistTmp + 1, sizeof(unsigned int **) * iNumRemainingColumnsToMove);
						memcpy(trace2[ikEnd] + 1, traceTmp + 1, sizeof(unsigned int **) * iNumRemainingColumnsToMove);
						memcpy(operation2[ikEnd] + 1, operationTmp + 1, sizeof(char **) * iNumRemainingColumnsToMove);
					}

				} else {


					for (ik = 1; ik < imax; ik++) {
						ijGoal = imax - ik;

						for (ij = 1; ij <= ijGoal; ij++) {
							memcpy(edist2[ik - 1][ij], edist[ik][ij], sizeof(unsigned int) * (im + 1));
							memcpy(operation2[ik - 1][ij], operation[ik][ij], sizeof(char) * (im + 1));
							memcpy(trace2[ik - 1][ij], trace[ik][ij], sizeof(unsigned int) * (im + 1));
						}
					}

					for (iRegion = 0; iRegion <= iLastRegion; iRegion++) {
						patternRegion = &patternRegionZero[iRegion];
						iR = patternRegion->iR;
						if (patternRegion->bIsUnpaired) {
							for (ii = patternRegion->iL; ii <= iR; ii++) {
								memcpy(compEntries2[iRegion].ijLastPerRowPerMatrix[ii], compEntries[iRegion].ijLastPerRowPerMatrix[ii], imax * sizeof(int));
							}
						} else {
							memcpy(compEntries2[iRegion].ijLast, compEntries[iRegion].ijLast, imax * sizeof(int));
						}
						if (patternRegion->bIsArc) {
							memcpy(compEntries2[iRegion].ikFirst, compEntries[iRegion].ikFirst, imax * sizeof(int));
						}
					}
					bFollowedSLink = true;

				}

			} else {
				bVisitedSuffix[iSufArray[iSuffix]] = true;
				break;
			}

		}

	}

	free(edistTmp);
	free(traceTmp);
	free(operationTmp);
	uiFree3DMatrixX(edist, im, imax);
	cFree3DMatrixX(operation, im, imax);
	uiFree3DMatrixX(trace, im, imax);

	uiFree3DMatrixX(edist2, im, imax);
	cFree3DMatrixX(operation2, im, imax);
	uiFree3DMatrixX(trace2, im, imax);

	free(cSeqOperations);
	free(cArcOperations);

	free(iRoundComputedPerRow);
	free(iRoundComputedPerRowInit);

	free(iEmptyArray);
//printf("=============================\n");
	freeComputedEntriesMatrix(compEntries, patternStrandDirection, true);
	freeComputedEntriesMatrix(compEntries2, patternStrandDirection, true);

	//printf("uiNumMatches=%d\n", uiNumMatches);

	//threadVars->uiNumMatchesPthread = uiNumMatches;
	pthread_mutex_lock(&mutex->mutexMatchCount);
	deliverMatches(threadVars);
	pthread_mutex_unlock(&mutex->mutexMatchCount);

	return EXIT_SUCCESS;
}

void *alignESABasedThreaded(void * tmanageIn) {
	ThreadVars * threadVars;
	threadVars = (ThreadVars *) tmanageIn;

	int startIdx = threadVars->startIdx;
	int endIdx   = threadVars->endIdx;

	PatternStrandDirection *patternStrandDirection = threadVars->patternStrandDirection;
	Pattern *pattern = patternStrandDirection->pattern;

	AffixArray *affixArray = threadVars->affixArray;

	int ii, iRemoveControlIndex = 0, iSuffix, im, iReadingDepth;
	int imin = 0, imax = 0, iLcp = 0, iMinLcp;
	int iDepthEnded, iLcpCheck;

	unsigned int ***edist, ***trace;
	unsigned int uiNumMatches = 0,
			uiNumMatchesPrev = 0;
	int iLength = affixArray->length;
	char *cSeqOperations, *cArcOperations;
	int *ijLastZero;

	bool bSuccess;
	LR LR;

	char ***operation;
	char *cu_str, *seqs;

	bool bComputed = false;

	cu_str = (char *) patternStrandDirection->structure;

	im = pattern->iLength;

	imax = im + pattern->uiIndels;
	imin = im > pattern->uiIndels ? im - pattern->uiIndels : 1;
	cSeqOperations = (char *) calloc(imax + 1, sizeof(char));
	cArcOperations = (char *) calloc(imax + 1, sizeof(char));

	edist = uiNew3DMatrix(im, imax);
	operation = cNew3DMatrix(im, imax);
	trace = uiNew3DMatrix(im, imax);

	initializeRowsAndColumns(edist, operation, trace, pattern->cost, cu_str, im, imax, patternStrandDirection->patternStructures);

	ijLastZero = (int *) calloc(imax, sizeof(int));

	ComputedEntries *compEntries = newComputedEntriesMatrix(patternStrandDirection, false);
	/*Setting iFilterOverlapsControl(checkpoints) for no-overlaps option*/
	for (ii = 0; ii < FILTERFREQ; ii++){
		if(threadVars->iFilterOverlapsControl[ii] > startIdx){
			iRemoveControlIndex = ii;
			break;
		}
	}

	iDepthEnded = im;

	for (iSuffix = startIdx; iSuffix <= endIdx; iSuffix++) {
		if (bComputed) {
			iLcpCheck = iDepthEnded + pattern->uiIndels;

			iMinLcp = iLcp = lcpvalue(iSuffix);

			if (uiNumMatches == uiNumMatchesPrev) {
				while (iLcp >= iLcpCheck && ++iSuffix < iLength) {
					if ((iLcp = lcpvalue(iSuffix)) < iMinLcp)
						iMinLcp = iLcp;
				}
			}
			iLcp = iMinLcp;

		} else {
			iLcp = 0;
		}

		seqs = (char *) affixArray->multiSeq->convSequences + affixArray->esa->xarray[iSuffix];
		iReadingDepth = 0;

		while (*(seqs + iReadingDepth) != $ && ++iReadingDepth < imax);

		if (iReadingDepth < imin) {
			bComputed = false;
			iDepthEnded = 0;
			continue;
		}

		uiNumMatchesPrev = uiNumMatches;
		bComputed = true;
		bSuccess = false;

		alignESA_LA_NEW(edist,
				operation,
				trace,
				patternStrandDirection,
				compEntries,
				ijLastZero,
				seqs,
				iReadingDepth,
				imax,
				affixArray,
				iLcp,
				iDepthEnded,
				true,
				&bSuccess,
				&LR);
		iDepthEnded = LR.iR;

		if (!bSuccess) continue;

		iLcpCheck = iDepthEnded + pattern->uiIndels;
		do {

			//pthread_mutex_lock(&mutexPattern);
			processAlignment_MT(edist,
					operation,
					trace,
					patternStrandDirection,
					affixArray->esa->xarray[iSuffix],
					cSeqOperations,
					cArcOperations,
					imin,
					iReadingDepth,
					affixArray,
					threadVars); //&(tmanage->uiNumMatchesPattern[tmanage->pattern->iId]));

			if (threadVars->searchParam->bFilterOverlaps && iSuffix > threadVars->iFilterOverlapsControl[iRemoveControlIndex]) {
				iRemoveControlIndex++;

				pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
				deliverMatches(threadVars);
				if (threadVars->patternStrandDirection->bForwardStrand) {
					sortMatchesByPos(threadVars->threadsManager->matchesPerPatternFor[pattern->iId], 0,
							threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] - 1);
					threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternFor[pattern->iId],
							threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId]);
				} else {
					sortMatchesByPos(threadVars->threadsManager->matchesPerPatternRev[pattern->iId], 0,
							threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] - 1);
					threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternRev[pattern->iId],
							threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId]);
				}
				pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);
			}

		} while (++iSuffix <= endIdx && lcpvalue(iSuffix) >= iLcpCheck);
		iSuffix--;
	}

	uiFree3DMatrix(edist, im, imax);
	cFree3DMatrix(operation, im, imax);
	uiFree3DMatrix(trace, im, imax);

	free(cSeqOperations);
	free(cArcOperations);

	free(ijLastZero);

	freeComputedEntriesMatrix(compEntries, patternStrandDirection, false);

	pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
	deliverMatches(threadVars);
	pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);

	return EXIT_SUCCESS;
}

void * alignOnlineScanMainThreaded(void *tmanageIn) {
	ThreadVars *threadVars;
	threadVars = (ThreadVars *) tmanageIn;

	AffixArray *affixArray = threadVars->affixArray;

	PatternStrandDirection *patternStrandDirection = threadVars->patternStrandDirection;
	Pattern *pattern = patternStrandDirection->pattern;

	int ii, im, ik, icount;
	int imax = 0, imin = 0;
	int iAlignmentLength;
	int iStartPos;
	int iEndPos;

	unsigned int ***edist, ***trace, **uiPPEdist, **uiPPTrace;

	char ***operation, **cPPOperation;
	char *cu_str, *seqs, *cSeqOperations, *cArcOperations;

	cu_str = (char *) patternStrandDirection->structure;

	im = pattern->iLength;

	imax = im + pattern->uiIndels;
	imin = im > pattern->uiIndels ? im - pattern->uiIndels : 1;

	edist     = uiNew3DMatrixX(im, imax);
	operation = cNew3DMatrixX(im, imax);
	trace     = uiNew3DMatrixX(im, imax);

	iStartPos = threadVars->iStartPos;
	iEndPos   = threadVars->iEndPos - imin + 1;
	//iSubstrLength = threadVars->iEndPos - iStartPos + 1;
	iAlignmentLength = imax;

	cSeqOperations = (char *) calloc(imax + 1, sizeof(char));
	cArcOperations = (char *) calloc(imax + 1, sizeof(char));

	initializeRowsAndColumnsRegions_S(edist, operation, trace, pattern->cost, cu_str, im, imax, patternStrandDirection->patternStructures);

	uiPPEdist    = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	uiPPTrace    = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	cPPOperation = (char **) calloc(imax + 1, sizeof(char *));

	int iRemoveControlIndex = 0;
	for (ii = 0; ii < FILTERFREQ; ii++){
		if(threadVars->iFilterOverlapsControl[ii] > iStartPos){
			iRemoveControlIndex = ii;
			break;
		}
	}

	bool bFirstWindow = true;

	for (ii = iStartPos; ii <= iEndPos; ii++) {

		if (bFirstWindow) {

			seqs             = affixArray->multiSeq->convSequences + ii;
			iAlignmentLength = imax;
			icount           = 0;

			while (*(seqs + ++icount) != $ && icount < iAlignmentLength);

			if (icount < iAlignmentLength) {
				if (icount >= imin) {
					iAlignmentLength = icount;
				} else {
					ii += icount;
					continue; //continue with next sequence after $
				}
			}

			alignESA_S(edist,
					operation,
					trace,
					patternStrandDirection,
					seqs,
					iAlignmentLength,
					affixArray,
					0);

			processAlignment_S_MT(edist,
					operation,
					trace,
					patternStrandDirection,
					ii, //iStartPos, //uiSeqStartPos,
					cSeqOperations,
					cArcOperations,
					imin,
					iAlignmentLength,
					affixArray,
					threadVars);

			bFirstWindow = false;
		} else {

			seqs = affixArray->multiSeq->convSequences + ii;

			while (ii + iAlignmentLength - 1 > threadVars->iEndPos)
				iAlignmentLength--;
			if (*(seqs + iAlignmentLength - 1) == $)
				iAlignmentLength--;

			if (iAlignmentLength < imin && *(seqs + iAlignmentLength) == $) {
				bFirstWindow = true;
				ii += iAlignmentLength;
				continue;
			}

			// Move columns
			memcpy(uiPPEdist, edist[0], sizeof(unsigned int **) * imax); // Backup first matrices
			memcpy(uiPPTrace, trace[0], sizeof(unsigned int **) * imax);
			memcpy(cPPOperation, operation[0], sizeof(char **) * imax);

			for (ik = 1; ik < imax; ik++) {
				memcpy(edist[ik - 1] + 1, edist[ik] + 1, sizeof(unsigned int **) * (imax - ik)); // Copy all columns
				memcpy(trace[ik - 1] + 1, trace[ik] + 1, sizeof(unsigned int **) * (imax - ik));
				memcpy(operation[ik - 1] + 1, operation[ik] + 1, sizeof(char **) * (imax - ik));

				edist[ik][imax - ik] = uiPPEdist[ik]; // Copy a column of backed-up matrix to the last column of another matrix
				trace[ik][imax - ik] = uiPPTrace[ik];
				operation[ik][imax - ik] = cPPOperation[ik];

				S(edist, ik, 0, imax - ik) = S(edist, ik, 0, imax - ik - 1) + pattern->cost->iCostDeletion;
			}
			// End Move columns

			alignLastColumn_S(edist,
					operation,
					trace,
					patternStrandDirection,
					seqs,
					iAlignmentLength,
					affixArray);

			processAlignment_S_MT(edist,
					operation,
					trace,
					patternStrandDirection,
					ii, //iStartPos + ij,//uiSeqStartPos + ij,
					cSeqOperations,
					cArcOperations,
					imin,
					iAlignmentLength,
					affixArray,
					threadVars);

			if (threadVars->searchParam->bFilterOverlaps && ii > threadVars->iFilterOverlapsControl[iRemoveControlIndex]) {
				iRemoveControlIndex++;

				pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
				deliverMatches(threadVars);
				if (threadVars->patternStrandDirection->bForwardStrand) {
					sortMatchesByPos(threadVars->threadsManager->matchesPerPatternFor[pattern->iId], 0,
							threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] - 1);
					threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternFor[pattern->iId],
							threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId]);
				} else {
					sortMatchesByPos(threadVars->threadsManager->matchesPerPatternRev[pattern->iId], 0,
							threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] - 1);
					threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternRev[pattern->iId],
							threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId]);
				}
				pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);
			}
		}

	}

	uiFree3DMatrixX(edist, im, imax);
	cFree3DMatrixX(operation, im, imax);
	uiFree3DMatrixX(trace, im, imax);

	free(cSeqOperations);
	free(cArcOperations);

	free(uiPPEdist);
	free(uiPPTrace);
	free(cPPOperation);

	//threadVars->uiNumMatchesPthread = uiNumMatches;
	pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
	deliverMatches(threadVars);
	pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);
	//printf("E\n");
	return EXIT_SUCCESS;

}

void * alignOnlineScanMainLA_NEWThreaded(void * tmanageIn) {
	ThreadVars * threadVars;
	threadVars = (ThreadVars *) tmanageIn;

	AffixArray *affixArray = threadVars->affixArray;

	PatternStrandDirection *patternStrandDirection = threadVars->patternStrandDirection;
	Pattern *pattern = patternStrandDirection->pattern;

	int ii, ij, im, ik, ikStart, ikEnd=0, iNumRemainingColumnsToMove, icount;
	int imax = 0, imin = 0;
	int iAlignmentLength, iRegion, iLastRegion=0, iLastRegionCheck;
	int iStartPos;
	int iEndPos;

	bool bSuccess;
	LR LR; // defined for compatibility only

	unsigned int ***edist, ***trace, **uiPPEdist, **uiPPTrace;
	PatternRegion *patternRegion;

	im = pattern->iLength;

	imax = im + pattern->uiIndels;
	imin = im > pattern->uiIndels ? im - pattern->uiIndels : 1;

	int iNumRegions = patternStrandDirection->patternStructures->iNumRegions;

	char ***operation, **cPPOperation;
	char *cu_str, *seqs, *cSeqOperations, *cArcOperations;

	cu_str = (char *) patternStrandDirection->structure;

	edist     = uiNew3DMatrixX(im, imax);
	operation = cNew3DMatrixX(im, imax);
	trace     = uiNew3DMatrixX(im, imax);

	iStartPos = threadVars->iStartPos;
	iEndPos   = threadVars->iEndPos - imin + 1;
	iAlignmentLength = imax;

	cSeqOperations = (char *) calloc(imax + 1, sizeof(char));
	cArcOperations = (char *) calloc(imax + 1, sizeof(char));

	initializeRowsAndColumnsRegions_S(edist, operation, trace, pattern->cost, cu_str, im, imax, patternStrandDirection->patternStructures);

	uiPPEdist    = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	uiPPTrace    = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	cPPOperation = (char **) calloc(imax + 1, sizeof(char *));

	ComputedEntries *compEntries = newComputedEntriesMatrix(patternStrandDirection, false);

	int *iRoundComputedPerRow     = (int *) calloc(imax + 1, sizeof(int));
	int *iRoundComputedPerRowInit = (int *) calloc(imax + 1, sizeof(int));
	int *iEmptyArray = (int *) calloc(imax + 1, sizeof(int));

	PatternRegion *patternRegionZero = patternStrandDirection->patternStructures->patternRegion;

	int iRemoveControlIndex = 0;
	for (ii = 0; ii < FILTERFREQ; ii++){
		if(threadVars->iFilterOverlapsControl[ii] > iStartPos){
			iRemoveControlIndex = ii;
			break;
		}
	}

	bool bFirstWindow = true;

	for (ii = iStartPos; ii <= iEndPos; ii++) {

		if (bFirstWindow) {

			seqs             = affixArray->multiSeq->convSequences + ii;
			iAlignmentLength = imax;
			icount           = 0;

			while (*(seqs + ++icount) != $ && icount < iAlignmentLength);

			if (icount < iAlignmentLength) {
				if (icount >= imin) {
					iAlignmentLength = icount;
				} else {
					ii += icount;
					continue; //continue with next sequence after $
				}
			}

			for (iRegion = 0; iRegion < iNumRegions; iRegion++) {
				patternRegion = &patternStrandDirection->patternStructures->patternRegion[iRegion];
				if (patternRegion->bIsUnpaired) {
					for (ij = patternRegion->iL; ij <= patternRegion->iR; ij++) {
						memcpy(compEntries[iRegion].ijLastPerRowPerMatrix[ij], iEmptyArray, imax * sizeof(int));
					}
				} else {
					memcpy(compEntries[iRegion].ijLast, iEmptyArray, imax * sizeof(int));
				}
				if (patternRegion->bIsArc) {
					memcpy(compEntries[iRegion].ikFirst, iEmptyArray, imax * sizeof(int));
				}
			}

			memcpy(iRoundComputedPerRow, iRoundComputedPerRowInit, (imax + 1) * sizeof(int));

			bSuccess = false;
			iLastRegionCheck = iLastRegion = alignESA_LA_NEW_S(edist,
					operation,
					trace,
					patternStrandDirection,
					compEntries,
					iEmptyArray,
					seqs,
					iAlignmentLength,
					iAlignmentLength,
					affixArray,
					0, //iLcp2,
					0, //iDepthEnded2,
					0, //iLastRegion2
					iRoundComputedPerRow,
					&bSuccess,
					&LR,
					false);

			ikEnd = patternRegionZero[iLastRegionCheck].ikBigEnd;

			if (bSuccess){
				processAlignment_S_MT(edist,
						operation,
						trace,
						patternStrandDirection,
						ii,
						cSeqOperations,
						cArcOperations,
						imin,
						iAlignmentLength,
						affixArray,
						threadVars);
			}

			bFirstWindow = false;

		} else {

			seqs = affixArray->multiSeq->convSequences + ii;

			while (ii + iAlignmentLength - 1 > threadVars->iEndPos)
				iAlignmentLength--;
			if (*(seqs + iAlignmentLength - 1) == $)
				iAlignmentLength--;

			if (iAlignmentLength < imin && *(seqs + iAlignmentLength) == $) {
				bFirstWindow = true;
				ii += iAlignmentLength;
				continue;
			}

			if (patternRegionZero[iLastRegionCheck].ikBigEnd > ikEnd)
				ikEnd = patternRegionZero[iLastRegionCheck].ikBigEnd;

			ikStart = 1;
			ik = ikStart - 1;
			memcpy(uiPPEdist, edist[ik], sizeof(unsigned int **) * (imax - ik)); // Backup first matrices
			memcpy(uiPPTrace, trace[ik], sizeof(unsigned int **) * (imax - ik));
			memcpy(cPPOperation, operation[ik], sizeof(char **) * (imax - ik));

			ij = imax - ikStart;
			for (ik = ikStart; ik <= ikEnd; ik++) {
				memcpy(edist[ik - 1] + 1, edist[ik] + 1, sizeof(unsigned int **) * (imax - ik)); // Copy all columns
				memcpy(trace[ik - 1] + 1, trace[ik] + 1, sizeof(unsigned int **) * (imax - ik));
				memcpy(operation[ik - 1] + 1, operation[ik] + 1, sizeof(char **) * (imax - ik));

				edist[ik][imax - ik] = uiPPEdist[ij]; // Copy a column of backed-up matrix to the last column of another matrix
				trace[ik][imax - ik] = uiPPTrace[ij];
				operation[ik][imax - ik] = cPPOperation[ij];

				ij--;
				// End Move columns
				S(edist, ik, 0, imax - ik) = S(edist, ik, 0, imax - ik - 1) + pattern->cost->iCostDeletion;
			}
			iNumRemainingColumnsToMove = imax - ikEnd - 1;
			if (iNumRemainingColumnsToMove > 0) {
				memcpy(edist[ikEnd] + 1, uiPPEdist + 1, sizeof(unsigned int **) * iNumRemainingColumnsToMove);
				memcpy(trace[ikEnd] + 1, uiPPTrace + 1, sizeof(unsigned int **) * iNumRemainingColumnsToMove);
				memcpy(operation[ikEnd] + 1, cPPOperation + 1, sizeof(char **) * iNumRemainingColumnsToMove);
			}

			bSuccess = false;

			iLastRegion = alignESA_LA_NEW_LastColumnS1(edist,
					operation,
					trace,
					patternStrandDirection,
					compEntries,
					iEmptyArray,
					seqs,
					iAlignmentLength,
					iAlignmentLength,
					affixArray,
					&LR,
					ii,
					iRoundComputedPerRow,
					&bSuccess,
					false,
					iLastRegion,
					&iLastRegionCheck);

			if (bSuccess) {
				processAlignment_S_MT(edist,
						operation,
						trace,
						patternStrandDirection,
						ii, //uiSeqStartPos + iPos,
						cSeqOperations,
						cArcOperations,
						imin,
						iAlignmentLength,
						affixArray,
						threadVars);

				if (threadVars->searchParam->bFilterOverlaps && ii > threadVars->iFilterOverlapsControl[iRemoveControlIndex]) {
					iRemoveControlIndex++;

					pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
					deliverMatches(threadVars);
					if (threadVars->patternStrandDirection->bForwardStrand) {
						sortMatchesByPos(threadVars->threadsManager->matchesPerPatternFor[pattern->iId], 0,
								threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] - 1);
						threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternFor[pattern->iId],
								threadVars->threadsManager->uiNumMatchesPerPatternFor[pattern->iId]);
					} else {
						sortMatchesByPos(threadVars->threadsManager->matchesPerPatternRev[pattern->iId], 0,
								threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] - 1);
						threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId] = removeDuplicateMatches(threadVars->threadsManager->matchesPerPatternRev[pattern->iId],
								threadVars->threadsManager->uiNumMatchesPerPatternRev[pattern->iId]);
					}
					pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);
				}
			}

			ikEnd--;
		}

	}

	uiFree3DMatrixX(edist, im, imax);
	cFree3DMatrixX(operation, im, imax);
	uiFree3DMatrixX(trace, im, imax);

	free(cSeqOperations);
	free(cArcOperations);

	free(uiPPEdist);
	free(uiPPTrace);
	free(cPPOperation);

	free(iRoundComputedPerRow);
	free(iRoundComputedPerRowInit);

	free(iEmptyArray);

	freeComputedEntriesMatrix(compEntries, patternStrandDirection, false);

	pthread_mutex_lock(&threadVars->mutex->mutexMatchCount);
	deliverMatches(threadVars);
	pthread_mutex_unlock(&threadVars->mutex->mutexMatchCount);
	return EXIT_SUCCESS;
}


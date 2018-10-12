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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "align.h"

void applySequenceBasedFilter(AffixArray *affixArray, PatternStrandDirection *patternStrandDirection, bool *bVisitedSuffix, int iStartIdx, int iEndIdx) {
	int ii, iiend, ii2, ij, ijend, iIndels, im, imax, imin, iMinDist;
	int iLcpCheck, iLcp=0, iMinLcp, iSuffix, iReadingDepth;
	int iCount=0, iCountSuccess=0;
	unsigned int **edist;
	bool bBreak, bShortLength;
	char *cSeqs;

	char *cu = (char *) patternStrandDirection->seq;

	bool **iupacTable = affixArray->alphabet->iupacTable;

	Pattern *pattern = patternStrandDirection->pattern;

	int iCostDeletion    = pattern->cost->iFCostDeletion;
	int iCostReplacement = pattern->cost->iFCostReplacement;
	int iLength          = affixArray->length;
	int iThreshold       = pattern->uiThreshold;
	int *iSufArray       = affixArray->esa->xarray;
	unsigned char *ucLcpTable = affixArray->esa->xlcp;

	iIndels = pattern->uiIndels;
	im   = pattern->iLength;
	imax = im + pattern->uiIndels;
	imin = im > pattern->uiIndels ? im - pattern->uiIndels : 1;

	bShortLength = im > 255 ? false : true;

	edist = (unsigned int **) calloc(im + 1, sizeof(unsigned int*));
	for (ii = 0; ii <= im; ii++) {
		if ((edist[ii] = (unsigned int *) calloc(imax + 1, sizeof(unsigned int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"edist\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	for (ij = 1; ij <= imax; ij++) {
		edist[0][ij] = edist[0][ij - 1] + iCostDeletion;
	}
	for (ii = 1; ii <= im; ii++) {
		edist[ii][0] = edist[ii - 1][0] + iCostDeletion;
	}

	iMinLcp = 0;
	bBreak  = false;

	for (iSuffix = iStartIdx; iSuffix <= iEndIdx; iSuffix++) {

		cSeqs = (char *) affixArray->multiSeq->convSequences + iSufArray[iSuffix];
		iReadingDepth = 0;

		if ((bShortLength && ucLcpTable[iSuffix] > imax) || (iLcp = lcpvalue(iSuffix)) > imax)
			iReadingDepth = imax;
		else
			while (*(cSeqs + iReadingDepth) != $ && ++iReadingDepth < imax);

		if (iReadingDepth < imin) {
			bVisitedSuffix[iSufArray[iSuffix]] = true;
			iMinLcp = 0;
			bBreak = true;
			iCount++;
			continue;
		}

		if (bBreak)
			iLcp = iMinLcp;

		if (iLcp > 0 && iLcp < imax) {
			iiend = iLcp > im ? im : iLcp;

			if ((ijend = iiend + iIndels) > iReadingDepth)
				ijend = iReadingDepth;

			for (ii = 1; ii <= iiend; ii++) {
				for (ij = iiend + 1; ij <= ijend; ij++) {
					iMinDist = edist[ii - 1][ij - 1] + (sigma(cu[ii - 1], cSeqs[ij - 1]) * iCostReplacement);
					if (iMinDist >= edist[ii - 1][ij] + iCostDeletion) {
						iMinDist = edist[ii - 1][ij] + iCostDeletion;
					} if (iMinDist >= edist[ii][ij - 1] + iCostDeletion) {
						iMinDist = edist[ii][ij - 1] + iCostDeletion;
					}
					edist[ii][ij] = iMinDist;
				}
			}

			bBreak = true;

			ii = iiend;
			ij = iLcp > iIndels ? iLcp - iIndels : 1;
			for (; ij <= ijend; ij++) {
				if (edist[ii][ij] <= iThreshold) {
					bBreak = false;
					break;
				}
			}
			if (bBreak) {
				if ((iSuffix + 1) < iLength) {
					iMinLcp = iLcp;
					iLcp = bShortLength ? ucLcpTable[iSuffix + 1] : lcpvalue(iSuffix + 1);
					iMinLcp = iMinLcp < iLcp ? iMinLcp : iLcp;
				}

				bVisitedSuffix[iSufArray[iSuffix]] = true;
				iCount++;
				continue;
			}
		}

		for (ii = iLcp + 1; ii <= im; ii++) {
			if ((ijend = ii + iIndels) > iReadingDepth)
				ijend = iReadingDepth;

			for (ii2 = 1; ii2 < ii; ii2++) {
				iMinDist = edist[ii2 - 1][ijend - 1] + (sigma(cu[ii2 - 1], cSeqs[ijend - 1]) * iCostReplacement);
				if (iMinDist >= edist[ii2 - 1][ijend] + iCostDeletion) {
					iMinDist = edist[ii2 - 1][ijend] + iCostDeletion;
				} if (iMinDist >= edist[ii2][ijend - 1] + iCostDeletion) {
					iMinDist = edist[ii2][ijend - 1] + iCostDeletion;
				}
				edist[ii2][ijend] = iMinDist;
			}

			for (ij = 1; ij <= ijend; ij++) {
				iMinDist = edist[ii - 1][ij - 1] + (sigma(cu[ii - 1], cSeqs[ij - 1]) * iCostReplacement);
				if (iMinDist >= edist[ii - 1][ij] + iCostDeletion) {
					iMinDist = edist[ii - 1][ij] + iCostDeletion;
				} if (iMinDist >= edist[ii][ij - 1] + iCostDeletion) {
					iMinDist = edist[ii][ij - 1] + iCostDeletion;
				}
				edist[ii][ij] = iMinDist;
			}

			bBreak = true;

			ij = ii > iIndels ? ii - iIndels : 1;
			for (; ij <= ijend; ij++) {
				if (edist[ii][ij] <= iThreshold) {
					bBreak = false;
					break;
				}
			}

			if (bBreak) {
				bVisitedSuffix[iSufArray[iSuffix]] = true;
				iCount++;

				if (++iSuffix >= iLength)
					break;
				iLcpCheck = ii + iIndels > imax ? imax : ii + iIndels;
				iLcp = bShortLength ? ucLcpTable[iSuffix] : lcpvalue(iSuffix);

				iMinLcp = ii < iLcp ? ii : iLcp;

				while (iLcp >= iLcpCheck) {
					bVisitedSuffix[iSufArray[iSuffix]] = true;
					iCount++;

					if (++iSuffix >= iLength)
						break;

					if ((iLcp = bShortLength ? ucLcpTable[iSuffix] : lcpvalue(iSuffix)) < iMinLcp)
						iMinLcp = iLcp;
				}
				iSuffix--;

				break;
			}
		}

		if (!bBreak) {
			iCountSuccess++;
			if (++iSuffix < iLength) {
				iLcpCheck = imax;
				iLcp      = bShortLength ? ucLcpTable[iSuffix] : lcpvalue(iSuffix);

				while (iLcp >= iLcpCheck) {
					iCountSuccess++;

					if (++iSuffix >= iLength)
						break;

					iLcp = bShortLength ? ucLcpTable[iSuffix] : lcpvalue(iSuffix);
				}
			}
			iSuffix--;
		}

	}

	for (ii = 0; ii <= im; ii++) {
		free(edist[ii]);
	}

}

int alignESA(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int **ijLastOfRegion, char *cv,
		const int in, const int inMax, AffixArray *affixArray, const int iMaxColOffset, int iDepthEnded, bool bFindFailingRow, bool *bSuccess, LR *LR) {

	int ii, ij, ik, iend, im, itemp, iColOffset, imin=0;
	int iLineToMerge1, iLineToMerge2;
	int ij2, ij2tmp;

	char ctemp;

	int iRegion, iL, iR;

	Pattern *pattern = patternStrandDirection->pattern;
	bool **compRules = patternStrandDirection->compRules;

	int iCostDeletion    = pattern->cost->iCostDeletion;
	int iCostReplacement = pattern->cost->iCostReplacement;
	int iCostArcRemoving = pattern->cost->iCostArcRemoving;
	int iCostArcAltering = pattern->cost->iCostArcAltering;
	int iCostArcBreak    = pattern->cost->iCostArcBreak;

	int iLeftMinus1, ikstart, ikend, icol, icol2, icols, istart2;

	int iIndels       = pattern->uiIndels;
	int iThreshold    = pattern->uiThreshold;
	int *iLeftMostPos = patternStrandDirection->patternStructures->iLeftMostPos;
	bool bBreak;
	bool **iupacTable = affixArray->alphabet->iupacTable;

	PatternStructures *patternStructures = patternStrandDirection->patternStructures;
	PatternRegion *patternRegion;
	char *cu = (char *) patternStrandDirection->seq;

	int *iLineIndex = patternStructures->iLineIndex;
	im = pattern->iLength;

	int endx, iCol;

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
		patternRegion = &patternStructures->patternRegion[iRegion];
		iL = patternRegion->iL;
		iR = patternRegion->iR;

		if (patternRegion->bIsUnpaired) {

			for (ik = 0; ik < in; ik++) {
				iend = in - ik;
				iColOffset = ik < iMaxColOffset ? iMaxColOffset - ik : 0;
				for (ij = 1 + iColOffset; ij <= iend; ij++) {
					for (ii = iL; ii <= iR; ii++) {

						imin = P(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
						ctemp = OpReplacement;

						if (imin >= (itemp = P(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpDeletion;
						} if (imin >= (itemp = P(edist, ik, ii, ij - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpInsertion;
						}

						P(edist, ik, ii, ij)     = imin;
						P(operation, ik, ii, ij) = ctemp;
					}
				}
			}

		} else if (patternRegion->bIsArc) {

			iColOffset = in > iMaxColOffset ? in - iMaxColOffset : 0;
			for (ij = 1; ij <= iColOffset; ij++) {
				iend = in - ij;
				for (ik = iMaxColOffset; ik <= iend; ik++) {
					iCol = ij;

					imin = P(edist, ik, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iR - 1], cv[iCol + ik - 1]) * iCostReplacement + iCostArcAltering;
					ctemp = OpArcAltering1;

					if (imin > (itemp = P(edist, ik + 1, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
						imin = itemp;
						ctemp = OpArcAltering2;
					}

					if (imin > (itemp = P(edist, ik, iLineIndex[iR - 1], iCol) + iCostArcRemoving)) {
						imin = itemp;
						ctemp = OpArcRemoving;
					}

					if (imin > (itemp = P(edist, ik, iR, iCol - 1) + iCostDeletion)) {
						imin = itemp;
						ctemp = OpBInsertion1;
					}

					if (imin > (itemp = P(edist, ik + 1, iR, iCol - 1) + iCostDeletion)) {
						imin = itemp;
						ctemp = OpBInsertion2;
					}

					if (iCol > 1) {
						if (imin >= (itemp = P(edist, ik + 1, iLineIndex[iR - 1], iCol - 2) + (sigma(cu[iR - 1], cv[iCol + ik - 1])
								+ sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[iCol + ik - 1]) == true ? 0 : 1))) {
							imin = itemp;
							if (complement(cv[ik], cv[iCol + ik - 1]))
								ctemp = OpArcBreaking;
							else
								ctemp = OpArcBreakingTrue;
						}
					}

					P(edist, ik, iR, iCol)     = imin;
					P(operation, ik, iR, iCol) = ctemp;
				}

			}

			if (iMaxColOffset > 0) {

				endx = iColOffset + 1;

				for (ij = 2; ij <= endx; ij++) {
					iCol = ij;

					for (ik = iMaxColOffset - 1; ik >= 0; ik--) {
						imin = P(edist, ik, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iR - 1], cv[iCol + ik - 1]) * iCostReplacement + iCostArcAltering;
						ctemp = OpArcAltering1;

						if (imin > (itemp = P(edist, ik + 1, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
							imin = itemp;
							ctemp = OpArcAltering2;
						}

						if (imin > (itemp = P(edist, ik, iLineIndex[iR - 1], iCol) + iCostArcRemoving)) {
							imin = itemp;
							ctemp = OpArcRemoving;
						}

						if (imin > (itemp = P(edist, ik, iR, iCol - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpBInsertion1;
						}

						if (imin > (itemp = P(edist, ik + 1, iR, iCol - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpBInsertion2;
						}

						if (iCol > 1) {
							if (imin >= (itemp = P(edist, ik + 1, iLineIndex[iR - 1], iCol - 2) + (sigma(cu[iR - 1], cv[iCol + ik - 1])
									+ sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[iCol + ik - 1]) == true ? 0 : 1))) {
								imin = itemp;
								if (complement(cv[ik], cv[iCol + ik - 1]))
									ctemp = OpArcBreaking;
								else
									ctemp = OpArcBreakingTrue;
							}
						}

						P(edist, ik, iR, iCol)     = imin;
						P(operation, ik, iR, iCol) = ctemp;
						iCol++;
					}

				}
			}



		} else {

			iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
			iLineToMerge2 = iR;

			for (ik = 0; ik < in; ik++) {
				iend = in - ik;
				iColOffset = ik < iMaxColOffset ? iMaxColOffset - ik : 1;
				for (ij = iColOffset; ij <= iend; ij++) {
					imin   = INT_MAX;
					ij2tmp = 0;
					for (ij2 = 0; ij2 <= ij; ij2++) {
						if (imin > (itemp = P(edist, ik, iLineToMerge1, ij2) + P(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
							imin   = itemp;
							ij2tmp = ij2;
						}
					}

					P(edist, ik, iLineIndex[iLineToMerge2], ij)     = imin;
					P(trace, ik, iLineToMerge2, ij)     = ij2tmp == 0 ? INT_MAX : ij2tmp;
				}
			}

		}

	}

	// Find first failing row (alignment cost > threshold) if there is one
	if (bFindFailingRow) {
		for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {

			patternRegion = &patternStructures->patternRegion[iRegion];
			iL = patternRegion->iL;
			iR = patternRegion->iR;

			if (patternRegion->bIsUnpaired) {

				if (patternRegion->bIsLeftBulge)
					iLeftMinus1 = iL - 1;
				else
					iLeftMinus1 = iLeftMostPos[iL] - 1;

				ikstart = iLeftMinus1 > iIndels ? iLeftMinus1 - iIndels : 0;
				ikend   = iLeftMinus1 + iIndels < in ? iLeftMinus1 + iIndels : in;
				icol   = iR - iLeftMinus1; //mid col of last line

				for (ii = iL; ii <= iR; ii++) {
					icol2 = ii - iLeftMinus1;

					bBreak = true;
					for (ik = ikstart; ik <= ikend; ik++) {
						icols = iIndels - abs(iLeftMinus1 - ik); /*number of columns to be checked - 1*/

						if ((iend = icol2 + icols) > in - ik)
							iend = in - ik;
						istart2 = icol2 > icols ? icol2 - icols : 1;

						for (ij = istart2; ij <= iend && bBreak; ij++) {
							if (P(edist, ik, ii, ij) <= iThreshold)
								bBreak = false;
						}

					}
					if (bBreak) {
						LR->iL = iLeftMostPos[iL];
						LR->iR = ii;
						return iRegion;
					}
				}

			} else if (patternRegion->bIsArc) {

				iLeftMinus1 = iL - 1;

				ikstart = iLeftMinus1 > iIndels ? iLeftMinus1 - iIndels : 0;
				ikend   = iLeftMinus1 + iIndels < in ? iLeftMinus1 + iIndels : in;
				icol   = iR - iLeftMinus1; //mid col of last line

				bBreak = true;
				for (ik = ikstart; ik <= ikend; ik++) {
					icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

					if ((iend = icol + icols) > in - ik)
						iend = in - ik;
					istart2 = icol > icols ? icol - icols : 1;

					for (ij = istart2; ij <= iend && bBreak; ij++) {
						if (P(edist, ik, iR, ij) <= iThreshold)
							bBreak = false;
					}
				}
				if (bBreak) {
					LR->iL = iL;
					LR->iR = iR;
					return iRegion;
				}

			} else {

				iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
				iLineToMerge2 = iR;

				iLeftMinus1 = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0] - 1;
				if (iLineToMerge1 == 0)
					iLeftMinus1++;
				ikstart = iLeftMinus1 > iIndels ? iLeftMinus1 - iIndels : 0;
				ikend   = iLeftMinus1 + iIndels < in ? iLeftMinus1 + iIndels : in - 1;
				icol = iLineToMerge2 - iLeftMinus1; //mid col for line

				bBreak = true;
				for (ik = ikstart; ik <= ikend; ik++) {
					icols = iIndels - abs(iLeftMinus1 - ik); // number of columns to be checked - 1

					if ((iend = icol + icols) > in - ik)
						iend = in - ik; // end of interval to be computed
					istart2 = icol > icols ? icol - icols : 1; // start of interval to be computed

					for (ij = istart2; ij <= iend && bBreak; ij++) {
						if (P(edist, ik, iLineIndex[iLineToMerge2], ij) <= iThreshold)
							bBreak = false;
					}
				}
				if (bBreak) {
					LR->iL = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0];
					LR->iR = iLineToMerge2;
					return iRegion;
				}

			}

		}
	}

	LR->iL = 1;
	LR->iR = im;
	*bSuccess = true;
	return patternStructures->iNumRegions - 1;
}

inline unsigned int xlcpExceptionValue(LcpException *xlcpException, const unsigned int k) {
	unsigned int low, high, mid;

	/* naive search */
	/*unsigned int i;
	for(i = 0; i < xlcpException->numExceptions; i++)
		if(xlcpException->index[i] == k)
			return xlcpException->value[i];*/

	low = 0;
	high = xlcpException->numExceptions - 1;

	while (low <= high) {
		mid = (low + high) >> 1;
		if (xlcpException->index[mid] < k)
			low = mid + 1;
		else if (xlcpException->index[mid] > k)
			high = mid - 1;
		else
			return xlcpException->value[mid];
	}

	return -1;
}

void getEditOperations(char *cSeqOperations, char *cArcOperations, char ***operation, unsigned int ***trace, int *iIndels, int *iDeletions, const int iMaxIndels,
		PatternStructures *pstr, unsigned int *uiIndex, unsigned int uik, unsigned int uii, unsigned int uij, bool bAllowBranch) {
	char cOp;
	unsigned int traceValue;

	if (*iIndels > iMaxIndels)
		return;

	if (P(trace, uik, uii, uij) != 0 && bAllowBranch) {
		cOp = 99;
	} else {
		cOp = P(operation, uik, uii, uij);
	}

	switch (cOp) {
	case OpReplacement:
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpInsertion:
		(*iIndels)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii, uij - 1, bAllowBranch);
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpDeletion:
		(*iIndels)++;
		(*iDeletions)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpArcBreaking:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '(';
		(*uiIndex)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii - 1, uij - 2, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = ')';
		(*uiIndex)++;
		break;
	case OpArcBreakingTrue:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii - 1, uij - 2, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpArcAltering1:
		(*iIndels)++;
		(*iDeletions)++;
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpArcAltering2:
		(*iIndels)++;
		(*iDeletions)++;
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpArcRemoving:
		(*iIndels) += 2;
		(*iDeletions) += 2;
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpBInsertion1:
		(*iIndels)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii, uij - 1, bAllowBranch);
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpBInsertion2:
		(*iIndels)++;
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii, uij - 1, bAllowBranch);
		break;
	case 99:
		if (P(trace, uik, uii, uij) == INT_MAX)
			traceValue = 0;
		else
			traceValue = P(trace, uik, uii, uij);
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, pstr->patternRegion[pstr->iPosRegion[uii]].iL - 1, traceValue, true);
		getEditOperations(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + traceValue, uii, uij - traceValue, false);
		break;
	}
}

int getSeqNumber(unsigned int *range, int length, unsigned int k) {
	int low = 0;
	int high = length - 1;
	int mid;
	int vMidLeft;
	int vMidRight;

	while (low <= high) {
		mid = (low + high) >> 1;
		vMidLeft = (mid == 0) ? 0 : (range[mid-1] + 1);
		vMidRight = range[mid];

		if (vMidRight < k)
			low = mid + 1;
		else if (vMidLeft > k)
			high = mid - 1;
		else
			return mid;
	}
	return -1;
}

/*optimal look-ahead*/
int alignESA_LA_NEW(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int *ijLastOfRegion,
		char *cv, const int in, const int inMax, AffixArray *affixArray, int iMaxColOffset, int iLastCompRow, bool bFindFailingRow, bool *bSuccess, LR *LR) {

	int ii, ii2, ij, ik, ik2, ikX, iend, im, itemp, iColOffset, imin=0;
	int iLineToMerge1, iLineToMerge2;
	int ij2, ij2tmp;
	int iLeftMinus1, icol2, istart2, icol, icols;
	int iRegion, iL, iR;
	int icount, ik2end;

	Pattern *pattern = patternStrandDirection->pattern;
	bool **compRules = patternStrandDirection->compRules;

	int iCostDeletion    = pattern->cost->iCostDeletion;
	int iCostReplacement = pattern->cost->iCostReplacement;
	int iCostArcRemoving = pattern->cost->iCostArcRemoving;
	int iCostArcAltering = pattern->cost->iCostArcAltering;
	int iCostArcBreak    = pattern->cost->iCostArcBreak;

	bool **iupacTable = affixArray->alphabet->iupacTable;
	char ctemp;

	int iIndels = pattern->uiIndels;
	int iThreshold = pattern->uiThreshold;

	PatternRegion *patternRegion, *currentPatternRegion;

	PatternStructures *patternStructures = patternStrandDirection->patternStructures;
	int *iLeftMostPos = patternStructures->iLeftMostPos;
	char *cu = (char *) patternStrandDirection->seq;

	int *iLineIndex = patternStructures->iLineIndex;
	im = pattern->iLength;

	int *ijLast, **ijLastPerRowPerMatrix, *ijGoal, **ijGoalPerRowPerMatrix, *ikFirst;
	int ikStart, ikEnd, iDepRegion, ijGoal_local, iDepRegionCheck = 0;
	bool bDepRegion, bBreak;

	int iLastRegion = patternStructures->iPosRegion[iLastCompRow];

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {

		currentPatternRegion = &patternStructures->patternRegion[iRegion];

		if (iLastRegion < iRegion || iMaxColOffset == 0) {
			if (currentPatternRegion->bIsUnpaired) {
				for (ii = currentPatternRegion->iL; ii <= currentPatternRegion->iR; ii++) {
					memcpy(compEntries[iRegion].ijLastPerRowPerMatrix[ii], ijLastOfRegion, inMax * sizeof(int));
				}
			} else {
				memcpy(compEntries[iRegion].ijLast, ijLastOfRegion, inMax * sizeof(int));
				if (currentPatternRegion->bIsArc) {
					memcpy(compEntries[iRegion].ikFirst, ijLastOfRegion, inMax * sizeof(int));
				}
			}
		} else {
			if (currentPatternRegion->bIsUnpaired) {
				ii = currentPatternRegion->iL;
				if (iLastRegion == iRegion) {
					for (; ii <= iLastCompRow; ii++) {
						for (ik = 0; ik < inMax; ik++) {
							iColOffset = ik < iMaxColOffset ? iMaxColOffset - ik : 0;
							if (compEntries[iRegion].ijLastPerRowPerMatrix[ii][ik] > iColOffset)
								compEntries[iRegion].ijLastPerRowPerMatrix[ii][ik] = iColOffset;
						}
					}
				}
				for (; ii <= currentPatternRegion->iR; ii++) {
					for (ik = 0; ik < inMax; ik++) {
						iColOffset = ik < iMaxColOffset ? iMaxColOffset - ik : 0;
						if (compEntries[iRegion].ijLastPerRowPerMatrix[ii][ik] > iColOffset)
							compEntries[iRegion].ijLastPerRowPerMatrix[ii][ik] = iColOffset;
					}
				}

			} else {
				for (ik = 0; ik < inMax; ik++) {
					iColOffset = ik < iMaxColOffset ? iMaxColOffset - ik : 0;
					if (compEntries[iRegion].ijLast[ik] > iColOffset)
						compEntries[iRegion].ijLast[ik] = iColOffset;
				}
			}
		}

		do {

			if ((iDepRegion = currentPatternRegion->iRegionDependsOnRegions[iDepRegionCheck]) != -1) {
				patternRegion = &patternStructures->patternRegion[iDepRegion];
				bDepRegion = true;
				iDepRegionCheck++;

				ijLast                = compEntries[iDepRegion].ijLast;
				ijLastPerRowPerMatrix = compEntries[iDepRegion].ijLastPerRowPerMatrix;
				ikFirst               = compEntries[iDepRegion].ikFirst;
			} else {
				patternRegion = currentPatternRegion;
				bDepRegion = false;
				iDepRegionCheck = 0;

				ijLast                = compEntries[iRegion].ijLast;
				ijLastPerRowPerMatrix = compEntries[iRegion].ijLastPerRowPerMatrix;
				ikFirst               = compEntries[iRegion].ikFirst;
			}

			iL      = patternRegion->iL;
			iR      = patternRegion->iR;
			ijGoal  = patternRegion->ijGoal;
			ikStart = patternRegion->ikStart;
			ikEnd   = patternRegion->ikEnd;

			if (patternRegion->bIsUnpaired) {

				if (bDepRegion) {
					ikStart = currentPatternRegion->ikStart;
					if ((ikEnd = currentPatternRegion->iR + iIndels - 1) >= in)
						ikEnd = in - 1;

					ijGoal = currentPatternRegion->ijGoal;

					for (ii = iL; ii <= iR; ii++) {
						for (ik = ikStart; ik <= ikEnd; ik++) {

							if ((ijGoal_local = ijGoal[ik]) > in - ik) {
								ijGoal_local = in - ik;
							}

							for (ij = 1 + ijLastPerRowPerMatrix[ii][ik]; ij <= ijGoal_local; ij++) {
								imin = P(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
								ctemp = OpReplacement;

								if (imin >= (itemp = P(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpDeletion;
								} if (imin >= (itemp = P(edist, ik, ii, ij - 1) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpInsertion;
								}
								P(edist, ik, ii, ij)     = imin;
								P(operation, ik, ii, ij) = ctemp;
							}
							ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
						}
					}

				} else {
					iLeftMinus1 = iLeftMostPos[iL] - 1;

					ijGoalPerRowPerMatrix = patternRegion->ijGoalPerRowPerMatrix;

					for (ii = iL; ii <= iR; ii++) {

						for (ik = ikStart; ik <= ikEnd; ik++) {

							if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik]) > in - ik) {
								ijGoal_local = in - ik;
							}

							ij = 1 + ijLastPerRowPerMatrix[ii][ik]
							                                   ;
							if (ii > iL && ijGoal_local >= ij) {
								ij2 = ijGoal_local;
								for (ii2 = iL; ii2 < ii; ii2++) {
									imin = P(edist, ik, iLineIndex[ii2 - 1], ij2 - 1) + (sigma(cu[ii2 - 1], cv[ij2 + ik - 1]) * iCostReplacement);
									ctemp = OpReplacement;

									if (imin >= (itemp = P(edist, ik, iLineIndex[ii2 - 1], ij2) + iCostDeletion)) {
										imin = itemp;
										ctemp = OpDeletion;
									} if (imin >= (itemp = P(edist, ik, ii2, ij2 - 1) + iCostDeletion)) {
										imin = itemp;
										ctemp = OpInsertion;
									}

									P(edist, ik, ii2, ij2)     = imin;
									P(operation, ik, ii2, ij2) = ctemp;
									ijLastPerRowPerMatrix[ii2][ik] = ijGoal_local;
								}
							}

							for (; ij <= ijGoal_local; ij++) {
								imin = P(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
								ctemp = OpReplacement;

								if (imin >= (itemp = P(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpDeletion;
								} if (imin >= (itemp = P(edist, ik, ii, ij - 1) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpInsertion;
								}
								P(edist, ik, ii, ij)     = imin;
								P(operation, ik, ii, ij) = ctemp;
								ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
							}

						}

						icol2 = ii - iLeftMinus1;
						bBreak = true;

						for (ik2 = ikStart; ik2 <= ikEnd; ik2++) {
							if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik2]) > in - ik2) {
								ijGoal_local = in - ik2;
							}

							icols = iIndels - abs(iLeftMinus1 - ik2);
							istart2 = icol2 > icols ? icol2 - icols : 0;

							for (ij = istart2; ij <= ijGoal_local; ij++) {
								if (P(edist, ik2, ii, ij) <= iThreshold) {
									bBreak = false;
									break;
								}
							}
							if (!bBreak) break;
						}
						if (bBreak) {
							LR->iL = iLeftMostPos[iL];
							LR->iR = ii;
							return iRegion;
						}
					}
				}

			} else if (patternRegion->bIsArc) {
				if (bDepRegion) {
					ik2end = currentPatternRegion->ikEnd;
					ikStart = currentPatternRegion->ikStart;
					ijGoal = currentPatternRegion->ijGoal;
				} else {
					ik2end = iL - 1;
				}

				for (ik2 = ikStart; ik2 <= ik2end; ik2++) {

					if ((ijGoal_local = ijGoal[ik2]) > in - ik2) {
						ijGoal_local = in - ik2;
					}

					iend = ijGoal_local - ijLast[ik2]; //number of steps to add

					for (icount = 0; icount < iend; icount++) {

						ikEnd = ijLast[ik2] + ik2;

						if (ijLast[ikEnd] > 0) {
							ikX = ikEnd;
							do {
								if ((ikX = ikFirst[ikX] - 1) < ikStart) {
									ikX = ikStart;
									break;
								}
							} while (ijLast[ikX] >= (ikEnd - ikX + 1));
							ij = ikEnd - ikX + 1;
							ikEnd = ikX;
						} else {
							ij = 1;
						}

						for (ik = ikEnd; ik >= ik2; ik--) {

							imin = P(edist, ik, iLineIndex[iR - 1], ij - 1) + sigma(cu[iR - 1], cv[ij + ik - 1]) * iCostReplacement + iCostArcAltering;
							ctemp = OpArcAltering1;

							if (imin > (itemp = P(edist, ik + 1, iLineIndex[iR - 1], ij - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
								imin = itemp;
								ctemp = OpArcAltering2;
							}

							if (imin > (itemp = P(edist, ik, iLineIndex[iR - 1], ij) + iCostArcRemoving)) {
								imin = itemp;
								ctemp = OpArcRemoving;
							}

							if (imin > (itemp = P(edist, ik, iR, ij - 1) + iCostDeletion)) {
								imin = itemp;
								ctemp = OpBInsertion1;
							}

							if (imin > (itemp = P(edist, ik + 1, iR, ij - 1) + iCostDeletion)) {
								imin = itemp;
								ctemp = OpBInsertion2;
							}

							if (ij > 1) {
								if (imin >= (itemp = P(edist, ik + 1, iLineIndex[iR - 1], ij - 2) + (sigma(cu[iR - 1], cv[ij + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[ij + ik - 1]) == true ? 0 : 1))) {
									imin = itemp;
									if (complement(cv[ik], cv[ij + ik - 1]))
										ctemp = OpArcBreaking;
									else
										ctemp = OpArcBreakingTrue;
								}
							}

							P(edist, ik, iR, ij)     = imin;
							P(operation, ik, iR, ij) = ctemp;

							ikFirst[ik] = ik2;
							ijLast[ik] = ij;

							ij++;
						}
					}

				}

				if (!bDepRegion) {
					iLeftMinus1 = iL - 1;
					icol   = iR - iLeftMinus1; //mid col of last line

					bBreak = true;
					for (ik = ikStart; ik <= ikEnd; ik++) {
						icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

						if ((iend = icol + icols) > in - ik)
							iend = in - ik;
						istart2 = icol > icols ? icol - icols : 0;

						for (ij = istart2; ij <= iend; ij++) {
							if (P(edist, ik, iR, ij) <= iThreshold) {
								bBreak = false;
								break;
							}
						}
						if (!bBreak) break;
					}
					if (bBreak) {
						LR->iL = iL;
						LR->iR = iR;
						return iRegion;
					}
				}

			} else {

				iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
				iLineToMerge2 = iR;

				if (bDepRegion) {
					ikStart = currentPatternRegion->ikStart;

					ijGoal = currentPatternRegion->ijGoal;

					if ((ikEnd = currentPatternRegion->ikEndFactual) >= in)
						ikEnd = in - 1;
				}

				for (ik = ikStart; ik <= ikEnd; ik++) {
					if ((ijGoal_local = ijGoal[ik]) > in - ik)
						ijGoal_local = in - ik;

					for (ij = 1 + ijLast[ik]; ij <= ijGoal_local; ij++) {
						imin   = INT_MAX;
						ij2tmp = 0;
						for (ij2 = 0; ij2 <= ij; ij2++) {
							if (imin > (itemp = P(edist, ik, iLineToMerge1, ij2) + P(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
								imin   = itemp;
								ij2tmp = ij2;
							}
						}
						P(edist, ik, iLineIndex[iLineToMerge2], ij) = imin;
						P(trace, ik, iLineToMerge2, ij) = ij2tmp == 0 ? INT_MAX : ij2tmp;
					}
					ijLast[ik] = ijGoal_local;
				}

				if (bDepRegion) {
					bBreak = false;
				} else {
					iLeftMinus1 = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0] - 1;
					if (iLineToMerge1 == 0)
						iLeftMinus1++;
					icol = iLineToMerge2 - iLeftMinus1; //mid col for line

					bBreak = true;

					for (ik = ikStart; ik <= ikEnd; ik++) {
						if ((ijGoal_local = ijGoal[ik]) > in - ik) {
							ijGoal_local = in - ik;
						}

						icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1
						istart2 = icol > icols ? icol - icols : 0;

						for (ij = istart2; ij <= ijGoal_local; ij++) {
							if (P(edist, ik, iLineIndex[iLineToMerge2], ij) <= iThreshold) {
								bBreak = false;
								break;
							}
						}
						if (!bBreak) break;
					}
				}

				if (bBreak) {
					LR->iL = iLeftMostPos[iLineToMerge1];
					LR->iR = iLineToMerge2;
					return iRegion;
				}

			}

		} while (bDepRegion);

	}

	LR->iL = 1;
	LR->iR = im;
	*bSuccess = true;
	return patternStructures->iNumRegions - 1;
}

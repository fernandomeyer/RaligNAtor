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

void getEditOperations_S(char *cSeqOperations, char *cArcOperations, char ***operation, unsigned int ***trace, int *iIndels, int *iDeletions, const int iMaxIndels,
		PatternStructures *pstr, unsigned int *uiIndex, unsigned int uik, unsigned int uii, unsigned int uij, bool bAllowBranch) {
	char cOp;
	unsigned int traceValue;

	if (*iIndels > iMaxIndels)
		return;

	if (S(trace, uik, uii, uij) != 0 && bAllowBranch) {
		cOp = 99;
	} else {
		cOp = S(operation, uik, uii, uij);
	}

	switch (cOp) {
	case OpReplacement:
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpInsertion:
		(*iIndels)++;
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii, uij - 1, bAllowBranch);
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpDeletion:
		(*iIndels)++;
		(*iDeletions)++;
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpArcBreaking:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '(';
		(*uiIndex)++;
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii - 1, uij - 2, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = ')';
		(*uiIndex)++;
		break;
	case OpArcBreakingTrue:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii - 1, uij - 2, true);
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
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij - 1, true);
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
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii - 1, uij - 1, true);
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
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii - 1, uij, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpBInsertion1:
		(*iIndels)++;
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, uii, uij - 1, bAllowBranch);
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpBInsertion2:
		(*iIndels)++;
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + 1, uii, uij - 1, bAllowBranch);
		break;
	case 99:
		if (S(trace, uik, uii, uij) == INT_MAX)
			traceValue = 0;
		else
			traceValue = S(trace, uik, uii, uij);
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik, pstr->patternRegion[pstr->iPosRegion[uii]].iL - 1, traceValue, true);
		getEditOperations_S(cSeqOperations, cArcOperations, operation, trace, iIndels, iDeletions, iMaxIndels, pstr, uiIndex, uik + traceValue, uii, uij - traceValue, false);
		break;
	}
}



int alignESA_S(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, char *cv,
		const int in, AffixArray *affixArray, const int iMaxColOffset) {

	int ii, ij, ik, iend, itemp, iColOffset, imin=0;
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

	bool **iupacTable = affixArray->alphabet->iupacTable;

	PatternStructures *patternStructures = patternStrandDirection->patternStructures;

	PatternRegion *patternRegion;
	char *cu = (char *) patternStrandDirection->seq;

	int *iLineIndex = patternStructures->iLineIndex;

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
						imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
						ctemp = OpReplacement;

						if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpDeletion;
						} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpInsertion;
						}

						S(edist, ik, ii, ij)     = imin;
						S(operation, ik, ii, ij) = ctemp;
					}
				}
			}

		} else if (patternRegion->bIsArc) {

			iColOffset = in > iMaxColOffset ? in - iMaxColOffset : 0;
			for (ij = 1; ij <= iColOffset; ij++) {
				iend = in - ij;
				for (ik = iMaxColOffset; ik <= iend; ik++) {

					iCol = ij;

					imin = S(edist, ik, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iR - 1], cv[iCol + ik - 1]) * iCostReplacement + iCostArcAltering;
					ctemp = OpArcAltering1;

					if (imin > (itemp = S(edist, ik + 1, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
						imin = itemp;
						ctemp = OpArcAltering2;
					}

					if (imin > (itemp = S(edist, ik, iLineIndex[iR - 1], iCol) + iCostArcRemoving)) {
						imin = itemp;
						ctemp = OpArcRemoving;
					}

					if (imin > (itemp = S(edist, ik, iR, iCol - 1) + iCostDeletion)) {
						imin = itemp;
						ctemp = OpBInsertion1;
					}

					if (imin > (itemp = S(edist, ik + 1, iR, iCol - 1) + iCostDeletion)) {
						imin = itemp;
						ctemp = OpBInsertion2;
					}

					if (iCol > 1) {
						if (imin >= (itemp = S(edist, ik + 1, iLineIndex[iR - 1], iCol - 2) + (sigma(cu[iR - 1], cv[iCol + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[iCol + ik - 1]) == true ? 0 : 1))) {
							imin = itemp;
							if (complement(cv[ik], cv[iCol + ik - 1]))
								ctemp = OpArcBreaking;
							else
								ctemp = OpArcBreakingTrue;
						}
					}

					S(edist, ik, iR, iCol)     = imin;
					S(operation, ik, iR, iCol) = ctemp;
				}

			}

			if (iMaxColOffset > 0) {

				endx = iColOffset + 1;

				for (ij = 2; ij <= endx; ij++) {
					iCol = ij;

					for (ik = iMaxColOffset - 1; ik >= 0; ik--) {

						imin = S(edist, ik, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iR - 1], cv[iCol + ik - 1]) * iCostReplacement + iCostArcAltering;
						ctemp = OpArcAltering1;

						if (imin > (itemp = S(edist, ik + 1, iLineIndex[iR - 1], iCol - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
							imin = itemp;
							ctemp = OpArcAltering2;
						}

						if (imin > (itemp = S(edist, ik, iLineIndex[iR - 1], iCol) + iCostArcRemoving)) {
							imin = itemp;
							ctemp = OpArcRemoving;
						}

						if (imin > (itemp = S(edist, ik, iR, iCol - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpBInsertion1;
						}

						if (imin > (itemp = S(edist, ik + 1, iR, iCol - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpBInsertion2;
						}

						if (iCol > 1) {
							if (imin >= (itemp = S(edist, ik + 1, iLineIndex[iR - 1], iCol - 2) + (sigma(cu[iR - 1], cv[iCol + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[iCol + ik - 1]) == true ? 0 : 1))) {
								imin = itemp;
								if (complement(cv[ik], cv[iCol + ik - 1]))
									ctemp = OpArcBreaking;
								else
									ctemp = OpArcBreakingTrue;
							}
						}

						S(edist, ik, iR, iCol)     = imin;
						S(operation, ik, iR, iCol) = ctemp;
						iCol++;
					}

				}
			}

		} else {

			iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
			iLineToMerge2 = iR;

			for (ik = 0; ik < in; ik++) {
				iend = in - ik;
				iColOffset = ik < iMaxColOffset ? iMaxColOffset - ik : 1; // TODO: 0 or 1?
				for (ij = iColOffset; ij <= iend; ij++) {
					imin   = INT_MAX;
					ij2tmp = 0;
					for (ij2 = 0; ij2 <= ij; ij2++) {
						if (imin > (itemp = S(edist, ik, iLineToMerge1, ij2) + S(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
							imin   = itemp;
							ij2tmp = ij2;
						}
					}

					S(edist, ik, iLineIndex[iLineToMerge2], ij) = imin;
					S(trace, ik, iLineToMerge2, ij)             = ij2tmp == 0 ? INT_MAX : ij2tmp;
				}
			}

		}

	}

	return patternStructures->iNumRegions - 1;
}

int alignLastColumn_S(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection,
		char *cv, const int in, AffixArray *affixArray) {

	int ii, ij, ik, itemp, imin=0;
	int iLineToMerge1, iLineToMerge2;
	int ij2, ij2tmp;

	Pattern *pattern = patternStrandDirection->pattern;
	bool **compRules = patternStrandDirection->compRules;

	int iCostDeletion    = pattern->cost->iCostDeletion;
	int iCostReplacement = pattern->cost->iCostReplacement;
	int iCostArcRemoving = pattern->cost->iCostArcRemoving;
	int iCostArcAltering = pattern->cost->iCostArcAltering;
	int iCostArcBreak    = pattern->cost->iCostArcBreak;

	bool **iupacTable = affixArray->alphabet->iupacTable;
	char ctemp;

	int iRegion, iL, iR;
	PatternRegion *patternRegion;

	PatternStructures *patternStructures = patternStrandDirection->patternStructures;
	char *cu = (char *) patternStrandDirection->seq;

	int *iLineIndex = patternStructures->iLineIndex;

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
		patternRegion = &patternStructures->patternRegion[iRegion];
		iL = patternRegion->iL;
		iR = patternRegion->iR;

		if (patternRegion->bIsUnpaired) {

			for (ik = 0; ik < in; ik++) {
				ij = in - ik;
				for (ii = iL; ii <= iR; ii++) {
					imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
					ctemp = OpReplacement;

					if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
						imin = itemp;
						ctemp = OpDeletion;
					} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
						imin = itemp;
						ctemp = OpInsertion;
					}

					S(edist, ik, ii, ij)     = imin;
					S(operation, ik, ii, ij) = ctemp;
				}
			}


		} else if (patternRegion->bIsArc) {

			ij = 1;
			for (ik = in - 1; ik >= 0; ik--) {

				imin = S(edist, ik, iLineIndex[iR - 1], ij - 1) + sigma(cu[iR - 1], cv[ij + ik - 1]) * iCostReplacement + iCostArcAltering;
				ctemp = OpArcAltering1;

				if (imin > (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
					imin = itemp;
					ctemp = OpArcAltering2;
				}

				if (imin > (itemp = S(edist, ik, iLineIndex[iR - 1], ij) + iCostArcRemoving)) {
					imin = itemp;
					ctemp = OpArcRemoving;
				}

				if (imin > (itemp = S(edist, ik, iR, ij - 1) + iCostDeletion)) {
					imin = itemp;
					ctemp = OpBInsertion1;
				}

				if (imin > (itemp = S(edist, ik + 1, iR, ij - 1) + iCostDeletion)) {
					imin = itemp;
					ctemp = OpBInsertion2;
				}

				if (ij > 1) {
					if (imin >= (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 2) + (sigma(cu[iR - 1], cv[ij + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[ij + ik - 1]) == true ? 0 : 1))) {
						imin = itemp;
						if (complement(cv[ik], cv[ij + ik - 1]))
							ctemp = OpArcBreaking;
						else
							ctemp = OpArcBreakingTrue;
					}
				}

				S(edist, ik, iR, ij)     = imin;
				S(operation, ik, iR, ij) = ctemp;

				ij++;
			}

		} else {

			iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
			iLineToMerge2 = iR;

			for (ik = 0; ik < in; ik++) {
				ij = in - ik;

				imin   = INT_MAX;
				ij2tmp = 0;
				for (ij2 = 0; ij2 <= ij; ij2++) {
					if (imin > (itemp = S(edist, ik, iLineToMerge1, ij2) + S(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
						imin   = itemp;
						ij2tmp = ij2;
					}
				}

				S(edist, ik, iLineIndex[iLineToMerge2], ij)     = imin;
				S(trace, ik, iLineToMerge2, ij)     = ij2tmp == 0 ? INT_MAX : ij2tmp;
			}

		}

	}

	return patternStructures->iNumRegions - 1;
}

int alignESA_LA_NEW_LastColumnS1(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int *iEmptyArray,
		char *cv, const int in, const int inMax, AffixArray *affixArray, LR *LR, int iCurrentRound, int *iRoundComputedPerRow, bool *bSuccess, bool bInOutAlignment, int iLastRegion, int *iLastRegionCheck) {

	int ii, ii2, ij, ik, ik2, ikX, iend, im, itemp, imin=0, iRow;
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

	PatternRegion *patternRegion, *currentPatternRegion, *patternRegionZero;

	PatternStructures *patternStructures = patternStrandDirection->patternStructures;
	int *iLeftMostPos = patternStructures->iLeftMostPos;
	char *cu = (char *) patternStrandDirection->seq;

	int *iLineIndex = patternStructures->iLineIndex;
	im = pattern->iLength;

	int *ijLast, **ijLastPerRowPerMatrix, *ijGoal, *ijGoal_target, **ijGoalPerRowPerMatrix, *ikFirst;
	int ikStart, ikStart_target, ikEnd, ikEnd_target, iDepRegion, ijGoal_local, iDepRegionCheck = 0, iRounds;
	int **iLargestRegionIndex, iTargetRegionIndex;
	bool bDepRegion, bBreak;

	if (bInOutAlignment) {
		patternRegionZero = patternStructures->patternRegionInsideOut;
		iLargestRegionIndex = patternStructures->iLargestRegionIndexInsideOut;
	} else {
		patternRegionZero = patternStructures->patternRegion;
		iLargestRegionIndex = patternStructures->iLargestRegionIndex;
	}

	for (iRegion = 0; iRegion <= iLastRegion; iRegion++) {
		patternRegion = &patternRegionZero[iRegion];

		if (patternRegion->bIsUnpaired) {
			for (ii = patternRegion->iL; ii <= patternRegion->iR; ii++) {
				iRounds = iCurrentRound - iRoundComputedPerRow[ii];

				if (iRounds < inMax) {
					memmove(compEntries[iRegion].ijLastPerRowPerMatrix[ii], compEntries[iRegion].ijLastPerRowPerMatrix[ii] + iRounds, (inMax - iRounds) * sizeof(int));
					memcpy(compEntries[iRegion].ijLastPerRowPerMatrix[ii] + (inMax - iRounds), iEmptyArray, iRounds * sizeof(int));
				} else {
					memcpy(compEntries[iRegion].ijLastPerRowPerMatrix[ii], iEmptyArray, inMax * sizeof(int));
				}
			}
		} else {
			iRounds = iCurrentRound - iRoundComputedPerRow[patternRegion->bIsArc ? patternRegion->iR : iLineIndex[patternRegion->iR]];

			if (iRounds < inMax) {
				memmove(compEntries[iRegion].ijLast, compEntries[iRegion].ijLast + iRounds, (inMax - iRounds) * sizeof(int));
				memcpy(compEntries[iRegion].ijLast + (inMax - iRounds), iEmptyArray, iRounds * sizeof(int));
			} else {
				memcpy(compEntries[iRegion].ijLast, iEmptyArray, inMax * sizeof(int));
			}

			if (patternRegion->bIsArc) {
				if (iRounds < inMax) {
					for (ik = iRounds; ik < inMax; ik++) {
						if ((compEntries[iRegion].ikFirst[ik - iRounds] = compEntries[iRegion].ikFirst[ik] - iRounds) < 0)
							compEntries[iRegion].ikFirst[ik - iRounds] = 0;
					}
				}
			}
		}


		ijLast                = compEntries[iRegion].ijLast;
		ijLastPerRowPerMatrix = compEntries[iRegion].ijLastPerRowPerMatrix;
		ikFirst               = compEntries[iRegion].ikFirst;
		iTargetRegionIndex    = iLargestRegionIndex[iLastRegion][iRegion];

		iL      = patternRegion->iL;
		iR      = patternRegion->iR;
		ijGoal  = patternRegion->ijGoal;
		ikStart = patternRegion->ikStart;
		ikEnd   = patternRegion->ikEnd > in ? in : patternRegion->ikEnd;

		if (patternRegion->bIsUnpaired) {

			iLeftMinus1 = iLeftMostPos[iL] - 1;
			ijGoalPerRowPerMatrix = patternRegion->ijGoalPerRowPerMatrix;

			ikStart_target = patternRegionZero[iTargetRegionIndex].ikStart;
			if ((ikEnd_target = patternRegionZero[iTargetRegionIndex].iR + iIndels - 1) >= in)
				ikEnd_target = in - 1;

			ijGoal_target = patternRegionZero[iTargetRegionIndex].ijGoal;

			for (ii = iL; ii <= iR; ii++) {
				iRoundComputedPerRow[ii] = iCurrentRound;

				for (ik = ikStart_target; ik <= ikEnd_target; ik++) {

					if ((ijGoal_local = ijGoal_target[ik]) > in - ik) {
						ijGoal_local = in - ik;
					}

					for (ij = 1 + ijLastPerRowPerMatrix[ii][ik]; ij <= ijGoal_local; ij++) {
						imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
						ctemp = OpReplacement;

						if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpDeletion;
						} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpInsertion;
						}
						S(edist, ik, ii, ij)     = imin;
						S(operation, ik, ii, ij) = ctemp;
					}
					ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
				}

				icol2 = ii - iLeftMinus1;
				bBreak = true;

				for (ik2 = ikStart; ik2 <= ikEnd; ik2++) {
					if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik2]) > in - ik2) {
						ijGoal_local = in - ik2;
					}

					icols = iIndels - abs(iLeftMinus1 - ik2); //number of columns to be checked - 1
					istart2 = icol2 > icols ? icol2 - icols : 0;

					for (ij = istart2; ij <= ijGoal_local; ij++) {
						if (S(edist, ik2, ii, ij) <= iThreshold) {
							bBreak = false;
							break;
						}
					}
					if (!bBreak) break;
				}

				if (bBreak) {
					LR->iL = iLeftMostPos[iL];
					LR->iR = ii;

					*iLastRegionCheck = iLastRegion;
					return iRegion;
				}
			}


		} else if (patternRegion->bIsArc) {

			iRoundComputedPerRow[iR] = iCurrentRound;

			ik2end         = patternRegionZero[iTargetRegionIndex].ikEnd;
			ikStart_target = patternRegionZero[iTargetRegionIndex].ikStart;
			ijGoal         = patternRegionZero[iTargetRegionIndex].ijGoal;

			for (ik2 = ikStart_target; ik2 <= ik2end; ik2++) {

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

						imin = S(edist, ik, iLineIndex[iR - 1], ij - 1) + sigma(cu[iR - 1], cv[ij + ik - 1]) * iCostReplacement + iCostArcAltering;
						ctemp = OpArcAltering1;

						if (imin > (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
							imin = itemp;
							ctemp = OpArcAltering2;
						}

						if (imin > (itemp = S(edist, ik, iLineIndex[iR - 1], ij) + iCostArcRemoving)) {
							imin = itemp;
							ctemp = OpArcRemoving;
						}

						if (imin > (itemp = S(edist, ik, iR, ij - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpBInsertion1;
						}

						if (imin > (itemp = S(edist, ik + 1, iR, ij - 1) + iCostDeletion)) {
							imin = itemp;
							ctemp = OpBInsertion2;
						}

						if (ij > 1) {
							if (imin >= (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 2) + (sigma(cu[iR - 1], cv[ij + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[ij + ik - 1]) == true ? 0 : 1))) {
								imin = itemp;
								if (complement(cv[ik], cv[ij + ik - 1]))
									ctemp = OpArcBreaking;
								else
									ctemp = OpArcBreakingTrue;
							}
						}

						S(edist, ik, iR, ij)     = imin;
						S(operation, ik, iR, ij) = ctemp;

						ikFirst[ik] = ik2;
						ijLast[ik] = ij;

						ij++;
					}
				}

			}

			iLeftMinus1 = iL - 1;
			icol   = iR - iLeftMinus1; //mid col of last line

			bBreak = true;
			ikStart = patternRegion->ikStart;
			ikEnd   = patternRegion->ikEnd > in ? in : patternRegion->ikEnd;
			for (ik = ikStart; ik <= ikEnd; ik++) {
				icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

				if ((iend = icol + icols) > in - ik)
					iend = in - ik;
				istart2 = icol > icols ? icol - icols : 0;

				for (ij = istart2; ij <= iend; ij++) {
					if (S(edist, ik, iR, ij) <= iThreshold) {
						bBreak = false;
						break;
					}
				}
				if (!bBreak) break;
			}
			if (bBreak) {
				LR->iL = iL;
				LR->iR = iR;
				*iLastRegionCheck = iLastRegion;
				return iRegion;
			}

		} else {

			iRoundComputedPerRow[iLineIndex[iR]] = iCurrentRound;

			iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
			iLineToMerge2 = iR;

			ikStart_target = patternRegionZero[iTargetRegionIndex].ikStart;
			ijGoal_target  = patternRegionZero[iTargetRegionIndex].ijGoal;

			if ((ikEnd = patternRegionZero[iTargetRegionIndex].ikEndFactual) >= in)
				ikEnd = in - 1;

			for (ik = ikStart_target; ik <= ikEnd; ik++) {

				if ((ijGoal_local = ijGoal_target[ik]) > in - ik)
					ijGoal_local = in - ik;

				for (ij = 1 + ijLast[ik]; ij <= ijGoal_local; ij++) {
					imin   = INT_MAX;
					ij2tmp = 0;
					for (ij2 = 0; ij2 <= ij; ij2++) {
						if (imin > (itemp = S(edist, ik, iLineToMerge1, ij2) + S(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
							imin   = itemp;
							ij2tmp = ij2;
						}
					}
					S(edist, ik, iLineIndex[iLineToMerge2], ij) = imin;
					S(trace, ik, iLineToMerge2, ij) = ij2tmp == 0 ? INT_MAX : ij2tmp;
				}
				ijLast[ik] = ijGoal_local;
			}

			iLeftMinus1 = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0] - 1;
			if (iLineToMerge1 == 0)
				iLeftMinus1++;
			icol = iLineToMerge2 - iLeftMinus1; //mid col for line

			bBreak = true;
			if (ikEnd < 0) ikEnd = 0;
			for (ik = ikStart; ik <= ikEnd; ik++) {
				if ((ijGoal_local = ijGoal[ik]) > in - ik) {
					ijGoal_local = in - ik;
				}

				icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1
				istart2 = icol > icols ? icol - icols : 0;

				for (ij = istart2; ij <= ijGoal_local; ij++) {
					if (S(edist, ik, iLineIndex[iLineToMerge2], ij) <= iThreshold) {
						bBreak = false;
						break;
					}
				}
				if (!bBreak) break;
			}


			if (bBreak) {
				LR->iL = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0];
				LR->iR = iLineToMerge2;
				*iLastRegionCheck = iLastRegion;
				return iRegion;
			}

		}

	}

	for (iRegion = iLastRegion + 1; iRegion < patternStructures->iNumRegions; iRegion++) {

		currentPatternRegion = &patternRegionZero[iRegion];

		if (currentPatternRegion->bIsUnpaired) {
			for (ii = currentPatternRegion->iL; ii <= currentPatternRegion->iR; ii++) {
				iRounds = iCurrentRound - iRoundComputedPerRow[ii];

				if (iRounds < inMax && iRoundComputedPerRow[ii] != 0) {
					memmove(compEntries[iRegion].ijLastPerRowPerMatrix[ii], compEntries[iRegion].ijLastPerRowPerMatrix[ii] + iRounds, (inMax - iRounds) * sizeof(int));
					memcpy(compEntries[iRegion].ijLastPerRowPerMatrix[ii] + (inMax - iRounds), iEmptyArray, iRounds * sizeof(int));
				} else {
					memcpy(compEntries[iRegion].ijLastPerRowPerMatrix[ii], iEmptyArray, inMax * sizeof(int));
				}
			}
		} else {
			iRow = currentPatternRegion->bIsArc ? currentPatternRegion->iR : iLineIndex[currentPatternRegion->iR];
			iRounds = iCurrentRound - iRoundComputedPerRow[iRow];

			if (iRounds < inMax && iRoundComputedPerRow[iRow] != 0) {
				memmove(compEntries[iRegion].ijLast, compEntries[iRegion].ijLast + iRounds, (inMax - iRounds) * sizeof(int));
				memcpy(compEntries[iRegion].ijLast + (inMax - iRounds), iEmptyArray, iRounds * sizeof(int));
			} else {
				memcpy(compEntries[iRegion].ijLast, iEmptyArray, inMax * sizeof(int));
			}

			if (currentPatternRegion->bIsArc) {
				if (iRounds < inMax) {
					for (ik = iRounds; ik < inMax; ik++) {
						if ((compEntries[iRegion].ikFirst[ik - iRounds] = compEntries[iRegion].ikFirst[ik] - iRounds) < 0)
							compEntries[iRegion].ikFirst[ik - iRounds] = 0;
					}
					memcpy(compEntries[iRegion].ikFirst + (inMax - iRounds), iEmptyArray, iRounds * sizeof(int));
				}
			}
		}

		do {

			if ((iDepRegion = currentPatternRegion->iRegionDependsOnRegions[iDepRegionCheck]) != -1) {
				patternRegion = &patternRegionZero[iDepRegion];
				bDepRegion    = true;
				iDepRegionCheck++;

				ijLast                = compEntries[iDepRegion].ijLast;
				ijLastPerRowPerMatrix = compEntries[iDepRegion].ijLastPerRowPerMatrix;
				ikFirst               = compEntries[iDepRegion].ikFirst;
			} else {
				patternRegion   = currentPatternRegion;
				bDepRegion      = false;
				iDepRegionCheck = 0;

				ijLast                = compEntries[iRegion].ijLast;
				ijLastPerRowPerMatrix = compEntries[iRegion].ijLastPerRowPerMatrix;
				ikFirst               = compEntries[iRegion].ikFirst;
			}

			iL      = patternRegion->iL;
			iR      = patternRegion->iR;
			ijGoal  = patternRegion->ijGoal;
			ikStart = patternRegion->ikStart;
			ikEnd   = patternRegion->ikEnd > in ? in : patternRegion->ikEnd;

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
								imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
								ctemp = OpReplacement;

								if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpDeletion;
								} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpInsertion;
								}
								S(edist, ik, ii, ij)     = imin;
								S(operation, ik, ii, ij) = ctemp;
							}
							ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
						}
					}

				} else {

					iLeftMinus1 = iLeftMostPos[iL] - 1;

					ijGoalPerRowPerMatrix = patternRegion->ijGoalPerRowPerMatrix;

					for (ii = iL; ii <= iR; ii++) {
						iRoundComputedPerRow[ii] = iCurrentRound;

						for (ik = ikStart; ik <= ikEnd; ik++) {

							if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik]) > in - ik) {
								ijGoal_local = in - ik;
							}

							ij = 1 + ijLastPerRowPerMatrix[ii][ik];

							if (ii > iL && ijGoal_local >= ij) {
								ij2 = ijGoal_local;
								for (ii2 = iL; ii2 < ii; ii2++) {
									imin = S(edist, ik, iLineIndex[ii2 - 1], ij2 - 1) + (sigma(cu[ii2 - 1], cv[ij2 + ik - 1]) * iCostReplacement);
									ctemp = OpReplacement;

									if (imin >= (itemp = S(edist, ik, iLineIndex[ii2 - 1], ij2) + iCostDeletion)) {
										imin = itemp;
										ctemp = OpDeletion;
									} if (imin >= (itemp = S(edist, ik, ii2, ij2 - 1) + iCostDeletion)) {
										imin = itemp;
										ctemp = OpInsertion;
									}

									S(edist, ik, ii2, ij2)     = imin;
									S(operation, ik, ii2, ij2) = ctemp;
									ijLastPerRowPerMatrix[ii2][ik] = ijGoal_local;
								}
							}

							for (; ij <= ijGoal_local; ij++) {
								imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
								ctemp = OpReplacement;

								if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpDeletion;
								} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpInsertion;
								}
								S(edist, ik, ii, ij)     = imin;
								S(operation, ik, ii, ij) = ctemp;

								ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
							}
						}

						icol2 = ii - iLeftMinus1;
						bBreak = true;

						for (ik2 = ikStart; ik2 <= ikEnd; ik2++) {
							if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik2]) > in - ik2) {
								ijGoal_local = in - ik2;
							}

							icols = iIndels - abs(iLeftMinus1 - ik2); //number of columns to be checked - 1
							istart2 = icol2 > icols ? icol2 - icols : 0;

							for (ij = istart2; ij <= ijGoal_local; ij++) {
								if (S(edist, ik2, ii, ij) <= iThreshold) {
									bBreak = false;
									break;
								}
							}
							if (!bBreak) break;
						}

						if (bBreak) {
							LR->iL = iLeftMostPos[iL];
							LR->iR = ii;
							*iLastRegionCheck = iRegion;
							return iRegion;
						}
					}
				}

			} else if (patternRegion->bIsArc) {

				iRoundComputedPerRow[iR] = iCurrentRound;

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

							imin = S(edist, ik, iLineIndex[iR - 1], ij - 1) + sigma(cu[iR - 1], cv[ij + ik - 1]) * iCostReplacement + iCostArcAltering;
							ctemp = OpArcAltering1;

							if (imin > (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
								imin = itemp;
								ctemp = OpArcAltering2;
							}

							if (imin > (itemp = S(edist, ik, iLineIndex[iR - 1], ij) + iCostArcRemoving)) {
								imin = itemp;
								ctemp = OpArcRemoving;
							}

							if (imin > (itemp = S(edist, ik, iR, ij - 1) + iCostDeletion)) {
								imin = itemp;
								ctemp = OpBInsertion1;
							}

							if (imin > (itemp = S(edist, ik + 1, iR, ij - 1) + iCostDeletion)) {
								imin = itemp;
								ctemp = OpBInsertion2;
							}

							if (ij > 1) {
								if (imin >= (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 2) + (sigma(cu[iR - 1], cv[ij + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[ij + ik - 1]) == true ? 0 : 1))) {
									imin = itemp;
									if (complement(cv[ik], cv[ij + ik - 1]))
										ctemp = OpArcBreaking;
									else
										ctemp = OpArcBreakingTrue;
								}
							}

							S(edist, ik, iR, ij)     = imin;
							S(operation, ik, iR, ij) = ctemp;

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
					ikStart = patternRegion->ikStart;
					ikEnd   = patternRegion->ikEnd > in ? in : patternRegion->ikEnd;
					for (ik = ikStart; ik <= ikEnd; ik++) {
						icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

						if ((iend = icol + icols) > in - ik)
							iend = in - ik;
						istart2 = icol > icols ? icol - icols : 0;

						for (ij = istart2; ij <= iend; ij++) {
							if (S(edist, ik, iR, ij) <= iThreshold) {
								bBreak = false;
								break;
							}
						}
						if (!bBreak) break;
					}
					if (bBreak) {
						LR->iL = iL;
						LR->iR = iR;
						*iLastRegionCheck = iRegion;
						return iRegion;
					}
				}

			} else {

				iRoundComputedPerRow[iLineIndex[iR]] = iCurrentRound;

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
							if (imin > (itemp = S(edist, ik, iLineToMerge1, ij2) + S(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
								imin   = itemp;
								ij2tmp = ij2;
							}
						}
						S(edist, ik, iLineIndex[iLineToMerge2], ij) = imin;
						S(trace, ik, iLineToMerge2, ij) = ij2tmp == 0 ? INT_MAX : ij2tmp;
					}
					ijLast[ik] = ijGoal_local;
				}

				if (bDepRegion) {
					bBreak = false; //x
				} else {

					iLeftMinus1 = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0] - 1;
					if (iLineToMerge1 == 0)
						iLeftMinus1++;
					icol = iLineToMerge2 - iLeftMinus1; //mid col for line

					bBreak = true;
					if (ikEnd < 0) ikEnd = 0;
					for (ik = ikStart; ik <= ikEnd; ik++) {
						if ((ijGoal_local = ijGoal[ik]) > in - ik) {
							ijGoal_local = in - ik;
						}

						icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1
						istart2 = icol > icols ? icol - icols : 0;

						for (ij = istart2; ij <= ijGoal_local; ij++) {
							if (S(edist, ik, iLineIndex[iLineToMerge2], ij) <= iThreshold) {
								bBreak = false;
								break;
							}
						}
						if (!bBreak) break;
					}
				}

				if (bBreak) {
					LR->iL = iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0];
					LR->iR = iLineToMerge2;
					*iLastRegionCheck = iRegion;
					return iRegion;
				}

			}

		} while (bDepRegion);

	}

	LR->iL = 1;
	LR->iR = im;
	*bSuccess = true;
	*iLastRegionCheck = patternStructures->iNumRegions - 1;
	return patternStructures->iNumRegions - 1;
}

int alignESA_LA_NEW_S(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int *ijLastOfRegion,
		char *cv, const int in, const int inMax, AffixArray *affixArray, int iMaxColOffset, int iLastCompRow, int iLastRegion, int *iRoundComputedPerRow, bool *bSuccess, LR *LR, bool bInOutAlignment) {

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

	PatternRegion *patternRegion, *currentPatternRegion, *patternRegionZero;

	PatternStructures *patternStructures = patternStrandDirection->patternStructures;
	int *iLeftMostPos = patternStructures->iLeftMostPos;
	char *cu = (char *) patternStrandDirection->seq;

	int *iLineIndex = patternStructures->iLineIndex;
	im = pattern->iLength;

	int *ijLast, **ijLastPerRowPerMatrix, *ijGoal, **ijGoalPerRowPerMatrix, *ikFirst;
	int ikStart, ikEnd, iDepRegion, ijGoal_local, iDepRegionCheck = 0;
	bool bDepRegion, bBreak;

	if (bInOutAlignment) {
		patternRegionZero = patternStructures->patternRegionInsideOut;
	} else {
		patternRegionZero = patternStructures->patternRegion;
	}

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {

		currentPatternRegion = &patternRegionZero[iRegion];

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
				patternRegion = &patternRegionZero[iDepRegion];
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
			ikEnd   = patternRegion->ikEnd > in ? in : patternRegion->ikEnd;

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
								imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
								ctemp = OpReplacement;

								if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpDeletion;
								} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpInsertion;
								}
								S(edist, ik, ii, ij)     = imin;
								S(operation, ik, ii, ij) = ctemp;
							}
							ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
						}
					}

				} else {

					iLeftMinus1 = iLeftMostPos[iL] - 1;

					ijGoalPerRowPerMatrix = patternRegion->ijGoalPerRowPerMatrix;

					for (ii = iL; ii <= iR; ii++) {
						iRoundComputedPerRow[ii] = 0;

						for (ik = ikStart; ik <= ikEnd; ik++) {

							if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik]) > in - ik) {
								ijGoal_local = in - ik;
							}

							ij = 1 + ijLastPerRowPerMatrix[ii][ik];

							if (ii > iL && ijGoal_local >= ij) {
								ij2 = ijGoal_local;
								for (ii2 = iL; ii2 < ii; ii2++) {
									imin = S(edist, ik, iLineIndex[ii2 - 1], ij2 - 1) + (sigma(cu[ii2 - 1], cv[ij2 + ik - 1]) * iCostReplacement);
									ctemp = OpReplacement;

									if (imin >= (itemp = S(edist, ik, iLineIndex[ii2 - 1], ij2) + iCostDeletion)) {
										imin = itemp;
										ctemp = OpDeletion;
									} if (imin >= (itemp = S(edist, ik, ii2, ij2 - 1) + iCostDeletion)) {
										imin = itemp;
										ctemp = OpInsertion;
									}

									S(edist, ik, ii2, ij2)     = imin;
									S(operation, ik, ii2, ij2) = ctemp;
									ijLastPerRowPerMatrix[ii2][ik] = ijGoal_local;
								}
							}

							for (; ij <= ijGoal_local; ij++) {
								imin = S(edist, ik, iLineIndex[ii - 1], ij - 1) + (sigma(cu[ii - 1], cv[ij + ik - 1]) * iCostReplacement);
								ctemp = OpReplacement;

								if (imin >= (itemp = S(edist, ik, iLineIndex[ii - 1], ij) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpDeletion;
								} if (imin >= (itemp = S(edist, ik, ii, ij - 1) + iCostDeletion)) {
									imin = itemp;
									ctemp = OpInsertion;
								}
								S(edist, ik, ii, ij)     = imin;
								S(operation, ik, ii, ij) = ctemp;
								ijLastPerRowPerMatrix[ii][ik] = ijGoal_local;
							}
						}

						icol2 = ii - iLeftMinus1;
						bBreak = true;

						for (ik2 = ikStart; ik2 <= ikEnd; ik2++) {
							if ((ijGoal_local = ijGoalPerRowPerMatrix[ii][ik2]) > in - ik2) {
								ijGoal_local = in - ik2;
							}

							icols = iIndels - abs(iLeftMinus1 - ik2); //number of columns to be checked - 1
							istart2 = icol2 > icols ? icol2 - icols : 0;

							for (ij = istart2; ij <= ijGoal_local; ij++) {
								if (S(edist, ik2, ii, ij) <= iThreshold) {
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

				iRoundComputedPerRow[iR] = 0;

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

							imin = S(edist, ik, iLineIndex[iR - 1], ij - 1) + sigma(cu[iR - 1], cv[ij + ik - 1]) * iCostReplacement + iCostArcAltering;
							ctemp = OpArcAltering1;

							if (imin > (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 1) + sigma(cu[iL - 1], cv[ik]) * iCostReplacement + iCostArcAltering)) {
								imin = itemp;
								ctemp = OpArcAltering2;
							}

							if (imin > (itemp = S(edist, ik, iLineIndex[iR - 1], ij) + iCostArcRemoving)) {
								imin = itemp;
								ctemp = OpArcRemoving;
							}

							if (imin > (itemp = S(edist, ik, iR, ij - 1) + iCostDeletion)) {
								imin = itemp;
								ctemp = OpBInsertion1;
							}

							if (imin > (itemp = S(edist, ik + 1, iR, ij - 1) + iCostDeletion)) {
								imin = itemp;
								ctemp = OpBInsertion2;
							}

							if (ij > 1) {
								if (imin >= (itemp = S(edist, ik + 1, iLineIndex[iR - 1], ij - 2) + (sigma(cu[iR - 1], cv[ij + ik - 1]) + sigma(cu[iL - 1], cv[ik])) * iCostReplacement + iCostArcBreak * (complement(cv[ik], cv[ij + ik - 1]) == true ? 0 : 1))) {
									imin = itemp;
									if (complement(cv[ik], cv[ij + ik - 1]))
										ctemp = OpArcBreaking;
									else
										ctemp = OpArcBreakingTrue;
								}
							}

							S(edist, ik, iR, ij)     = imin;
							S(operation, ik, iR, ij) = ctemp;

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
					ikStart = patternRegion->ikStart;
					ikEnd   = patternRegion->ikEnd > in ? in : patternRegion->ikEnd;
					for (ik = ikStart; ik <= ikEnd; ik++) {
						icols = iIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

						if ((iend = icol + icols) > in - ik)
							iend = in - ik;
						istart2 = icol > icols ? icol - icols : 0;

						for (ij = istart2; ij <= iend; ij++) {
							if (S(edist, ik, iR, ij) <= iThreshold) {
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

				iRoundComputedPerRow[iLineIndex[iR]] = 0;

				iLineToMerge1 = iLineIndex[patternRegion->iLineToMerge1];
				iLineToMerge2 = iR;

				if (bDepRegion) {
					ikStart = currentPatternRegion->ikStart;
					ijGoal  = currentPatternRegion->ijGoal;

					if ((ikEnd = currentPatternRegion->ikEndFactual) >= in) {
						ikEnd = in - 1;
					}
				}

				for (ik = ikStart; ik <= ikEnd; ik++) {
					if ((ijGoal_local = ijGoal[ik]) > in - ik) {
						ijGoal_local = in - ik;
					}

					for (ij = 1 + ijLast[ik]; ij <= ijGoal_local; ij++) {
						imin   = INT_MAX;
						ij2tmp = 0;
						for (ij2 = 0; ij2 <= ij; ij2++) {
							if (imin > (itemp = S(edist, ik, iLineToMerge1, ij2) + S(edist, ik + ij2, iLineToMerge2, ij - ij2))) {
								imin   = itemp;
								ij2tmp = ij2;
							}
						}

						S(edist, ik, iLineIndex[iLineToMerge2], ij) = imin;
						S(trace, ik, iLineToMerge2, ij) = ij2tmp == 0 ? INT_MAX : ij2tmp;
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
							if (S(edist, ik, iLineIndex[iLineToMerge2], ij) <= iThreshold) {
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

	*bSuccess = true;
	LR->iL = 1;
	LR->iR = im;
	return patternStructures->iNumRegions - 1;
}

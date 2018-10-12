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

void getEditOperationsGlobal_S(char *cSeqOperations, char *cArcOperations, char ***operation, unsigned int ***trace,
		PatternStructures *pstr, unsigned int *uiIndex, unsigned int uik, unsigned int uii, unsigned int uij, bool bAllowBranch) {
	char cOp;
	unsigned int traceValue;

	if (S(trace, uik, uii, uij) != 0 && bAllowBranch) {
		cOp = 99;
	} else {
		cOp = S(operation, uik, uii, uij);
	}

	switch (cOp) {
	case OpReplacement:
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpInsertion:
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, uii, uij - 1, bAllowBranch);
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpDeletion:
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, uii - 1, uij, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpArcBreaking:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '(';
		(*uiIndex)++;
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik + 1, uii - 1, uij - 2, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = ')';
		(*uiIndex)++;
		break;
	case OpArcBreakingTrue:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik + 1, uii - 1, uij - 2, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpArcAltering1:
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpArcAltering2:
		cSeqOperations[*uiIndex] = OpReplacement;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik + 1, uii - 1, uij - 1, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpArcRemoving:
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, uii - 1, uij, true);
		cSeqOperations[*uiIndex] = OpDeletion;
		cArcOperations[*uiIndex] = '-';
		(*uiIndex)++;
		break;
	case OpBInsertion1:
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, uii, uij - 1, bAllowBranch);
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		break;
	case OpBInsertion2:
		cSeqOperations[*uiIndex] = OpInsertion;
		cArcOperations[*uiIndex] = '.';
		(*uiIndex)++;
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik + 1, uii, uij - 1, bAllowBranch);
		break;
	case 99:
		if (S(trace, uik, uii, uij) == INT_MAX)
			traceValue = 0;
		else
			traceValue = S(trace, uik, uii, uij);
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik, pstr->patternRegion[pstr->iPosRegion[uii]].iL - 1, traceValue, true);
		getEditOperationsGlobal_S(cSeqOperations, cArcOperations, operation, trace, pstr, uiIndex, uik + traceValue, uii, uij - traceValue, false);
		break;
	}
}

unsigned int alignGlobal(PatternStrandDirection *patternStrandDirection, AffixArray *affixArray, SearchParam *searchParam) {
	int ii, im;
	int imax = 0, iMaxOperations;
	int iSeqLength;

	unsigned int ***edist, ***trace, **uiPPEdist, **uiPPTrace;
	unsigned int uiNumSeqs = affixArray->multiSeq->numSeqs;
	unsigned int uiSeqStartPos; // position of the first sequence character in the concatenated sequence
	unsigned int uiNumOperations;

	Pattern *pattern = patternStrandDirection->pattern;

	char ***operation, **cPPOperation;
	char *cu_str, *seqs;

	Match match;

	cu_str = (char *) patternStrandDirection->structure;

	im   = pattern->iLength;

	imax = affixArray->multiSeq->seqEndPos[0];

	for (ii = 1; ii < uiNumSeqs; ii++) {
		uiSeqStartPos = affixArray->multiSeq->seqEndPos[ii - 1] + 1;
		iSeqLength = affixArray->multiSeq->seqEndPos[ii] - affixArray->multiSeq->seqEndPos[ii - 1] - 1;

		if (iSeqLength > imax)
			imax = iSeqLength;
	}

	iMaxOperations = (imax > im) ? (imax + 1) * 2 : (im + 1) * 2;
	match.cSeqOperations = (char *) calloc(iMaxOperations, sizeof(char));
	match.cArcOperations = (char *) calloc(iMaxOperations, sizeof(char));

	edist     = uiNew3DMatrixX(im, imax);
	operation = cNew3DMatrixX(im, imax);
	trace     = uiNew3DMatrixX(im, imax);

	initializeRowsAndColumnsRegions_S(edist, operation, trace, pattern->cost, cu_str, im, imax, patternStrandDirection->patternStructures);

	uiPPEdist    = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	uiPPTrace    = (unsigned int **) calloc(imax + 1, sizeof(unsigned int *));
	cPPOperation = (char **) calloc(imax + 1, sizeof(char *));

	int iRemoveOverlapsInterval = affixArray->length / FILTERFREQ;
	int iRemoveControl[FILTERFREQ], iRemoveControlIndex = 0;

	iRemoveControl[0] = iRemoveOverlapsInterval;
	for (ii = 1; ii < FILTERFREQ; ii++) {
		iRemoveControl[ii] = iRemoveControl[ii - 1] + iRemoveOverlapsInterval;
	}

	for (ii = 0; ii < uiNumSeqs; ii++) {
		if (ii == 0) {
			uiSeqStartPos = 0;
			seqs = affixArray->multiSeq->convSequences;
			iSeqLength = affixArray->multiSeq->seqEndPos[0];
		} else {
			uiSeqStartPos = affixArray->multiSeq->seqEndPos[ii - 1] + 1;
			seqs = affixArray->multiSeq->convSequences + uiSeqStartPos;
			iSeqLength = affixArray->multiSeq->seqEndPos[ii] - affixArray->multiSeq->seqEndPos[ii - 1] - 1;
		}

		alignESA_S(edist,
				operation,
				trace,
				patternStrandDirection,
				seqs,
				iSeqLength,
				affixArray,
				0);

		uiNumOperations = 0;

		getEditOperationsGlobal_S(match.cSeqOperations, match.cArcOperations, operation, trace, patternStrandDirection->patternStructures, &uiNumOperations, 0, im, iSeqLength, true);

		match.cSeqOperations[uiNumOperations] = 0;
		match.patternStrandDirection = patternStrandDirection;
		match.iSeqId     = ii;
		match.uiPos      = uiSeqStartPos;
		match.uiEndPos   = iSeqLength;
		match.iCost           = S(edist, 0, patternStrandDirection->patternStructures->iLineIndex[im], iSeqLength);
		match.iScore          = pattern->weight - match.iCost;
		match.bForwardStrand  = patternStrandDirection->bForwardStrand;
		match.uiNumOperations = uiNumOperations;

		printf(">%s\n", affixArray->multiSeq->seqDescription[match.iSeqId]);
		printf("\n Pattern desc.: %s, Strand = %s, Edist = %d\n\n", match.patternStrandDirection->pattern->desc, match.patternStrandDirection->bForwardStrand ? "forward" : "reverse", match.iCost);
		printAlignment(&match, affixArray->multiSeq->convSequences + match.uiPos, 0, affixArray);
		printf("\n");

		if (searchParam->bShowProgress && uiSeqStartPos > iRemoveControl[iRemoveControlIndex]) {
			iRemoveControlIndex++;
			fprintf(stderr, "*%.1f%% done\n", ((float) iRemoveControlIndex / (float) FILTERFREQ) * 100.0);
			fflush(stderr);
		}

	}

	uiFree3DMatrixX(edist, im, imax);
	cFree3DMatrixX(operation, im, imax);
	uiFree3DMatrixX(trace, im, imax);

	free(match.cSeqOperations);
	free(match.cArcOperations);

	free(uiPPEdist);
	free(uiPPTrace);
	free(cPPOperation);

	return 0;
}


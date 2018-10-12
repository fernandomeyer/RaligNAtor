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
#include <limits.h>
#include <string.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "align.h"
#include "redblack/redblack.h"

bool adjustNumIndelsAccordingToEdist (unsigned int* uiIndels, unsigned int* uiThreshold, Cost *cost) {
	unsigned int uiTestNumIndels;
	bool bStatus = false;

	if (*uiThreshold == -1)
		*uiThreshold = 0;

	uiTestNumIndels = *uiThreshold / cost->iCostDeletion;

	if (*uiIndels == -1)
		*uiIndels = uiTestNumIndels;
	if (*uiIndels > uiTestNumIndels) {
		*uiIndels = uiTestNumIndels;
		bStatus = true;
	}

	return bStatus;
}

void init(AffixArray *affixArray, SearchParam *searchParam,
		bool *argSuf, bool *argLcp, bool *argSufinv, bool *argLnk, bool *argAflk, bool *argSufr, bool *argLcpr, bool *argAflkr,
		Alphabet *alphabet, char *argAlphabetFile, const int predefAlphabet) {

	*argSuf   = false;
	*argSufinv = false;
	*argLcp   = false;
	*argAflk  = false;
	*argSufr  = false;
	*argLcpr  = false;
	*argAflkr = false;

	switch(searchParam->cVariant) {
		case V_LESA:
			*argSuf = true;
			*argLcp = true;
			break;
		case V_ESA:
			*argSuf = true;
			*argLcp = true;
			break;
		case V_LGSLINK:
			*argSuf = true;
			*argLcp = true;
			*argSufinv = true;
			break;
		case V_LSLINK:
			*argSuf = true;
			*argLcp = true;
			*argSufinv = true;
			break;
		case V_GSLINK:
			*argSuf = true;
			*argLcp = true;
			*argSufinv = true;
			break;
		case V_SLINK:
			*argSuf = true;
			*argLcp = true;
			*argSufinv = true;
			break;
		default:
			break;
	}

	affixArray->alphabet = alphabet;

	if((affixArray->multiSeq = (MultiSeq *) malloc(sizeof(MultiSeq))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"affixArray->multiSeq\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	affixArray->multiSeq->convSequences = NULL;
	affixArray->multiSeq->convSequencesMmapped = false;

	if (*argSuf) {
		if((affixArray->esa = (EArray *) malloc(sizeof(EArray))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"esa\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		affixArray->esa->xarray = NULL;
		affixArray->esa->xlcp = NULL;
		affixArray->esa->xlcpException = NULL;
		affixArray->esa->affixLink = NULL;
	} else {
		affixArray->esa  = NULL;
	}

	if (*argSufr) {
		if((affixArray->erpa = (EArray *) malloc(sizeof(EArray))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"erpa\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		affixArray->erpa->xarray = NULL;
		affixArray->erpa->xlcp = NULL;
		affixArray->erpa->xlcpException = NULL;
		affixArray->erpa->affixLink = NULL;
	} else {
		affixArray->erpa = NULL;
	}

	if (alphabet != NULL) {
		alphabet->conversionTablePattern   = NULL;
		alphabet->conversionTableTargetSeq = NULL;
		if (searchParam->cVariant == V_SCAN || searchParam->cVariant == V_SCANLA || searchParam->cVariant == V_STD
				|| searchParam->cVariant == V_STDARRAY || searchParam->cVariant == V_GLOBAL) {
			if (argAlphabetFile != NULL) {
				loadAlphabetFile(argAlphabetFile, alphabet);
			} else {
				if (predefAlphabet > 2)
					setPredefinedAlphabet(0, alphabet);
				else
					setPredefinedAlphabet(predefAlphabet, alphabet);
			}
		}
	}

	if (searchParam->chainParam->isactive || searchParam->bPrintMatchesBySeq || searchParam->bFilterOverlaps)
		searchParam->bStoreInRB = true;
	else
		searchParam->bStoreInRB = false;

	if (searchParam->bPrintMatchesBySeq || searchParam->cPrintMatchesByScore > 0)
		searchParam->bPrintOutMatches = true;
}

void freeAll(AffixArray *affixArray, MultiPattern *multiPattern) {
	unsigned int ui;

	if (affixArray->multiSeq != NULL) {
		if (affixArray->multiSeq->convSequences != NULL) {
			if (affixArray->multiSeq->convSequencesMmapped)
				munmap(affixArray->multiSeq->convSequences, affixArray->length * sizeof(unsigned char));
			else
				free(affixArray->multiSeq->convSequences);
		}
		if (affixArray->multiSeq->seqDescLength != NULL)
			free(affixArray->multiSeq->seqDescLength);
		if (affixArray->multiSeq->seqDescription != NULL) {
			for (ui = 0; ui < affixArray->multiSeq->numSeqs; ui++) {
				free(affixArray->multiSeq->seqDescription[ui]);
			}
			free(affixArray->multiSeq->seqDescription);
		}
		if (affixArray->multiSeq->seqEndPos != NULL)
			free(affixArray->multiSeq->seqEndPos);
		if (affixArray->multiSeq->sequences != NULL) {
			if (affixArray->multiSeq->sequencesMmapped)
				munmap(affixArray->multiSeq->sequences, affixArray->length * sizeof(unsigned char));
			else
				free(affixArray->multiSeq->sequences);
		}
		free(affixArray->multiSeq);
	}

	if (affixArray->alphabet != NULL) {
		if (affixArray->alphabet->classRepresentative)
			free(affixArray->alphabet->classRepresentative);
		if (affixArray->alphabet->classSize)
			free(affixArray->alphabet->classSize);
		if (affixArray->alphabet->eqClass) {
			for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
				free(affixArray->alphabet->eqClass[ui]);
			}
			free(affixArray->alphabet->eqClass);
		}
		if (affixArray->alphabet->isWildCard)
			free(affixArray->alphabet->isWildCard);

		if (affixArray->alphabet->iupacTable != NULL) {
			for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
				free(affixArray->alphabet->iupacTable[ui]);
			}
			free(affixArray->alphabet->iupacTable);
		}

		if (affixArray->alphabet->conversionTablePattern != NULL) {
			affixArray->alphabet->conversionTablePattern = NULL;
			free(affixArray->alphabet->conversionTablePattern);
		}
		if (affixArray->alphabet->conversionTableTargetSeq != NULL) {
			affixArray->alphabet->conversionTableTargetSeq = NULL;
			free(affixArray->alphabet->conversionTableTargetSeq);
		}
	}

	if (affixArray->esa != NULL) {
		if (affixArray->esa->affixLink != NULL)
			munmap(affixArray->esa->affixLink, affixArray->length * sizeof(int));
		if (affixArray->esa->xarray != NULL)
			munmap(affixArray->esa->xarray, affixArray->length * sizeof(int));
		if (affixArray->esa->xlcp != NULL)
			munmap(affixArray->esa->xlcp, affixArray->length * sizeof(unsigned char));
		if (affixArray->esa->xlcpException != NULL) {
			if (affixArray->esa->xlcpException->numExceptions > 0) {
				free(affixArray->esa->xlcpException->index);
				free(affixArray->esa->xlcpException->value);
			}
			free(affixArray->esa->xlcpException);
		}
		free(affixArray->esa);
	}

	if (affixArray->erpa != NULL) {
		if (affixArray->erpa->affixLink != NULL)
			munmap(affixArray->erpa->affixLink, affixArray->length * sizeof(int));
		if (affixArray->erpa->xarray != NULL)
			munmap(affixArray->erpa->xarray, affixArray->length * sizeof(int));
		if (affixArray->erpa->xlcp != NULL)
			munmap(affixArray->erpa->xlcp, affixArray->length * sizeof(unsigned char));
		if (affixArray->erpa->xlcpException != NULL) {
			if (affixArray->erpa->xlcpException->numExceptions > 0) {
				free(affixArray->erpa->xlcpException->index);
				free(affixArray->erpa->xlcpException->value);
			}
			free(affixArray->erpa->xlcpException);
		}
		free(affixArray->erpa);
	}

	if (multiPattern != NULL) {

		if (multiPattern->pattern != NULL) {
			if (multiPattern->pattern[0].forwardStrand != NULL) {
				for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
					free(multiPattern->pattern[0].forwardStrand->compRules[ui]);
				}
				free(multiPattern->pattern[0].forwardStrand->compRules);
			}
			if (multiPattern->pattern[0].reverseStrand != NULL) {
				for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
					free(multiPattern->pattern[0].reverseStrand->compRules[ui]);
				}
				free(multiPattern->pattern[0].reverseStrand->compRules);
			}

			for (ui = 0; ui < multiPattern->iNumPatterns; ui++) {
				if (multiPattern->pattern[ui].desc != NULL)
					free(multiPattern->pattern[ui].desc);
				if (multiPattern->pattern[ui].forwardStrand != NULL) {
					if (multiPattern->pattern[ui].forwardStrand->seq != NULL)
						free(multiPattern->pattern[ui].forwardStrand->seq);
					if (multiPattern->pattern[ui].forwardStrand->structure != NULL)
						free(multiPattern->pattern[ui].forwardStrand->structure);
				}
				if (multiPattern->pattern[ui].reverseStrand != NULL) {
					if (multiPattern->pattern[ui].reverseStrand->seq != NULL)
						free(multiPattern->pattern[ui].reverseStrand->seq);
					if (multiPattern->pattern[ui].reverseStrand->structure != NULL)
						free(multiPattern->pattern[ui].reverseStrand->structure);
				}
			}
			free(multiPattern->pattern);
		}
		free(multiPattern);
	}
}

Cost* getUnitCosts() {
	Cost *cost = (Cost *) malloc(sizeof(Cost));

	cost->iCostReplacement = 1;
	cost->iCostDeletion    = 1;
	cost->iCostArcBreak    = 1;
	cost->iCostArcAltering = 1;
	cost->iCostArcRemoving = 2;

	return cost;
}

unsigned int*** uiNew3DMatrixWithEntryCopies(unsigned int uim, unsigned int uin, char *cu_str) {
	unsigned int uii, uik;
	unsigned int ***matrix = (unsigned int ***) calloc(uin + 1, sizeof(unsigned int**)); // uin + 1 matrices

	int i, iOpeningBracket = -1;
	bool *bInitPositions = (bool *) calloc(uim, sizeof(bool));

	// Process pattern to determine the matrix lines that will be copies of lines 0
	for (i = 0; i < uim; i++) {
		if (cu_str[i] == '(') {
			iOpeningBracket = i;
		} else if (cu_str[i] == '.') {
			bInitPositions[iOpeningBracket] = true;
		} else if (cu_str[i] == ')') {
			bInitPositions[iOpeningBracket] = true;
			break;
		}
	}

	for (uik = 0; uik <= uin; uik++) {
		matrix[uik] = (unsigned int **) calloc(uim + 1, sizeof(unsigned int*));  //uim + 1 lines
		if (uik == 0)
			matrix[uik][0] = (unsigned int *) calloc(uin + 1, sizeof(unsigned int));
		else
			matrix[uik][0] = matrix[0][0];

		for (uii = 1; uii <= uim; uii++) {
			if (bInitPositions[uii - 1] == true)
				matrix[uik][uii] = matrix[0][0];
			else {
				matrix[uik][uii] = (unsigned int *) calloc(uin + 1 - uik, sizeof(unsigned int)); // uin + 1 - k columns
			}
		}
	}

	free(bInitPositions);

	return matrix;
}

char*** cNew3DMatrixWithEntryCopies(unsigned int uim, unsigned int uin, char *cu_str) {
	unsigned int uii, uik;
	char ***matrix = (char ***) calloc(uin + 1, sizeof(char**)); // uin + 1 matrices

	int i, iOpeningBracket = -1;
	bool *bInitPositions = (bool *) calloc(uim, sizeof(bool));

	//Preprocess pattern
	for (i = 0; i < uim; i++) {
		if (cu_str[i] == '(') {
			iOpeningBracket = i;
		} else if (cu_str[i] == '.') {
			bInitPositions[iOpeningBracket] = true;
		} else if (cu_str[i] == ')') {
			bInitPositions[iOpeningBracket] = true;
			break;
		}
	}

	for (uik = 0; uik <= uin; uik++) {
		matrix[uik] = (char **) calloc(uim + 1, sizeof(char*));  //uim + 1 lines
		if (uik == 0)
			matrix[uik][0] = (char *) calloc(uin + 1, sizeof(char));
		else
			matrix[uik][0] = matrix[0][0];

		for (uii = 1; uii <= uim; uii++) {
			if (bInitPositions[uii - 1] == true)
				matrix[uik][uii] = matrix[0][0];
			else
				matrix[uik][uii] = (char *) calloc(uin + 1 - uik, sizeof(char)); // uin + 1 - k columns
		}
	}

	free(bInitPositions);

	return matrix;
}

#if KIJ == 1

unsigned int*** uiNew3DMatrix(unsigned int uim, unsigned int uin) {
	unsigned int uik, uii;
	unsigned int ***matrix = (unsigned int ***) calloc(uin + 1, sizeof(unsigned int**));

	for (uik = 0; uik <= uin; uik++) {
		matrix[uik] = (unsigned int **) calloc(uim + 1, sizeof(unsigned int*));  //uim + 1 lines
		for (uii = 0; uii <= uim; uii++) {
			matrix[uik][uii] = (unsigned int *) calloc(uin + 1 - uik, sizeof(unsigned int)); // uin + 1 - k columns
		}
	}

	return matrix;
}

char*** cNew3DMatrix(unsigned int uim, unsigned int uin) {
	unsigned int uik, uii;
	char ***matrix = (char ***) calloc(uin + 1, sizeof(char**));

	for (uik = 0; uik <= uin; uik++) {
		matrix[uik] = (char **) calloc(uim + 1, sizeof(char*));  //uim + 1 lines
		for (uii = 0; uii <= uim; uii++) {
			matrix[uik][uii] = (char *) calloc(uin + 1 - uik, sizeof(char)); // uin + 1 - k columns
		}
	}

	return matrix;
}

#else

unsigned int*** uiNew3DMatrix(unsigned int uim, unsigned int uin) {
	unsigned int uik, uii;
	unsigned int ***matrix = (unsigned int ***) calloc(uim + 1, sizeof(unsigned int**));

	for (uii = 0; uii <= uim; uii++) {
		matrix[uii] = (unsigned int **) calloc(uin + 1, sizeof(unsigned int*));
		for (uik = 0; uik <= uin; uik++) {
			matrix[uii][uik] = (unsigned int *) calloc(uin + 1 - uik, sizeof(unsigned int));
		}
	}

	return matrix;
}

int uiFree3DMatrix(unsigned int ***matrix, unsigned int uim, unsigned int uin) {
	unsigned int uik, uii;

	for (uii = 0; uii <= uim; uii++) {
		for (uik = 0; uik <= uin; uik++) {
			free(matrix[uii][uik]);
		}
		free(matrix[uii]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

char*** cNew3DMatrix(unsigned int uim, unsigned int uin) {
	unsigned int uik, uii;
	char ***matrix = (char ***) calloc(uim + 1, sizeof(char**));

	for (uii = 0; uii <= uim; uii++) {
		matrix[uii] = (char **) calloc(uin + 1, sizeof(char*));
		for (uik = 0; uik <= uin; uik++) {
			matrix[uii][uik] = (char *) calloc(uin + 1 - uik, sizeof(char));
		}
	}

	return matrix;
}

int cFree3DMatrix(char ***matrix, unsigned int uim, unsigned int uin) {
	unsigned int uik, uii;

	for (uii = 0; uii <= uim; uii++) {
		for (uik = 0; uik <= uin; uik++) {
			free(matrix[uii][uik]);
		}
		free(matrix[uii]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

unsigned int*** uiNew3DMatrixX(unsigned int uim, unsigned int uin) {
	unsigned int uik, uij, uiend;
	unsigned int ***matrix = (unsigned int ***) calloc(uin + 1, sizeof(unsigned int**));

	for (uik = 0; uik <= uin; uik++) {
		uiend = uin + 1 - uik;

		matrix[uik] = (unsigned int **) calloc(uiend, sizeof(unsigned int*));  //uim + 1 - k columns
		for (uij = 0; uij < uiend; uij++) {
			matrix[uik][uij] = (unsigned int *) calloc(uim + 1, sizeof(unsigned int)); // uim + 1 lines
		}
	}

	return matrix;
}

int uiFree3DMatrixX(unsigned int ***matrix, unsigned int uim, unsigned int uin) {
	unsigned int uik, uij, uiend;

	for (uik = 0; uik <= uin; uik++) {
		uiend = uin + 1 - uik;
		for (uij = 0; uij < uiend; uij++) {
			free(matrix[uik][uij]);
		}
		free(matrix[uik]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

char*** cNew3DMatrixX(unsigned int uim, unsigned int uin) {
	unsigned int uik, uij, uiend;
	char ***matrix = (char ***) calloc(uin + 1, sizeof(char**));

	for (uik = 0; uik <= uin; uik++) {
		uiend = uin + 1 - uik;

		matrix[uik] = (char **) calloc(uiend, sizeof(char*));  //uim + 1 - k columns
		for (uij = 0; uij < uiend; uij++) {
			matrix[uik][uij] = (char *) calloc(uim + 1, sizeof(char)); // uim + 1 lines
		}
	}

	return matrix;
}

int cFree3DMatrixX(char ***matrix, unsigned int uim, unsigned int uin) {
	unsigned int uik, uij, uiend;

	for (uik = 0; uik <= uin; uik++) {
		uiend = uin + 1 - uik;
		for (uij = 0; uij < uiend; uij++) {
			free(matrix[uik][uij]);
		}
		free(matrix[uik]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

#endif

int free3DMatrix_OLD(unsigned int ***matrix, unsigned int lengthX, unsigned int lengthY, unsigned int lengthZ) {
	unsigned int x, y;

	for (x = 0; x < lengthX; x++) {
		for (y = 0; y < lengthY; y++) {
			free(matrix[x][y]);
		}
		free(matrix[x]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

int free3DMatrixC_OLD(char ***matrix, unsigned int lengthX, unsigned int lengthY, unsigned int lengthZ) {
	unsigned int x, y;

	for (x = 0; x < lengthX; x++) {
		for (y = 0; y < lengthY; y++) {
			free(matrix[x][y]);
		}
		free(matrix[x]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

ComputedEntries *newComputedEntriesMatrix(PatternStrandDirection *patternStrandDirection, bool bInsideOut) {
	int ii, iRegion;
	int iNumRegions = patternStrandDirection->patternStructures->iNumRegions;
	int im = patternStrandDirection->pattern->iLength;
	int in = im + patternStrandDirection->pattern->uiIndels;

	PatternRegion *patternRegion;
	ComputedEntries *compEntries  = (ComputedEntries *) calloc(iNumRegions, sizeof(ComputedEntries));

	if (bInsideOut) {
		patternRegion = patternStrandDirection->patternStructures->patternRegionInsideOut;
	} else {
		patternRegion = patternStrandDirection->patternStructures->patternRegion;
	}

	for (iRegion = 0; iRegion < iNumRegions; iRegion++) {
		if (patternRegion[iRegion].bIsUnpaired) {
			compEntries[iRegion].ijLast  = NULL;
			compEntries[iRegion].ijLastPerRowPerMatrix = (int **) calloc(im + 1, sizeof(int *));
			for (ii = 1; ii <= im; ii++) {
				compEntries[iRegion].ijLastPerRowPerMatrix[ii] = (int *) calloc(in, sizeof(int));
			}
		} else {
			compEntries[iRegion].ijLastPerRowPerMatrix = NULL;
			compEntries[iRegion].ijLast = (int *) calloc(in, sizeof(int));
		}
		if (patternRegion[iRegion].bIsArc) {
			compEntries[iRegion].ikFirst = (int *) calloc(in, sizeof(int));
		}
	}

	return compEntries;
}

int *freeComputedEntriesMatrix(ComputedEntries *compEntries, PatternStrandDirection *patternStrandDirection, bool bInsideOut) {
	int ii, iRegion;
	int iNumRegions = patternStrandDirection->patternStructures->iNumRegions;
	int im = patternStrandDirection->pattern->iLength;

	PatternRegion *patternRegion;

	if (bInsideOut) {
		patternRegion = patternStrandDirection->patternStructures->patternRegionInsideOut;
	} else {
		patternRegion = patternStrandDirection->patternStructures->patternRegion;
	}

	for (iRegion = 0; iRegion < iNumRegions; iRegion++) {
		if (patternRegion[iRegion].bIsUnpaired) {
			for (ii = 1; ii <= im; ii++) {
				free(compEntries[iRegion].ijLastPerRowPerMatrix[ii]);
			}
			free(compEntries[iRegion].ijLastPerRowPerMatrix);
		} else {
			free(compEntries[iRegion].ijLast);
		}
		if (patternRegion[iRegion].bIsArc) {
			free(compEntries[iRegion].ikFirst);
		}
	}

	return EXIT_SUCCESS;
}

#define copyPatternRegion(patternRegion1, i1, patternRegion2, i2) \
		patternRegion1[i1].bArcAfterDot  = patternRegion2[i2].bArcAfterDot; \
		patternRegion1[i1].bDotBeforeArc = patternRegion2[i2].bDotBeforeArc; \
		patternRegion1[i1].bIsArc        = patternRegion2[i2].bIsArc; \
		patternRegion1[i1].bIsLeftBulge  = patternRegion2[i2].bIsLeftBulge; \
		patternRegion1[i1].bIsUnpaired      = patternRegion2[i2].bIsUnpaired; \
		patternRegion1[i1].iLineToMerge1 = patternRegion2[i2].iLineToMerge1; \
		patternRegion1[i1].iL            = patternRegion2[i2].iL; \
		patternRegion1[i1].iR            = patternRegion2[i2].iR; \

PatternStructures* processPattern(char *cu_str, unsigned int uim, unsigned int uiIndels) {
	int i, i2, i3;

	int iTop = -1;
	int iRegion = -1;

	int           *iStack  = (int *) calloc(uim, sizeof(int));
	PatternRegion *prStack = (PatternRegion *) calloc(uim, sizeof(PatternRegion));

	PatternStructures *patternStructures      = (PatternStructures *) malloc(sizeof(PatternStructures));;
	patternStructures->patternRegion          = (PatternRegion *) calloc(uim, sizeof(PatternRegion));
	patternStructures->patternRegionInsideOut = (PatternRegion *) calloc(uim, sizeof(PatternRegion));
	patternStructures->iRegionRegion          = (int *) calloc(uim, sizeof(int));
	patternStructures->iRegionRegionInv       = (int *) calloc(uim, sizeof(int));
	patternStructures->iRegion                = (int *) calloc(uim, sizeof(int));
	patternStructures->iPosRegion             = (int *) calloc(uim + 1, sizeof(int));
	patternStructures->iPosRegionInsideOut    = (int *) calloc(uim + 1, sizeof(int));
	patternStructures->iLineIndex             = (int *) calloc(uim + 1, sizeof(int));
	patternStructures->iLeftMostPos           = (int *) calloc(uim + 1, sizeof(int));

	PatternRegion *patternRegion          = patternStructures->patternRegion;
	PatternRegion *patternRegionInsideOut = patternStructures->patternRegionInsideOut;

	int *iLineIndex = patternStructures->iLineIndex;

	bool bJustOpened = true;

	iStack[iTop = 0] = 0;
	for (i = 0; i < uim; i++) {
		if (cu_str[i] == '(') {
			patternStructures->iLeftMostPos[i + 1] = i + 1;

			//push
			iStack[++iTop] = i + 1;
		} else if (cu_str[i] == ')') {
			if (iStack[iTop] - 2 >= 0 && (cu_str[iStack[iTop] - 2] == '.' || cu_str[iStack[iTop] - 2] == ')')) {
				patternStructures->iLeftMostPos[i + 1] = patternStructures->iLeftMostPos[iStack[iTop] - 1];
			} else {
				patternStructures->iLeftMostPos[i + 1] = iStack[iTop];
			}

			//pop
			iTop--;
		} else {
			do {
				patternStructures->iLeftMostPos[i + 1] = iStack[iTop] + 1;
			} while (++i < uim && cu_str[i] == '.');
			i--;
		}
	}

	iTop = -1;
	iRegion = -1;

	for (i = 0; i < uim; i++) {
		if (cu_str[i] == '(') {
			//push
			iStack[++iTop] = i + 1;
			bJustOpened = true;
		} else if (cu_str[i] == ')') {
			//pop
			++iRegion;
			patternRegion[iRegion].iL = iStack[iTop];
			patternRegion[iRegion].iR = i + 1;
			patternRegion[iRegion].bIsArc = true;
			patternRegion[iRegion].bIsUnpaired = false;

			if (patternRegion[iRegion].iL > 1 && cu_str[patternRegion[iRegion].iL - 2] != '(') {
				patternRegion[iRegion].bDotBeforeArc = true;

				++iRegion;

				patternRegion[iRegion].iL = patternStructures->iLeftMostPos[iStack[iTop] - 1];
				patternRegion[iRegion].iR = i + 1;
				patternRegion[iRegion].bIsArc = false;
				patternRegion[iRegion].bIsUnpaired = false;
				patternRegion[iRegion].iLineToMerge1 = patternRegion[iRegion - 1].iL - 1;

			} else {
				patternRegion[iRegion].bDotBeforeArc = false;
			}

			iTop--;
			bJustOpened = false;
		} else {

			iRegion++;
			i2 = i;
			do {

				if (cu_str[i2 + 1] != '.') {
					patternRegion[iRegion].iL = i + 1;
					patternRegion[iRegion].iR = i2 + 1;
					patternRegion[iRegion].bIsArc = false;
					patternRegion[iRegion].bIsUnpaired = true;
					patternRegion[iRegion].bIsLeftBulge = bJustOpened;
					break;
				}
				i2++;
			} while (i2 < uim);
			i = i2;

			if (i + 1 < uim && cu_str[i2 + 1] == '(')
				patternRegion[iRegion].bArcAfterDot = true;
			else
				patternRegion[iRegion].bArcAfterDot = false;
		}
	}
	patternStructures->iNumRegions = iRegion + 1;

	iTop = -1;
	iRegion = -1;
	for (i = 0; i < patternStructures->iNumRegions; i++) {
		if (!patternRegion[i].bIsArc) {
			if (patternRegion[i].bIsUnpaired) {
				if (i > 0 && (patternRegion[i - 1].bIsArc || (!patternRegion[i - 1].bIsArc && !patternRegion[i - 1].bIsUnpaired))) {
					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, patternRegion, i)
				} else {
					++iTop;
					copyPatternRegion(prStack, iTop, patternRegion, i)
				}
			} else {
				++iRegion;
				copyPatternRegion(patternRegionInsideOut, iRegion, patternRegion, i)
			}
		} else {
			if (iTop >= 0) {
				if (prStack[iTop].iR == patternRegion[i].iR - 1) {
					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, prStack, iTop)
					iTop--;

					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, patternRegion, i)

				} else if (prStack[iTop].iL == patternRegion[i].iR + 1 || prStack[iTop].iR == patternRegion[i].iL - 1) {
					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, patternRegion, i)

					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, prStack, iTop)
					iTop--;
				} else {
					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, patternRegion, i)
				}

				if (iTop >= 0 && !prStack[iTop].bIsArc && prStack[iTop].iR == patternRegion[i].iL - 1) {
					++iRegion;
					copyPatternRegion(patternRegionInsideOut, iRegion, prStack, iTop)
					iTop--;
				}
			} else {
				++iRegion;
				copyPatternRegion(patternRegionInsideOut, iRegion, patternRegion, i)
			}
		}
	}
	if (iTop == 0) {
		iRegion++;
		copyPatternRegion(patternRegionInsideOut, iRegion, prStack, iTop)
		iTop--;
	}

	for (i = 0; i < patternStructures->iNumRegions; i++) {
		if (patternRegion[i].bIsArc) {
			patternStructures->iPosRegion[patternRegion[i].iL] = i;
			patternStructures->iPosRegion[patternRegion[i].iR] = i;
		} else if (patternRegion[i].bIsUnpaired) {
			for (i2 = patternRegion[i].iL; i2 <= patternRegion[i].iR; i2++) {
				patternStructures->iPosRegion[i2] = i;
			}
		}

		if (patternRegionInsideOut[i].bIsArc) {
			patternStructures->iPosRegionInsideOut[patternRegionInsideOut[i].iL] = i;
			patternStructures->iPosRegionInsideOut[patternRegionInsideOut[i].iR] = i;
		} else if (patternRegionInsideOut[i].bIsUnpaired) {
			for (i2 = patternRegionInsideOut[i].iL; i2 <= patternRegionInsideOut[i].iR; i2++) {
				patternStructures->iPosRegionInsideOut[i2] = i;
			}
		}
	}
	for (i = 0; i < patternStructures->iNumRegions; i++) {
		if (!patternRegionInsideOut[i].bIsArc && !patternRegionInsideOut[i].bIsUnpaired) {
			for (i2 = 0; i2 < patternStructures->iNumRegions; i2++) {
				if (patternRegionInsideOut[i].iL == patternRegion[i2].iL && patternRegionInsideOut[i].iR == patternRegion[i2].iR) {
					patternStructures->iRegionRegion[i] = i2;
				}
			}
		} else {
			patternStructures->iRegionRegion[i] = patternStructures->iPosRegion[patternRegionInsideOut[i].iL];
		}
		patternStructures->iRegion[i] = i;
	}
	for (i = 0; i < patternStructures->iNumRegions; i++) {
		patternStructures->iRegionRegionInv[patternStructures->iRegionRegion[i]] = i;
	}

	// Set iLineIndex
	for (i = 0; i <= uim; i++) {
		iLineIndex[i] = i;
	}
	for (i = 1; i < uim; i++) {
		if (cu_str[i - 1] == '(') {
			iLineIndex[i] = 0;
		} else if (cu_str[i - 1] == ')' &&
				patternStructures->patternRegion[patternStructures->iPosRegion[i]].bDotBeforeArc) {
			iLineIndex[i] = patternStructures->patternRegion[patternStructures->iPosRegion[i]].iL;
		}
	}
	if (cu_str[uim - 1] == ')' && patternStructures->patternRegion[patternStructures->iPosRegion[uim]].iL > 1)
		iLineIndex[uim] = patternStructures->patternRegion[patternStructures->iPosRegion[uim]].iL;

	for (i = 0; i < patternStructures->iNumRegions; i++) {
		if (!patternRegion[i].bIsArc && !patternRegion[i].bIsUnpaired) {
			patternStructures->iPosRegion[iLineIndex[patternRegion[i].iR]] = i;
		}
		if (!patternRegionInsideOut[i].bIsArc && !patternRegionInsideOut[i].bIsArc) {
			patternStructures->iPosRegionInsideOut[iLineIndex[patternRegionInsideOut[i].iR]] = i;
		}
	}

	free(iStack);

	for (i = 0; i < patternStructures->iNumRegions; i++) {
		patternRegion[i].iRegionDependsOnRegions = (int *) calloc(patternStructures->iNumRegions + 1, sizeof(int));
		iTop = -1;
		i2 = i;
		do {
			if (patternRegion[i2].bIsUnpaired) {
				if (patternRegion[i2].iL == 1 || cu_str[patternRegion[i2].iL - 2] == '(') {
					break;
				} else {
					iTop++;
					if (!patternRegion[i2 - 1].bIsUnpaired && !patternRegion[i2 - 1].bIsArc) {
						patternRegion[i].iRegionDependsOnRegions[iTop] = --i2;
					} else {
						patternRegion[i].iRegionDependsOnRegions[iTop] = i2 = patternStructures->iPosRegion[patternRegion[i2].iL - 1];
					}
				}
			} else if (patternRegion[i2].bIsArc) {
				iTop++;
				if ((i2 > 0) && !patternRegion[i2 - 1].bIsUnpaired && !patternRegion[i2 - 1].bIsArc) {
					if (patternRegion[i2 - 1].iL > patternRegion[i2].iL && patternRegion[i2 - 1].iR < patternRegion[i2].iR) //outer arc of a multi-loop, i2-1 inner region
						patternRegion[i].iRegionDependsOnRegions[iTop] = --i2;
					else {
						iTop--;
						break;
					}
				} else {
					if ((patternRegion[i2].iR - patternRegion[i2].iL) > 1) {
						patternRegion[i].iRegionDependsOnRegions[iTop] = i2 = patternStructures->iPosRegion[patternRegion[i2].iR - 1];
					} else {
						iTop--;
						break;
					}
				}

			} else {
				iTop++;
				patternRegion[i].iRegionDependsOnRegions[iTop] = i2 - 1;

				i3 = 0;
				while (patternRegion[i2 - 1].iRegionDependsOnRegions[i3] != -1) {
					iTop++;
					patternRegion[i].iRegionDependsOnRegions[iTop] = patternRegion[i2 - 1].iRegionDependsOnRegions[i3];
					i3++;
				}

				i2 = patternStructures->iPosRegion[iLineIndex[patternRegion[i2].iLineToMerge1]];

				iTop++;
				patternRegion[i].iRegionDependsOnRegions[iTop] = i2;

				i3 = 0;
				while (patternRegion[i2].iRegionDependsOnRegions[i3] != -1) {
					iTop++;
					patternRegion[i].iRegionDependsOnRegions[iTop] = patternRegion[i2].iRegionDependsOnRegions[i3];
					i3++;
				}
				break;
			}
		} while (true);
		patternRegion[i].iRegionDependsOnRegions[iTop + 1] = -1;
	}

	int itmp, *ptmp;
	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
		iTop = 0;
		while (patternRegion[iRegion].iRegionDependsOnRegions[iTop] != -1) {
			iTop++;
		}
		i2 = iTop / 2;
		for(i3 = 0; i3 < i2; i3++) {
			itmp = patternRegion[iRegion].iRegionDependsOnRegions[i3];
			patternRegion[iRegion].iRegionDependsOnRegions[i3] = patternRegion[iRegion].iRegionDependsOnRegions[iTop - 1 - i3];
			patternRegion[iRegion].iRegionDependsOnRegions[iTop - 1 - i3] = itmp;
		}
	}

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
		patternRegionInsideOut[iRegion].iRegionDependsOnRegions = (int *) calloc(patternStructures->iNumRegions + 1, sizeof(int));

		i = 0;
		ptmp = patternRegion[patternStructures->iRegionRegion[iRegion]].iRegionDependsOnRegions;
		while (ptmp[i] != -1) {
			patternRegionInsideOut[iRegion].iRegionDependsOnRegions[i] = patternStructures->iRegionRegionInv[ptmp[i]];
			i++;
		}
		patternRegionInsideOut[iRegion].iRegionDependsOnRegions[i] = -1;
	}

	patternStructures->iLargestRegionIndex          = (int **) calloc(patternStructures->iNumRegions, sizeof(int *));
	patternStructures->iLargestRegionIndexInsideOut = (int **) calloc(patternStructures->iNumRegions, sizeof(int *));
	int ii, ii2, ij, ij2;
	bool bNotFound;

	for (ii = 0; ii < patternStructures->iNumRegions; ii++) {
		patternStructures->iLargestRegionIndex[ii]          = (int *) calloc(patternStructures->iNumRegions, sizeof(int));
		patternStructures->iLargestRegionIndexInsideOut[ii] = (int *) calloc(patternStructures->iNumRegions, sizeof(int));
		for (ij = 0; ij < ii; ij++) {
			bNotFound = true;
			ii2 = ii;
			do {
				ij2 = 0;
				while (patternRegion[ii2].iRegionDependsOnRegions[ij2] != -1) {
					if (patternRegion[ii2].iRegionDependsOnRegions[ij2] == ij) {
						patternStructures->iLargestRegionIndex[ii][ij] = ii2;
						bNotFound = false;
						break;
					}
					ij2++;
				}
				ii2--;
			} while (bNotFound && ii2 >= 0);
			if (bNotFound && ii2 < 0) patternStructures->iLargestRegionIndex[ii][ij] = ij;
		}
		patternStructures->iLargestRegionIndex[ii][ii] = ii;
		for (ij = 0; ij < ii; ij++) {
			bNotFound = true;
			ii2 = ii;
			do {
				ij2 = 0;
				while (patternRegionInsideOut[ii2].iRegionDependsOnRegions[ij2] != -1) {
					if (patternRegionInsideOut[ii2].iRegionDependsOnRegions[ij2] == ij) {
						patternStructures->iLargestRegionIndexInsideOut[ii][ij] = ii2;
						bNotFound = false;
						break;
					}
					ij2++;
				}
				ii2--;
			} while (bNotFound && ii2 >= 0);
			if (bNotFound && ii2 < 0) patternStructures->iLargestRegionIndexInsideOut[ii][ij] = ij;
		}
		patternStructures->iLargestRegionIndexInsideOut[ii][ii] = ii;
	}


	int ik, iL, iR, iLeftMinus1, icol, icol2, icols, ikEnd, icount;
	int in = uim + uiIndels;
	int iLineToMerge1, iLineToMerge2;
	PatternRegion *patternRegionZero;
	for (icount = 0; icount < 2; icount++) {
		if (icount == 0)
			patternRegionZero = patternStructures->patternRegion;
		else
			patternRegionZero = patternStructures->patternRegionInsideOut;

		int ikBigStart = in;
		int ikBigEnd   = -1;

		for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
			patternRegion = &patternRegionZero[iRegion];
			iL = patternRegion->iL;
			iR = patternRegion->iR;
			patternRegion->ijGoal     = (int *) calloc(in, sizeof(int));
			patternRegion->ijBigGoal  = (int *) calloc(in, sizeof(int));
			if (patternRegion->bIsUnpaired) {
				patternRegion->ijGoalPerRowPerMatrix = (int **) calloc(uim + 1, sizeof(int*));
				for (ii = 0; ii <= uim; ii++) {
					patternRegion->ijGoalPerRowPerMatrix[ii] = (int *) calloc(in, sizeof(int));
				}
			}

			if (patternRegion->bIsUnpaired) {

				iLeftMinus1            = patternStructures->iLeftMostPos[iL] - 1;
				patternRegion->ikStart = iLeftMinus1 > uiIndels ? iLeftMinus1 - uiIndels : 0;
				patternRegion->ikEnd   = iLeftMinus1 + uiIndels < in ? iLeftMinus1 + uiIndels : in - 1;
				icol                   = iR - iLeftMinus1; //mid col of last line

				if (patternRegion->iRegionDependsOnRegions[0] != -1)
					ikEnd = iR + uiIndels - 1;
				else
					ikEnd = patternRegion->ikEnd;

				patternRegion->ikEndFactual = ikEnd;

				if (ikBigStart > patternRegion->ikStart)
					ikBigStart = patternRegion->ikBigStart = patternRegion->ikStart;
				else
					patternRegion->ikBigStart = ikBigStart;

				if (ikBigEnd < patternRegion->ikEndFactual)
					ikBigEnd = patternRegion->ikBigEnd = patternRegion->ikEndFactual;
				else
					patternRegion->ikBigEnd = ikBigEnd;

				for (ii = iL; ii <= iR; ii++) {
					icol2 = ii - iLeftMinus1;

					for (ik = patternRegion->ikStart; ik <= ikEnd; ik++) {
						icols = uiIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

						if ((patternRegion->ijGoal[ik] = patternRegion->ijGoalPerRowPerMatrix[ii][ik] = icol2 + icols) > in - ik) {
							patternRegion->ijGoalPerRowPerMatrix[ii][ik] = in - ik;
							patternRegion->ijGoal[ik] = in - ik;
						}
						if (patternRegion->ijGoalPerRowPerMatrix[ii][ik] < 0) {
							patternRegion->ijGoalPerRowPerMatrix[ii][ik] = 0;
							patternRegion->ijGoal[ik] = 0;
						}
					}
				}

			} else if (patternRegion->bIsArc) {

				iLeftMinus1            = iL - 1;
				patternRegion->ikStart = iLeftMinus1 > uiIndels ? iLeftMinus1 - uiIndels : 0;
				patternRegion->ikEnd   = iLeftMinus1 + uiIndels < in ? iLeftMinus1 + uiIndels : in;
				icol                   = iR - iLeftMinus1; //mid col of last line

				patternRegion->ikEndFactual = iR + uiIndels - 1;
				if (patternRegion->ikEndFactual >= in)
					patternRegion->ikEndFactual = in - 1;

				if (ikBigStart > patternRegion->ikStart)
					ikBigStart = patternRegion->ikBigStart = patternRegion->ikStart;
				else
					patternRegion->ikBigStart = ikBigStart;

				if (ikBigEnd < patternRegion->ikEndFactual)
					ikBigEnd = patternRegion->ikBigEnd = patternRegion->ikEndFactual;
				else
					patternRegion->ikBigEnd = ikBigEnd;

				for (ik = patternRegion->ikStart; ik < in; ik++) {
					icols = uiIndels - abs(iLeftMinus1 - ik); //number of columns to be checked - 1

					if ((patternRegion->ijGoal[ik] = icol + icols) > in - ik)
						patternRegion->ijGoal[ik] = in - ik;
					if (patternRegion->ijGoal[ik] < 0)
						patternRegion->ijGoal[ik] = 0;
				}

			} else {

				iLineToMerge1 = patternRegion->iLineToMerge1;
				iLineToMerge2 = iR;

				iLeftMinus1 = patternStructures->iLeftMostPos[iLineToMerge1 > 0 ? iLineToMerge1 : 0] - 1;
				if (iLineToMerge1 == 0)
					iLeftMinus1++;
				patternRegion->ikStart = iLeftMinus1 > uiIndels ? iLeftMinus1 - uiIndels : 0;
				patternRegion->ikEnd   = iLeftMinus1 + uiIndels < in ? iLeftMinus1 + uiIndels : in - 1;
				icol = iLineToMerge2 - iLeftMinus1; //mid col for line

				patternRegion->ikEndFactual = iR + uiIndels - 1;
				if (patternRegion->ikEndFactual >= in)
					patternRegion->ikEndFactual = in - 1;

				if (ikBigStart > patternRegion->ikStart)
					ikBigStart = patternRegion->ikBigStart = patternRegion->ikStart;
				else
					patternRegion->ikBigStart = ikBigStart;

				if (ikBigEnd < patternRegion->ikEndFactual)
					ikBigEnd = patternRegion->ikBigEnd = patternRegion->ikEndFactual;
				else
					patternRegion->ikBigEnd = ikBigEnd;

				for (ik = patternRegion->ikStart; ik < in; ik++) {
					icols = uiIndels - abs(iLeftMinus1 - ik); // number of columns to be checked - 1

					if ((patternRegion->ijGoal[ik] = icol + icols) > in - ik)
						patternRegion->ijGoal[ik] = in - ik; // end of interval to be computed
					if (patternRegion->ijGoal[ik] < 0)
						patternRegion->ijGoal[ik] = 0;
				}

			}
		}
	}

	for (icount = 0; icount < 2; icount++) {
		if (icount == 0)
			patternRegionZero = patternStructures->patternRegion;
		else
			patternRegionZero = patternStructures->patternRegionInsideOut;

		int *ijBigGoal = (int *) calloc(in, sizeof(int));

		for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
			patternRegion = &patternRegionZero[iRegion];

			for (ik = 0; ik < in; ik++) {
				if (ijBigGoal[ik] < patternRegion->ijGoal[ik])
					ijBigGoal[ik] = patternRegion->ijGoal[ik];
				patternRegion->ijBigGoal[ik] = ijBigGoal[ik];
			}
		}
		free(ijBigGoal);
	}

	patternStructures->iLCC = (int **) calloc(uim + 1, sizeof(int *));
	for (i = 0; i <= uim; i++) {
		patternStructures->iLCC[i] = (int *) calloc(uim + 1, sizeof(int));
	}

	return patternStructures;
}

void initializeRowsAndColumnsRegions_S(unsigned int ***edist, char ***operation, unsigned int ***trace, Cost *cost,
		char *cu_str, int im, int in, PatternStructures *patternStructures) {
	int ii, ij, ik/*, iCountArcs*/, iend;
	int iL, iR, iRegion;

	int *iLineIndex = patternStructures->iLineIndex;
	PatternRegion *patternRegion /*= patternStructures->patternRegion*/;

	// Initialize first rows
	S(edist, 0, 0, 0) = 0;
	for (ij = 1; ij <= in; ij++) {
		S(edist, 0, 0, ij)     = S(edist, 0, 0, ij - 1) + cost->iCostDeletion;
		S(operation, 0, 0, ij) = OpInsertion; //OpDeletion;
	}
	for (ik = 1; ik <= in; ik++) {
		iend = in - ik;
		for (ij = 1; ij <= iend; ij++) {
			S(edist, ik, 0, ij) = S(edist, 0, 0, ij);
			S(operation, ik, 0, ij) = S(operation, 0, 0, ij);
		}
	}
	// End Initialize first rows

	//iCountArcs = 0;

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
		patternRegion = &patternStructures->patternRegion[iRegion];
		iL = patternRegion->iL;
		iR = patternRegion->iR;

		//printf("iL=%d iR=%d\n", iL, iR);

		if (patternRegion->bIsUnpaired) {

			for (ii = iL; ii <= iR; ii++) {
				S(edist, 0, ii, 0)     = S(edist, 0, iLineIndex[ii - 1], 0) + cost->iCostDeletion;
				S(operation, 0, ii, 0) = OpDeletion;
			}

		} else if (patternRegion->bIsArc) {

			//iCountArcs++;

			//P(edist, 0, iR, 0)     = iCountArcs * cost->iCostArcRemoving + iR + 1 - iL - 2 * iCountArcs;
			S(edist, 0, iR, 0)     = S(edist, 0, iLineIndex[iR - 1], 0) + cost->iCostArcRemoving;
			S(operation, 0, iR, 0) = OpArcRemoving;



			for (ik = 0; ik <= in; ik++) {
				iend = in - ik;
				for (ij = 1; ij <= iend; ij++) {
					S(operation, ik, iL, ij) = OpInsertion;
				}
			}



		} else {

			S(edist, 0, iLineIndex[iR], 0) = S(edist, 0, iLineIndex[patternRegion->iLineToMerge1], 0) + S(edist, 0, iR, 0);
			S(trace, 0, iR, 0)             = INT_MAX;

		}

	}

	for (ik = 1; ik <= in; ik++) {
		for (ii = 0; ii <= im; ii++) {
			S(edist, ik, ii, 0)     = S(edist, 0, ii, 0);
			S(operation, ik, ii, 0) = S(operation, 0, ii, 0);
		}
	}

}

void initializeRowsAndColumns(unsigned int ***edist, char ***operation, unsigned int ***trace, Cost *cost,
		char *cu_str, int im, int in, PatternStructures *patternStructures) {
	int ii, ij, ik/*, iCountArcs, iend*/;
	int iL, iR, iRegion;

	int *iLineIndex = patternStructures->iLineIndex;
	PatternRegion *patternRegion;

	// Initialize first rows
	P(edist, 0, 0, 0) = 0;
	for (ij = 1; ij <= in; ij++) {
		P(edist, 0, 0, ij)     = P(edist, 0, 0, ij - 1) + cost->iCostDeletion;
		P(operation, 0, 0, ij) = OpInsertion; //OpDeletion;
	}
	for (ik = 1; ik <= in; ik++) {
		memcpy(P1(edist, ik, 0), P1(edist, 0, 0), sizeof(unsigned int) * (in - ik + 1));
		memcpy(P1(operation, ik, 0), P1(operation, 0, 0), sizeof(char) * (in - ik + 1));
	}
	// End Initialize first rows

	//iCountArcs = 0;

	for (iRegion = 0; iRegion < patternStructures->iNumRegions; iRegion++) {
		patternRegion = &patternStructures->patternRegion[iRegion];
		iL = patternRegion->iL;
		iR = patternRegion->iR;

		if (patternRegion->bIsUnpaired) {

			for (ii = iL; ii <= iR; ii++) {
				P(edist, 0, ii, 0)     = P(edist, 0, iLineIndex[ii - 1], 0) + cost->iCostDeletion;
				P(operation, 0, ii, 0) = OpDeletion;
			}

		} else if (patternRegion->bIsArc) {

			//iCountArcs++;

			//P(edist, 0, iR, 0)     = iCountArcs * cost->iCostArcRemoving + iR + 1 - iL - 2 * iCountArcs;
			P(edist, 0, iR, 0)     = P(edist, 0, iLineIndex[iR - 1], 0) + cost->iCostArcRemoving;
			P(operation, 0, iR, 0) = OpArcRemoving;

			for (ik = 0; ik <= in; ik++) {
				memcpy(P1(operation, ik, iL), P1(operation, 0, 0), sizeof(char) * (in - ik + 1));
			}
			/*for (ik = 0; ik <= in; ik++) {
				iend = in - ik;
				for (ij = 1; ij <= iend; ij++) {
					P(operation, ik, iL, ij) = OpInsertion;
				}
			}*/

		} else {

			P(edist, 0, iLineIndex[iR], 0) = P(edist, 0, iLineIndex[patternRegion->iLineToMerge1], 0) + P(edist, 0, iR, 0);
			P(trace, 0, iR, 0)             = INT_MAX;

		}

	}

	for (ik = 1; ik <= in; ik++) {
		for (ii = 0; ii <= im; ii++) {
			P(edist, ik, ii, 0)     = P(edist, 0, ii, 0);
			P(operation, ik, ii, 0) = P(operation, 0, ii, 0);
		}
	}

}


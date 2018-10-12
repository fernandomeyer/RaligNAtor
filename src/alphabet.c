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
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "align.h"

void setPredefinedAlphabet(int type, Alphabet *alphabet) {
	int i;
	char *cClasses[] = {"Aa", "Cc", "Gg", "UuTt",
			"RAG",
			"YCU",
			"SCG",
			"WAU",
			"KGU",
			"MAC",
			"BCGU",
			"DAGU",
			"HACU",
			"VACG",
			"NACGU",
			"rAG",
			"yCU",
			"sCG",
			"wAU",
			"kGU",
			"mAC",
			"bCGU",
			"dAGU",
			"hACU",
			"vACG",
			"nACGU",
			"NnRrYySsWwKkMmBbDdHhVv"};

	/********* DNA or RNA *********/
	if(type == 0 || type == 1) {
		alphabet->numClasses = 27;
		if ((alphabet->eqClass    = (char **) calloc(alphabet->numClasses, sizeof(char*))) == NULL) {
			fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((alphabet->classSize  = (int *) calloc(alphabet->numClasses, sizeof(int))) == NULL) {
			fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((alphabet->classRepresentative = (unsigned char *) calloc(alphabet->numClasses, sizeof(unsigned char))) == NULL) {
			fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((alphabet->isWildCard = (bool *) calloc(alphabet->numClasses, sizeof(bool))) == NULL) {
			fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		for(i = 0; i < alphabet->numClasses; i++) {
			alphabet->classSize[i] = strlen(cClasses[i]);

			if (i > 3)
				alphabet->isWildCard[i] = true;
			else
				alphabet->isWildCard[i] = false;

			if ((alphabet->eqClass[i] = (char *) calloc(alphabet->classSize[i] + 1, sizeof(char))) == NULL) {
				fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}

			strcpy(alphabet->eqClass[i], cClasses[i]);

			alphabet->classRepresentative[i] = cClasses[i][0];
		}
		if (type == 0)
			alphabet->classRepresentative[3] = 'T';

		alphabet->isWildCard[alphabet->numClasses - 1] = false;
	}
}

int setAlphabetConversionTable(Alphabet* alphabet) {
	int i, j;
	unsigned char *tablePattern, *tableTargetSeq;

	if (alphabet == NULL) {
		fprintf(stderr, "Error. Alphabet was not initialized. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	if((tablePattern = (unsigned char *) calloc(UCHAR_MAX, sizeof(unsigned char))) == NULL) {
		fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	if((tableTargetSeq = (unsigned char *) calloc(UCHAR_MAX, sizeof(unsigned char))) == NULL) {
		fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	for (i = 0; i < UCHAR_MAX; i++) {
		tablePattern[i]   = 0;
		tableTargetSeq[i] = 0;
	}

	//Check if the classes' representative is unique, except the last class (which can have any representative)
	for (i = 0; i < alphabet->numClasses - 1; i++) {
		for (j = 0; j < alphabet->numClasses - 1; j++) {
			if (i != j && alphabet->classRepresentative[i] == alphabet->classRepresentative[j]) {
				fprintf(stderr,"Character \"%c\" can represent only one character class.\n", alphabet->eqClass[i][j]);
				return 1;
			}
		}
	}

	for (i = 0; i < alphabet->numClasses - 1; i++) {
		if (!alphabet->isWildCard[i]) {
			for(j = 0; j < alphabet->classSize[i]; j++) {
				tablePattern[(int) alphabet->eqClass[i][j]] = i + 1; // Characters in the 1st eqClass receive value (code) 1, in the 2nd eqClass value 2...
			}
		}
	}
	for (i = 0; i < alphabet->numClasses - 1; i++) {
		if (alphabet->isWildCard[i]) {
			for(j = 0; j < alphabet->classSize[i]; j++) {
				if (tablePattern[(int) alphabet->eqClass[i][j]] == 0 && alphabet->eqClass[i][j] != alphabet->classRepresentative[i]) {
					fprintf(stderr,"Character \"%c\" belongs to a wildcard character class, therefore it must also be defined in a non-wildcard class.\n", alphabet->eqClass[i][j]);
					return 1;
				}
				if (tablePattern[alphabet->classRepresentative[i]] != 0) {
					fprintf(stderr,"Character \"%c\" cannot simultaneously belong to a non-wildcard character class and be the representative of a wildcard class.\n", alphabet->classRepresentative[i]);
					return 1;
				}
			}
			tablePattern[(int) alphabet->classRepresentative[i]] = i + 1;
		}
	}

	for (i = 0; i < alphabet->numClasses; i++) {
		if (!alphabet->isWildCard[i]) {
			for(j = 0; j < alphabet->classSize[i]; j++) {
				tableTargetSeq[(int) alphabet->eqClass[i][j]] = i + 1;
			}
		}
	}

	tablePattern[$] = tableTargetSeq[$] = $;

	alphabet->conversionTablePattern   = tablePattern;
	alphabet->conversionTableTargetSeq = tableTargetSeq;

	return EXIT_SUCCESS;
}

bool convertToAlphabet(char *seq, char *convSeq, int length, bool isPattern, Alphabet* alphabet) {
	int iS;
	unsigned char *tablePattern, *tableTargetSeq;
	unsigned char ucCodePattern, ucCodeTargetSeq;

	if (alphabet == NULL) {
		fprintf(stderr, "Error. Alphabet was not initialized. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	if (alphabet->conversionTablePattern == NULL) {
		fprintf(stderr, "Error. \"alphabet->conversionTable\" was not initialized. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	tablePattern   = alphabet->conversionTablePattern;
	tableTargetSeq = alphabet->conversionTableTargetSeq;

	for(iS = 0; iS < length; iS++) {
		if (isPattern) {
			ucCodePattern   = tablePattern[(unsigned char) seq[iS]];

			if (ucCodePattern == 0) {
				fprintf(stderr,"Character \"%c\" is not defined in the alphabet for patterns.\n", seq[iS]);
				return 1;
			}
			convSeq[iS] = ucCodePattern;
		} else {
			ucCodeTargetSeq = tableTargetSeq[(unsigned char) seq[iS]];

			if (ucCodeTargetSeq == 0) {
				fprintf(stderr,"Character \"%c\" is not defined in the alphabet for the target sequences.\n", seq[iS]);
				return 1;
			}
			convSeq[iS] = ucCodeTargetSeq;
		}
	}

	return EXIT_SUCCESS;
}

/* Converts each value of an alphabetically transformed array to the
 * value of the representative character of the corresponding class. */
bool convertToRepresentative(char *seq, int length, Alphabet *alphabet) {
	int i;

	for(i = 0; i < length; i++) {
		if(seq[i] == $)
			seq[i] = $;
		else
			seq[i] = alphabet->classRepresentative[seq[i] - 1];
	}

	return 0;
}

bool** setIupacTable (Alphabet *alphabet) {
	bool **table;
	unsigned int i, j;
	unsigned char *conversionTable;
	unsigned char uiClassCode;

	if (alphabet == NULL) {
		fprintf(stderr, "Error. Alphabet was not initialized - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	if (alphabet->conversionTablePattern == NULL) {
		fprintf(stderr,"\"alphabet->conversionTable\" not set. - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	if((table = (bool **) calloc(alphabet->numClasses + 1, sizeof(bool*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"iupacTable\". - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	alphabet->iupacTable = table;
	conversionTable = alphabet->conversionTablePattern;

	for (i = 0; i <= alphabet->numClasses; i++) {
		if ((table[i] = (bool *) calloc(alphabet->numClasses + 1, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"iupacTable\". - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	}
	for (i = 1; i <= alphabet->numClasses; i++) {
		if (!alphabet->isWildCard[i - 1]) {
			table[i][i] = true;
		} else {
			for (j = 0; j < alphabet->classSize[i - 1]; j++) {
				uiClassCode = conversionTable[(unsigned int) alphabet->eqClass[i - 1][j]];

				if (!alphabet->isWildCard[uiClassCode - 1]) {
					table[i][uiClassCode] = true;
				}

			}
		}
	}

	return table;
}

bool** loadReverseComplementarityRules(bool **bCompCheck, char *cArrayComplementaryCode, AffixArray *affixArray) {
	int i, j;
	bool **bRevCompCheck;
	Alphabet *alphabet = affixArray->alphabet;

	if ((bRevCompCheck = (bool **) calloc(alphabet->numClasses, sizeof(bool*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"bRevCompCheck\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	for (i = 0; i < alphabet->numClasses; i++) {
		if ((bRevCompCheck[i] = (bool *) calloc(alphabet->numClasses, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"bRevCompCheck_i\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	}

	for (i = 1; i < affixArray->alphabet->numClasses; i++) {
		for (j = 1; j < affixArray->alphabet->numClasses; j++) {
			bRevCompCheck[cArrayComplementaryCode[i] - 1][cArrayComplementaryCode[j] - 1] = bCompCheck[i - 1][j - 1];
		}
	}

	return bRevCompCheck;
}

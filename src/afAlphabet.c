/*Copyright (C) 2010  Fernando Meyer

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

#include "af.h"

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

bool convertToAlphabet(unsigned char *seq, unsigned char *convSeq, int length, bool isPattern, Alphabet* alphabet) {
	int i, j, iS;
	unsigned char *tablePattern, *tableTargetSeq;
	unsigned char ucCodeTargetSeq;

	if (alphabet == NULL) {
		fprintf(stderr, "Error. Alphabet was not initialized. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	if((tablePattern = (unsigned char *) calloc(256, sizeof(unsigned char))) == NULL) {
		fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	if((tableTargetSeq = (unsigned char *) calloc(256, sizeof(unsigned char))) == NULL) {
		fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	for (i = 0; i < 256; i++) {
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

	for(iS = 0; iS < length; iS++) {
		ucCodeTargetSeq = tableTargetSeq[(unsigned char) seq[iS]];

		if(ucCodeTargetSeq == 0) {
			fprintf(stderr,"Character \"%c\" is not defined in the alphabet for the target sequence.\n", seq[iS]);
			return 1;
		}

		convSeq[iS] = ucCodeTargetSeq;
	}

	free(tablePattern);
	free(tableTargetSeq);

	return EXIT_SUCCESS;
}

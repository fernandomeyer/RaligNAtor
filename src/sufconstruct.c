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
#include <unistd.h>
#include <getopt.h>

#include <limits.h>

#include "af.h"

int main(int argc, char *argv[]) {
	AffixArray affixArray;
	Alphabet   alphabet;

	char *argFastaFile = NULL;
    char *argAlphabetFile  = NULL;
    int predefAlphabet = 0;
	bool argSuf     = 0;
	bool argLcp     = 0;
	bool argSufinv  = 0;

	alphabet.bPredefined = true;
	bool argSaveTransfSeq = 1;
	bool argShowTimes = 0;
	char *argIndex    = NULL;
	int argOScreen  = 0;

	int opt, longIndex = 0;
	extern int optind;
    static const char *optString = "z:drabs:xci";

	static const struct option longOpts[] = {
		{"alph", required_argument, NULL, 'z'},
		{"dna", no_argument, NULL, 'd'},
		{"rna", no_argument, NULL, 'r'},
		{"lesa", no_argument, NULL, 'a'},
		{"lgslink", no_argument, NULL, 'b'},

		{"s", required_argument, NULL, 's'},
		{"x", no_argument, NULL, 'x'},
		{"c", no_argument, NULL, 'c'},
		{"time", no_argument, NULL, 'i'},
		{NULL, no_argument, NULL, 0}
	};

	while((opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex)) != -1) {
	    switch(opt) {
			case 'z':
				argAlphabetFile = optarg;
				alphabet.bPredefined = false;
				break;
			case 'd':
				predefAlphabet = 0;
				break;
			case 'r':
				predefAlphabet = 1;
				break;
			case 'a':
				argSuf = 1;
				argLcp = 1;
				break;
			case 'b':
				argSuf    = 1;
				argLcp    = 1;
				argSufinv = 1;
				break;

			case 's':
				argIndex = optarg;
				break;
			case 'x':
				argSaveTransfSeq = 0;
				break;
			case 'c':
				argOScreen = 1;
				break;
			case 'i':
				argShowTimes = 1;
				break;
			case '?':
				exit(1);
		}
	}

    if(argc < optind || argc == 1) {
    	fprintf(stderr,"Usage: %s [OPTION]\n", argv[0]);
	    fprintf(stderr,"  <file>        Target FASTA file\n");
	    fprintf(stderr,"  -alph <file>	Use alphabet defined in file\n");
	    fprintf(stderr,"  -dna          Use DNA alphabet {A, C, G, T} and IUPAC wildcards (default)\n");
		fprintf(stderr,"  -rna          Use RNA alphabet {A, C, G, U} and IUPAC wildcards\n");
	    fprintf(stderr,"\n");
	    fprintf(stderr,"  -lesa         Construct index for LESAAlign (tables suf and lcp)\n");
	    fprintf(stderr,"  -lgslink      Construct index for LGSlinkAlign and LESAAlign (tables suf, lcp, and suf^-1)\n");
	    fprintf(stderr,"\n");
	    fprintf(stderr,"  -s <index>    Save constructed structures to given index name\n");
	    fprintf(stderr,"  -x            Do not save alphabetically transformed sequence\n");
	    fprintf(stderr,"  -c            Output constructed structures to screen\n");
	    fprintf(stderr,"  -time         Display elapsed times\n");
	    return 1;
    }
	if(optind < argc)
		argFastaFile = argv[optind];

    init(&affixArray, &alphabet, argAlphabetFile, predefAlphabet, &argSuf, &argLcp, &argSufinv, true);

    if(argFastaFile == NULL) {
		fprintf(stderr,"Please specify a FASTA file.\n");
		return 1;
    } else if(loadFastaFile(&affixArray, argFastaFile, true))
		return 1;

	if(affixArray.multiSeq == NULL)
		return 1;

	if(argFastaFile != NULL)
		printf("FASTA file:          %s\n", argFastaFile);
	printf("Number of sequences: %i\n", affixArray.multiSeq->numSeqs);
	printf("Total length:        %i\n", affixArray.length - affixArray.multiSeq->numSeqs);

	if(convertToAlphabet(affixArray.multiSeq->sequences,
						affixArray.multiSeq->convSequences,
						affixArray.length,
						false,
						&alphabet))
		return 1;

	if(argIndex != NULL && saveAffFiles(&affixArray, argIndex, argSaveTransfSeq))
		return 1;

	free(affixArray.multiSeq->sequences);
	affixArray.multiSeq->sequences = NULL;

	if(affixArray.multiSeq != NULL && affixArray.esa->xarray == NULL && argSuf) {
		computeAffixArray(&affixArray, argSuf, argLcp, argSufinv,
    						argOScreen, argIndex, argShowTimes);
	} else {
		printf("Please select the tables to be constructed.\n");
		exit(EXIT_FAILURE);
	}

	if(argOScreen) {
		printf("\n"); printEArray(&affixArray);
	}

	freeAll(&affixArray, NULL);

	return 0;
}

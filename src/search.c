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
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <pthread.h>

#include <unistd.h>

#include "align.h"
#include "alignThread.h"

int main(int argc, char *argv[]) {
	AffixArray   affixArray;
	Alphabet     alphabet;
	MultiPattern *multiPattern = NULL;
	bool         **complementarityRules = NULL;
	SearchParam  searchParam;
	Chainparam chainParam;
	searchParam.chainParam = &chainParam;

	int opt, longIndex = 0;

	int iNumThreads    = 0;
	int iNumPartitions = 0;

	bool argSuf   = false;
	bool argLcp   = false;
	bool argLnk   = false;
	bool argAflk  = false;
	bool argSufr  = false;
	bool argLcpr  = false;
	bool argAflkr = false;
	bool argSufinv = false;
	bool bIndexBased = false;

	int  iPredefAlphabet = 0;
	char *argAlphabetFile = NULL;
	char *argLoadFile = NULL;
	char *argStruct   = NULL;
	char *argCompFile = NULL;

	searchParam.uiThreshold          = -1;
	searchParam.uiIndels             = -1;
	searchParam.cVariant             = V_SCANLA;
	searchParam.bSearchForwardString = false;
	searchParam.bSearchReverseString = false;
	searchParam.bBED                 = false;
	searchParam.bReportAllMatches    = false;
	searchParam.bPrintOutMatches     = true;
	searchParam.bPrintMatchesBySeq   = false;
	searchParam.cPrintMatchesByScore = 0; //0=no sorting, 1=descending, 2=ascending
	searchParam.bPrintTable          = false;
	searchParam.uiMinDiffMatches     = 0;
	searchParam.bFilterOverlaps      = false;
	searchParam.bShowProgress        = false;
	searchParam.bUseSequenceBasedFilter = true; //Use filter for lgslink variant

	searchParam.cost.iCostArcAltering  = -1;
	searchParam.cost.iCostArcBreak     = -1;
	searchParam.cost.iCostArcRemoving  = -1;
	searchParam.cost.iCostDeletion     = -1;
	searchParam.cost.iCostReplacement  = -1;
	searchParam.cost.iFCostDeletion    = -1;
	searchParam.cost.iFCostReplacement = -1;

	chainParam.reportfile   = NULL;
	chainParam.isactive     = false;
	chainParam.islocal      = false;
	chainParam.globalparm   = NULL;
	chainParam.localparm    = NULL;
	chainParam.maxgap       = 0;
	chainParam.weightfactor = 1.0;
	chainParam.show         = false;
	chainParam.show2        = false;
	chainParam.minchainlength = 0;
	chainParam.minchainscore = 0;
	chainParam.topkscoring   = UINT_MAX;

	extern int optind;
	static const char *optString = "h:kpqs:fbc:n:0123456789zmglw:x:e:odiu:v:@!:%:^:&:*:=:a:()-[]t:j:_$:#:";

	static const struct option longOpts[] = {
			{"alph", required_argument, NULL, 'h'},
			{"dna", no_argument, NULL, 'k'},
			{"rna", no_argument, NULL, 'q'},
			//{"protein", no_argument, NULL, 'p'},
			{"pat", required_argument, NULL, 's'},
			{"for", no_argument, NULL, 'f'},
			{"rev", no_argument, NULL, 'b'},
			{"comp", required_argument, NULL, 'c'},
			{"match", required_argument, NULL, 'n'},
			{"bed", no_argument, NULL, '@'},
			{"allm", no_argument, NULL, 'm'},
			{"byseq", no_argument, NULL, 'd'},
			{"byscore", no_argument, NULL, ')'},
			{"byscorea", no_argument, NULL, '-'},
			{"table", no_argument, NULL, '('},
			{"silent", no_argument, NULL, 'i'},
			{"no-overlaps", no_argument, NULL, '['},
			//{"progress", no_argument, NULL, '+'},
			/* insert the option -thread to parse the number of threads*/
			{"threads", required_argument, NULL, 't'},
			{"parts", required_argument, NULL, 'j'},

			{"lesa", no_argument, NULL, V_LESA},
			{"lgslink", no_argument, NULL, V_LGSLINK},
			{"lgslink_nof", no_argument, NULL, '_'},
			{"scan", no_argument, NULL, V_SCAN},
			{"lscan", no_argument, NULL, V_SCANLA},
			{"aligngl", no_argument, NULL, V_GLOBAL},

			{"replacement", required_argument, NULL, '!'},
			{"freplacement", required_argument, NULL, '$'},
			{"deletion", required_argument, NULL, '%'},
			{"fdeletion", required_argument, NULL, '#'},
			{"arc-breaking", required_argument, NULL, '^'},
			{"arc-altering", required_argument, NULL, '&'},
			{"arc-removing", required_argument, NULL, '*'},
			{"cost", required_argument, NULL, '='},
			{"indels", required_argument, NULL, 'a'},

			{"global", no_argument, NULL, 'g'},
			{"local", no_argument, NULL, 'l'},
			{"allglobal", no_argument, NULL, 'y'},
			{"wf", required_argument, NULL, 'w'},
			{"maxgap", required_argument, NULL, 'x'},
			{"minscore", required_argument, NULL, 'u'},
			{"minlen", required_argument, NULL, 'e'},
			{"top", required_argument, NULL, 'v'},
			{"show", no_argument, NULL, 'o'},
			{"show2", no_argument, NULL, ']'},
			{NULL, no_argument, NULL, 0}
	};

	while((opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex)) != -1) {
		switch(opt) {
		case 'h':
			argAlphabetFile = optarg;
			break;
		case 'k':
			iPredefAlphabet = 0;
			break;
		case 'q':
			iPredefAlphabet = 1;
			break;
		case 's':
			argStruct = optarg;
			break;
		case 'f':
			searchParam.bSearchForwardString = true;
			break;
		case 'b':
			searchParam.bSearchReverseString = true;
			break;
		case 'c':
			argCompFile = optarg;
			break;
		case 'n':
			searchParam.uiMinDiffMatches = atoi(optarg);
			break;
		case '@':
			searchParam.bBED = true;
			break;
		case 'm':
			searchParam.bReportAllMatches = true;
			break;
		case 'd':
			searchParam.bPrintMatchesBySeq = true;
			break;
		case ')':
			searchParam.cPrintMatchesByScore = 1;
			break;
		case '-':
			searchParam.cPrintMatchesByScore = 2;
			break;
		case '(':
			searchParam.bPrintTable = true;
			break;
		case 'i':
			searchParam.bPrintOutMatches = false;
			break;
		case '[':
			searchParam.bFilterOverlaps = true;
			break;
		/*case for parsing number of thread and convert the input string to an integer*/
		case 't':
		        iNumThreads = atoi(optarg);
			break;
		case 'j':
				iNumPartitions = atoi(optarg);
			break;

		case '!':
			searchParam.cost.iCostReplacement = atoi(optarg);
			break;
		case '%':
			searchParam.cost.iCostDeletion = atoi(optarg);
			break;
		case '^':
			searchParam.cost.iCostArcBreak = atoi(optarg);
			break;
		case '&':
			searchParam.cost.iCostArcAltering = atoi(optarg);
			break;
		case '*':
			searchParam.cost.iCostArcRemoving = atoi(optarg);
			break;
		case '=':
			searchParam.uiThreshold = atoi(optarg);
			break;
		case 'a':
			searchParam.uiIndels = atoi(optarg);
			break;

		case '$':
			searchParam.cost.iFCostReplacement = atoi(optarg);
			break;
		case '#':
			searchParam.cost.iFCostDeletion = atoi(optarg);
			break;

		case V_LESA:
			searchParam.cVariant = V_LESA;
			bIndexBased = true;
			break;
		case V_LGSLINK:
			searchParam.cVariant = V_LGSLINK;
			bIndexBased = true;
			break;
		case '_':
			searchParam.cVariant = V_LGSLINK;
			bIndexBased = true;
			searchParam.bUseSequenceBasedFilter = false;
			break;
		case V_SCAN:
			searchParam.cVariant = V_SCAN;
			break;
		case V_SCANLA:
			searchParam.cVariant = V_SCANLA;
			break;
		case V_GLOBAL:
			searchParam.cVariant = V_GLOBAL;
			break;

		case 'g':
			chainParam.islocal  = false;
			chainParam.isactive = true;
			break;
		case 'l':
			chainParam.islocal  = true;
			chainParam.isactive = true;
			break;
		case 'y':
			chainParam.isactive = true;
			chainParam.globalparm = (char *) malloc(sizeof(char) * 10);
			strcpy(chainParam.globalparm, "all");
			break;
		case 'w':
			chainParam.weightfactor = strtod(optarg, NULL);
			break;
		case 'x':
			chainParam.maxgap = atoi(optarg);
			break;
		case 'u':
			chainParam.minchainscore = atoi(optarg);
			break;
		case 'e':
			chainParam.minchainlength = atoi(optarg);
			break;
		case 'v':
			chainParam.topkscoring = atoi(optarg);
			break;
		case 'o':
			chainParam.show = true;
			break;
		case ']':
			chainParam.show2 = true;
			break;
		case '?':
			exit(1);
		}
	}

	if(argc < optind || argc == 1) {
		fprintf(stderr,"Usage: %s [OPTION]\n", argv[0]);
		fprintf(stderr,"  <data>                    Index name or FASTA file\n");
		fprintf(stderr,"  -alph <file>              Use alphabet defined by file (option applies only to FASTA file)\n");
		fprintf(stderr,"  -dna                      Use DNA alphabet {A, C, G, T} and IUPAC wildcards (default)\n");
		fprintf(stderr,"  -rna                      Use RNA alphabet {A, C, G, U} and IUPAC wildcards\n");
		fprintf(stderr,"  -pat <file>               Structural pattern(s) to search for\n");
		fprintf(stderr,"  -for                      Search in the forward sequence (default)\n");
		fprintf(stderr,"  -rev                      Search in the reverse complement sequence. For searching in the forward sequence as well, combine it with -for\n");
		fprintf(stderr,"  -comp <file>              Load base pair complementarity rules from file\n");
		fprintf(stderr,"  -byseq                    Sort matches by sequence and matching position\n");
		fprintf(stderr,"  -byscore                  Sort matches of the same pattern by descending score\n");
		fprintf(stderr,"  -byscorea                 Sort matches of the same pattern by ascending score\n");
		fprintf(stderr,"  -table                    Print matches in table format\n");
		fprintf(stderr,"  -no-overlaps              Filter out low-scoring overlapping matches of the same pattern\n");
		fprintf(stderr,"  -silent                   Do not output matches\n");
		//fprintf(stderr,"  -progress                 Show progress message for each ~5%% processed data (DISABLED IN THIS VERSION)\n");
		/*Usage message for -pthread*/
		fprintf(stderr,"  -threads <x>              Use x CPU threads. Use x = 0 for automatic selection (Linux only) (default = 0)\n");
		fprintf(stderr,"  -parts <x>                Partition data in x parts for multithreaded search. Use x = 0 for automatic selection (Linux only) (default = 0)\n");

		fprintf(stderr,"\nOperation costs and thresholds. These override parameters set in the patterns file\n");
		fprintf(stderr,"  -replacement <cost>       Cost of a base mismatch (default = 1)\n");
		fprintf(stderr,"  -deletion <cost>          Cost of base deletion/insertion (default = 1)\n");
		fprintf(stderr,"  -arc-breaking <cost>      Cost of an arc-breaking (default = 1)\n");
		fprintf(stderr,"  -arc-altering <cost>      Cost of an arc-altering (default = 1)\n");
		fprintf(stderr,"  -arc-removing <cost>      Cost of an arc-removing (default = 2)\n");
		fprintf(stderr,"  -cost <x>                 Allow edit distance <= x (default = 0)\n");
		fprintf(stderr,"  -indels <x>               Allow number of indels <= x (default = allowed edit distance / cost of one indel)\n");

		fprintf(stderr,"\nIndex-based algorithmic variants*\n");
		fprintf(stderr,"  -lgslink                  Uses early-stop acceleration, enhanced suffix array, and generalized suffix links\n");
		fprintf(stderr,"  -lgslink_nof              Variant lgslink with disabled sequence-based filter\n");
		fprintf(stderr,"  -lesa                     Uses early-stop acceleration and enhanced suffix array\n");
		fprintf(stderr,"*lgslink requires tables suf, lcp, and sufinv. lesa requires only suf and lcp.\n");

		fprintf(stderr,"\nOnline algorithmic variants\n");
		fprintf(stderr,"  -scan                     Slides a window over the target sequence reusing matrix entries\n");
		fprintf(stderr,"  -lscan                    Scanning variant with early-stop acceleration\n");
		fprintf(stderr,"  -aligngl                  Aligns globally reporting the best alignment (no pattern matching) - SINGLE THREADED ONLY\n");

		fprintf(stderr,"                          \nChaining options\n");
		fprintf(stderr,"  -global                   Perform global chaining\n");
		fprintf(stderr,"  -local                    Perform local chaining\n");
		fprintf(stderr,"  -wf <wf>                  Apply weight factor > 0.0 to fragments\n");
		fprintf(stderr,"  -maxgap <width>           Allow chain gaps with up to the specified width\n");
		fprintf(stderr,"  -minscore <score>         Report only chains with at least the specified score\n");
		fprintf(stderr,"  -minlen <length>          Report only chains with number of fragments >= length\n");
		fprintf(stderr,"  -top <#>                  Report only top # scoring chains of each sequence\n");
		fprintf(stderr,"  -allglobal                Report for each sequence all global chains satisfying above criteria\n");
		fprintf(stderr,"  -show                     Show chains in the report\n");
		fprintf(stderr,"  -show2                    Print complete sequences and omit all other matching information\n");
		return EXIT_FAILURE;
	}
	if(optind < argc)
		argLoadFile = argv[optind];

	init(&affixArray, &searchParam,
			&argSuf, &argLcp, &argSufinv, &argLnk, &argAflk, &argSufr, &argLcpr, &argAflkr,
			&alphabet, argAlphabetFile, iPredefAlphabet);

	if(argLoadFile == NULL) {
		fprintf(stderr,"Please specify an index or FASTA file.\n");
		return EXIT_FAILURE;
	}

	if(bIndexBased) {
		if (loadAffFiles(&affixArray, argLoadFile,
				argSuf, argLcp, argSufinv, argLnk, argAflk,
				argSufr, argLcpr, argAflkr)) {
			fprintf(stderr,"Error loading index files.\n");
			return EXIT_FAILURE;
		}
	} else if(loadFastaFile(&affixArray, argLoadFile, false)) {
		fprintf(stderr,"Error loading FASTA file.\n");
		return EXIT_FAILURE;
	}

	setAlphabetConversionTable(&alphabet);
	setIupacTable(&alphabet);

	if(affixArray.multiSeq == NULL) {
		fprintf(stderr,"ERROR.\n");
		return EXIT_FAILURE;
	}

	if (affixArray.multiSeq->convSequences == NULL) {
		allocConvSequences(&affixArray);

		if (!chainParam.show2) fprintf(stderr, "\n%cPerforming alphabet conversion... ", LINESYMBOL);
		fflush(stderr);

		if (convertToAlphabet(affixArray.multiSeq->sequences,
							affixArray.multiSeq->convSequences,
							affixArray.length,
							false,
							affixArray.alphabet))
			return 1;

		if (!chainParam.show2) fprintf(stderr, "done\n");
	}

	if (argStruct != NULL) {
		multiPattern = loadPatternFile(argStruct, &searchParam, &affixArray);

		if(argCompFile != NULL)
			complementarityRules = loadComplementarityFile(argCompFile, &affixArray);
		else
			complementarityRules = loadDefaultComplementarityRules(&affixArray, iPredefAlphabet);

		setSearchPatterns(multiPattern, complementarityRules, &searchParam, &affixArray);
	}

	if (!chainParam.show2) fprintf(stderr, "%cNumber of sequences: %i\n", LINESYMBOL, affixArray.multiSeq->numSeqs);
	if (!chainParam.show2) fprintf(stderr, "%cTotal length:        %i\n", LINESYMBOL, affixArray.length - affixArray.multiSeq->numSeqs);

#if	__gnu_linux__

	if (iNumThreads == 0 || iNumThreads > sysconf(_SC_NPROCESSORS_ONLN))
		iNumThreads = sysconf(_SC_NPROCESSORS_ONLN);
	if (iNumPartitions == 0)
		iNumPartitions = sysconf(_SC_NPROCESSORS_ONLN);

#else
	if (iNumThreads == 0)
		iNumThreads = 1;
#endif

	/*inserted instead of aligne to call the threadMain with is our main for threading*/
	threadMain(multiPattern, &affixArray, complementarityRules, &searchParam, iNumThreads, iNumPartitions);
	/*commented out by Clemens and Ole*/
	//align(multiPattern, &affixArray, complementarityRules, &searchParam);

	freeAll(&affixArray, multiPattern);

	pthread_exit(NULL);
	//return EXIT_SUCCESS;
}


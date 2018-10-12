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
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h> //includes function getpagesize()
#include <fcntl.h> //includes function open()
#include <sys/mman.h>

#include "align.h"

int loadFastaFile(AffixArray *affixArray, char *fileName, bool bAllocConvSeq) {
	FILE *fp;
	char temp, *info;
	int i;

	if((fp=fopen(fileName,"r")) == NULL) {
		fprintf(stderr,"Error opening file \"%s\". If this is an index, please select the desired tables.\n", fileName);
		return 1;
	}

	if (fseek(fp, 0L, SEEK_END)) {
		fprintf(stderr,"Error reading file \"%s\".\n", fileName);
		return EXIT_FAILURE;
	}

	i=ftell(fp);
	if (i == 0) {
		fprintf(stderr, "File \"%s\" is empty.\n", fileName);
    	return EXIT_FAILURE;
	}

	rewind(fp);

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr, "Error. \"multiSeq\" was not initialized - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	affixArray->multiSeq->convSequences  = NULL;
	if ((affixArray->multiSeq->sequences = (char *) calloc(BUFFER1, sizeof(char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	if ((affixArray->multiSeq->seqEndPos      = (unsigned int *) calloc(BUFFER3, sizeof(unsigned int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	if ((affixArray->multiSeq->seqDescription = (char **) calloc(BUFFER3, sizeof(char*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqDescription\" - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}
	if ((affixArray->multiSeq->seqDescLength  = (int *) calloc(BUFFER3, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqDescLength\" - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	affixArray->multiSeq->sequencesMmapped = false;

	affixArray->length = affixArray->multiSeq->numSeqs = 0;
	temp = getc(fp);
	if(temp != '>') {
		fprintf(stderr,"The file does not have a valid fasta format.\n");
		return 1;
	}
	do {
		if(temp == '>') {
			if(affixArray->length > 0) {
				affixArray->multiSeq->sequences[affixArray->length++] = $;
				if(affixArray->length % BUFFER1 == 0) {
					if((affixArray->multiSeq->sequences = realloc(affixArray->multiSeq->sequences, (affixArray->length+BUFFER1)*sizeof(char))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
						return 1;
					}
				}
				affixArray->multiSeq->seqEndPos[affixArray->multiSeq->numSeqs-1] = affixArray->length - 1;
			}
			i = 0;
			if((info = (char *) calloc(BUFFER3, sizeof(char))) == NULL){
				fprintf(stderr,"Memory allocation failed for \"info\" - %s %d.\n", __FILE__, __LINE__);
				return 1;
			}
			while((temp = getc(fp)) != EOF && temp != '\r' && temp != '\n') {
				info[i] = temp;
				info[i+1] = 0;
				if(((++i)+1) % BUFFER3 == 0 && (info = realloc(info, (i+BUFFER3+1)*sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"info\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
			}
			if(affixArray->multiSeq->numSeqs % BUFFER3 == 0) {
				if((affixArray->multiSeq->seqEndPos = realloc(affixArray->multiSeq->seqEndPos, (affixArray->multiSeq->numSeqs+BUFFER3)*sizeof(int))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
				if((affixArray->multiSeq->seqDescription  = realloc(affixArray->multiSeq->seqDescription, (affixArray->multiSeq->numSeqs+BUFFER3)*sizeof(char*))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"seqInfo\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
				if((affixArray->multiSeq->seqDescLength = realloc(affixArray->multiSeq->seqDescLength, (affixArray->multiSeq->numSeqs+BUFFER3)*sizeof(int))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"seqDescLength\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
			}
			affixArray->multiSeq->seqDescription[affixArray->multiSeq->numSeqs] = info;
			affixArray->multiSeq->seqDescLength[affixArray->multiSeq->numSeqs] = i + 1; //1 extra space for 0
			affixArray->multiSeq->numSeqs++;
		}
		do {
			if(temp != '\r' && temp != '\n') {
				affixArray->multiSeq->sequences[affixArray->length++] = temp;
				if(affixArray->length % BUFFER1 == 0){
					if((affixArray->multiSeq->sequences = realloc(affixArray->multiSeq->sequences, (affixArray->length+BUFFER1)*sizeof(unsigned char))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
						return 1;
					}
				}
			}
		} while((temp = getc(fp)) != EOF && temp != '>');
	} while(temp != EOF);
	affixArray->multiSeq->sequences[affixArray->length++] = $;
	affixArray->multiSeq->seqEndPos[affixArray->multiSeq->numSeqs - 1] = affixArray->length - 1;

	fclose(fp);

	return EXIT_SUCCESS;
}

bool loadAffFiles(AffixArray *affixArray, char *fileName,
					bool argSuf, bool argLcp, bool argSufinv, bool argLnk, bool argAfsuf,
					bool argRpref, bool argRlcp, bool argAfrpref) {
	FILE *fp;
	int i, fd;

	int nameLength = strlen(fileName);
	char  fileNameCpy[1000];
	strcpy(fileNameCpy, fileName);

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr, "Error. \"multiSeq\" was not initialized - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	/********* base file *********/
	strcat(fileNameCpy, ".base");
	if((fp = fopen(fileNameCpy, "rb")) == NULL) {
	   fprintf(stderr,"Error opening file %s.\n", fileNameCpy);
	   return 1;
	}

	if(!fread(&affixArray->length, sizeof(int), 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s.\n", fileNameCpy);
		return 1;
	}
	if(!fread(&affixArray->multiSeq->numSeqs, sizeof(int), 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s.\n", fileNameCpy);
		return 1;
	}

	if((affixArray->multiSeq->seqEndPos = (unsigned int *) calloc(affixArray->multiSeq->numSeqs + 1, sizeof(unsigned int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
			return 1;
		}
	if(!fread(affixArray->multiSeq->seqEndPos, sizeof(int), affixArray->multiSeq->numSeqs + 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s.\n", fileNameCpy);
		return 1;
	}

	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* alphabet file *********/
	strcat(fileNameCpy, ".alph");
	if((fp = fopen(fileNameCpy, "rb")) == NULL) {
		fprintf(stderr,"Error opening file %s.\n", fileNameCpy);
		return 1;
	}

	if (affixArray->alphabet == NULL) {
		affixArray->alphabet = (Alphabet *) malloc(sizeof(Alphabet));
	}

	if(!fread(&affixArray->alphabet->numClasses, sizeof(int), 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
		return 1;
	}

	if((affixArray->alphabet->classSize = (int *) calloc(affixArray->alphabet->numClasses, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->classSize\".\n");
	if((affixArray->alphabet->classRepresentative = (unsigned char *) calloc(affixArray->alphabet->numClasses, sizeof(unsigned char))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->classRepresentative\".\n");
	if((affixArray->alphabet->eqClass = (char **) calloc(affixArray->alphabet->numClasses, sizeof(char*))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->eqClass\".\n");
	if((affixArray->alphabet->isWildCard = (bool *) calloc(affixArray->alphabet->numClasses, sizeof(bool))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->isWildCard\".\n");

	for(i = 0; i < affixArray->alphabet->numClasses; i++) {
		if(!fread(&affixArray->alphabet->classSize[i], sizeof(int), 1, fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
		if(!fread(&affixArray->alphabet->classRepresentative[i], sizeof(char), 1, fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}

		if (affixArray->alphabet->classSize[i] > 0) {
			if((affixArray->alphabet->eqClass[i] = (char *) calloc(affixArray->alphabet->classSize[i], sizeof(char*))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"alphabet->eqClass\".\n");
				return 1;
			}
			if(!fread(affixArray->alphabet->eqClass[i], sizeof(char), affixArray->alphabet->classSize[i], fp) != 0) {
				fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
				return 1;
			}
		} else {
			affixArray->alphabet->eqClass[i] = NULL;
		}
		if(!fread(&affixArray->alphabet->isWildCard[i], sizeof(bool), 1, fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
	}

	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* description file *********/
	strcat(fileNameCpy, ".des");
	if((fp = fopen(fileNameCpy, "rb")) == NULL) {
	   fprintf(stderr,"Error opening file %s.\n", fileNameCpy);
	   return 1;
	}

	if((affixArray->multiSeq->seqDescLength = (int *) calloc(affixArray->multiSeq->numSeqs, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"seqDescLength\".\n");
	if(!fread(affixArray->multiSeq->seqDescLength, sizeof(int), affixArray->multiSeq->numSeqs, fp) != 0) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
		return 1;
	}

	if((affixArray->multiSeq->seqDescription = (char **) calloc(affixArray->multiSeq->numSeqs, sizeof(char*))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"seqDescription\".\n");
	for(i=0; i < affixArray->multiSeq->numSeqs; i++){
		if((affixArray->multiSeq->seqDescription[i] = (char *) calloc(affixArray->multiSeq->seqDescLength[i], sizeof(char))) == NULL)
			fprintf(stderr,"Memory allocation failed for \"seqDescription\".\n");
		if(!fread(affixArray->multiSeq->seqDescription[i], sizeof(char), affixArray->multiSeq->seqDescLength[i], fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
	}

	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* sequences file *********/
	strcat(fileNameCpy, ".seq");
	if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
		return 1;
	}
	affixArray->multiSeq->sequences = mmap(NULL, affixArray->length * sizeof(unsigned char), PROT_READ, MAP_SHARED, fd, 0);
	affixArray->multiSeq->sequencesMmapped = true;

	fileNameCpy[nameLength] = 0;
	//fclose(fd);

	/********* alphabetically transformed sequences file *********/
	strcat(fileNameCpy, ".tseq");
	if((fd = open(fileNameCpy, O_RDONLY)) != -1) {
		affixArray->multiSeq->convSequences = mmap(NULL, affixArray->length * sizeof(unsigned char), PROT_READ, MAP_SHARED, fd, 0);
		affixArray->multiSeq->convSequencesMmapped = true;

		fileNameCpy[nameLength] = 0;
		//fclose(fd);
	} else {
		affixArray->multiSeq->convSequences = NULL;
		affixArray->multiSeq->convSequencesMmapped = false;
	}

	bool loadEArray(EArray *earray, int length, bool isEsa, char *fileName,
					bool argXarray, bool argXlcp, bool argXsufinv, bool argXlnk, bool argAfx);
	if(loadEArray(affixArray->esa, affixArray->length, 1, fileName,
					argSuf, argLcp, argSufinv, argLnk, argAfsuf)) return 1;
	if(loadEArray(affixArray->erpa, affixArray->length, 0, fileName,
    				argRpref, argRlcp, false, false, argAfrpref)) return 1;

	return 0;
}

bool loadEArray(EArray *earray, int length, bool isEsa, char *fileName,
				bool argXarray, bool argXlcp, bool argXsufinv, bool argXlnk, bool argAfx) {
	int i, j, fd, *dataInt32;

	int nameLength = strlen(fileName);
	char  fileNameCpy[100];
	strcpy(fileNameCpy, fileName);

	/********* suf | rpref file *********/
	if(argXarray) {
		if(isEsa)
			strcat(fileNameCpy, ".suf");
		else
			strcat(fileNameCpy, ".sufr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->xarray = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	/********* lcp | rlcp file *********/
	if(argXlcp) {
		if(isEsa)
			strcat(fileNameCpy, ".lcp");
		else
			strcat(fileNameCpy, ".lcpr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->xlcp = mmap(NULL, length * sizeof(unsigned char), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);

		/********* lcp | rlcpe exception file *********/

		if(isEsa)
			strcat(fileNameCpy, ".lcpe");
		else
			strcat(fileNameCpy, ".lcper");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
			earray->xlcpException = (LcpException *) malloc(sizeof(LcpException));
			earray->xlcpException->numExceptions = 0;
		} else {
			dataInt32 = mmap(NULL, sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
			if((earray->xlcpException = (LcpException *) malloc(sizeof(LcpException))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"lcpException\".\n");
			earray->xlcpException->numExceptions = dataInt32[0];

			dataInt32 = mmap(NULL, earray->xlcpException->numExceptions * sizeof(int) * 2 + sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
			if((earray->xlcpException->index = (int *) calloc(earray->xlcpException->numExceptions, sizeof(int))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"lcpException->index\".\n");
			if((earray->xlcpException->value = (int *) calloc(earray->xlcpException->numExceptions, sizeof(int))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"lcpException->value\".\n");

			j = 1;
			for(i=0; i < earray->xlcpException->numExceptions; i++) {
				earray->xlcpException->index[i] = dataInt32[j++];
				earray->xlcpException->value[i] = dataInt32[j++];
			}
			//close(fd);
		}
		fileNameCpy[nameLength] = 0;

	}

	/********* afsuf | afsufr file *********/
	if(argAfx) {
		if(isEsa)
			strcat(fileNameCpy, ".aflk");
		else
			strcat(fileNameCpy, ".aflkr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->affixLink = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	/********* sufinv | sufinvr file *********/
	if(argXsufinv) {
		if(isEsa)
			strcat(fileNameCpy, ".sufinv");
		else
			strcat(fileNameCpy, ".sufinvr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
			fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
			return 1;
		}

		earray->xsufinv = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	/********* lnk | lnkr file *********/
	if(argXlnk) {
		if(isEsa)
			strcat(fileNameCpy, ".lnk");
		else
			strcat(fileNameCpy, ".lnkr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
			fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
			return 1;
		}

		earray->xlnk = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		//fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	return 0;
}

MultiPattern* loadPatternFile(char *fileName, SearchParam *searchParam, AffixArray *affixArray) {
	MultiPattern *multiPattern;
	char *p, cNumber[12], maxLength=12, parameters[100], *pParameters;
	int iWhichPipeSymbol;

	int i, j, k, x, pos, fd;
	int iBrackets = 0;
	Pattern *pattern;
	bool bParams;
	struct stat sb;

	Cost *cost = getUnitCosts();

	if((fd = open(fileName, O_RDONLY)) == -1) {
	   fprintf(stderr,"Error opening file \"%s\".\n", fileName);
	   exit(1);
	}
	if (fstat(fd, &sb) == -1) {
		fprintf(stderr,"fstat\n");
		exit(1);
	}
	p = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	if((multiPattern = (MultiPattern *) malloc(sizeof(MultiPattern))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"multiPattern\".\n");
		exit(1);
	}

	if((pattern = (Pattern *) calloc(BUFFER2, sizeof(Pattern))) == NULL){
		fprintf(stderr,"Memory allocation failed for \"patternP\".\n");
		exit(1);
	}

	multiPattern->bUsesStartPosParam = false;

	pos = i = 0;
	do {
		if(p[pos] != '>') {
			fprintf(stderr, "Format of \"%s\" is invalid. Delete blank lines and make sure that the first line of a pattern begins with \">\".\n", fileName);
			exit(1);
		}

		pattern[i].iId = i;

		/******************* Read description *******************/
		if((pattern[i].desc = (char *) calloc(BUFFER2, sizeof(char))) == NULL){
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].desc\".\n",i);
			exit(1);
		}

		for (j = 0; j < 100; ++j) {
			parameters[j] = 0;
		}
		for (j = 0; j < maxLength; ++j) {
			cNumber[j] = 0;
		}

		k = j = iWhichPipeSymbol = 0;
		bParams = 0;
		while(++pos < sb.st_size && p[pos] != '\r' && p[pos] != '\n') {
			if(p[pos] == '|') {
				bParams = 1;
			}

			if (!bParams) {
				pattern[i].desc[j++] = p[pos];

				if((j + 1) % BUFFER2 == 0 && (pattern[i].desc = realloc(pattern[i].desc, (j + BUFFER2 + 1) * sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"patternP[%d].desc\".\n", i);
					exit(1);
				}
			} else {
				parameters[k++] = p[pos];
			}
		}

		if((pattern[i].desc = realloc(pattern[i].desc, (j + 1) * sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].desc\".\n", i);
			exit(1);
		}
		pattern[i].desc[j] = 0;

		if((pattern[i].cost = (Cost *) malloc(sizeof(Cost))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"pattern[i].cost\".\n");
			exit(1);
		}

		if (searchParam->cost.iCostReplacement != -1) { //overrides value in file
			pattern[i].cost->iCostReplacement = searchParam->cost.iCostReplacement;
		} else if ((pParameters = strstr(parameters, "replacement")) != NULL || (pParameters = strstr(parameters, "rep")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].cost->iCostReplacement = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].cost->iCostReplacement = cost->iCostReplacement;
		}

		if (searchParam->cost.iCostDeletion != -1) { //overrides value in file
			pattern[i].cost->iCostDeletion = searchParam->cost.iCostDeletion;
		} else if ((pParameters = strstr(parameters, "deletion")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].cost->iCostDeletion = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].cost->iCostDeletion = cost->iCostDeletion;
		}


		if (searchParam->cost.iFCostDeletion == -1)
			pattern[i].cost->iFCostDeletion = pattern[i].cost->iCostDeletion;
		else
			pattern[i].cost->iFCostDeletion = searchParam->cost.iFCostDeletion;

		if (searchParam->cost.iFCostReplacement == -1)
			pattern[i].cost->iFCostReplacement = pattern[i].cost->iCostReplacement;
		else
			pattern[i].cost->iFCostReplacement = searchParam->cost.iFCostReplacement;


		if (searchParam->cost.iCostArcBreak != -1) { //overrides value in file
			pattern[i].cost->iCostArcBreak = searchParam->cost.iCostArcBreak;
		} else if ((pParameters = strstr(parameters, "arc-breaking")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].cost->iCostArcBreak = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].cost->iCostArcBreak = cost->iCostArcBreak;
		}

		if (searchParam->cost.iCostArcAltering != -1) { //overrides value in file
			pattern[i].cost->iCostArcAltering = searchParam->cost.iCostArcAltering;
		} else if ((pParameters = strstr(parameters, "arc-altering")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].cost->iCostArcAltering = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].cost->iCostArcAltering = cost->iCostArcAltering;
		}

		if (searchParam->cost.iCostArcRemoving != -1) { //overrides value in file
			pattern[i].cost->iCostArcRemoving = searchParam->cost.iCostArcRemoving;
		} else if ((pParameters = strstr(parameters, "arc-removing")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].cost->iCostArcRemoving = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].cost->iCostArcRemoving = cost->iCostArcRemoving;
		}

		if (searchParam->uiThreshold != -1) { //overrides value in file
			pattern[i].uiThreshold = searchParam->uiThreshold;
		} else if ((pParameters = strstr(parameters, "cost")) != NULL || (pParameters = strstr(parameters, "edist")) != NULL || (pParameters = strstr(parameters, "threshold")) != NULL || (pParameters = strstr(parameters, "thr")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].uiThreshold = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].uiThreshold = 0;
		}
		if (searchParam->uiIndels != -1) { //overrides value in file
			pattern[i].uiIndels = searchParam->uiIndels;
		} else if ((pParameters = strstr(parameters, "indel")) != NULL) { //also works with parameter "indels"
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].uiIndels = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].uiIndels = pattern[i].uiThreshold / cost->iCostDeletion;
		}
		if ((pParameters = strstr(parameters, "startpos")) != NULL) {
			multiPattern->bUsesStartPosParam = true;
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].startpos = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		}
		if ((pParameters = strstr(parameters, "weight")) != NULL) {
			k = 0;
			while (pParameters[0] != '|' && pParameters[0] != '\0') {
				if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
					cNumber[k++] = pParameters[0];
				}
				pParameters++;
			}
			if (k > 0) {
				pattern[i].weight = atoi(cNumber);
				for (x = 0; x <= k; x++) {
					cNumber[x] = 0;
				}
			}
		} else {
			pattern[i].weight = 0;
		}

		if (adjustNumIndelsAccordingToEdist (&pattern[i].uiIndels, &pattern[i].uiThreshold, pattern[i].cost))
			fprintf(stderr, "Warning. Adjusting the number of allowed indels of pattern %s to %d - the cost of the specified number of indels exceeded the cost threshold enforced by the allowed edit distance.\n", pattern[i].desc, pattern[i].uiIndels);

		/******************* Read sequence *******************/
		if((pattern[i].forwardStrand = (PatternStrandDirection *) malloc(sizeof(PatternStrandDirection))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"pattern[%d].forwardStrand\".\n", i);
			exit(1);
		}
		pattern[i].reverseStrand = NULL;

		if((pattern[i].forwardStrand->seq = (char *) calloc(BUFFER2, sizeof(char))) == NULL){
			fprintf(stderr,"Memory allocation failed for \"pattern[%d].forwardStrand->seq\".\n", i);
			exit(1);
		}
		j = 0;
		while(++pos < sb.st_size && p[pos] != '\r' && p[pos] != '\n') {
			if(p[pos] != ' ') {
				pattern[i].forwardStrand->seq[j++] = p[pos];
				if((j + 1) % BUFFER2 == 0 && (pattern[i].forwardStrand->seq = realloc(pattern[i].forwardStrand->seq, (j + BUFFER2 + 1) * sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"patternP[%d].forwardStrand->seq\".\n", i);
					exit(1);
				}
			}
		}
		if((pattern[i].forwardStrand->seq = realloc(pattern[i].forwardStrand->seq, (j + 1) * sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"pattern[%d].forwardStrand->seq\".\n", i);
			exit(1);
		}
		pattern[i].forwardStrand->seq[j] = 0;
		pattern[i].iLength = j;

		/******************* Read structure *******************/
		if((pattern[i].forwardStrand->structure = (char *) calloc(BUFFER2, sizeof(char))) == NULL){
			fprintf(stderr,"Memory allocation failed for \"pattern[%d].forwardStrand->structure\".\n", i);
			exit(1);
		}
		j = 0;
		while(++pos < sb.st_size && p[pos] != '\r' && p[pos] != '\n') {
			if(p[pos] != ' ') {
				if(p[pos] == '(') {
					iBrackets++;
				} else if(p[pos] == ')')
					iBrackets--;
				if(iBrackets < 0) {
					fprintf(stderr,"Structure of pattern \"%s\" is invalid.\n", pattern[i].desc);
					exit(1);
				}
				pattern[i].forwardStrand->structure[j++] = p[pos];
				if((j + 1) % BUFFER2 == 0 && (pattern[i].forwardStrand->structure = realloc(pattern[i].forwardStrand->structure, (j + BUFFER2 + 1) * sizeof(unsigned char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"pattern[%d].forwardStrand->structure\".\n", i);
					exit(1);
				}
			}
		}
		pos++;
		pattern[i].forwardStrand->structure[j] = 0;

		if(iBrackets != 0) {
			fprintf(stderr,"Structure of pattern \"%s\" is invalid.\n", pattern[i].desc);
			exit(1);
		}
		if(strlen((char*) pattern[i].forwardStrand->structure) != pattern[i].iLength) {
			fprintf(stderr,"Sequence and structure string lengths of pattern \"%s\" do not match.\n", pattern[i].desc);
			exit(1);
		}
		if((pattern[i].forwardStrand->structure = realloc(pattern[i].forwardStrand->structure, (j + 1) * sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"pattern[%d].forwardStrand->structure\".\n", i);
			exit(1);
		}

		if (pattern[i].uiIndels > (pattern[i].iLength / 2)) {
			fprintf(stderr, "Warning. The number of allowed indels %d for pattern %s corresponds to more than 50%% of the pattern length. This can lead to a large number of matches.\n", pattern[i].uiIndels, pattern[i].desc);
		}

		if(((++i)+1) % BUFFER2 == 0 && (pattern = realloc(pattern, (i + BUFFER2 + 1) * sizeof(Pattern))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"pattern\".\n");
			exit(1);
		}
	} while(pos < sb.st_size);

	free(cost);

	multiPattern->iNumPatterns = i;
	multiPattern->pattern      = pattern;

	return multiPattern;
}

Alphabet* loadAlphabetFile(char *fileName, Alphabet *alphabet) {
	FILE *fp;
	int i, j, l;
	int iNonMatchingClassSize = 0;
	bool bFoundRep, bFoundNonMatching;
	char temp;
	char *eqClassNonMatching = NULL;
	unsigned char ucClassRepresentativeNonMatching = 0;

	if((fp=fopen(fileName,"r")) == NULL) {
		fprintf(stderr,"Error opening alphabet file %s.\n", fileName);
		exit(1);
	}

	if (fseek(fp, 0L, SEEK_END)) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileName, __FILE__, __LINE__);
		exit(1);
	}

	i=ftell(fp);
	if (i == 0) {
		fprintf(stderr, "File %s is empty - %s %d.\n", fileName, __FILE__, __LINE__);
    	exit(1);
	}

	rewind(fp);

	if ((alphabet->eqClass    = (char **) calloc(127, sizeof(char*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"eqClass\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
	if ((alphabet->classSize  = (int *) calloc(127, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"classSize\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
	if ((alphabet->isWildCard = (bool *) calloc(127, sizeof(bool))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"isWildCard\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
	if ((alphabet->classRepresentative = (unsigned char *) calloc(127, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"classRepresentative\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	i = 0;
	temp = getc(fp);
	do {
		if(temp == '\r' || temp == '\n') {
			temp = getc(fp);
		}
		else if(temp != EOF) {
			j = bFoundRep = 0;
			bFoundNonMatching = false;

			if(temp == '*') {
				alphabet->isWildCard[i] = 1;
				temp = getc(fp);
			} else if(temp == '!') {
				if (iNonMatchingClassSize > 0) {
					fprintf(stderr,"Class of non-matching characters in the target sequence can only be defined once.\n");
					exit(1);
				}

				bFoundNonMatching = true;
				temp = getc(fp);
			} else {
				alphabet->isWildCard[i] = 0;
			}

			if (!bFoundNonMatching) {
				if ((alphabet->eqClass[i] = (char *) calloc(UCHAR_MAX, sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
					exit(1);
				}
			} else {
				if ((eqClassNonMatching = (char *) calloc(UCHAR_MAX, sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
					exit(1);
				}
			}

			if (temp != EOF && temp != '\r' && temp != '\n') {
				do {
					if(temp != ' ') {
						if(bFoundRep && j > 0) {
							if (bFoundNonMatching)
								ucClassRepresentativeNonMatching = temp;
							else
								alphabet->classRepresentative[i] = temp;
						} else {
							if (bFoundNonMatching)
								eqClassNonMatching[j++] = temp;
							else
								alphabet->eqClass[i][j++] = temp;
						}
					} else {
						bFoundRep = 1;
					}
				} while((temp = getc(fp)) != EOF && temp != '\r' && temp != '\n');
			}
			if(j > 0) {
				if (!bFoundRep) {
					if (bFoundNonMatching)
						ucClassRepresentativeNonMatching = eqClassNonMatching[0];
					else
						alphabet->classRepresentative[i] = alphabet->eqClass[i][0];
				}

				// Look for representative as an element of eqClass
				bFoundRep = false;
				if (bFoundNonMatching) {
					for (l = 0; l < j; l++) {
						if (eqClassNonMatching[l] == ucClassRepresentativeNonMatching) {
							bFoundRep = true;
							break;
						}
					}
					if (!bFoundRep) // If representative is not found, add it
						eqClassNonMatching[j++] = ucClassRepresentativeNonMatching;

					iNonMatchingClassSize = j;
				}
				else  {
					for (l = 0; l < j; l++) {
						if (alphabet->eqClass[i][l] == alphabet->classRepresentative[i]) {
							bFoundRep = true;
							break;
						}
					}
					if (!bFoundRep) // If representative is not found, add it
						alphabet->eqClass[i][j++] = alphabet->classRepresentative[i];

					alphabet->classSize[i++] = j;
				}

			}
		}
	} while(temp != EOF);

	if (iNonMatchingClassSize == 0) {
		alphabet->eqClass[i]             = NULL;
		alphabet->classRepresentative[i] = ' ';
		alphabet->classSize[i]           = 0;
		alphabet->isWildCard[i]          = 0;
	} else {
		alphabet->eqClass[i]             = eqClassNonMatching;
		alphabet->classRepresentative[i] = ucClassRepresentativeNonMatching;
		alphabet->classSize[i]           = iNonMatchingClassSize;
		alphabet->isWildCard[i]          = 0;
	}

	alphabet->numClasses = i + 1;

	fclose(fp);

	return alphabet;
}

void allocConvSequences (AffixArray *affixArray) {
	if ((affixArray->multiSeq->convSequences = (char *) calloc(affixArray->length, sizeof(char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"convSequences\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
}

void printfMatrices_S_up_toLastRegion(unsigned int ***edist, char ***operation, char *cu, char *cv, AffixArray *affixArray, int iLastRegion, PatternStrandDirection *patternStrandDirection) {
	int ii, ij, ik;
	for (ik = patternStrandDirection->patternStructures->patternRegion[iLastRegion].ikBigStart; ik <= patternStrandDirection->patternStructures->patternRegion[iLastRegion].ikBigEnd; ik++) {
		fprintf(stderr, "k=%d\n", ik);

		if (cv != NULL) {
			fprintf(stderr, "			");
			for (ij = 0; ij < patternStrandDirection->patternStructures->patternRegion[iLastRegion].ijBigGoal[ik]; ij++) {
				fprintf(stderr, "%c	", c(cv[ij + ik]));
			}
			fprintf(stderr, "\n");
		}
		ii = 0;
		for (ij = 0; ij <= patternStrandDirection->patternStructures->patternRegion[iLastRegion].ijBigGoal[ik]; ij++) {
			if (ij == 0) {
				fprintf(stderr, "	");
				if (ii == 0)
					fprintf(stderr, "	");
			}
			fprintf(stderr, "%d	", S(edist, ik, ii, ij) == INT_MAX ? 999 : S(edist, ik, ii, ij));
		}
		fprintf(stderr, "\n");

		for (ii = 0; ii <= patternStrandDirection->patternStructures->patternRegion[iLastRegion].iR; ii++) {

			if (ii > 0) fprintf(stderr, "%d	%c", ii, c(cu[ii - 1]));
			for (ij = 0; ij <= patternStrandDirection->patternStructures->patternRegion[iLastRegion].ijBigGoal[ik]; ij++) {
				if (ij == 0) {
					fprintf(stderr, "	");
					if (ii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", S(edist, ik, ii, ij) == INT_MAX ? 999 : S(edist, ik, ii, ij));
			}
			fprintf(stderr, "\n");
		} //return;
		fprintf(stderr, "\n");
	}
}

void printMatrices_A(unsigned int *edist, char *operation, unsigned int ***position, unsigned int *Kfactor, unsigned int Jfactor, char *cu, char *cv, unsigned int uim, unsigned int uin, AffixArray *affixArray) {
	unsigned int uii, uij, uik, uitemp;

	for (uik = 0; uik < uin; uik++) {
		fprintf(stderr, "k=%d\n", uik);

		for (uii = 0; uii <= uim; uii++) {
			uitemp = uin - uik;

			if (cv != NULL && uii == 0) {
				fprintf(stderr, "			");
				for (uij = 0; uij < uitemp; uij++) {
					fprintf(stderr, "%c	", c(cv[uij + uik]));
				}
				fprintf(stderr, "\n");
			}

			if (uii > 0) fprintf(stderr, "%d	%c", uii, c(cu[uii - 1]));
			for (uij = 0; uij <= uitemp; uij++) {
				if (uij == 0) {
					fprintf(stderr, "	");
					if (uii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", A(edist, Kfactor, Jfactor, uik, uii, uij) == INT_MAX ? 999 : A(edist, Kfactor, Jfactor, uik, uii, uij));
			}
			fprintf(stderr, "\n");
		} //return;
		fprintf(stderr, "\n");

		for (uii = 0; uii <= uim; uii++) {
			uitemp = uin - uik;

			if (cv != NULL && uii == 0) {
				fprintf(stderr, "			");
				for (uij = 0; uij < uitemp; uij++) {
					fprintf(stderr, "%c	", c(cv[uij + uik]));
				}
				fprintf(stderr, "\n");
			}

			if (uii > 0) fprintf(stderr, "%d	%c", uii, c(cu[uii - 1]));
			for (uij = 0; uij <= uitemp; uij++) {
				if (uij == 0) {
					fprintf(stderr, "	");
					if (uii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", A(operation, Kfactor, Jfactor, uik, uii, uij));
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
}

void printMatrices_S(unsigned int ***edist, char ***operation, char *cu, char *cv, unsigned int uim, unsigned int uin, AffixArray *affixArray) {
	unsigned int uii, uij, uik, uitemp;

	for (uik = 0; uik < uin; uik++) {
		fprintf(stderr, "k=%d\n", uik);

		for (uii = 0; uii <= uim; uii++) {
			uitemp = uin - uik;

			if (cv != NULL && uii == 0) {
				fprintf(stderr, "			");
				for (uij = 0; uij < uitemp; uij++) {
					fprintf(stderr, "%c	", c(cv[uij + uik]));
				}
				fprintf(stderr, "\n");
			}

			if (uii > 0) fprintf(stderr, "%d	%c", uii, c(cu[uii - 1]));
			for (uij = 0; uij <= uitemp; uij++) {
				if (uij == 0) {
					fprintf(stderr, "	");
					if (uii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", S(edist, uik, uii, uij) == INT_MAX ? 999 : S(edist, uik, uii, uij));
			}
			fprintf(stderr, "\n");
		} //return;
		fprintf(stderr, "\n");

		/*for (uii = 0; uii <= uim; uii++) {
			uitemp = uin - uik;

			if (cv != NULL && uii == 0) {
				fprintf(stderr, "			");
				for (uij = 0; uij < uitemp; uij++) {
					fprintf(stderr, "%c	", c(cv[uij + uik]));
				}
				fprintf(stderr, "\n");
			}

			if (uii > 0) fprintf(stderr, "%d	%c", uii, c(cu[uii - 1]));
			for (uij = 0; uij <= uitemp; uij++) {
				if (uij == 0) {
					fprintf(stderr, "	");
					if (uii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", S(operation, uik, uii, uij));
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");*/
		//if (uik > 4) return;
	}
}

void printMatrices(unsigned int ***edist, char ***operation, char *cu, char *cv, unsigned int uim, unsigned int uin, AffixArray *affixArray) {
	unsigned int uii, uij, uik, uitemp;

	for (uik = 0; uik < uin; uik++) {
		fprintf(stderr, "k=%d\n", uik);

		for (uii = 0; uii <= uim; uii++) {
			uitemp = uin - uik;

			if (cv != NULL && uii == 0) {
				fprintf(stderr, "			");
				for (uij = 0; uij < uitemp; uij++) {
					fprintf(stderr, "%c	", c(cv[uij + uik]));
				}
				fprintf(stderr, "\n");
			}

			if (uii > 0) fprintf(stderr, "%d	%c", uii, c(cu[uii - 1]));
			for (uij = 0; uij <= uitemp; uij++) {
				if (uij == 0) {
					fprintf(stderr, "	");
					if (uii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", P(edist, uik, uii, uij) == INT_MAX ? 999 : P(edist, uik, uii, uij));
			}
			fprintf(stderr, "\n");
		} //return;
		fprintf(stderr, "\n");

		for (uii = 0; uii <= uim; uii++) {
			uitemp = uin - uik;

			if (cv != NULL && uii == 0) {
				fprintf(stderr, "			");
				for (uij = 0; uij < uitemp; uij++) {
					fprintf(stderr, "%c	", c(cv[uij + uik]));
				}
				fprintf(stderr, "\n");
			}

			if (uii > 0) fprintf(stderr, "%d	%c", uii, c(cu[uii - 1]));
			for (uij = 0; uij <= uitemp; uij++) {
				if (uij == 0) {
					fprintf(stderr, "	");
					if (uii == 0)
						fprintf(stderr, "	");
				}
				fprintf(stderr, "%d	", P(operation, uik, uii, uij));
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
}

bool** loadComplementarityFile(char *fileName, AffixArray *affixArray) {
	char *p, *c;
	int i, j, k, pos, fd;
	struct stat sb;
	bool **bCompCheck;

	if((fd = open(fileName, O_RDONLY)) == -1) {
	   fprintf(stderr,"Error opening file \"%s\".\n", fileName);
	   exit(1);
	}
	if (fstat(fd, &sb) == -1) {
		fprintf(stderr,"fstat\n");
		exit(1);
	}
	p = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	if ((bCompCheck = (bool **) calloc(affixArray->alphabet->numClasses, sizeof(bool*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"bCompCheck\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		if ((bCompCheck[i] = (bool *) calloc(affixArray->alphabet->numClasses, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"bCompCheck_i\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			bCompCheck[i][j] = 0;
		}
	}

	pos = 0;
	do {
		if(p[pos] != '\r' && p[pos] != '\n') {
			c = (char *) calloc(BUFFER2, sizeof(char));

			j = 0;
			do {
				c[j++] = p[pos++];
			} while (pos < sb.st_size && p[pos] != '\r' && p[pos] != '\n');

			if(convertToAlphabet(c, c, j, true, affixArray->alphabet))
				exit(1);

			for (k = 1; k < j; k++) {
				bCompCheck[c[0] - 1][c[k] - 1] = 1;
			}

			free(c);
		}

		pos++;
	} while(pos < sb.st_size);

	/* If base X is complement of base Y, make vice-versa true */
	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			if(bCompCheck[i][j]) {
				bCompCheck[j][i] = 1;
			}
		}
	}

	if (munmap(p, sb.st_size) == -1)
		fprintf(stderr,"Error unmapping file \"%s\".\n", fileName);

	close(fd);

	return bCompCheck;
}

bool** loadDefaultComplementarityRules(AffixArray *affixArray, const int predefAlphabet) {
	int i, j, k, numComps=/*3*/2, length=2;
	//char c[3][2] = {{'A', 'U'}, {'C', 'G'}, {'G', 'U'}}, cc[length];
	char c[2][2] = {{'A', 'U'}, {'C', 'G'}}, cc[length];

	bool **bCompCheck = (bool **) calloc(affixArray->alphabet->numClasses, sizeof(bool*));

	if (predefAlphabet == 0) {
		c[0][1] = 'T';
		fprintf(stderr, "\n%cUsing Watson-Crick complementarity rules A-T, T-A, C-G, G-C\n", LINESYMBOL);
	} else {
		fprintf(stderr, "\n%cUsing Watson-Crick complementarity rules A-U, U-A, C-G, G-C\n", LINESYMBOL);
	}

	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		bCompCheck[i] = (bool *) calloc(affixArray->alphabet->numClasses, sizeof(bool));
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			bCompCheck[i][j] = 0;
		}
	}

	for (i = 0; i < numComps; i++) {
		if(convertToAlphabet(c[i], cc, 2, true, affixArray->alphabet))
			exit(1);

		for (k = 1; k < length; k++) {
			bCompCheck[cc[0] - 1][cc[k] - 1] = 1;
		}
	}

	/* If base X is complement of base Y, make vice-versa true */
	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			if(bCompCheck[i][j]) {
				bCompCheck[j][i] = 1;
			}
		}
	}

	return bCompCheck;
}

//Reverse arbitrary string and return a new pointer
char* reverseStringNewPointer(char *seq, int length) {
	int i, iLengthBy2 = length / 2;
	char *newSeq = (char *) malloc((length + 1) * sizeof(char));

	for(i = 0; i < iLengthBy2; i++) {
		newSeq[i] = seq[length - 1 - i];
		newSeq[length - 1 - i] = seq[i];
	}
	if (length % 2 != 0)
		newSeq[iLengthBy2] = seq[iLengthBy2];
	newSeq[length] = 0;

	return newSeq;
}

char getComplementaryIupacSymbol (char cSymbol, bool bIsDNA) {
	switch (toupper(cSymbol)) {
		case 'R':
			return 'Y';
		case 'Y':
			return 'R';
		case 'S':
			return 'S';
		case 'W':
			return 'W';
		case 'K':
			return 'M';
		case 'M':
			return 'K';
		case 'B':
			return 'V';
		case 'V':
			return 'B';
		case 'D':
			return 'H';
		case 'H':
			return 'D';
		case 'N':
			return 'N';
		case 'A':
			return bIsDNA ? 'T' : 'U';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'U':
			return 'A';
		default:
			fprintf(stderr, "Warning. \"%c\" is not an IUPAC character and will not be mapped to a complement symbol for searching the reverse complement strand.\n", cSymbol);
			return cSymbol;
	}
}

void setSearchPatterns(MultiPattern *multiPattern, bool **compRules, SearchParam *searchParam, AffixArray *affixArray) {
	unsigned int i, j, iStrLength, iCountArcs;

	bool argSearchForward = searchParam->bSearchForwardString;
	bool argSearchReverse = searchParam->bSearchReverseString;
	bool **revComplementarityRules = NULL;
	Pattern *pattern = multiPattern->pattern;
	Alphabet *alphabet = affixArray->alphabet;

	char cComplementarySymbol;
	char cCodeComplementarySymbol;
	bool bIsDNA = false;
	char *cArrayComplementaryCode;

	if ((!argSearchForward && !argSearchReverse) || argSearchForward) {
		argSearchForward = searchParam->bSearchForwardString = true; // search in the forward seq (default)
		for (i = 0; i < multiPattern->iNumPatterns; i++) {
			pattern[i].forwardStrand->patternStructures = processPattern(pattern[i].forwardStrand->structure, pattern[i].iLength, pattern[i].uiIndels);
		}
	}

	// If no pattern weight (score) was set...
	for (i = 0; i < multiPattern->iNumPatterns; i++) {
		iStrLength = pattern[i].iLength;
		iCountArcs = 0;
		for (j = 0; j < iStrLength; j++) {
			if (pattern[i].forwardStrand->structure[j] == '(')
				iCountArcs++;
		}
		if (pattern[i].weight == 0) {
			pattern[i].weight = iStrLength * pattern->cost->iCostReplacement + iCountArcs * pattern->cost->iCostArcRemoving;
		} else if (pattern[i].weight < (iStrLength * pattern->cost->iCostReplacement + iCountArcs * pattern->cost->iCostArcRemoving)) {
			fprintf(stderr, "Warning. The low weight of pattern %s can lead to matches with negative score and unexpected chains.\n", pattern[i].desc);
		}
	}

	if (argSearchReverse) {
        // Begin computation of cArrayComplementaryCode
		if ((cArrayComplementaryCode = (char *) calloc(alphabet->numClasses + 1, sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		for (i = 1; i <= alphabet->numClasses; i++) {
			if (toupper(alphabet->classRepresentative[i - 1]) == 'T') {
				bIsDNA = true;
				break;
			}
		}

		for (i = 1; i < alphabet->numClasses; i++) {
			cComplementarySymbol = getComplementaryIupacSymbol(alphabet->classRepresentative[i - 1], bIsDNA);

			cCodeComplementarySymbol = 0;
			for (j = 1; j < alphabet->numClasses; j++) {
				if (toupper(alphabet->classRepresentative[j - 1]) == toupper(cComplementarySymbol)) {
					cCodeComplementarySymbol = j;
					break;
				}
			}
			if (cCodeComplementarySymbol == 0) {
				fprintf(stderr, "Warning. Complement IUPAC character class \"%c\" of character \"%c\" cannot be found. \"%c\" will not be mapped to a complement symbol for searching the reverse complement strand.\n", cComplementarySymbol, alphabet->classRepresentative[i], alphabet->classRepresentative[i]);
				cArrayComplementaryCode[i] = i;
			} else {
				cArrayComplementaryCode[i] = cCodeComplementarySymbol;
			}
		}
		// End computation of cArrayComplementaryCode

		revComplementarityRules = loadReverseComplementarityRules(compRules, cArrayComplementaryCode, affixArray);

		for (i = 0; i < multiPattern->iNumPatterns; i++) {
			iStrLength = pattern[i].iLength;

			if((pattern[i].reverseStrand = (PatternStrandDirection *) malloc(sizeof(PatternStrandDirection))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"pattern[%d].reverseStrand\".\n", i);
				exit(1);
			}

			pattern[i].reverseStrand->seq       = reverseStringNewPointer(pattern[i].forwardStrand->seq, iStrLength);
			pattern[i].reverseStrand->structure = reverseStringNewPointer(pattern[i].forwardStrand->structure, iStrLength);

			if(convertToAlphabet(pattern[i].reverseStrand->seq,
					pattern[i].reverseStrand->seq,
					iStrLength,
					true,
					affixArray->alphabet))
				exit(1);

			for (j = 0; j < iStrLength; j++) {
				pattern[i].reverseStrand->seq[j] = cArrayComplementaryCode[(int) pattern[i].reverseStrand->seq[j]];

				if (pattern[i].reverseStrand->structure[j] == ')') {
					pattern[i].reverseStrand->structure[j] = '(';
				} else if (pattern[i].reverseStrand->structure[j] == '(') {
					pattern[i].reverseStrand->structure[j] = ')';
				}
			}

			pattern[i].reverseStrand->patternStructures = processPattern(pattern[i].reverseStrand->structure, pattern[i].iLength, pattern[i].uiIndels);
		}
	}

	if (argSearchForward) {
		for (i = 0; i < multiPattern->iNumPatterns; i++) {
			if(convertToAlphabet(pattern[i].forwardStrand->seq,
					pattern[i].forwardStrand->seq,
					pattern[i].iLength,
					true,
					affixArray->alphabet))
				exit(1);
		}
	} else {
		for (i = 0; i < multiPattern->iNumPatterns; i++) {
			free(pattern[i].forwardStrand);
			pattern[i].forwardStrand = NULL;
		}
	}

	for (i = 0; i < multiPattern->iNumPatterns; i++) {
		if (pattern[i].forwardStrand != NULL) {
			pattern[i].forwardStrand->pattern = &pattern[i];
			pattern[i].forwardStrand->bForwardStrand = true;
			pattern[i].forwardStrand->compRules = compRules;
		}
		if (pattern[i].reverseStrand != NULL) {
			pattern[i].reverseStrand->pattern = &pattern[i];
			pattern[i].reverseStrand->bForwardStrand = false;
			pattern[i].reverseStrand->compRules = revComplementarityRules;
		}
	}

}

void printAlignment(Match *match, char *cv, unsigned int uik, AffixArray *affixArray) {
	unsigned int uitemp, uitemp2;
	bool **iupacTable = affixArray->alphabet->iupacTable;
	int iSeqPos;

	char *cSeqOperions           = match->cSeqOperations;
	char *cArcOperations         = match->cArcOperations;
	unsigned int uiNumOperations = match->uiNumOperations;

	char *cu;
	char *cu_str;

	cu     = match->patternStrandDirection->seq;
	cu_str = match->patternStrandDirection->structure;

	printf("Pattern ");
	// Print structure
	uitemp2 = 0;
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperions[uitemp] == OpInsertion)
			printf("%c", '-');
		else if (cSeqOperions[uitemp] == OpDeletion || cSeqOperions[uitemp] == OpReplacement)
			printf("%c", cu_str[uitemp2++]);
		else
			printf("ERROR");
	}
	printf("\n");

	printf("        ");
	// Print pattern
	uitemp2 = 0;
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperions[uitemp] == OpInsertion)
			printf("%c", '-');
		else if (cSeqOperions[uitemp] == OpDeletion || cSeqOperions[uitemp] == OpReplacement)
			printf("%c", c(cu[uitemp2++]));
		else
			printf("ERROR");
	}
	printf("\n");

	printf("        ");
	// Print additional line
	uitemp2 = 0;
	iSeqPos = uik;
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperions[uitemp] == OpInsertion || cSeqOperions[uitemp] == OpDeletion) {
			printf("%c", ' ');
		} else if (cSeqOperions[uitemp] == OpReplacement) {
			printf("%c", cu[uitemp2] == cv[iSeqPos] ? '|' : (iupacTable[(unsigned int) cu[uitemp2]][(unsigned int) cv[iSeqPos]] ? '+' : ' '));
		} else {
			printf("ERROR");
		}
		if (cSeqOperions[uitemp] == OpDeletion || cSeqOperions[uitemp] == OpReplacement)
			uitemp2++;
		if (cSeqOperions[uitemp] == OpInsertion || cSeqOperions[uitemp] == OpReplacement)
			iSeqPos++;
	}
	printf("\n");

	printf("Target  ");
	iSeqPos = uik;
	// Print database sequence
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		if (cSeqOperions[uitemp] == OpDeletion)
			printf("%c", '-');
		else if (cSeqOperions[uitemp] == OpInsertion || cSeqOperions[uitemp] == OpReplacement)
			printf("%c", c(cv[iSeqPos++]));
		else
			printf("ERROR");
	}
	printf("\n");

	printf("        ");
	// Print database structure
	for (uitemp = 0; uitemp < uiNumOperations; uitemp++) {
		printf("%c", cArcOperations[uitemp]);
	}
	printf("\n");
}


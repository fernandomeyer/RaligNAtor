/*Copyright (C) 2011  Fernando Meyer

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
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h> //includes function getpagesize()
#include <fcntl.h> //includes function open()
#include <sys/mman.h>

#include "af.h"

/* Initializes affixArray. Parameter bConstruct must be set to true if affixArray will be used
 * for affix array construction */
void init(AffixArray *affixArray, Alphabet *alphabet, char *argAlphabetFile, int predefAlphabet,
			bool *argSuf, bool *argLcp, bool *argSufinv, bool bConstruct) {

	/* Consistency checking */
	if(*argLcp || *argSufinv)
		*argSuf = 1;

	affixArray->alphabet = alphabet;

	if((affixArray->multiSeq = (MultiSeq *) malloc(sizeof(MultiSeq))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"affixArray->multiSeq\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	affixArray->multiSeq->convSequences = NULL;
	affixArray->multiSeq->convSequencesMmapped = false;

	if (*argSuf || bConstruct) {
		if((affixArray->esa = (EArray *) malloc(sizeof(EArray))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"esa\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		if((affixArray->erpa = (EArray *) malloc(sizeof(EArray))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"erpa\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}

		affixArray->esa->xarray = NULL;
		affixArray->esa->xlcp = NULL;
		affixArray->esa->xlcpException = NULL;
		affixArray->esa->xskp = NULL;
		affixArray->esa->xsufinv = NULL;
		affixArray->esa->affixLink = NULL;

		affixArray->erpa->xarray = NULL;
		affixArray->erpa->xlcp = NULL;
		affixArray->erpa->xlcpException = NULL;
		affixArray->erpa->xskp = NULL;
		affixArray->erpa->xsufinv = NULL;
		affixArray->erpa->affixLink = NULL;

	} else {
		affixArray->esa  = NULL;
		affixArray->erpa = NULL;
    }

	if (alphabet != NULL) {
		if (predefAlphabet > -1) {
			if (argAlphabetFile != NULL) {
				loadAlphabetFile(argAlphabetFile, alphabet);
			} else {
				if (predefAlphabet != 0 && predefAlphabet != 1 && predefAlphabet != 2)
					setPredefinedAlphabet(0, alphabet);
				else
					setPredefinedAlphabet(predefAlphabet, alphabet);
			}
		}
	}

}

void freeAll(AffixArray *affixArray, bool **complementarityRules) {
	unsigned int ui;
	//char *iupac_characters  = "RYSWKMBDHVN*";

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

		if (complementarityRules != NULL) {
			for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
					free(complementarityRules[ui]);
			}
			free(complementarityRules);
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

		if (affixArray->esa->xskp != NULL)
			munmap(affixArray->esa->xskp, affixArray->length * sizeof(int));
		if (affixArray->esa->xsufinv != NULL)
			munmap(affixArray->esa->xsufinv, affixArray->length * sizeof(int));
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

		if (affixArray->erpa->xskp != NULL)
			munmap(affixArray->erpa->xskp, affixArray->length * sizeof(int));
		if (affixArray->erpa->xsufinv != NULL)
			munmap(affixArray->erpa->xsufinv, affixArray->length * sizeof(int));
		free(affixArray->erpa);
	}
}

bool loadFastaFile(AffixArray *affixArray, char *fileName, bool bAllocConvSeq) {
	FILE *fp;
	char temp, *info;
	int i;

	if((fp=fopen(fileName,"r")) == NULL) {
		fprintf(stderr,"Error opening file \"%s\". If this is an index, please select the desired tables.\n", fileName);
		return 1;
	}

	if (fseek(fp, 0L, SEEK_END)) {
		fprintf(stderr,"Error reading file \"%s\".\n", fileName);
		return 1;
	}

	i=ftell(fp);
	if (i == 0) {
		fprintf(stderr, "File \"%s\" is empty.\n", fileName);
    	return 1;
	}

	rewind(fp);

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr, "Error. \"multiSeq\" was not initialized - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	affixArray->multiSeq->convSequences  = NULL;
	if ((affixArray->multiSeq->sequences = (unsigned char *) calloc(BUFFER1, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}
	if ((affixArray->multiSeq->seqEndPos      = (int *) calloc(BUFFER3, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}
	if ((affixArray->multiSeq->seqDescription = (char **) calloc(BUFFER3, sizeof(char*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqDescription\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}
	if ((affixArray->multiSeq->seqDescLength  = (int *) calloc(BUFFER3, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqDescLength\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
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
	affixArray->multiSeq->seqEndPos[affixArray->multiSeq->numSeqs-1] = affixArray->length - 1;

	fclose(fp);

	if (bAllocConvSeq && (affixArray->multiSeq->convSequences = (unsigned char *) calloc(affixArray->length, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"convSequences\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	return 0;
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

bool saveAffFiles(AffixArray *affixArray, char *fileName, bool bSaveTransfSeq) {
	FILE *fp;
	int i;

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr,"Please load a file.\n");
		return 1;
	}

	int nameLength = strlen(fileName);
	char fileNameCpy[nameLength + 7];
	strcpy(fileNameCpy, fileName);

	/********* base file *********/
	strcat(fileNameCpy, ".base");
	if((fp = fopen(fileNameCpy, "wb")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
	   return 1;
	}

	if(fwrite(&affixArray->length, sizeof(int), 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	if(fwrite(&affixArray->multiSeq->numSeqs, sizeof(int), 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	if(fwrite(affixArray->multiSeq->seqEndPos, sizeof(int), affixArray->multiSeq->numSeqs + 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* alphabet file *********/
	if(affixArray->alphabet->numClasses > 0) {
		strcat(fileNameCpy, ".alph");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}

		if(fwrite(&affixArray->alphabet->numClasses, sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}

		for(i = 0; i < affixArray->alphabet->numClasses; i++) {
			if(fwrite(&affixArray->alphabet->classSize[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if(fwrite(&affixArray->alphabet->classRepresentative[i], sizeof(char), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if (affixArray->alphabet->classSize[i] > 0) {
				if(fwrite(affixArray->alphabet->eqClass[i], sizeof(char), affixArray->alphabet->classSize[i], fp) < 1) {
					fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
					return 1;
				}
			}
			if(fwrite(&affixArray->alphabet->isWildCard[i], sizeof(bool), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		}

		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* description file *********/
	strcat(fileNameCpy, ".des");
	if((fp = fopen(fileNameCpy, "wb")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
	   return 1;
	}
	if(fwrite(affixArray->multiSeq->seqDescLength, sizeof(int), affixArray->multiSeq->numSeqs, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	for(i=0; i < affixArray->multiSeq->numSeqs; i++)
		if(fwrite(affixArray->multiSeq->seqDescription[i], sizeof(char), affixArray->multiSeq->seqDescLength[i], fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* sequences file *********/
	strcat(fileNameCpy, ".seq");
	if((fp = fopen(fileNameCpy, "wb")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
	   return 1;
	}
	/*if(fwrite(affixArray->multiSeq->seqEndPos, sizeof(int), affixArray->multiSeq->numSeqs + 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}*/
	if(fwrite(affixArray->multiSeq->sequences, sizeof(char), affixArray->length, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* alphabetically transformed sequences file *********/
	if (bSaveTransfSeq) {
		strcat(fileNameCpy, ".tseq");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(affixArray->multiSeq->convSequences, sizeof(char), affixArray->length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	return 0;
}

bool saveEArray(EArray* earray, int length, bool isEsa, char *fileName) {
	FILE *fp;
	int i;

	int nameLength = strlen(fileName);
	char fileNameCpy[nameLength + 7];
	strcpy(fileNameCpy, fileName);

	/********* suf | rpref file *********/
	if(earray->xarray != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".suf");
		else
			strcat(fileNameCpy, ".sufr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->xarray, sizeof(int), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* lcp | rlcp file *********/
	if(earray->xlcp != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".lcp");
		else
			strcat(fileNameCpy, ".lcpr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->xlcp, sizeof(char), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}


	/********* lcp | rlcp exception file *********/
	if(earray->xlcp != NULL && earray->xlcpException->numExceptions > 0) {
		if(isEsa)
			strcat(fileNameCpy, ".lcpe");
		else
			strcat(fileNameCpy, ".lcper");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}

		if(fwrite(&earray->xlcpException->numExceptions, sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		for(i=0; i < earray->xlcpException->numExceptions; i++) {
			if(fwrite(&earray->xlcpException->index[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if(fwrite(&earray->xlcpException->value[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* sufinv file *********/
	if(earray->xsufinv != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".sufinv");
		else
			strcat(fileNameCpy, ".sufinvr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		if(fwrite(earray->xsufinv, sizeof(int), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	return 0;
}

void allocConvSequences (AffixArray *affixArray) {
	if ((affixArray->multiSeq->convSequences = (unsigned char *) calloc(affixArray->length, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"convSequences\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
}

bool printEArray(AffixArray *affixArray) {
	int i;
	unsigned char charValue;
	unsigned char *seqTmp;
	EArray *earray;

	if(affixArray->length == 0) {
		fprintf(stderr,"Please load a file.\n");
		return 1;
	}

	earray = affixArray->esa;
	if(earray->xarray == NULL) return 0;

	printf("i	suf[i]	");
	if(earray->xlcp != NULL)
		printf("lcp[i]	");
	if(earray->xskp != NULL)
		printf("skp[i]	");
	if(earray->xsufinv != NULL)
		printf("sufinv[i]");
	if(earray->affixLink != NULL)
		printf("aflk[i]");
	printf("Ssuf[i]\n");

	for(i = 0; i < affixArray->length; i++) {
		printf("%d	%d	",i, earray->xarray[i]);
		if(earray->xlcp != NULL)
			printf("%d	", xlcpvalue(i));
		if(earray->xskp != NULL)
			printf("%d	", earray->xskp[i]);
		if(earray->xsufinv != NULL)
			printf("%d	", earray->xsufinv[i]);
		if(earray->affixLink != NULL)
			earray->affixLink[i] != INT_MAX ? printf("%d	", earray->affixLink[i]) : printf("	");

		seqTmp = affixArray->multiSeq->convSequences + earray->xarray[i];
		do {
			charValue = *seqTmp++;
			if(charValue == 127)
				printf("$");
			else
				printf("%c", affixArray->alphabet->classRepresentative[charValue - 1]);
		} while (charValue != $);
		printf("\n");
	}

	return 0;
}

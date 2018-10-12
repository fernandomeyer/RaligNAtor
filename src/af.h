/*Copyright (C) 2012  Fernando Meyer

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

#ifndef ESA_H_
#define ESA_H_ 1

#include "stdbool.h"

#define BUFFER1 1000
#define BUFFER2 100
#define BUFFER3 10
#define $ 127
#define LINESYMBOL '!'

typedef struct {
	int  numSeqs;
	unsigned char *sequences;
	unsigned char *convSequences; //alphabetically converted sequences
	bool sequencesMmapped;
	bool convSequencesMmapped;
	int  *seqEndPos; //ending position of each sequence + '$' in the concatenated array
	char **seqDescription;
	int  *seqDescLength;
} MultiSeq;

typedef struct {
	int *index;
	int *value;
	int numExceptions;
} LcpException;

typedef struct {
	bool bPredefined;
	char **eqClass;
	int  *classSize;
	unsigned char *classRepresentative;
	int  numClasses;
	bool *isWildCard;
	bool **iupacTable;
} Alphabet;

typedef struct {
	int *xarray; 					//suf or rpref
	unsigned char *xlcp; 		//lcp or rlcp
	LcpException  *xlcpException; 	//lcpException or rlcpException
	int *xskp; 					//skp or rskp
	int *affixLink;
	int *xsufinv;
} EArray;

typedef struct {
	MultiSeq *multiSeq;
	Alphabet *alphabet;
	int      length;
	EArray   *esa;
	EArray   *erpa;
} AffixArray;


typedef struct {
	int i;
	//Int32 j; //save space by not storing j
	int lcp;
} LInterval;

typedef struct {
	bool *isPairing;
	bool *isLeftRightExt;
	unsigned int length;
} StructureInfo;

typedef struct {
	unsigned char extension[100];
	unsigned char toRight;
} StringExtension;

#define lcpvalueX(data, index) (data->xlcp[index] == 255 ? \
		xlcpExceptionValue(data->xlcpException, index) : data->xlcp[index])

#define lcpvalue(index) (affixArray->esa->xlcp[index] == 255 ? \
		xlcpExceptionValue(affixArray->esa->xlcpException, index) : affixArray->esa->xlcp[index])

#define rlcpvalue(index) (affixArray->erpa->xlcp[index] == 255 ? \
		xlcpExceptionValue(affixArray->erpa->xlcpException, index) : affixArray->erpa->xlcp[index])

#define xlcpvalue(index) (earray->xlcp[index] == 255 ? \
		xlcpExceptionValue(earray->xlcpException, index) : earray->xlcp[index])

#define c(A) affixArray->alphabet->classRepresentative[(A) - 1]

//afCons.c
void affXSuf(EArray*, unsigned char*, int);
int  xlcpExceptionValue(LcpException*, int);
void affXLcp(EArray*, Alphabet*, unsigned char*, int);

void computeAffixArray(AffixArray*, bool, bool, bool, bool, char*, bool);
void computeEArray(EArray*, Alphabet*, unsigned char*, int, int, bool, bool, bool);

//afAlphabet.c
void  setPredefinedAlphabet(int, Alphabet*);
bool convertToAlphabet(unsigned char*, unsigned char*, int, bool, Alphabet*);
bool convertToRepresentative(char*, int, Alphabet*);
bool** setIupacTable (Alphabet*);

//afStreamHandling.c
void init(AffixArray*, Alphabet*, char*, int, bool*, bool*, bool*, bool);
bool loadAffFiles(AffixArray*, char*, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool);
void freeAll(AffixArray*, bool**);
bool loadFastaFile(AffixArray*, char*, bool);
Alphabet* loadAlphabetFile(char*, Alphabet*);
bool saveAffFiles(AffixArray*, char*, bool);
bool saveEArray(EArray*, int, bool, char*);
bool saveAflk(EArray*, int, bool, char*);
void allocConvSequences (AffixArray*);
bool printEArray(AffixArray*);

#endif /*ESA_H_*/

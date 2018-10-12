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

#ifndef ALIGN_H_
#define ALIGN_H_ 1

#include "stdbool.h"
#include "chaining/chain2dim.h"

#define V_ESA '0'
#define V_LESA '1'
#define V_LGSLINK '2'
#define V_LSLINK '3'
#define V_GSLINK '4'
#define V_SLINK '5'
#define V_SCAN '6'
#define V_SCANLA '7'
#define V_STD '8'
#define V_STDARRAY '9'
#define V_GLOBAL 'z'

#ifndef GLOBALCOUNT
#define GLOBALCOUNT 0
unsigned long int globalcount;
#endif

#define sigma(A, B) ((iupacTable[(unsigned int) (A)][(unsigned int) (B)]) ? 0 : 1)
//#define sigma(A, B) ((A) == (B) ? 0 : 1)

#define complement(A, B) compRules[(A) - 1][(B) - 1]

#define KIJ 0
#define P(MATRIX, A, B, C) MATRIX[B][A][C]
#define P1(MATRIX, A, B) MATRIX[B][A]
//#define P(MATRIX, A, B, C) MATRIX[A][B][C]
//#define P1(MATRIX, A, B) MATRIX[A][B]

#define S(MATRIX, A, B, C) MATRIX[A][C][B]

#define A(ARRAY, Kfactor, Jfactor, K, I, J) ARRAY[(Kfactor[K]) + (Jfactor) * (J) + (I)]
//#define A(ARRAY, Kfactor, Jfactor, K, I, J) ARRAY[P(position, K, I, J)]

#define BUFFER1 1000
#define BUFFER2 100
#define BUFFER3 10
#define $ 127
#define LINESYMBOL '!'

#define FILTERFREQ 21

#define c(A) affixArray->alphabet->classRepresentative[(A) - 1]

#define OpReplacement  1
#define OpInsertion    2
#define OpDeletion     3
#define OpArcBreaking  4
#define OpArcAltering1 5
#define OpArcAltering2 6
#define OpArcRemoving  7
#define OpBInsertion1  8
#define OpBInsertion2  9
#define OpArcBreakingTrue 10

#define lcpvalueX(data, index) (data->xlcp[index] == 255 ? \
		xlcpExceptionValue(data->xlcpException, index) : data->xlcp[index])

#define lcpvalue(index) (affixArray->esa->xlcp[index] == 255 ? \
		xlcpExceptionValue(affixArray->esa->xlcpException, index) : affixArray->esa->xlcp[index])

typedef struct rbtree rbtree;
typedef struct Match Match;
typedef struct Pattern Pattern;

typedef struct {
	int iL;
	int iR;
} LR;

typedef struct {
	int iCostReplacement;
	int iCostDeletion;
	int iCostArcBreak;
	int iCostArcAltering;
	int iCostArcRemoving;

	int iFCostReplacement; //cost for filtering phase
	int iFCostDeletion; //cost for filtering phase
} Cost;

typedef struct {
	int iL;
	int iR;
	bool bIsArc;
	bool bIsUnpaired;
	bool bDotBeforeArc;  //true if arc is preceded by an unpaired position
	bool bArcAfterDot; //true if unpaired region precedes an arc
	bool bIsLeftBulge; //is a left bulge or first unpaired region of a multi-loop
	//bool bIsMerge;
	int iLineToMerge1;

	int *iRegionDependsOnRegions;
	//int *iKLastColumnComputed; // Last computed column of matrix k
	int ikStart;
	int ikEnd;
	int ikEndFactual; //last truly computed matrix
	int ikBigStart;
	int ikBigEnd;

	int *ijGoal; // last column to be computed for matrix k
	int *ijBigGoal; // left-most column to be computed for matrix k of all regions up to this region regardless of dependency
	int **ijGoalPerRowPerMatrix; // last column to be computed for matrix k at row i
} PatternRegion;

typedef struct {
	int *ijLast;
	int **ijLastPerRowPerMatrix;
	int *ikFirst; //matrix k where a computed arc alignment begins
} ComputedEntries;

typedef struct {
	PatternRegion *patternRegion;
	PatternRegion *patternRegionInsideOut;
	int *iPosRegion;  // the region a pattern position belongs to
	int *iPosRegionInsideOut;  // the region a pattern position belongs to
	int *iRegionRegion; //[iPosRegionInsideOut] = iPosRegion
	int *iRegionRegionInv; //[iPosRegion] = iPosRegionInsideOut
	int *iRegion; //[iPosRegion] = iPosRegion
	int iNumRegions;
	int *iLineIndex;  // line indices for reading / writing
	int *iLeftMostPos; // iLeftMostPos stores the left position of the pattern that is already aligned/computed up to the moment when the given right position is aligned
						// It is only defined for right alignment extensions of the pattern

	int **iLCC; // Last computed column x of row i of matrix k, i.e. iLCC[k][i] = x
	int **iLargestRegionIndex;
	int **iLargestRegionIndexInsideOut;
} PatternStructures;

typedef struct {
	bool bForwardStrand; // current searched strand
	char *seq;
	char *structure;
	bool **compRules;
	PatternStructures *patternStructures;
	Pattern *pattern;
	Match **matches;
} PatternStrandDirection;

struct Pattern {
	char *desc;
	PatternStrandDirection *forwardStrand;
	PatternStrandDirection *reverseStrand;

	int iId;
	int iLength;
	int  weight;
	int  startpos;
	int  instance;
	unsigned int uiThreshold;
	unsigned int uiIndels;
	unsigned int uiMinDiffMatches;
	Cost *cost;
};

typedef struct {
	Pattern *pattern;
	int     iNumPatterns;
	bool bUsesStartPosParam; // local chaining
} MultiPattern;

typedef struct {
	int  numSeqs;
	char *sequences;
	char *convSequences; //alphabetically converted sequences
	bool sequencesMmapped;
	bool convSequencesMmapped;
	unsigned int  *seqEndPos; //ending position of each sequence + '$' in the concatenated array
	char **seqDescription;
	int  *seqDescLength;
} MultiSeq;

typedef struct {
	int *index;
	int *value;
	int numExceptions;
} LcpException;

typedef struct {
	char **eqClass;
	int  *classSize;
	unsigned char *classRepresentative;
	unsigned char * conversionTablePattern;
	unsigned char * conversionTableTargetSeq;
	int  numClasses;
	bool *isWildCard;
	bool **iupacTable;
} Alphabet;

typedef struct {
	int  *xarray; //suf or rpref
	unsigned char *xlcp; //lcp or rlcp
	//int  *xskp; //skp or rskp
	int  *xsufinv;
	int  *xlnk;
	int           *affixLink;
	LcpException  *xlcpException; //lcpException or rlcpException
} EArray;

typedef struct {
	MultiSeq *multiSeq;
	Alphabet *alphabet;
	int      length;
	EArray   *esa;
	EArray   *erpa;
} AffixArray;

struct Match {
	PatternStrandDirection *patternStrandDirection;
	int iSeqId;
	int iCost;
	int iScore;
	unsigned int uiPos;
	unsigned int uiEndPos;
	unsigned int uiNumOperations; // number of edit operations
	char *cSeqOperations;
	char *cArcOperations;
	bool bForwardStrand;
};

typedef struct {
	unsigned int startpos[2];
	unsigned int endpos[2];
	int weight;
	Match *match;
} ChainFragment;

typedef struct {
	int  seqId;
	int  chainLength;
	int  chainScore;
	ChainFragment *fragmentInfo;
	bool forwardSeq;
	bool bStoredInReverseOrder;
} Chain;

typedef struct {
	bool isactive;
	char *reportfile;
	bool islocal;
	double weightfactor;
	GtChain2Dimpostype maxgap;
	char *localparm, *globalparm, *origopts;
	bool show;
	bool show2;
	unsigned int minchainlength;
	unsigned int minchainscore;
	unsigned int topkscoring;
} Chainparam;

typedef struct {
	unsigned int uiThreshold;
	unsigned int uiIndels;
	unsigned int uiMinDiffMatches;
	Cost cost;
	char cVariant;
	char cPrintMatchesByScore; //0=no sorting, 1=descending, 2=ascending
	bool bSearchForwardString;
	bool bSearchReverseString;
	bool bUseSequenceBasedFilter;
	bool bPrintOutMatches;
	bool bPrintMatchesBySeq;
	bool bPrintTable;
	bool bBED;
	bool bReportAllMatches;
	bool bIncludeSeqDesc;
	bool bStoreInRB; // store match in red-black tree? chainparam->isactive || bPrintMatchesBySeq
	bool bFilterOverlaps;
	bool bShowProgress;
	Chainparam *chainParam;
	unsigned int *uiNumMatchesSeq;
	unsigned int *uiNumMatchesReverseSeq;
} SearchParam;

// init.c
void init(AffixArray *affixArray, SearchParam *searchParam,
		bool *argSuf, bool *argLcp, bool *argSufinv, bool *argLnk, bool *argAflk, bool *argSufr, bool *argLcpr, bool *argAflkr,
		Alphabet *alphabet, char *argAlphabetFile, const int predefAlphabet);
bool adjustNumIndelsAccordingToEdist (unsigned int* uiIndels, unsigned int* uiThreshold, Cost *cost);
void freeAll(AffixArray *affixArray, MultiPattern *multiPattern);
Cost* getUnitCosts();

unsigned int*** uiNew3DMatrix(unsigned int uim, unsigned int uin);
unsigned int*** uiNew3DMatrixX(unsigned int uim, unsigned int uin);
int uiFree3DMatrixX(unsigned int ***matrix, unsigned int uim, unsigned int uin);

int uiFree3DMatrix(unsigned int ***matrix, unsigned int uim, unsigned int uin);
unsigned int*** uiNew3DMatrixWithEntryCopies(unsigned int uim, unsigned int uin, char *cu_str);
char*** cNew3DMatrixWithEntryCopies(unsigned int uim, unsigned int uin, char *cu_str);

char*** cNew3DMatrix(unsigned int uim, unsigned int uin);
char*** cNew3DMatrixX(unsigned int uim, unsigned int uin);
int cFree3DMatrixX(char ***matrix, unsigned int uim, unsigned int uin);

int cFree3DMatrix(char ***matrix, unsigned int uim, unsigned int uin);
int free3DMatrix(unsigned int ***matrix, unsigned int lengthX, unsigned int lengthY, unsigned int lengthZ);
int free3DMatrixC(char ***matrix, unsigned int lengthX, unsigned int lengthY, unsigned int lengthZ);
ComputedEntries *newComputedEntriesMatrix(PatternStrandDirection *patternStrandDirection, bool bInsideOut);
int *freeComputedEntriesMatrix(ComputedEntries *compEntries, PatternStrandDirection *patternStrandDirection, bool bInsideOut);
PatternStructures* processPattern(char *cu_str, unsigned int uim, unsigned int uiIndels);
void initializeRowsAndColumns(unsigned int ***edist, char ***operation, unsigned int ***trace, Cost *cost,
		char *cu_str, int im, int in, PatternStructures *patternStructures);
void initializeRowsAndColumnsRegions_S(unsigned int ***edist, char ***operation, unsigned int ***trace, Cost *cost,
		char *cu_str, int im, int in, PatternStructures *patternStructures);
void initializeRowsAndColumnsRegions_STEST(unsigned int ***edist, char ***operation, unsigned int ***trace, Cost *cost,
		char *cu_str, int im, int in, PatternStructures *patternStructures);

// alignSLink.c
void getEditOperations_S(char *cSeqOperations, char *cArcOperations, char ***operation, unsigned int ***trace, int *iIndels, int *iDeletions, const int iMaxIndels,
		PatternStructures *pstr, unsigned int *uiIndex, unsigned int uik, unsigned int uii, unsigned int uij, bool bAllowBranch);
int alignESA_LA_NEW_S(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int *ijLastOfRegion,
		char *cv, const int in, const int inMax, AffixArray *affixArray, int iMaxColOffset, int iLastCompRow, int iLastRegion, int *iRoundComputedPerRegion, bool *bSuccess, LR *LR, bool bInOutAlignment);
int alignESA_S(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, char *cv,
		const int in, AffixArray *affixArray, const int iMaxColOffset);
int alignLastColumn_S(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection,
		char *cv, const int in, AffixArray *affixArray);
int alignESA_LA_NEW_LastColumnS1(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int *iEmptyArray,
		char *cv, const int in, const int inMax, AffixArray *affixArray, LR *LR, int iCurrentRound, int *iRoundComputedPerRegion, bool *bSuccess, bool bInOutAlignment, int iLastRegion, int *iLastRegionCheck);

// align.c
void applySequenceBasedFilter(AffixArray *affixArray, PatternStrandDirection *patternStrandDirection, bool *bVisitedSuffix, int iStartIdx, int iEndIdx);
int alignESA_LA_NEW(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, ComputedEntries *compEntries, int *ijLastOfRegion,
		char *cv, const int in, const int inMax, AffixArray *affixArray, const int iMaxColOffset, int iDepthEnded, bool bFindFailingRow, bool *bSuccess, LR *LR);
int alignESA(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *PatternStrandDirection,	ComputedEntries *compEntries, int **ijLastOfRegion, char *cv,
		const int in, const int inMax, AffixArray *affixArray, const int iMaxColOffset, int iDepthEnded, bool bCheckFailingRow, bool *bSuccess, LR *LR);
inline unsigned int xlcpExceptionValue(LcpException *xlcpException, const unsigned int k);
int getSeqNumber(unsigned int *range, int length, unsigned int k);

// alignOnline.c
void getEditOperations(char *cSeqOperations, char *cArcOperations, char ***operation, unsigned int ***trace, int *iIndels, int *iDeletions, const int iMaxIndels,
		PatternStructures *pstr, unsigned int *uiIndex, unsigned int uik, unsigned int uii, unsigned int uij, bool bAllowBranch);
unsigned int alignGlobal(PatternStrandDirection *patternStrandDirection, AffixArray *affixArray, SearchParam *searchParam);

// alphabet.c
void setPredefinedAlphabet(int type, Alphabet *alphabet);
int setAlphabetConversionTable(Alphabet* alphabet);
bool** setIupacTable (Alphabet *alphabet);
bool convertToAlphabet(char *seq, char *convSeq, int length, bool isPattern, Alphabet* alphabet);
bool** loadReverseComplementarityRules(bool **bCompCheck, char *cArrayComplementaryCode, AffixArray *affixArray);

// streamHandling.c
int loadFastaFile(AffixArray *affixArray, char *fileName, bool bAllocConvSeq);
bool loadAffFiles(AffixArray*, char*, bool, bool, bool, bool, bool, bool, bool, bool);
MultiPattern* loadPatternFile(char *fileName, SearchParam *searchParam, AffixArray *affixArray);
Alphabet* loadAlphabetFile(char *fileName, Alphabet *alphabet);
void allocConvSequences (AffixArray *affixArray);
bool** loadComplementarityFile(char *fileName, AffixArray *affixArray);
bool** loadDefaultComplementarityRules(AffixArray *affixArray, const int predefAlphabet);
bool** getReverseComplementarityTable(bool **bCompCheck, AffixArray *affixArray);
void setSearchPatterns(MultiPattern *multiPattern, bool **compRules, SearchParam *searchParam, AffixArray *affixArray);
void printAlignment(Match *match, char *cv, unsigned int uik, AffixArray *affixArray);
void printMatrices(unsigned int ***edist, char ***operation, char *cu, char *cv, unsigned int uim, unsigned int uin, AffixArray *affixArray);
void printMatrices_S(unsigned int ***edist, char ***operation, char *cu, char *cv, unsigned int uim, unsigned int uin, AffixArray *affixArray);
void printMatrices_A(unsigned int *edist, char *operation, unsigned int ***position, unsigned int *Kfactor, unsigned int Jfactor, char *cu, char *cv, unsigned int uim, unsigned int uin, AffixArray *affixArray);
void printfMatrices_S_up_toLastRegion(unsigned int ***edist, char ***operation, char *cu, char *cv, AffixArray *affixArray, int iLastRegion, PatternStrandDirection *patternStrandDirection);

// resultsProcessing.c
void chainAll(void *vrb, void *vrbCR, unsigned int *uiNumMatchesSeq, SearchParam *searchParam,
		MultiPattern *multiPattern, AffixArray *affixArray, bool bFreeMatches, bool bForwardSeq);
void printMatch(Match *match, AffixArray *affixArray);
void printSortedMatches(void *vrbF, void *vrbR, AffixArray *affixArray, SearchParam *searchParam);
int compareRB(const void *pa, const void *pb, const void *config);
int compareRBCR(const void *pa, const void *pb, const void *config);
void reportChainingResults(void *vrb, Chainparam *chainparam, AffixArray *affixArray);
void printMatchInOneLine(Match *match, AffixArray *affixArray);
void sortMatchesByPos(Match **match, int left, int right);
unsigned int removeDuplicateMatches(Match **match, const int iMatches);

#endif /*ALIGN_H_*/


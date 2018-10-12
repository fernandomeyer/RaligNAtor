/*Copyright (C) 2012  Ole Eigenbrod, Clemens Hoebsch

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

/*
*In this header file we declare 2 structs to manage the pthreads and the pointers for the Red-Black tree.
*The parameters must be stored to call the function alignThreaded.
*The threadmanage struct needs to store information for every thread that will be called.
*Additionally in the headerfiles there are definitions of global variables
*for the mutexes (pthread) and the thread counter.
*initThread is responsible for the initialization of the threadManage structs.
*The function initRBPointers initializes the Red-Black tree.
*CallChaining function will perform chaining(not modified) if the parameters for that are set.
*The function freeTmanage releases the space that was reserved for the threadManage structs.
*The function insertMatchesRB removes overlapping matches and inserts the unique ones into the rb-tree.
*The function matchPrinter is responsible for printing out all matches.
*The function alignESABasedThreaded computes matches for a given pattern in a given index partition with the LESA variant.
*The functions alignOnlineScanMainThreaded, alignOnlineScanMainLA_NEWThreaded compute matches for a given pattern using the
*SCAN and LSCAN variants, respectively.
*The function generalSLinkAlignNEWThreaded computes matches for a given pattern using LGSLINK.
*The function threadMain is the main method in this module and calls the other functions, where
*alignThreaded is the most important one.
*alignThreaded is responsible for aligning the patterns to the sequences.
*/

#ifndef ALIGNTHREAD_H_
#define ALIGNTHREAD_H_ 1

pthread_mutex_t mutexPrintSTDOUT;
pthread_mutex_t mutexPrintSTDERROR;
pthread_mutex_t mutexThreadCounter;

pthread_cond_t condThreadCount;

typedef struct ThreadsManager ThreadsManager;

typedef struct {
	pthread_mutex_t mutexTree;
	pthread_mutex_t mutexMatchCount;
} Mutex;

typedef struct {
	rbtree *rbF;
	rbtree *rbR;
	rbtree *rbCR;
} RBPointers;

typedef struct {
	int iThreadID;
	AffixArray *affixArray;
	SearchParam *searchParam;
	RBPointers *rbPointers;
	PatternStrandDirection *patternStrandDirection;

	int startIdx;
	int endIdx;
	int iStartPos;
	int iEndPos;

	int *iFilterOverlapsControl;
	bool *bVisitedSuffix;

	Mutex *mutex;
	Match **matches;
	unsigned int uiNumMatchesPthread;

	ThreadsManager *threadsManager;
} ThreadVars;

struct ThreadsManager {
	int iThreadCounter;
	int iNumThreads;

	int iNumJobs;

	Match ***matchesPerPatternFor; //[patID][arrayPos][match]
	Match ***matchesPerPatternRev; //[patID][arrayPos][match]
	unsigned int *uiNumMatchesPerPatternFor; //[patID]
	unsigned int *uiNumMatchesPerPatternRev; //[patID]

	ThreadVars *threadVars;
};

int *initCheckpoints(AffixArray *affixArray);
void chain(SearchParam *searchParam, RBPointers *rbPointers, MultiPattern *multiPattern, AffixArray *affixArray);
RBPointers *initRBPointers(SearchParam *searchParam, AffixArray *affixArray);
void * alignThreaded(void * tmanageIn); 
ThreadVars *initThreads(ThreadsManager *threadsManager, int numThreads, int iNumPartitions,
		MultiPattern *multiPattern, SearchParam *searchParam,
		AffixArray *affixArray, RBPointers *rbPointers, int nPerPart);
void freeThreads(ThreadVars *tmanage, int numPatterns, int iNumPartitions);
void threadMain(MultiPattern *multiPattern, AffixArray *affixArray, bool **compRules, SearchParam *searchParam, int numThreads, int partifactor);
void *generalSLinkAlignNEWThreaded(void *tmanageIn);
void *alignESABasedThreaded(void *params);
void *alignOnlineScanMainThreaded(void *tmanageIn);
void *alignOnlineScanMainLA_NEWThreaded(void *tmanageIn);
void insertMatchesRB(ThreadVars *tmanage);
void processAlignment_MT(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, const int iAbsoluteMatchPos,
		char *cSeqOperations, char *cArcOperations, const int iMinLength, const int iMaxLength, AffixArray *affixArray, ThreadVars *threadVars);
void processAlignment_S_MT(unsigned int ***edist, char ***operation, unsigned int ***trace, PatternStrandDirection *patternStrandDirection, const int iAbsoluteMatchPos,
		char *cSeqOperations, char *cArcOperations, const int iMinLength, const int iMaxLength, AffixArray *affixArray, ThreadVars *threadVars);


#endif /*ALIGNTHREAD_H_*/


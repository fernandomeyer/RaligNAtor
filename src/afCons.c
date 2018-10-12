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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>

#include "af.h"
#include "divsufsort/config.h"
#include "divsufsort/divsufsort.h"
#include "divsufsort/lfs.h"

void computeSufinvTable(AffixArray *affixArray, bool argShowTimes);

void computeAffixArray(AffixArray *affixArray, bool argSuf, bool argLcp, bool argSufinv,
					bool argOScreen, char *fileName, bool argShowTimes) {

	/* Compute Enhanced Suffix Array */
	computeEArray(affixArray->esa,
				affixArray->alphabet,
				affixArray->multiSeq->convSequences,
				affixArray->multiSeq->numSeqs,
				affixArray->length,
				argSuf, argLcp,
				argShowTimes);

	if (argSufinv)
		computeSufinvTable(affixArray, argShowTimes);

	if (fileName != NULL)
		if (saveEArray(affixArray->esa, affixArray->length, 1, fileName)){
			fprintf(stderr, "Error occurred while storing the suffix array");
			exit(1);
		}

	/* Free unnecessary structures */
	if(!argOScreen) {
		if (affixArray->esa->xarray != NULL) {
			free(affixArray->esa->xarray);
			affixArray->esa->xarray = NULL;
		}
		if(affixArray->esa->xlcp != NULL) {
			free(affixArray->esa->xlcp);
			affixArray->esa->xlcp = NULL;
		}
		if (affixArray->esa->xlcpException != NULL) {
			free(affixArray->esa->xlcpException);
			affixArray->esa->xlcpException = NULL;
		}
		if(affixArray->esa->xsufinv != NULL) {
			free(affixArray->esa->xsufinv);
			affixArray->esa->xsufinv = NULL;
		}
	}
}

void computeEArray(EArray *earray, Alphabet *alphabet, unsigned char *seq, int numSeqs, int length,
					bool argXarray, bool argXlcp, bool argShowTimes) {
	struct timeval start, end;

	if(argXarray) {
		printf("\nComputing suf... "); fflush(stdout);
		if (argShowTimes) gettimeofday(&start, NULL);
		affXSuf(earray, seq, length);
		if (argShowTimes) gettimeofday(&end, NULL);
		printf("done\n");
		if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
	}

	if(argXlcp) {
		printf("\nComputing lcp... "); fflush(stdout);
		if (argShowTimes) gettimeofday(&start, NULL);
		affXLcp(earray, alphabet, seq, length);
		if (argShowTimes) gettimeofday(&end, NULL);
		printf("done\n");
		if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
	}
}

void computeSufinvTable(AffixArray *affixArray, bool argShowTimes) {
	if(affixArray->esa->xarray == NULL) {
		fprintf(stderr,"lnk computation requires table suf.\n");
		exit(1);
	}
	struct timeval start, end;
	int i, iLength, *sufinv, *suf;

	if (argShowTimes) gettimeofday(&start, NULL);

	printf("\nComputing suf^-1... "); fflush(stdout);

	iLength = affixArray->length;
	suf = affixArray->esa->xarray;

	if((affixArray->esa->xsufinv = (int *) calloc(iLength, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xsufinv\" - %s %d.\n", __FILE__, __LINE__);
	sufinv = affixArray->esa->xsufinv;

	// Compute inverse suffix array
	for (i = 0; i < iLength; i++) {
		sufinv[suf[i]] = i;
	}

	if (argShowTimes) gettimeofday(&end, NULL);
	printf("done\n");
	if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
}

//Reverse sequences terminating with $
void reverseSequences(unsigned char *seq, int length) {
	int i;
	char temp;

	for(i = 0; i < length / 2; i++) {
		temp = seq[i];
		seq[i] = seq[length - 2 - i];
		seq[length - 2 - i] = temp;
	}
}

//Reverse arbitrary string and return a new pointer
unsigned char* reverseStringNewPointer(unsigned char *seq, unsigned int length) {
	unsigned int i, lengthBy2 = length / 2;
	unsigned char *newSeq = (unsigned char *) malloc((length + 1) * sizeof(unsigned char));

	for(i = 0; i < lengthBy2; i++) {
		newSeq[i] = seq[length - 1 - i];
		newSeq[length - 1 - i] = seq[i];
	}
	if (length % 2 != 0)
		newSeq[lengthBy2] = seq[lengthBy2];
	newSeq[length] = 0;

	return newSeq;
}

#define swap(x, a, b) { temp = x[a]; \
                     x[a] = x[b]; x[b] = temp; }

void sortLcpException(int *index, int *value, int left, int right) {
	int i, last, temp;
	if(left >= right)
		return;
	//i = (left + right)/2;
	i = (left + right) >> 1;
	swap(index, left, i);
	swap(value, left, i);
	last = left;
	for(i = left+1; i <= right; i++){
		if(index[i] < index[left]) {
			++last;
			swap(index, last, i);
			swap(value, last, i);
		}
	}
	swap(index, left, last);
	swap(value, left, last);
	sortLcpException(index, value, left, last-1);
	sortLcpException(index, value, last+1, right);
}

void affXSuf(EArray *earray, unsigned char *seq, int length) {
	if((earray->xarray = (int *) calloc(length, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xarray\".\n");

	if(divsufsort(seq, earray->xarray, (saidx_t)length) != 0) {
		fprintf(stderr, "Cannot allocate memory.\n");
		exit(1);
	}
}

int xlcpExceptionValue(LcpException *xlcpException, int k) {
	int low, high, mid;

	/* naive search */
	/*unsigned int i;
	for(i = 0; i < xlcpException->numExceptions; i++)
		if(xlcpException->index[i] == k)
			return xlcpException->value[i];*/

	low = 0;
	high = xlcpException->numExceptions - 1;

	while (low <= high) {
		mid = (low + high) >> 1;
		if (xlcpException->index[mid] < k)
			low = mid + 1;
		else if (xlcpException->index[mid] > k)
			high = mid - 1;
		else
			return xlcpException->value[mid];
	}

	return -1;
}

void affXLcp(EArray *earray, Alphabet *alphabet, unsigned char *seq, int length) {
	if(earray->xarray == NULL) {
		fprintf(stderr,"lcp (lcpr) computation requires table suf (sufr).\n");
		exit(1);
	}

	int i, j, k, h=0;
	int *rank = (int *) calloc(length, sizeof(int));

	if((earray->xlcp = (unsigned char *) calloc(length, sizeof(unsigned char))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xlcp\".\n");
	if((earray->xlcpException = (LcpException *) malloc(sizeof(LcpException))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xlcpException\".\n");

	earray->xlcpException->index = (int *) calloc(BUFFER2, sizeof(int));
	earray->xlcpException->value = (int *) calloc(BUFFER2, sizeof(int));
	earray->xlcpException->numExceptions = 0;

	for(i = 0; i < length; i++)
		rank[earray->xarray[i]] = i;

	for(i = 0; i < length; i++){
		k = rank[i];
		if(k == 0)
			earray->xlcp[k] = 0;
		else {
			j = earray->xarray[k - 1];

			while(seq[i + h] != $ && seq[j + h] != $ && seq[i + h] == seq[j + h] && seq[i + h] <  alphabet->numClasses && seq[j + h] < alphabet->numClasses)
				h++;

			if(h > 254) {
				earray->xlcp[k] = 255;

				if(earray->xlcpException->numExceptions % BUFFER2 == 0) {
					if((earray->xlcpException->index = realloc(earray->xlcpException->index, (earray->xlcpException->numExceptions + BUFFER2)*sizeof(int))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"xlcpException.index\".\n");
						exit(1);
					}
					if((earray->xlcpException->value = realloc(earray->xlcpException->value, (earray->xlcpException->numExceptions + BUFFER2)*sizeof(int))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"xlcpException.value\".\n");
						exit(1);
					}
				}
				earray->xlcpException->value[earray->xlcpException->numExceptions] = h;
				earray->xlcpException->index[earray->xlcpException->numExceptions] = k;
				earray->xlcpException->numExceptions++;
			} else
				earray->xlcp[k] = h;
		}
		if(h > 0)
			h--;
	}

	free(rank);

	if(earray->xlcpException->numExceptions > 0)
		sortLcpException(earray->xlcpException->index, earray->xlcpException->value, 0, earray->xlcpException->numExceptions - 1);
}

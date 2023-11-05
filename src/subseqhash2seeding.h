/*
  subseqhash2 seeding.

  By: Ke, Xiang@PSU
  Last edited: 10/6/2023
*/

#include "seeding.h"
#include "util.h"
#include <cstdint>
#include <inttypes.h>
#include <algorithm>
#include <vector>

#ifndef _SUBSEQHASH2SEEDING_H
#define _SUBSEQHASH2SEEDING_H

struct DPCell{
    int64_t f_max, f_min;
    int g_max, g_min;
};

#define MAXK 64
#define MAXD 31

class subseqhash2seeding: private seeding
{

private:

	int64_t A[MAXK][ALPHABETSIZE][MAXD];
	int B1[MAXK][ALPHABETSIZE][MAXD];
	int B2[MAXK][ALPHABETSIZE][MAXD];

	int64_t revA[MAXK][ALPHABETSIZE][MAXD];
	int revB1[MAXK][ALPHABETSIZE][MAXD];
	int revB2[MAXK][ALPHABETSIZE][MAXD];

	int64_t A3[MAXK][ALPHABETSIZE];

	int combine1[MAXK][ALPHABETSIZE];
	int combine2[MAXK][ALPHABETSIZE];
	int combine3[MAXK][ALPHABETSIZE];

	int C1[MAXK][ALPHABETSIZE];
	int C2[MAXK][ALPHABETSIZE];
	int C3[MAXK][ALPHABETSIZE]; 

	int num_valid = 0;
	int valid[MAXK] = {0};

	int d;
	int dim1, dim2, dim3;
	int chunk_size = 500;
	int64_t threshold;

	int dpIndex(int d1, int d2, int d3, int d4);
	int hIndex(int d2, int d3, int d4);
	void DP(std::string s, size_t start, size_t end, DPCell* dp, int* h);
	void revDP(std::string s, size_t start, size_t end, DPCell* revdp, int* revh);
	void combine(std::string s, size_t start, size_t end, DPCell* dp, DPCell* revdp, std::vector<std::vector<seed>>& seeds);

public:

	subseqhash2seeding(int n1, int k1, int d1, int subsample, int64_t threshold1 = ((int64_t)1<<63)): seeding(n1, k1)
	{
			d = d1;
	    threshold = threshold1;
	    num_valid = subsample;
	    dim3 = d;
	    dim2 = (k+1) * dim3;
	    dim1 = (n+1) * dim2;
	}

	void init(const char* table_filename);
	void getSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, int* h, int* revh,
			     std::vector<std::vector<seed>>& seeds);
	void writeSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, int* h, int* revh,
		     std::vector<std::vector<seed>>& seeds, std::vector<FILE*> fout);

	int getChunkSize()
	{
		return chunk_size;
	}	

	int getNumPerWindow()
	{
		return num_valid;
	}
};
#endif
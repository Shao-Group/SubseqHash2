/*
  subseqhash2 SIMD seeding, MAXD = 32.

  By: Xiang@PSU
  Last edited: 12/5/2023
*/

#include "seeding.h"
#include "util.h"
#include <bits/stdc++.h>
#include <x86intrin.h>

#ifndef _SUBSEQHASH2SIMD_32_H
#define _SUBSEQHASH2SIMD_32_H

#define MAXK 64
#define MAXD 32

struct DPCell{
	  short f_max[MAXD], g_max[MAXD];
    short f_min[MAXD], g_min[MAXD];
};


class subseqhash2SIMDseeding: private seeding
{

private:

	alignas(64) short A[MAXK][ALPHABETSIZE][MAXD];
	alignas(64) short B1[MAXK][ALPHABETSIZE][MAXD];
	alignas(64) short B2[MAXK][ALPHABETSIZE][MAXD];

	alignas(64) short revA[MAXK][ALPHABETSIZE][MAXD];
	alignas(64) short revB1[MAXK][ALPHABETSIZE][MAXD];
	alignas(64) short revB2[MAXK][ALPHABETSIZE][MAXD];

	short A3[MAXK][ALPHABETSIZE];

	short combine1[MAXK][ALPHABETSIZE];
	short combine2[MAXK][ALPHABETSIZE];
	short combine3[MAXK][ALPHABETSIZE];

	short C1[MAXK][ALPHABETSIZE];
	short C2[MAXK][ALPHABETSIZE];
	short C3[MAXK][ALPHABETSIZE]; 

	int num_valid = 0;
	int valid[MAXK] = {0};
	alignas(64) short shift[MAXD+1][MAXD];
	alignas(64) short comp[MAXD+1][MAXD];

	short d;
	int dim1, dim2;
	int chunk_size = 512;
	int64_t threshold;

	int dpIndex(int d1, int d2, int d3);
	int hIndex(int d2, int d3, int d4);
	void DP(std::string s, size_t start, size_t end, DPCell* dp, short* h);
	void revDP(std::string s, size_t start, size_t end, DPCell* revdp, short* h);
	void combine(std::string s, size_t start, size_t end, DPCell* dp, 
		DPCell* revdp, std::vector<std::vector<seed>>& seeds);

public:

	subseqhash2SIMDseeding(int n1, int k1, int d1, int subsample, int64_t threshold1 = ((int64_t)1<<63)): seeding(n1, k1)
	{
			d = d1;
	    threshold = threshold1;
	    num_valid = subsample;

	    dim2 = (k+1);
	    dim1 = (n+1) * dim2;

	    for(int i = 1; i <= d; i++)
	    {
	    	for(int j = 0; j < d; j++)
	    		shift[i][j] = ((i+j)%d);
	    	for(int j = d; j < MAXD; j++)
	    		shift[i][j] = j;
	    }

	    for(int i = 0; i < d; i++)
	    {
	    	for(int j = 0; j < d; j++)
	    		comp[i][j] = ((i - j + d)%d);
	    	for(int j = d; j < MAXD; j++)
	    		comp[i][j] = j;
	    }
	}

	void init(const char* table_filename);
	void getSubseq2Seeds(std::string s, DPCell* dp, short* h,
			DPCell* revdp, short* revh, std::vector<std::vector<seed>>& seeds);
	void writeSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, short* h, short* revh,
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
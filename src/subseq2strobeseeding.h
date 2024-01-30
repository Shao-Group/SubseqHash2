/*
  strobemer style subseqhash2 seeding, small kmer + subsequence

  By: Ke, Xiang@PSU
  Last edited: 10/07/2023
*/

#include "seeding.h"
#include "util.h"
#include <cstdint>
#include <inttypes.h>
#include <algorithm>
#include <vector>

#ifndef _SUBSEQ2STROBESEEDING_H
#define _SUBSEQ2STROBESEEDING_H

struct DPCell{
    int64_t f_max, f_min;
    int g_max, g_min;
};

#define MAXK 64
#define MAXD 31

class subseq2strobeseeding: public seeding
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
    int w, prek;

    int dpIndex(int d1, int d2, int d3, int d4);
    int hIndex(int d2, int d3, int d4);
    void DP(std::string s, size_t start, size_t end, DPCell* dp, int* h);
    void revDP(std::string s, size_t start, size_t end, DPCell* revdp, int* revh);
    void combine(std::string s, size_t start, size_t end, DPCell* dp, DPCell* revdp, std::vector<std::vector<seed>>& seeds);
	
public:

    subseq2strobeseeding(int n1, int k1, int d1, int subsample, int w1, int prek1): seeding(n1, k1){
	d = d1;
	num_valid = subsample;
	dim3 = d;
	dim2 = (k+1) * dim3;
	dim1 = (n+1) * dim2;
	w = w1;
	prek = prek1;
    }

    subseq2strobeseeding(int n1, int k1, int d1, int subsample,
			 int w1, int prek1,
			 const char* table_filename):
	seeding(n1, k1), num_valid(subsample),
	d(d1), dim3(d), w(w1), prek(prek1){
	dim2 = (k+1) * dim3;
	dim1 = (n+1) * dim2;
	init(table_filename);
    }

    void init(const char* table_filename);
    void getSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, int* h, int* revh,
			 std::vector<std::vector<seed>>& seeds);
    double getSeeds(std::string& s, const size_t s_idx,
		    const char* output_dir, const int dir_len,
		    DPCell* dp, DPCell* revdp, int* h, int* revh);
	void writeSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, int* h, int* revh,
		     std::vector<std::vector<seed>>& seeds, std::vector<FILE*> fout);

    int getChunkSize() const{
	return chunk_size;
    }	

    int getNumPerWindow() const{
	return num_valid;
    }
};
#endif

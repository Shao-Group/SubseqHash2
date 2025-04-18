/*
  subseqhash1 seeding.

  By: Ke, Xiang@PSU
  Last edited: 05/16/2023
*/

#include "seeding.hpp"
#include "util.hpp"
#include <cstdint>
#include <inttypes.h>
#include <algorithm>
#include <vector>
#include <cmath>

#ifndef _SUBSEQHASH1SEEDING_H
#define _SUBSEQHASH1SEEDING_H

struct DPCell{
    double f_max, f_min;
    int g_max, g_min;
};

#define MAXK 64
#define MAXD 31

class subseqhash1seeding: private seeding
{

private:

	double A[MAXK][ALPHABETSIZE][MAXD];
	int B1[MAXK][ALPHABETSIZE][MAXD];
	int B2[MAXK][ALPHABETSIZE][MAXD];

	int C[MAXK][ALPHABETSIZE];

	int subsample;
	int d;
	int dim1, dim2;

	int dpIndex(int d1, int d2, int d3);

public:

	subseqhash1seeding(int n1, int k1, int d1): seeding(n1, k1)
	{
			d = d1;
	    dim2 = d;
	    dim1 = (k+1) * dim2;
	}
	void DP(std::string s, DPCell* dp, int* h, std::vector<seed>& seeds);
	void init(const char* table_filename);
};
#endif
/*
  minimap2-style minimizer seeding.

  By: Ke, Xiang@PSU
  Last edited: 06/15/2023
*/

#include "seeding.h"
#include "util.h"
#include <cstdint>
#include <inttypes.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <list>

#ifndef _MINIMIZERSEEDING_H
#define _MINIMIZERSEEDING_H

struct MMentry{
	uint64_t hash;
	size_t pos;
	kmer str;
};

class minimizerseeding: private seeding
{

private:

	inline uint64_t hash64(kmer s, uint64_t mask);
	int w;

public:

	minimizerseeding(int k1, int w1): seeding(k1 + w1 - 1, k1)
	{
		w = w1;
	}
	void get_minimizers(std::string s, std::vector<seed>& seeds);

};
#endif
/*
  minimap2-style minimizer seeding.

  By: Ke, Xiang@PSU
  Last edited: 09/25/2023
*/

#include "seeding.h"
#include "util.h"
#include <bits/stdc++.h>

#ifndef _MINIMIZERSEEDING_H
#define _MINIMIZERSEEDING_H

struct MMentry{
    uint64_t hash;
    int pos;
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
    //produce seeds for both s and revComp(s)
    //stored as 0-s_idx.mmseed and 1-s_idx.mmseed respectively
    double getSeeds(std::string& s, const size_t s_idx,
		    const char* ouput_dir, const int dir_len);
};
#endif

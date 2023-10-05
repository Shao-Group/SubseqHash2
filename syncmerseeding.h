/*
  minimap2-style syncmer seeding.

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

#ifndef _SYNCMERSEEDING_H
#define _SYNCMERSEEDING_H

struct MMentry{
    uint64_t hash;
    size_t pos;
    kmer str;
};

class syncmerseeding: private seeding
{

private:

    uint64_t murmur64(kmer key);
    int s;
    std::vector<int> pos;
public:

    syncmerseeding(int k1, int s1): seeding(k1, k1)
	{
	    s = s1;
	}
    void add(int p)
	{
	    pos.push_back(p);
	}

    void get_syncmers(std::string s, std::vector<seed>& seeds);
    //produce seeds for both s and revcomp(s),
    //stored as 0-s_idx.syncmerseed and 1-s_idx.syncmerseed respectively
    double getSeeds(std::string& s, const size_t s_idx,
		    const char* output_dir, const int dir_len);

};
#endif

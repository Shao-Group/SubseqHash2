/*
  Utility functions for SubseqHash2.

  By: Ke, Xiang@PSU
  Last edited: 05/14/2023
*/

#ifndef _UTIL_H
#define _UTIL_H 1

#include <cstring>
#include <map>
#include <vector>
#include <unordered_map>

typedef __uint128_t kmer;

#define ALPHABETSIZE 4
const char ALPHABET[ALPHABETSIZE] = {'A', 'C', 'G', 'T'};
static inline int alphabetIndex(char c)
{
    return 3 & ((c>>2) ^ (c>>1));
}

/*
  Encode the string representation of a k-mer.
*/
kmer encode(const char* s, const int k);

/*
  Decode an k-mer into its string representation.
  If str is not null, it is used to store the resulting string; 
  otherwise a new char array is allocated.
*/
char* decode(const kmer enc, const int k, char* str);


struct seed
{
    int64_t hashval;
    kmer str, str_rc;
    size_t st, ed;
    uint64_t index = 0;
};


/*
  Compute the reverse compliment of a given kmer s.
*/
kmer revComp(kmer s, int k);

struct seedinfo
{
    uint32_t read;
    seed* s;
};

struct seedmatch
{
    seed *s1, *s2;
};

typedef std::unordered_map<int64_t, std::vector<seedinfo>> ssh_index;

ssh_index* index_build(std::vector<seed>& seeds);
void index_get(ssh_index* ht, std::vector<seed>& seeds, std::vector<seedmatch>& matches);
void saveSeeds(const char* filename, int k, const std::vector<seed>& seeds);
void saveSeedsPosition(const char* filename, const std::vector<seed>& seeds);
void saveSeedsStrPosition(const char* filename, const std::vector<seed>& seeds);
void loadSeeds(const char* filename, std::vector<seed> &seeds);
void loadSeedsStr(const char* filename, std::vector<seed> &seeds, int k);

#endif //util.h

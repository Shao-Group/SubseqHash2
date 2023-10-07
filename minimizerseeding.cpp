#include "minimizerseeding.h"

inline uint64_t minimizerseeding::hash64(kmer s, uint64_t mask)
{
    uint64_t key = (uint64_t) s;
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

void minimizerseeding::get_minimizers(std::string s, std::vector<seed>& seeds)
{
    size_t len = s.length();
    std::list<MMentry> que;

    kmer mask = (1ULL<<(k<<1)) - 1;
    char cur[k];
    s.copy(cur, k, 0);

    kmer now = encode(cur, k);
    uint64_t v = hash64(now, mask);

    que.push_back((MMentry){v, 0, now});

    for(size_t i = 1; i < w; i++)
    {
        now = (now<<2) & mask;
        now |= alphabetIndex(s[i+k-1]);
        v = hash64(now, mask);

        while(!que.empty() && v <= que.back().hash)
            que.pop_back();

        que.push_back((MMentry){v, i, now});
    }

    for(size_t i = w; i <= len - k; i++)
    {       
        while(!que.empty() && que.front().pos < i - w)
            que.pop_front();

        if(seeds.empty() || que.front().pos != seeds.back().st)
        {
            seed tmp;

            tmp.st = que.front().pos;
            tmp.ed = que.front().pos + k - 1;
            tmp.hashval = tmp.str = que.front().str;

            seeds.push_back(tmp);
        }

        now = (now<<2) & mask;
        now |= alphabetIndex(s[i+k-1]);
        v = hash64(now, mask);

        while(!que.empty() && v <= que.back().hash)
            que.pop_back();

        que.push_back((MMentry){v, i, now});
    }        

    while(!que.empty() && que.front().pos <= len - k - w)
        que.pop_front();

    if(seeds.empty() || que.front().pos != seeds.back().st)
    {
        seed tmp;

        tmp.st = que.front().pos;
        tmp.ed = que.front().pos + k - 1;
        tmp.hashval = tmp.str = que.front().str;

        seeds.push_back(tmp);
    }
}

double minimizerseeding::getSeeds(std::string& s, const size_t s_idx,
				  const char* output_dir, const int dir_len){
    std::vector<seed> seeds;
    get_minimizers(s, seeds);

    char output_filename[500];
    sprintf(output_filename, "%.*s/%d-%zu.mmseed",
	    dir_len, output_dir, 0, s_idx);
    saveSeeds(output_filename, k, seeds);

    double density = (double) seeds.size();

    seeds.clear();
    get_minimizers(revComp(s), seeds);
    sprintf(output_filename, "%.*s/%d-%zu.mmseed",
	    dir_len, output_dir, 1, s_idx);
    saveSeeds(output_filename, k, seeds);

    return (density + seeds.size())/(s.length()<<1);
}

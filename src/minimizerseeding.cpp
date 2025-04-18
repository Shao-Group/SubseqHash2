#include "minimizerseeding.hpp"

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
    int len = s.length();
    std::list<MMentry> que;
    std::vector<MMentry> kmerlist;

    kmer mask = (1ULL<<(k<<1)) - 1;

    kmer now = 0;
    int l = 0;
    uint64_t v;

    for(int i = 0; i < len; i++)
    {
        if(s[i] == 'N')
        {
            l = 0;
            now = 0;
            continue;
        }

        now = ((now<<2) | alphabetIndex(s[i])) & mask;
        
        if(++l >= k)
        {
            v = hash64(now, mask);
            kmerlist.push_back((MMentry){v, i + 1 - k, now});
        }
    }

    l = 0;
    int sz = kmerlist.size();

    for(int i = 0; i <= len - k; i++)
    {
        if(l < sz && kmerlist[l].pos == i)
        {
            while(!que.empty() && kmerlist[l].hash <= que.back().hash)
                que.pop_back();
            que.push_back(kmerlist[l++]);
        }

        if(i < w - 1)
            continue;

        while(!que.empty() && que.front().pos <= i - w)
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

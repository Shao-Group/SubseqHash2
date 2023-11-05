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
    std::vector<MMentry> kmerlist;

    kmer mask = (1ULL<<(k<<1)) - 1;

    kmer now = 0;
    size_t l = 0;
    uint64_t v;

    for(size_t i = 0; i < len; i++)
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
            kmerlist.push_back((MMentry){v, i - k + 1, now});
        }
    }

    l = 0;
    int sz = kmerlist.size();

    for(size_t i = 0; i <= len - k; i++)
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

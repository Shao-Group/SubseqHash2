#include "syncmerseeding.h"


uint64_t syncmerseeding::murmur64(kmer key)
{    
    uint64_t h = (uint64_t) key;
    h+=1;
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdL;
    h ^= (h >> 33);
    h *= 0xc4ceb9fe1a85ec53L;
    h ^= (h >> 33);
    return h;
}

void syncmerseeding::get_syncmers(std::string str, std::vector<seed>& seeds)
{
    size_t len = str.length();
    std::list<MMentry> que;

    kmer mask = (1ULL<<(s<<1)) - 1;
    kmer mask_k = (1ULL<<(k<<1)) - 1;

    kmer now_k = 0;
    kmer now = 0;
    size_t l = 0;
    uint64_t v;

    for(size_t i = 0; i < len; i++)
    {
        if(str[i] == 'N')
        {
            l = 0;
            now = 0;
            continue;
        }

        now = ((now<<2) | alphabetIndex(str[i])) & mask;
        now_k = ((now_k<<2) | alphabetIndex(str[i])) & mask_k;
        
        if(++l >= s)
        {
            v = murmur64(now);
            while(!que.empty() && v <= que.back().hash)
                que.pop_back();
            que.push_back((MMentry){v, i-s+1, now});
        }

        if(l >= k)
        {        
            while(!que.empty() && que.front().pos <= i - k)
                que.pop_front();

            for(int j:pos)
                if(que.front().pos - (i-k+1) == j)
                {
                    seed tmp;

                    tmp.st = i - k + 1;
                    tmp.ed = i;
                    tmp.hashval = now_k;

                    seeds.push_back(tmp);
                }

        }
    }
}

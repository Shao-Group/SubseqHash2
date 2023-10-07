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
    char cur[k];

    str.copy(cur, s, 0);
    kmer now = encode(cur, s);
    uint64_t v = murmur64(now);

    str.copy(cur, k, 0);
    kmer strval = encode(cur, k);

    que.push_back((MMentry){v, 0, now});

    int num = k - s + 1;

    for(size_t i = 1; i < num; i++)
    {
        now = (now<<2) & mask;
        now |= alphabetIndex(str[i+s-1]);
        v = murmur64(now);

        while(!que.empty() && v <= que.back().hash)
            que.pop_back();

        que.push_back((MMentry){v, i, now});
    }

    for(size_t i = k-s+1; i <= len - s; i++)
    {       
        while(!que.empty() && que.front().pos < i - num)
            que.pop_front();

        for(int j:pos)
            if(que.front().pos - (i-num) == j)
            {
                seed tmp;

                tmp.st = i - num;
                tmp.ed = tmp.st + k - 1;
                tmp.hashval = strval;
		tmp.str = strval;

                seeds.push_back(tmp);
		break;
            }

        now = (now<<2) & mask;
        now |= alphabetIndex(str[i+s-1]);
        v = murmur64(now);

        strval = (strval<<2) & mask_k;
        strval |= alphabetIndex(str[i+s-1]);

        // char* tmp1 = (char*)malloc(sizeof(char) *k), *tmp2= (char*)malloc(sizeof(char) *k);
        // decode(now, s, tmp1);
        // decode(strval, k, tmp2);
        // for(int z = 0; z<s; z++)
        //     printf("%c", tmp1[z]);
        // printf(" ");
        // for(int z = 0; z<k; z++)
        //     printf("%c", tmp2[z]);
        // printf("\n");
        // printf("%s\n", str.substr(i - num + 1, k).c_str());

        while(!que.empty() && v <= que.back().hash)
            que.pop_back();

        que.push_back((MMentry){v, i, now});
    }        

    while(!que.empty() && que.front().pos < len - k)
        que.pop_front();

    for(int j:pos)
        if(que.front().pos - (len - k) == j)
        {
            seed tmp;

            tmp.st = len - k;
            tmp.ed = len - 1;
            tmp.hashval = strval;
	    tmp.str = strval;

            seeds.push_back(tmp);
	    break;
        }
}

double syncmerseeding::getSeeds(std::string& s, const size_t s_idx,
				const char* output_dir, const int dir_len){
    std::vector<seed> seeds;
    get_syncmers(s, seeds);
    
    char output_filename[500];
    sprintf(output_filename, "%.*s/%d-%zu.syncmerseed",
	    dir_len, output_dir, 0, s_idx);
    saveSeeds(output_filename, k, seeds);

    double density = (double) seeds.size();

    seeds.clear();
    get_syncmers(revComp(s), seeds);
    sprintf(output_filename, "%.*s/%d-%zu.syncmerseed",
	    dir_len, output_dir, 1, s_idx);
    saveSeeds(output_filename, k, seeds);

    return (density + seeds.size())/(s.length()<<1);
}

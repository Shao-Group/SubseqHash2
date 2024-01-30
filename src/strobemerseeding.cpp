#include "strobemerseeding.h"

void strobemerseeding::get_strobemers(std::string str, std::vector<std::vector<seed>>& seeds)
{
    strobes_vector randstrobes;

    for(int i = 0; i < seednum; i++)
    {
        if(w == 2)
        {
            randstrobes = seq_to_randstrobes2(w, k, w_min, w_max, str, 0, randseeds[i]);

            for(auto &t : randstrobes) // iterate over the strobemer tuples
            {
                seed tmp;
                tmp.st = std::get<2>(t);
                tmp.ed = std::get<3>(t) + k - 1;
                tmp.index = std::get<3>(t) - tmp.st;

                std::string strobe = str.substr(tmp.st, k) + str.substr(std::get<3>(t), k);
                char cur[k*2];

                strobe.copy(cur, k*2, 0);
                tmp.str = encode(cur, k*2);
                tmp.hashval = (uint64_t)(tmp.str>>64) ^ ((uint64_t)tmp.str << 1);
                seeds[i].push_back(tmp);
            }
        }
        else if(w == 3)
        {
            randstrobes = seq_to_randstrobes3(w, k, w_min, w_max, str, 0, randseeds[i]);

            for(auto &t : randstrobes) // iterate over the strobemer tuples
            {
                seed tmp;
                tmp.st = std::get<2>(t);
                tmp.ed = std::get<4>(t) + k - 1;
                tmp.index = std::get<3>(t) - tmp.st + ((std::get<4>(t) - tmp.st)<<10);

                std::string strobe = str.substr(tmp.st, k) + str.substr(std::get<3>(t), k) + str.substr(std::get<4>(t), k);
                char cur[k*3];

                strobe.copy(cur, k*3, 0);
                tmp.str = encode(cur, k*3);
                tmp.hashval = (uint64_t)(tmp.str>>64) ^ ((uint64_t)tmp.str << 1);
                
                seeds[i].push_back(tmp);
            }       
        }
    }
}

double strobemerseeding::getSeeds(std::string& s, const size_t s_idx,
				  const char* output_dir, const int dir_len){
    std::vector<std::vector<seed>> seeds(seednum);
    get_strobemers(s, seeds);
    int seed_len = k*w;

    double density = 0.0;
    char output_filename[500];
    for(int i=0; i<seednum; ++i){
	sprintf(output_filename, "%.*s/%d-%zu.strobemerseed",
		dir_len, output_dir, i<<1, s_idx);
	saveSeeds(output_filename, seed_len, seeds[i]);
	density += seeds[i].size();
	seeds[i].clear();
    }

    get_strobemers(revComp(s), seeds);
    for(int i=0; i<seednum; ++i){
	sprintf(output_filename, "%.*s/%d-%zu.strobemerseed",
		dir_len, output_dir, (i<<1)+1, s_idx);
	saveSeeds(output_filename, seed_len, seeds[i]);
	density += seeds[i].size();
    }

    return density/(s.length()*seednum*2);
}

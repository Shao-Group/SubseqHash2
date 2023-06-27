#include "strobemerseeding.h"

void strobemerseeding::get_strobemers(std::string str, std::vector<seed>& seeds)
{
    strobes_vector randstrobes;

    if(w == 2)
    {
        randstrobes = seq_to_randstrobes2(w, k, w_min, w_max, str, 0);

        for(auto &t : randstrobes) // iterate over the strobemer tuples
        {
            seed tmp;
            tmp.st = std::get<2>(t);
            tmp.ed = std::get<3>(t) + k - 1;
            tmp.index = std::get<3>(t) - tmp.st;

            std::string strobe = str.substr(tmp.st, k) + str.substr(std::get<3>(t), k);
            char cur[k*2];

            strobe.copy(cur, k*2, 0);
            tmp.hashval = encode(cur, k*2);
            
            seeds.push_back(tmp);
        }
    }
    else if(w == 3)
    {
        randstrobes = seq_to_randstrobes3(w, k, w_min, w_max, str, 0);

        for(auto &t : randstrobes) // iterate over the strobemer tuples
        {
            seed tmp;
            tmp.st = std::get<2>(t);
            tmp.ed = std::get<4>(t) + k - 1;
            tmp.index = std::get<3>(t) - tmp.st + ((std::get<4>(t) - tmp.st)<<10);

            std::string strobe = str.substr(tmp.st, k) + str.substr(std::get<3>(t), k) + str.substr(std::get<4>(t), k);
            char cur[k*3];

            strobe.copy(cur, k*3, 0);
            tmp.hashval = encode(cur, k*3);
            
            seeds.push_back(tmp);
        }       
    }
}

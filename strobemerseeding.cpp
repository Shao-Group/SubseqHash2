#include "strobemerseeding.h"

void strobemerseeding::get_strobemers(std::string str, std::vector<seed>& seeds)
{
    strobes_vector randstrobes;
    randstrobes = seq_to_randstrobes2(w, k, w_min, w_max, str, 0);

    for(auto &t : randstrobes) // iterate over the strobemer tuples
    {
        seed tmp;
        tmp.st = std::get<2>(t);
        tmp.ed = std::get<3>(t) + k - 1;
        tmp.index = std::get<3>(t);

        std::string strobe = str.substr(tmp.st, k) + str.substr(tmp.index, k);
        char cur[k*2];

        strobe.copy(cur, k*2, 0);
        tmp.hashval = encode(cur, k*2);
        
        seeds.push_back(tmp);
    }
}

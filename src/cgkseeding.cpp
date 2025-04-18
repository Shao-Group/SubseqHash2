#include "cgkseeding.h"

void cgkseeding::get_cgk(std::string s, std::vector<seed>& seeds)
{
    int len = s.length();
    for(int i = 0; i + k <= len; i++)
    {
        seed tmp;
        tmp.hashval = tmp.str = cgk_embedding(s, i);
        tmp.st = i;
        tmp.ed = i + k - 1;
        seeds.push_back(tmp);
    }
}

std::string cgkseeding::cgk (const std::string& input) 
{
    std::string output = "";
    for (int i = 0, j = 0; j < 3 * k; j++) 
    {
        char s = i < k ? input.at(i) : 'N';

        i = i + *(p[dict[s]] + j);
        output = output + s;
    }
    return output;
};

inline int64_t cgkseeding::cgk_embedding (const std::string& input, int st) 
{
    kmer output = 0;
    for (int i = st, j = 0; j < 3 * k; j++) 
    {
        char s = i < st + k ? input.at(i) : 'N';

        i = i + *(p[dict[s]] + j);
        output = output * 5 + dict[s];
    }
    return kmerHash{}(output);
};

inline int64_t cgkseeding::smoothq (std::string& input, int st)
{
    int64_t output = 0;
    int i = st, partdigit = 0;

    for (int j = 0; partdigit < hash_lsh.size(); ++j)
    {
        char s = i < st + k ? input[i] : 'N';
        //only record the bits used in LSH
        while (j == hash_lsh[partdigit])
        {
            output = output * 5 + dict[s];
            partdigit += 1;
        }
        i = i + *(p[dict[s]] + j);
    }
    return output;
}

double cgkseeding::getSeeds(std::string& s, const size_t s_idx,
			    const char* output_dir, const int dir_len){
    std::vector<seed> seeds;
    get_cgk(s, seeds);

    char output_filename[500];
    sprintf(output_filename, "%.*s/%d-%zu.cgkseed",
	    dir_len, output_dir, 0, s_idx);
    saveSeeds(output_filename, 3*k, seeds);

    double density = (double) seeds.size();

    seeds.clear();
    get_cgk(revComp(s), seeds);
    sprintf(output_filename, "%.*s/%d-%zu.cgkseed",
	    dir_len, output_dir, 1, s_idx);
    saveSeeds(output_filename, 3*k, seeds);

    return (density + seeds.size())/(s.length()<<1);
}

/*
  CGK embedding seeding.

  By: Xiang@PSU
  Last edited: 04/18/2025
*/

#include "seeding.hpp"
#include "util.hpp"
#include <bits/stdc++.h>

#ifndef _CGKSEEDING_H
#define _CGKSEEDING_H


class cgkseeding: private seeding
{

private:
    std::vector<int> generate_random_sampling_positions(int size_cgk, int size_smoothq)
    {
        std::random_device rd;
        std::mt19937 generator(rd());
        std::vector<int> hash_lsh;
        hash_lsh.reserve(size_cgk);
        for (int j = 0; j < size_cgk; j++)
            hash_lsh.push_back(j);
        shuffle(hash_lsh.begin(), hash_lsh.end(), generator);
        hash_lsh.resize(size_smoothq);

        sort(hash_lsh.begin(), hash_lsh.end());
        return hash_lsh;
    };

    int** generate_random_binary(int size_cgk) 
    {
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<int> distribution(0, 1);	

        const int num_char = 5;
        int **p = new int *[num_char];
        for (int t = 0; t < num_char; t++) 
        {
            p[t] = new int [size_cgk];
            for (int d = 0; d < size_cgk; d++)
                p[t][d] = distribution(generator);
        }
        return p;
    };

public:
    int ** p;
    std::vector<int> hash_lsh;

    std::unordered_map<char, int> dict = { {'A', 0 }, { 'C', 1 }, { 'G', 2 }, { 'T', 3 }, { 'N', 4 } };

    cgkseeding(int k1, int r=0): seeding(k1, k1) {
        p = generate_random_binary(3 * k1);
    };
    
    inline int64_t cgk_embedding (const std::string& input, int st); 
    void get_cgk(std::string s, std::vector<seed>& seeds);
    inline int64_t smoothq (std::string& input, int st);
    std::string cgk (const std::string& input);
    //produce seeds for both s and revComp(s)
    //stored as 0-s_idx.cgkseed and 1-s_idx.cgkseed respectively
    double getSeeds(std::string& s, const size_t s_idx,
		    const char* ouput_dir, const int dir_len);

};
#endif

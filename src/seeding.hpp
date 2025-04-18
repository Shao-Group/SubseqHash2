/*
  Base class for seeding method.

  By: Ke, Xiang@PSU
  Last edited: 05/16/2023
*/

#ifndef _SEEDING_H
#define _SEEDING_H

#include <stdexcept>

class seeding{
public:
    int k;
    int n;

    seeding(int n1, int k1); //kmer size, window size
    
    template<typename...Resources>
    double getSeeds(std::string& s, const char* output_dir, const int dir_len,
		    Resources&&...args); //throws unimplemented error
};


#endif

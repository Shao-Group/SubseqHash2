#include "seeding.hpp"

seeding::seeding(int n1, int k1)
{
	n = n1;
	k = k1;
}

template<typename...Resources>
double seeding::getSeeds(std::string& s,
			 const char* output_dir, const int dir_len,
			 Resources&&... args){
    throw std::logic_error("getSeeds is not implemented");
}

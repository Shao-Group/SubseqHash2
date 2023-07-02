/*
  strobemer seeding.

  By: Xiang@PSU
  Last edited: 06/21/2023
*/

#include "seeding.h"
#include "strobemer/index.hpp"
#include "util.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <random>

#ifndef _STROBEMER_H
#define _STROBEMER_H

typedef std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> strobematch;
typedef std::vector<strobematch> strobes_vector;

class strobemerseeding: private seeding
{

private:

	int w, w_min, w_max;
	int seednum;
	std::vector<uint64_t> randseeds;
public:

	strobemerseeding(int k1, int w1, int w_min1, int w_max1, int repeat): seeding(2 * (w_max1-w_min) + 3 * k1, k1)
	{
		w = w1; w_min = w_min1; w_max = w_max1; seednum = repeat;

		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<int64_t> distribution(0, (uint64_t)-1);

		for(int i = 0; i < seednum; i++)
			randseeds.push_back(distribution(generator));
	}

	void get_strobemers(std::string s, std::vector<std::vector<seed>>& seeds);

};
#endif
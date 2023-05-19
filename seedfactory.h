/*
  seeding threads.

  By: Ke, Xiang@PSU
  Last edited: 05/18/2023
*/


#include "subseqhash2seeding.h"

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <functional>

#ifndef _SEEDFACTORY_H
#define _SEEDFACTORY_H

#define NUMTHREADS 15

struct Read
{
    std::string seq;
    size_t idx;

    Read(std::string&& s, size_t i): seq(move(s)), idx(i) {};
    Read(Read&& o): seq(move(o.seq)), idx(std::exchange(o.idx, 0)) {};
};

class seedFactory
{
	//static_assert(std::is_base_of<seeding, seedType>::value, "SeedType must derive from base seeding class");

    const char* output_dir;
    const int dir_len;

    std::queue<Read> jobs;
    std::vector<std::thread> minions;
    bool done;
    std::mutex door, out_door;
    std::condition_variable trumpet;
    double total_density;
    int num_jobs;

    subseqhash2seeding myseeding;
    int n, k, d, dim1;

    void getSubseqSeeds(const Read &r, DPCell* dp, DPCell* revdp, int* h, int* revh);
    void atWork(int x);

public:
    seedFactory(const char* output_dir, int dir_len, int n1, int k1, int d1, const char* tablefile)	
    :output_dir(output_dir), dir_len(dir_len), done(false), total_density(0.0), num_jobs(0), n(n1), k(k1), d(d1), myseeding(n1, k1, d1)
	{
	    dim1 = (n1+1) * (k1+1) * d1;
	    myseeding.init(tablefile);
	    
		minions.reserve(NUMTHREADS);
		for(int i=0; i<NUMTHREADS; ++i)
		    minions.emplace_back(std::bind(&seedFactory::atWork, this, i));
		
	}

    ~seedFactory()
    {
		std::unique_lock<std::mutex> lock(door);
		done = true;
		lock.unlock();
		trumpet.notify_all();

		for(auto& x : minions)
		    x.join();

		printf("Density: %.8f\n", total_density/num_jobs);
	}
    void addJob(std::string&& r, size_t idx);
};
#endif
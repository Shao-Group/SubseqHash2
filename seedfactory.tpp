#include "seedfactory.h"
/*    
void seedFactory::getSubseqSeeds(const Read &r, DPCell* dp, DPCell* revdp, int* h, int* revh)
{
	int num_valid = myseeding.getNumPerWindow();

	std::vector<std::vector<seed>> seeds(num_valid, std::vector<seed>(0));
	myseeding.getSubseq2Seeds(r.seq, dp, revdp, h, revh, seeds);

	double density = 0.0;
	char output_filename[200];

	for(int i = 0; i < num_valid; i++)
	{
	    sprintf(output_filename, "%.*s/%d-%zu.subseqseed2",
		    dir_len, output_dir, i, r.idx);
	    saveSeeds(output_filename, k, seeds[i]);

	    density += seeds[i].size();
	}
	
	const std::lock_guard<std::mutex> lock(out_door);
	total_density += density/(r.seq.length()*num_valid);
}
*/

template <typename SeedType>
template <typename... Resources>
void seedFactory<SeedType>::atWork(int id, Resources&&... args)
{
    /*
    DPCell* dp = (DPCell*) malloc(sizeof *dp * myseeding.getChunkSize() * dim1);
    DPCell* revdp = (DPCell*) malloc(sizeof *revdp * myseeding.getChunkSize() * dim1);
    int* h = (int*) malloc(sizeof *h * dim1);
    int* revh = (int*) malloc(sizeof *revh * dim1);
    */
    
    std::unique_lock<std::mutex> lock(door, std::defer_lock);
    while(true)
    {
	lock.lock();
	while(!done && jobs.empty())
	    trumpet.wait(lock);
	    
	if(!jobs.empty())
	{
	    Read r = std::move(jobs.front());

	    jobs.pop();

	    lock.unlock();
	    //getSubseqSeeds(r, dp, revdp, h, revh);
	    double density = myseeding.getSeeds(r.seq, r.idx,
	    	   	     			output_dir, dir_len,
						std::forward<Resources>(args)...);
	    {
		std::lock_guard<std::mutex> out_lock(out_door);
		total_density += density;
	    }
	}
	else
	{
	    /*
	    free(dp);
	    free(revdp);
	    free(h);
	    free(revh);
	    */
	    return;
	}
    }
}

template <typename SeedType>
void seedFactory<SeedType>::addJob(std::string&& r, size_t idx)
{
    std::unique_lock<std::mutex> lock(door);
    ++ num_jobs;
    jobs.emplace(move(r), idx);
    trumpet.notify_one();
}

template <typename SeedType>
template <typename... Resources>
void seedFactory<SeedType>::addMinions(int id, Resources&&... args){
    minions.emplace_back(std::bind(&seedFactory<SeedType>::atWork<Resources&&...>,
				   this, id, std::forward<Resources>(args)...));
}

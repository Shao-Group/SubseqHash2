#include "../seedfactory.h"
#include "../subseq2minstrobeseeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char * argv[]){    
    if(argc != 9){
	fprintf(stderr, "usage: genSubseq2minstrobeSeeds.out readFile n k d randTableFile w pre_k mm_w\n");
	return 1;
    }

    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    int d = atoi(argv[4]);
    int w = atoi(argv[6]);
    int pre_k = atoi(argv[7]);
    int mm_w = atoi(argv[8]);

    int dim1 = (n+1) * (k+1) * d;

    struct stat test_table;
    if(stat(argv[5], &test_table) != 0){//table file exists
	fprintf(stderr, "cannot read table file %s\n", argv[5]);
	return 1;
    }

    char output_dir[500];
    int dir_len = sprintf(output_dir, "%s-subseq2minstrobeSeeds-n%d-k%d-d%d-w%d-pre_k%d-mm_w%d/",
			  argv[1], n, k, d, w, pre_k, mm_w);
    
    ifstream fin(argv[1]);

    string read;
    size_t read_idx = 0;

    //resources for minions
    DPCell* dp[NUMTHREADS];
    DPCell* revdp[NUMTHREADS];
    int* h[NUMTHREADS];
    int* revh[NUMTHREADS];
    
    {
	seedFactory<subseq2minstrobeseeding> factory(output_dir, dir_len,
						  n, k, d, k, w, pre_k, mm_w,
						  argv[5]);
	const subseq2minstrobeseeding& myseeding = factory.getMySeeding();
	int chunk_size = myseeding.getChunkSize();

	mkdir(output_dir, 0744);

	for(int i=0; i<NUMTHREADS; ++i){
	    dp[i] = (DPCell*) malloc(sizeof *dp[i] * chunk_size * dim1);
	    revdp[i] = (DPCell*) malloc(sizeof *revdp[i] * chunk_size * dim1);
	    h[i] = (int*) malloc(sizeof *h[i] * dim1);
	    revh[i] = (int*) malloc(sizeof *revh[i] * dim1);
	    factory.addMinions(i, dp[i], revdp[i], h[i], revh[i]);
	}

	while(fin.get()=='>'){
	    //skip the header
	    fin.ignore(numeric_limits<streamsize>::max(), '\n');
	    getline(fin, read);
	    ++ read_idx;
	    factory.addJob(move(read), read_idx);
	    //skip the mapping info
	    fin.ignore(numeric_limits<streamsize>::max(), '\n');
	}
    }

    for(int i=0; i<NUMTHREADS; ++i){
	free(dp[i]);
	free(revdp[i]);
	free(h[i]);
	free(revh[i]);
    }

    return 0;
}

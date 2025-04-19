#include "../src/seedfactory.hpp"
#include "../src/minimizerseeding.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char * argv[]){    
    if(argc != 4){
	fprintf(stderr, "usage: genMinimizerSeeds.out readFile w k\n");
	return 1;
    }

    int w = atoi(argv[2]);
    int k = atoi(argv[3]);

    char output_dir[500];
    int dir_len = sprintf(output_dir, "%s-minimizerSeeds-w%d-k%d/",
			  argv[1], w, k);
    
    ifstream fin(argv[1]);

    string read;
    size_t read_idx = 0;
    
    {
	seedFactory<minimizerseeding> factory(output_dir, dir_len, k, w);

	mkdir(output_dir, 0744);

	for(int i=0; i<NUMTHREADS; ++i){
	    factory.addMinions(i);
	}

	while(fin.get()=='>'){
	    //skip the header
	    fin.ignore(numeric_limits<streamsize>::max(), '\n');
	    getline(fin, read);
	    ++ read_idx;
	    factory.addJob(move(read), read_idx);
	    //skip the mapping info
	    //fin.ignore(numeric_limits<streamsize>::max(), '\n');
	}
    }

    return 0;
}

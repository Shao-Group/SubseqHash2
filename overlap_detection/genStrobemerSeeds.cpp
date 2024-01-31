#include "../src/seedfactory.h"
#include "../src/strobemerseeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char * argv[]){    
    if(argc != 7){
	fprintf(stderr, "usage: genStrobemerSeeds.out readFile k w w_min w_max seednum\n");
	return 1;
    }

    int k = atoi(argv[2]);
    int w = atoi(argv[3]);
    int w_min = atoi(argv[4]);
    int w_max = atoi(argv[5]);
    int seednum = atoi(argv[6]);

    char output_dir[500];
    int dir_len = sprintf(output_dir, "%s-strobemerSeeds-k%d-w%d-w_min%d-w_max%d/",
			  argv[1], k, w, w_min, w_max);
    
    ifstream fin(argv[1]);

    string read;
    size_t read_idx = 0;

    
    {
	seedFactory<strobemerseeding> factory(output_dir, dir_len,
						k, w, w_min, w_max, seednum);

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

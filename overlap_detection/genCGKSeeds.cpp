#include "../src/seedfactory.h"
#include "../src/cgkseeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char * argv[]){    
    if(argc != 4){
	fprintf(stderr, "usage: genCGKSeeds.out readFile k r\n");
	return 1;
    }

    int k = atoi(argv[2]);
    int r = atoi(argv[3]);

    char output_dir[500];
    int dir_len = sprintf(output_dir, "%s-CGKSeeds-k%d-r$d/",
			  argv[1], k, r);
    
    ifstream fin(argv[1]);

    string read;
    size_t read_idx = 0;
    
    {
	seedFactory<cgkseeding> factory(output_dir, dir_len, k);

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

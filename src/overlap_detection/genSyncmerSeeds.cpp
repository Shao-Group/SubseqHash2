#include "../seedfactory.h"
#include "../syncmerseeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char * argv[]){    
    if(argc != 4){
	fprintf(stderr, "usage: genSyncmerSeeds.out readFile k s\n");
	return 1;
    }

    int k = atoi(argv[2]);
    int s = atoi(argv[3]);

    char output_dir[500];
    int dir_len = sprintf(output_dir, "%s-syncmerSeeds-k%d-s%d/",
			  argv[1], k, s);
    
    ifstream fin(argv[1]);

    string read;
    size_t read_idx = 0;
    
    {
	seedFactory<syncmerseeding> factory(output_dir, dir_len, k, s);
	syncmerseeding& myseeding = factory.getMySeeding();
	//closed syncmer
	myseeding.add(0);
	myseeding.add(k-s);

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
	    fin.ignore(numeric_limits<streamsize>::max(), '\n');
	}
    }

    return 0;
}

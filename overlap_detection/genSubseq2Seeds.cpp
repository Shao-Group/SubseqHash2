#include "../seedfactory.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char * argv[]){    
    if(argc != 6){
	fprintf(stderr, "usage: genSubseq2Seeds.out readFile n k d randTableFile\n");
	return 1;
    }

    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    int d = atoi(argv[4]);

    int dim1 = (n+1) * (k+1) * d;

    struct stat test_table;
    if(stat(argv[5], &test_table) != 0){//table file exists
	fprintf(stderr, "cannot read table file %s\n", argv[5]);
	return 1;
    }

    char output_dir[500];
    int dir_len = sprintf(output_dir, "%s-locSeeds-n%d-k%d-d%d/",
			  argv[1], n, k, d);
    
    ifstream fin(argv[1]);

    string read;
    size_t read_idx = 0;
    
    {
	seedFactory factory(output_dir, dir_len,
			    n, k, d, n, argv[5]);

	mkdir(output_dir, 0744);

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

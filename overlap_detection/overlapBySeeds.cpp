/*
  Given a set of seed files (readable by loadSeedsUnordered), output pairs of reads
  with the number of unique seeds they share.
  
  By: Ke@PSU
  Last edited: 09/19/2023
*/

#include "../src/util.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

/*
  An upper diagonal matrix without the main diagonal,
  valid indices for Table::access(i, j) are 1<=i<j<=n
*/
class Table{
    size_t n;
    unsigned int* arr;

public:
    Table(size_t n);
    ~Table();
    unsigned int& access(size_t i, size_t j);
    void saveNoneZeroEntries(const char* filename);
};

int main(int argc, const char * argv[])    
{   
    if(argc != 6){
	fprintf(stderr, "usage: overlapBySeeds.out seedsDir k repeatId numReads fileExt\nrepeatId should be between 0 and floor((k-1)/2), seeds from both batches repeatId*2 and repeatId*2+1 are used (self symmetric if repeadId*2+1==k)\n");
	return 1;
    }

    int k = atoi(argv[2]);
    int r = atoi(argv[3])<<1;
    int is_subseq2 = strcmp(argv[5], "subseqseed2");
    if(is_subseq2 == 0 && (r >= k || r < 0)){
	fprintf(stderr, "r=%d is not within valid range [0, %d]\n", r>>1, (k-1)>>1);
	return 1;
    }
    int n = atoi(argv[4]);

    int r2 = r + 1;

    char filename[500];
    int i = strlen(argv[1]);
    memcpy(filename, argv[1], i);
    if(filename[i-1] != '/'){
	filename[i] = '/';
	++i;
    }
    
    unordered_map<kmer, vector<size_t>, kmerHash> all_seeds;
    int j;
    
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+i, "%d-%d.%s", r, j, argv[5]);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r, j, argv[5]);
	    break;
	}
	loadSeedsUnordered(filename, j, all_seeds);
	
	if(is_subseq2 == 0 && r2 >= k){
	    continue;
	}else{
            sprintf(filename+i, "%d-%d.%s", r2, j, argv[5]);
            if(stat(filename, &test_file) != 0){//seed file does not exist
                fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r2, j, argv[5]);
                break;
            }
            loadSeedsUnordered(filename, j, all_seeds);
        }
    }

    sprintf(filename+i, "overlap-n%d-r%d.all-pair", n, r>>1);

    Table share_ct(n);

    int a, b, c;
    
    for(auto& seed : all_seeds){
	c = seed.second.size();
	for(i=0; i<c; ++i){
	    a = seed.second[i];
	    for(j=i+1; j<c; ++j){
		b = seed.second[j];
		++ share_ct.access(a, b);
		/*
		if(share_ct.access(a, b) == threshold){
		    fprintf(fout, "%d %d\n", a, b);
		}
		*/
	    }
	}
    }

    share_ct.saveNoneZeroEntries(filename);
    
    return 0;
}

Table::Table(size_t n): n(n){
    size_t size = (n*(n-1))>>1;
    arr = new unsigned int[size];
    memset(arr, 0, sizeof *arr * size);
}

Table::~Table(){
    delete[] arr;
}

unsigned int& Table::access(size_t i, size_t j){
    return arr[((((n<<1)-i)*(i-1))>>1)+j-i-1];
}

void Table::saveNoneZeroEntries(const char* filename){
    FILE* fout = fopen(filename, "w");

    size_t i, j, k, size=(n*(n-1))>>1;
    for(k=0, i=1, j=2; k<size; ++k){
	if(arr[k] > 0){
	    fprintf(fout, "%zu %zu %u\n", i, j, arr[k]);
	}
	j += 1;
	if(j > n){
	    i += 1;
	    j = i + 1;
	}
    }
    
    fclose(fout);
}

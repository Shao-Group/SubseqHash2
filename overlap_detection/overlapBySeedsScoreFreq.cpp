/*
  Given a set of seed files (readable by loadSeedsUnordered), output pairs of reads
  with the number of unique seeds they share;
  also output the min score seed match (and its frequency) between the pair and the max score seed match for investigating a good score/frequency threshold.
  
  By: Ke@PSU
  Last edited: 11/14/2023
*/

#include "../util.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

struct AnchorInfo{
    unsigned int num_uniq_seeds;
    int64_t min_score, max_score;
    size_t min_freq, max_freq; //frequencies of the min/max_score seed matches
};

/*
  An upper diagonal matrix without the main diagonal,
  valid indices for Table::access(i, j) are 1<=i<j<=n
*/
class Table{
    size_t n;
    AnchorInfo* arr;

public:
    Table(size_t n);
    ~Table();
    AnchorInfo& access(size_t i, size_t j);
    void saveNoneZeroEntries(const char* filename);
};

struct SeedInfo{
    int64_t score;
    size_t frequency;
    vector<size_t> read_ids;

    SeedInfo(int64_t sc):score(sc),frequency(0),read_ids(){}
};

void loadSeedsWithScoreUnordered(const char* filename, const size_t read_id,
				 unordered_map<kmer, SeedInfo, kmerHash> &all_seeds){
    FILE* fin = fopen(filename, "rb");
    size_t ret;
    int k;
    ret = fread(&k, sizeof(k), 1, fin);
    kmer s;
    int64_t score;
    while(true){
	ret = fread(&s, sizeof(s), 1, fin);
	if(ret != 1) break;
	ret = fread(&score, sizeof(score), 1, fin);
	//skip psi, st, index
	fseek(fin, sizeof(short)+(sizeof(uint64_t)<<1), SEEK_CUR);

	auto result = all_seeds.emplace(piecewise_construct,
					forward_as_tuple(s),
					forward_as_tuple(score));
	result.first->second.frequency += 1;
	
	if(result.first->second.frequency == 1){//first seen
	    result.first->second.read_ids.push_back(read_id);
	}else if(result.second == false && result.first->second.read_ids.back()<read_id){//each seed is only added once for a read
	    result.first->second.read_ids.push_back(read_id);
	}
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}
				 

int main(int argc, const char * argv[])    
{   
    if(argc != 6){
	fprintf(stderr, "usage: overlapBySeedsScoreFreq.out seedsDir k repeatId numReads fileExt\nrepeatId should be between 0 and floor((k-1)/2), seeds from both batches repeatId*2 and repeatId*2+1 are used (self symmetric if repeadId*2+1==k)\n");
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
    
    unordered_map<kmer, SeedInfo, kmerHash> all_seeds;
    int j;
    
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+i, "%d-%d.%s", r, j, argv[5]);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r, j, argv[5]);
	    break;
	}
	loadSeedsWithScoreUnordered(filename, j, all_seeds);
	
	if(is_subseq2 == 0 && r2 >= k){
	    continue;
	}else{
            sprintf(filename+i, "%d-%d.%s", r2, j, argv[5]);
            if(stat(filename, &test_file) != 0){//seed file does not exist
                fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r2, j, argv[5]);
                break;
            }
            loadSeedsWithScoreUnordered(filename, j, all_seeds);
        }
    }

    sprintf(filename+i, "overlapScoreFreq-n%d-r%d.all-pair", n, r>>1);

    Table share_ct(n);

    int a, b, c;
    
    for(auto& seed : all_seeds){
	c = seed.second.read_ids.size();
	for(i=0; i<c; ++i){
	    a = seed.second.read_ids[i];
	    for(j=i+1; j<c; ++j){
		b = seed.second.read_ids[j];
		AnchorInfo& entry = share_ct.access(a, b);
		++ entry.num_uniq_seeds;
		if(entry.num_uniq_seeds == 1){//first seed match
		    entry.max_score = entry.min_score = seed.second.score;
		    entry.max_freq = entry.min_freq = seed.second.frequency;
		}else if(entry.min_score > seed.second.score){
		    entry.min_score = seed.second.score;
		    entry.min_freq = seed.second.frequency;
		}else if(entry.max_score < seed.second.score){
		    entry.max_score = seed.second.score;
		    entry.max_freq = seed.second.frequency;
		}
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
    arr = new AnchorInfo[size];
    memset(arr, 0, sizeof *arr * size);
}

Table::~Table(){
    delete[] arr;
}

AnchorInfo& Table::access(size_t i, size_t j){
    return arr[((((n<<1)-i)*(i-1))>>1)+j-i-1];
}

void Table::saveNoneZeroEntries(const char* filename){
    FILE* fout = fopen(filename, "w");

    size_t i, j, k, size=(n*(n-1))>>1;
    for(k=0, i=1, j=2; k<size; ++k){
	if(arr[k].num_uniq_seeds > 0){
	    fprintf(fout, "%zu %zu %u %ld %zu %ld %zu\n", i, j,
		    arr[k].num_uniq_seeds,
		    arr[k].max_score, arr[k].max_freq,
		    arr[k].min_score, arr[k].min_freq);
	}
	j += 1;
	if(j > n){
	    i += 1;
	    j = i + 1;
	}
    }
    
    fclose(fout);
}

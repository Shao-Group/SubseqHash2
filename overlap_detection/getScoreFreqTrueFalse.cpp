/*
  Given a set of seed files (readable by loadSeedsWithScore defined below),
  output frequency for different scores.
  
  By: Ke@PSU
  Last edited: 11/09/2023
*/

#include "../util.h"
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

struct Seed{
    int64_t score;
    size_t ct;
    vector<size_t> read_ids;
    Seed(int64_t x):score(x),ct(0),read_ids{}{}
};

void loadSeedsWithScore(const char* filename, const size_t read_id,
			unordered_map<kmer, Seed, kmerHash>& score_freq){
    FILE* fin = fopen(filename, "rb");
    size_t ret;
    int k;
    kmer s;
    int64_t score;
    ret = fread(&k, sizeof(k), 1, fin);
    while(true){
	ret = fread(&s, sizeof(s), 1, fin);
	if(ret != 1) break;
	//skip kmer
	//ret = fseek(fin, sizeof(kmer), SEEK_CUR);
	ret = fread(&score, sizeof(score), 1, fin);
	//if(ret != 1) break;
	//skip psi, st, index
	fseek(fin, sizeof(short)+(sizeof(uint64_t)<<1), SEEK_CUR);
	
	auto result = score_freq.emplace(piecewise_construct,
			                 forward_as_tuple(s),
					 forward_as_tuple(score));
	if(result.second == false){
	    if(result.first->second.read_ids.empty() ||
	       result.first->second.read_ids.back() < read_id){
		result.first->second.read_ids.push_back(read_id);
	    }
	}
	result.first -> second.ct += 1;
    }

    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}

int main(int argc, const char * argv[])    
{   
    if(argc != 7){
	fprintf(stderr, "usage: getScoreFreq.out seedsDir k repeatId numReads fileExt truePairFile\nrepeatId should be between 0 and floor((k-1)/2), seeds from both batches repeatId*2 and repeatId*2+1 are used (self symmetric if repeadId*2+1==k)\n");
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

    //load true pairs
    Table share_ct(n);
    
    FILE* fin = fopen(argv[6], "r");
    unsigned int x, y;
    int len;
    while(fscanf(fin, "%u %u %d\n", &x, &y, &len)==3){
	share_ct.access(x, y) = 1;
    }

    char filename[500];
    int i = strlen(argv[1]);
    memcpy(filename, argv[1], i);
    if(filename[i-1] != '/'){
	filename[i] = '/';
	++i;
    }

    unordered_map<kmer, Seed, kmerHash> score_freq;
    //unordered_map<kmer, vector<size_t>, kmerHash> all_seeds;
    int j;
    
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+i, "%d-%d.%s", r, j, argv[5]);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r, j, argv[5]);
	    break;
	}
	loadSeedsWithScore(filename, j, score_freq);
	//loadSeedsUnordered(filename, j, all_seeds);
	
	if(is_subseq2 == 0 && r2 >= k){
	    continue;
	}else{
            sprintf(filename+i, "%d-%d.%s", r2, j, argv[5]);
            if(stat(filename, &test_file) != 0){//seed file does not exist
                fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r2, j, argv[5]);
                break;
            }
	    loadSeedsWithScore(filename, j, score_freq);
            //loadSeedsUnordered(filename, j, all_seeds);
        }
    }

    sprintf(filename+i, "score_freq_tf-n%d-r%d.txt", n, r>>1);
    FILE* fout = fopen(filename, "w");
    
    char buf[k+1];

    size_t a, b, c, true_ct, false_ct;
    for(const auto& x : score_freq){
	true_ct = false_ct = 0;
	c = x.second.read_ids.size();
	for(i=0; i<c; ++i){
	    a = x.second.read_ids[i];
	    for(j=i+1; j<c; ++j){
		b = x.second.read_ids[j];
		if(share_ct.access(a, b)==1){
		    true_ct += 1;
		}else{
		    false_ct += 1;
		}
	    }
	}
	fprintf(fout, "%ld %.*s %zu %zu %zu\n", x.second.score,
		k, decode(x.first, k, buf), x.second.ct,
		true_ct, false_ct);
    }
    fclose(fout);
    
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

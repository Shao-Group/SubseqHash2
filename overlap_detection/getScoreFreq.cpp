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

void loadSeedsWithScore(const char* filename, map<int64_t, size_t>& score_freq){
    FILE* fin = fopen(filename, "rb");
    size_t ret;
    int k;
    int64_t score;
    ret = fread(&k, sizeof(k), 1, fin);
    while(true){
	//skip kmer
	ret = fseek(fin, sizeof(kmer), SEEK_CUR);
	ret = fread(&score, sizeof(score), 1, fin);
	if(ret != 1) break;
	//skip psi, st, index
	fseek(fin, sizeof(short)+(sizeof(uint64_t)<<1), SEEK_CUR);
	
	auto result = score_freq.emplace(score, 0);
	result.first -> second += 1;
    }

    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}

int main(int argc, const char * argv[])    
{   
    if(argc != 6){
	fprintf(stderr, "usage: getScoreFreq.out seedsDir k repeatId numReads fileExt\nrepeatId should be between 0 and floor((k-1)/2), seeds from both batches repeatId*2 and repeatId*2+1 are used (self symmetric if repeadId*2+1==k)\n");
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

    map<int64_t, size_t> score_freq;
    //unordered_map<kmer, vector<size_t>, kmerHash> all_seeds;
    int j;
    
    struct stat test_file;
    for(j=1; j<=n; j+=1){
	sprintf(filename+i, "%d-%d.%s", r, j, argv[5]);
	if(stat(filename, &test_file) != 0){//seed file does not exist
	    fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r, j, argv[5]);
	    break;
	}
	loadSeedsWithScore(filename, score_freq);
	//loadSeedsUnordered(filename, j, all_seeds);
	
	if(is_subseq2 == 0 && r2 >= k){
	    continue;
	}else{
            sprintf(filename+i, "%d-%d.%s", r2, j, argv[5]);
            if(stat(filename, &test_file) != 0){//seed file does not exist
                fprintf(stderr, "Stopped, cannot find file %d-%d.%s\n", r2, j, argv[5]);
                break;
            }
	    loadSeedsWithScore(filename, score_freq);
            //loadSeedsUnordered(filename, j, all_seeds);
        }
    }

    sprintf(filename+i, "score_freq-n%d-r%d.txt", n, r>>1);
    FILE* fout = fopen(filename, "w");
    
    for(const auto& x : score_freq){
	fprintf(fout, "%ld %zu\n", x.first, x.second);
    }
    fclose(fout);
    
    return 0;
}

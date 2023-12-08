/*
  Print out seeds on one read (readable by loadSeedsWithScore defined below),
  in the order they appear on the read.
  
  By: Ke@PSU
  Last edited: 11/27/2023
*/

#include "../util.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

struct Seed{
    int64_t score;
    size_t ct;
    vector<size_t> read_ids;
    Seed(int64_t x):score(x),ct(0),read_ids{}{}
};

void ioSeedsWithScore(const char* input_filename,
		      const char* output_filename){
    FILE* fin = fopen(input_filename, "rb");
    FILE* fout = fopen(output_filename, "w");
    size_t ret;
    int k;
    kmer s;
    int64_t score;
    short psi;
    uint64_t st;
    ret = fread(&k, sizeof(k), 1, fin);
    char buf[k+1];
    buf[k] = '\0';
    
    while(true){
	ret = fread(&s, sizeof(s), 1, fin);
	if(ret != 1) break;
	ret = fread(&score, sizeof(score), 1, fin);
	ret = fread(&psi, sizeof(psi), 1, fin);
	ret = fread(&st, sizeof(st), 1, fin);
	//skip index
	fseek(fin, sizeof(uint64_t), SEEK_CUR);

	
	fprintf(fout, "%lu %ld %hd %.*s\n", st, score, psi, k, decode(s, k, buf));
    }

    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", input_filename);
    }
    fclose(fin);
    fclose(fout);
}

int main(int argc, const char * argv[])    
{   
    if(argc != 3){
	fprintf(stderr, "usage: getSeedsOneRead.out seedFile outputFile\n");
	return 1;
    }
    
    struct stat test_file;
    if(stat(argv[1], &test_file) != 0){//seed file does not exist
	fprintf(stderr, "Stopped, cannot find file %s\n", argv[1]);
	return 1;
    }
    ioSeedsWithScore(argv[1], argv[2]);
    
    return 0;
}

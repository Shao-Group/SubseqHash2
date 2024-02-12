/*
For pairs of sequences at the same edit distance, compute their sets of subseq2seeds, compare each pair of seeds by hamming distance.

Input file is assumed to have 2*N lines for N pairs of sequences, each sequence is of length n. k seeds are generated for each sequence (where k is the parameter for subseqhash2). Output N lines, one for each pair; each line contains k numbers which are the hamming distances between the k pairs of seeds.

Output to stdout.

By: Ke @ PSU
Last edited: 02/12/2024
*/

#include "subseqhash2seeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>

using namespace std;

int calcHammingDist(kmer a, kmer b, int k){
    int r = 0;
    for(int i=0; i<k; ++i, a>>=2, b>>=2){
	if((a&0x3)!=(b&0x3)) r+=1;
    }
    return r;
}

int main(int argc, const char * argv[]){
    if(argc != 6){
	fprintf(stderr, "Usage: subseq2_hamming.out seqFile n k d randTableFile\n!!Output to stdout!!\n");
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
    
    subseqhash2seeding sub2(n, k, d, k);
    sub2.init(argv[5]);

    int chunk_size = sub2.getChunkSize();

    DPCell* dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
    DPCell* revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
    int* h = (int*) malloc(sizeof *h * dim1);
    int* revh = (int*) malloc(sizeof *revh * dim1);

    ifstream fin(argv[1]);

    string seq;

    while(getline(fin, seq))
    {
	vector<vector<seed>> seeds1(k, vector<seed>(0));
	sub2.getSubseq2Seeds(seq, dp, revdp, h, revh, seeds1);

	getline(fin, seq);
	vector<vector<seed>> seeds2(k, vector<seed>(0));
	sub2.getSubseq2Seeds(seq, dp, revdp, h, revh, seeds2);

	for(int i = 0; i < k; ++i){
	    kmer a = seeds1[i][0].str; //should be only one seed per repeat
	    kmer b = seeds2[i][0].str;

	    int hd = calcHammingDist(a, b, k);

	    printf("%d ", hd);
	}
	printf("\n");
    }

    free(dp);
    free(revdp);
    free(h);
    free(revh);

    return 0;
}

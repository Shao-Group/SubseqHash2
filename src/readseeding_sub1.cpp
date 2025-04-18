#include "subseqhash1seeding.hpp"
#include "util.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

DPCell* dp;
int* h;

FILE* fout;

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * d;

    clock_t start,end;

    subseqhash1seeding sub1(n, k, d);
    sub1.init(argv[5]);

	dp = (DPCell*) malloc(sizeof *dp * (dim1));
	h = (int*) malloc(sizeof *h * dim1);

	string path = "./readseed_sub1/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(subsample);
	fout = fopen(path.c_str(), "wb");


	ifstream refin(argv[6]);

	string info, refseq;
    uint64_t tmp = 0;
    kmer zkmer = 0;
    int num = 0;

    double anstime = 0;
    
	while(getline(refin, refseq))
	{
		getline(refin, info);
		vector<seed> seeds;
		
    	start = clock();
		sub1.DP(refseq, dp, h, seeds);
		end = clock();
	    anstime += (double)(end-start)/CLOCKS_PER_SEC;

	    uint64_t pos[2];//st, index
	    for(auto s : seeds)
	    {
			fwrite(&(s.hashval), sizeof(int64_t), 1, fout);
			fwrite(&(s.str), sizeof(kmer), 1, fout);
			pos[0] = s.st;
			pos[1] = s.index;
			fwrite(pos, sizeof(uint64_t), 2, fout);
	    }			

	    fwrite(&(tmp), sizeof(uint64_t), 1, fout);
		fwrite(&(zkmer), sizeof(kmer), 1, fout);
		pos[0] = 0;
		pos[1] = 0;
		fwrite(pos, sizeof(uint64_t), 2, fout);


		int lens = refseq.length();
		reverse(refseq.begin(), refseq.end());
		for(int i = 0; i < lens; i++)
			refseq[i] = ALPHABET[3-alphabetIndex(refseq[i])];

		seeds.clear();
		sub1.DP(refseq, dp, h, seeds);

	    for(auto s : seeds)
	    {
			fwrite(&(s.hashval), sizeof(int64_t), 1, fout);
			fwrite(&(s.str), sizeof(kmer), 1, fout);
			pos[0] = s.st;
			pos[1] = s.index;
			fwrite(pos, sizeof(uint64_t), 2, fout);
	    }			

	    fwrite(&(tmp), sizeof(uint64_t), 1, fout);
		fwrite(&(zkmer), sizeof(kmer), 1, fout);
		pos[0] = 0;
		pos[1] = 0;
		fwrite(pos, sizeof(uint64_t), 2, fout);

		num++;
	}
	//printf("%.2lf\n", anstime);

    return 0;
}
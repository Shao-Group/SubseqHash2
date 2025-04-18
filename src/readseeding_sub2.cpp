#include "subseq2simdseeding.hpp"
#include "util.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

const int maxd = 32;
int n,k,d, dim1;

DPCell* dp, *revdp;
short* h, *revh;

FILE* fout[50];

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * maxd;

    subseqhash2SIMDseeding sub2(n, k, d, subsample);
    sub2.init(argv[5]);

    clock_t start,end;


	int chunk_size = sub2.getChunkSize();

	dp = (DPCell*) aligned_alloc(32, sizeof *dp * chunk_size * (n+1) * (k+1));
	revdp = (DPCell*) aligned_alloc(32, sizeof *dp * chunk_size * (n+1) * (k+1));
	h = (short*) aligned_alloc(64, sizeof(short) * (dim1));
	revh = (short*) aligned_alloc(64, sizeof(short) * (dim1));

	ifstream refin(argv[6]);
	string info, refseq;

	for(int i = 0; i < subsample; i++)
	{		
		string path = "./readseed/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
		//saveSeedsPosotion(path.c_str(), seeds[i]);
		fout[i] = fopen(path.c_str(), "wb");
	}


    uint64_t tmp = 0;
    kmer zkmer = 0;
    int num = 0;
    double anstime = 0;

	while(getline(refin, refseq))
	{
		getline(refin, info);
		vector<vector<seed>> seeds(subsample, vector<seed>(0));

		//cout<<refseq<<endl;

    	start = clock();
		sub2.getSubseq2Seeds(refseq, dp, h, revdp, revh, seeds);
		end = clock();


	    anstime += (double)(end-start)/CLOCKS_PER_SEC;

		for(int i = 0; i < subsample; i++)
		{		
		    uint64_t pos[2];//st, index
		    for(auto s : seeds[i])
		    {
				fwrite(&(s.hashval), sizeof(int64_t), 1, fout[i]);
				fwrite(&(s.str), sizeof(kmer), 1, fout[i]);
				pos[0] = s.st;
				pos[1] = s.index;
				fwrite(pos, sizeof(uint64_t), 2, fout[i]);
		    }			

		    fwrite(&(tmp), sizeof(uint64_t), 1, fout[i]);
		    fwrite(&(zkmer), sizeof(kmer), 1, fout[i]);
			pos[0] = 0;
			pos[1] = 0;
			fwrite(pos, sizeof(uint64_t), 2, fout[i]);
		}

		num++;
	}

	//printf("%.2lf\n", anstime);
    return 0;
}
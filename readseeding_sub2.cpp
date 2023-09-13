#include "subseqhash2seeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

DPCell* dp, *revdp;
int* h, *revh;

FILE* fout[50];

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * d;

    subseqhash2seeding sub2(n, k, d, subsample);
    sub2.init(argv[5]);


	int chunk_size = sub2.getChunkSize();

	dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
	revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
	h = (int*) malloc(sizeof *h * dim1);
	revh = (int*) malloc(sizeof *revh * dim1);

	string species = argv[6];

	for(int i = 0; i < subsample; i++)
	{		
		string path = "./readsub2/" + species + "/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
		//saveSeedsPosotion(path.c_str(), seeds[i]);
		fout[i] = fopen(path.c_str(), "wb");
	}



	ifstream refin("./reads/" + species);
	string info, refseq;
    uint64_t tmp = 0;
    kmer zkmer = 0;
    int num = 0;

	while(getline(refin, refseq))
	{
		getline(refin, info);
		getline(refin, info);
		vector<vector<seed>> seeds(subsample, vector<seed>(0));

		//cout<<refseq<<endl;

		sub2.getSubseq2Seeds(refseq, dp, revdp, h, revh, seeds);

		for(int i = 0; i < subsample; i++)
		{		
			//cout<<i<<" "<<seeds[i].size()<<endl;
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
		    fwrite(&(zkmer), sizeof(uint64_t), 1, fout[i]);
			pos[0] = 0;
			pos[1] = 0;
			fwrite(pos, sizeof(uint64_t), 2, fout[i]);
		}

		num++;
	}

    return 0;
}
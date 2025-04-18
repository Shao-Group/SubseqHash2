#include "subseq2strobeseeding.hpp"
#include "util.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;
int prek;

DPCell* dp, *revdp;
int* h, *revh;

FILE* fout[50];

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);
    prek = atoi(argv[5]);

    dim1 = (n+1) * (k+1) * d;

    subseq2strobeseeding sub2(n, k, d, subsample, 1, prek);
    sub2.init(argv[6]);

    clock_t start,end;

	int chunk_size = sub2.getChunkSize();

	dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
	revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
	h = (int*) malloc(sizeof *h * dim1);
	revh = (int*) malloc(sizeof *revh * dim1);

	ifstream refin(argv[7]);
	string info, refseq;

	for(int i = 0; i < subsample; i++)
	{		
		string path = "./readseed_substro/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i) + "_" + to_string(prek);
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

    	start = clock();
		sub2.getSubseq2Seeds(refseq, dp, revdp, h, revh, seeds);
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

		int lens = refseq.length();
		reverse(refseq.begin(), refseq.end());
		for(int i = 0; i < lens; i++)
			refseq[i] = ALPHABET[3-alphabetIndex(refseq[i])];

		for(int i = 0; i < subsample; i++)
			seeds[i].clear();

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
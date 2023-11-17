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

	string path = argv[7];
	FILE* fout = fopen(path.c_str(), "w");



	ifstream refin(argv[6]);

	string refseq;

	while(getline(refin, refseq))
	{
		vector<vector<seed>> seeds(subsample, vector<seed>(0));

		sub2.getSubseq2Seeds(refseq, dp, revdp, h, revh, seeds);

		for(int i = 0; i < subsample; i++)
		{		
		    string name = "Order " + to_string(i+1) + '\n';

		    fputs(name.c_str(), fout);
		    
		    for(auto s : seeds[i])
		    {
		    	fprintf(fout, "(%d, %ld)\t", s.psi, s.hashval);

				char* seedstr = (char*)malloc(sizeof(char) *k);
				decode(s.str, k, seedstr);
				fputs(seedstr, fout);
				fputs("\t", fout);

				fputs("\n", fout);
		    }			
			fputs("\n", fout);
		}
	}

    return 0;
}
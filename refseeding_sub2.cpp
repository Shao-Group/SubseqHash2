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

	string species = argv[6];

	ifstream refin("./ref/" + species);
	string info, refseq;

	while(getline(refin, info))
	{
		getline(refin, refseq);
		vector<vector<seed>> seeds(subsample, vector<seed>(0));

		info = info.substr(1, info.find(' ') - 1);


		sub2.getSubseq2Seeds(refseq, dp, revdp, h, revh, seeds);


		for(int i = 0; i < subsample; i++)
		{		
			string path = "./refsub2/" + species + "/" + info + "_" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
			saveSeedsStrPosition(path.c_str(), seeds[i]);
		}
	}

    return 0;
}
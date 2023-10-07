#include "subseq2strobeseeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;
int prek;

DPCell* dp, *revdp;
int* h, *revh;

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


	int chunk_size = sub2.getChunkSize();

	dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
	revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
	h = (int*) malloc(sizeof *h * dim1);
	revh = (int*) malloc(sizeof *revh * dim1);

	string species = argv[7];

	ifstream refin("./ref/" + species);
	string info, refseq;

	while(getline(refin, info))
	{
		getline(refin, refseq);
		vector<vector<seed>> seeds(subsample, vector<seed>(0));
		vector<FILE*> fout;

		info = info.substr(1, info.find(' ') - 1);

		for(int i = 0; i < subsample; i++)
		{		
			string path = "./refsubstro/" + species + "/" + info + "_" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i) + "_" + to_string(prek);
			fout.push_back(fopen(path.c_str(), "wb"));
		}
		sub2.writeSubseq2Seeds(refseq, dp, revdp, h, revh, seeds, fout);
	}

    return 0;
}
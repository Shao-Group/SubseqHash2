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

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * maxd;

    subseqhash2SIMDseeding sub2(n, k, d, subsample);
    sub2.init(argv[5]);


	int chunk_size = sub2.getChunkSize();

	dp = (DPCell*) aligned_alloc(32, sizeof *dp * chunk_size * (n+1) * (k+1));
	revdp = (DPCell*) aligned_alloc(32, sizeof *dp * chunk_size * (n+1) * (k+1));
	h = (short*) aligned_alloc(64, sizeof(short) * (dim1));
	revh = (short*) aligned_alloc(64, sizeof(short) * (dim1));

	// string species = argv[6];

	ifstream refin(argv[6]);
	string info, refseq;

	getline(refin, refseq);
	getline(refin, refseq);
	vector<vector<seed>> seeds(subsample, vector<seed>(0));
	vector<FILE*> fout;


	for(int i = 0; i < subsample; i++)
	{		
		string path = "./refseed/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
		fout.push_back(fopen(path.c_str(), "wb"));
	}

	sub2.writeSubseq2Seeds(refseq, dp, revdp, h, revh, seeds, fout);


    return 0;
}
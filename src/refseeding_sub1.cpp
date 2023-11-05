#include "subseqhash1seeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

DPCell* dp;
int* h;

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * d;

    subseqhash1seeding sub1(n, k, d);
    sub1.init(argv[5]);


	dp = (DPCell*) malloc(sizeof *dp * (dim1));
	h = (int*) malloc(sizeof *h * dim1);

	string species = argv[6];

	ifstream refin("./ref/" + species);
	string info, refseq;

	while(getline(refin, info))
	{
		getline(refin, refseq);
		vector<seed> seeds;

		info = info.substr(1, info.find(' ') - 1);


		sub1.DP(refseq, dp, h, seeds);

		string path = "./refsub1/" + species + "/" + info + "_" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(subsample);
		saveSeedsPosition(path.c_str(), seeds);
	}

    return 0;
}
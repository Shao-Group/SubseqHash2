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


	ifstream refin(argv[6]);
	string refseq;

	getline(refin, refseq);
	getline(refin, refseq);

	vector<seed> seeds;
	sub1.DP(refseq, dp, h, seeds);

	string path = "./refseed_sub1/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(subsample);
	saveSeedsStrPosition(path.c_str(), seeds);
    return 0;
}
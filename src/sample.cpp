#include "subseqhash2seeding.hpp"
#include "util.hpp"
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

	ifstream refin(argv[6]);

	string refseq;

	while(getline(refin, refseq))
	{
		vector<vector<seed>> seeds(subsample, vector<seed>(0));
		getline(refin, refseq);

		sub2.getSubseq2Seeds(refseq, dp, revdp, h, revh, seeds);

		for(int i = 0; i < subsample; i++)
		{		
			printf("Order %d:", i+1);
		    for(auto s : seeds[i]) {
				printf("%ld %I64ld ", s.st, s.hashval);
				for(int i = 0; i < n; i++)
					if((s.index>>i) & 1)
						printf("%c", refseq[s.st + i]);
				printf("\n");
			}
		}
	}

    return 0;
}
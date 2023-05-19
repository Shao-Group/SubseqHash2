#include "subseqhash2seeding.h"
#include "seedfactory.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

int main(int argc, const char * argv[])
{    
    if(argc != 6)
    {
		printf("usage: genSubseq2Seeds.out readFile n k d randTableFile\n");
		return 1;
    }

    n = atoi(argv[2]);
    k = atoi(argv[3]);
    d = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * d;


	ifstream fin(argv[1]);

	string seq, seq2, tmp;
	int x;
	int num = 0;
    
    char output_dir[200];
    int dir_len = sprintf(output_dir, "./seeds");

    seedFactory factory("./seeds", dir_len, n, k, d, argv[5]);

	while(fin>>seq)
	{
		fin>>seq2;				

		int len = seq.length();

    	vector<int> align;
		for(int i = 0; i < len; i++)
		{
			fin>>x;
			align.push_back(x);
		}
	    factory.addJob(move(seq), num++);
	    factory.addJob(move(seq2), num++);
		num++;
	}	

    return 0;
}
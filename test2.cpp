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
    if(argc != 7)
    {
		printf("usage: test2.out readFile n k d subsample randTableFile\n");
		return 1;
    }

    n = atoi(argv[2]);
    k = atoi(argv[3]);
    d = atoi(argv[4]);
    int subsample = atoi(argv[5]);

    dim1 = (n+1) * (k+1) * d;


	ifstream fin(argv[1]);

	string seq, seq2, tmp;
	int x;
	int num = 0;
    
    char output_dir[200];
    int dir_len = sprintf(output_dir, "./seeds");

    seedFactory factory("./seeds", dir_len, n, k, d, subsample , argv[6]);

	while(fin>>tmp)
	{			
		fin>>tmp;
		fin>>seq;
		fin>>tmp;		
		fin>>tmp;
		fin>>seq2;				

		int len = seq.length();

cout<<seq.length()<<" "<<seq2.length()<<endl;
	    factory.addJob(move(seq), num++);
	    factory.addJob(move(seq2), num++);
		num++;
	}	

    return 0;
}
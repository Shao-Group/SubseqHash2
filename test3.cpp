#include "minimizerseeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 
const int maxlen = 1<<20;

vector<seedmatch> matches;

vector<bool> scover(maxlen, 0);
vector<bool> scover2(maxlen, 0);
vector<bool> mcover(maxlen, 0);
vector<bool> mcover2(maxlen, 0);
vector<bool> fscover(maxlen, 0);
vector<bool> fscover2(maxlen, 0);

void pseudo_match(string s, string t, minimizerseeding & seeding)
{
	int lens = s.length();
	int lent = t.length();


	// int seednum = sub2.getNumPerWindow();

	// vector<vector<seed>> seeds(seednum, vector<seed>(0));
	// vector<vector<seed>> seedt(seednum, vector<seed>(0));
	vector<seed> seeds;
	vector<seed> seedt;

    clock_t start,end;

    start = clock();

	seeding.get_minimizers(s, seeds);
	seeding.get_minimizers(t, seedt);
    end = clock();

    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;

    fill(scover.begin(), scover.end(), false);
    fill(scover2.begin(), scover2.end(), false);
    fill(mcover.begin(), mcover.end(), false);
    fill(mcover2.begin(), mcover2.end(), false);
    fill(fscover.begin(), fscover.end(), false);
    fill(fscover2.begin(), fscover2.end(), false);

	int totalmatches = 0;
	int totaltruematches = 0;

	int sz = seeds.size();
	for(int i = 0; i < sz; i++)
	{
		seed tmp;
		tmp.st = seeds[i].st;
		tmp.ed = seeds[i].ed;
		tmp.str = tmp.hashval = revComp(seeds[i].hashval, n);
		seeds.push_back(tmp);
	}

	sz = seedt.size();
	for(int i = 0; i < sz; i++)
	{
		seed tmp;
		tmp.st = seedt[i].st;
		tmp.ed = seedt[i].ed;
		tmp.str = tmp.hashval = revComp(seedt[i].hashval, n);
		seedt.push_back(tmp);
	}

	uint64_t x1;
	cout<<s.size()<<" "<<t.size()<<endl;
	cout<<seeds.size()<<" "<<seedt.size()<<endl;
	start = clock();
	ssh_index* ht = index_build(seeds);
	matches.clear();

	index_get(ht, seedt, matches);
	end = clock();
	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;


	int truematches = 0;
	cout<<matches.size()<<endl;
	for(auto m: matches)
		cout<<m.s1->st<<" "<<m.s2->st<<endl;
	ans[7] += double(seeds.size())/s.length() + double(seedt.size())/t.length();
	totalmatches += matches.size();
	totaltruematches += truematches;
}

int main(int argc, const char * argv[])
{    
    if(argc != 7)
    {
		printf("usage: genSubseq2Seeds.out readFile n k d subsample randTableFile\n");
		return 1;
    }

    n = atoi(argv[2]);
    k = atoi(argv[3]);
    d = atoi(argv[4]);
    int subsample = atoi(argv[5]);
    dim1 = (n+1) * (k+1) * d;
	
	minimizerseeding seeding(n, k);

    //subseqhash2seeding sub2(n, k, d, subsample);
    //sub2.init(argv[6]);

	//int chunk_size = sub2.getChunkSize();

	// dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
	// revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
	// h = (int*) malloc(sizeof *h * dim1);
	// revh = (int*) malloc(sizeof *revh * dim1);

	ifstream fin(argv[1]);

	string seq, seq2, tmp;
	int x;
	int num = 0;

	fin>>tmp;
	fin>>tmp;
	fin>>seq;
	fin>>tmp;
	fin>>tmp;
	fin>>seq2;				

	int len = seq.length();
	int lent = seq2.length();

	pseudo_match(seq, seq2, seeding);
	num++;
	printf("%d/%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", n, k, d, subsample,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (2*num), ans[4] / (2*num), ans[5] / (2*num), 
		ans[6] / (2*num), ans[7] / (2*num), ans[8], ans[9]);

    // const char* table_filename = argv[5];

    return 0;
}
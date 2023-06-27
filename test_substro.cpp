#include "subseq2strobeseeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;
int prek;

double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 
DPCell* dp, *revdp;
int* h, *revh;
const int maxlen = 1<<20;

vector<seedmatch> matches;

vector<bool> scover(maxlen, 0);
vector<bool> scover2(maxlen, 0);
vector<bool> mcover(maxlen, 0);
vector<bool> mcover2(maxlen, 0);
vector<bool> fscover(maxlen, 0);
vector<bool> fscover2(maxlen, 0);

void pseudo_match(string s, string t, vector<int> &align, subseq2strobeseeding & sub2)
{
	int lens = s.length();
	int lent = t.length();


	int seednum = sub2.getNumPerWindow();

	vector<vector<seed>> seeds(seednum, vector<seed>(0));
	vector<vector<seed>> seedt(seednum, vector<seed>(0));

    clock_t start,end;

    start = clock();
	sub2.getSubseq2Seeds(s, dp, revdp, h, revh, seeds);
	sub2.getSubseq2Seeds(t, dp, revdp, h, revh, seedt);
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

	uint64_t x1;

	for(int j = 0; j < seednum; j++)
	{	
    	start = clock();
		ssh_index* ht = index_build(seeds[j]);
		matches.clear();

		index_get(ht, seedt[j], matches);
    	end = clock();
    	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;


		int truematches = 0;
		for(seedmatch m: matches)
		{
			int tp = 0;

			for(int i = 0; i < n + prek; i++)
				if((m.s1->index>>i) & 1)
				{
					x1 = align[m.s1->st + i]; 
					if(x1 >= m.s2->st && x1 <= m.s2->ed)
						tp += (m.s2->index>>(x1-m.s2->st)) & 1;
				}

			if(2 * tp >= k)
			{
				truematches++;	
				for(int i = 0; i < n + prek; i++)
					if((m.s1->index>>i) & 1)
						scover[m.s1->st + i] = 1;
				for(int i = 0; i < n + prek; i++)
					if((m.s2->index>>i) & 1)
						scover2[m.s2->st + i] = 1;

				for(int i = m.s1->st; i <= m.s1->ed; i++)
					mcover[i] = 1;
				for(int i = m.s2->st; i <= m.s2->ed; i++)
					mcover2[i] = 1;
			}
			else
			{			
				for(int i = 0; i < n + prek; i++)
					if((m.s1->index>>i) & 1)
						fscover[m.s1->st + i] = 1;
				for(int i = 0; i < n + prek; i++)
					if((m.s2->index>>i) & 1)
						fscover2[m.s2->st + i] = 1;
			}
		}

		ans[7] += double(seeds[j].size())/s.length() + double(seedt[j].size())/t.length();
		totalmatches += matches.size();
		totaltruematches += truematches;
	}		

	int sc = 0, sc2 = 0;
	int fsc = 0, fsc2 = 0;
	int mc = 0, mc2 = 0;
	int island1 = 0, island2 = 0;
	int inv = 0;

	for(int i = 0; i < lens; i++)
	{
		sc += scover[i];
		fsc += fscover[i];
		mc += mcover[i];

		if(!scover[i])
			inv++;
		else if(inv)
		{
			island1 += inv * inv;
			inv = 0;
		}
	}

	if(inv)
	{
		island1 += inv * inv;
		inv = 0;
	}

	for(int i = 0; i < lent; i++)
	{
		sc2 += scover2[i];
		fsc2 += fscover2[i];
		mc2 += mcover2[i];

		if(!scover2[i])
			inv++;
		else if(inv)
		{
			island2 += inv * inv;
			inv = 0;
		}
	}
	if(inv)
	{
		island2 += inv * inv;
		inv = 0;
	}

	ans[0] += totalmatches;
	ans[1] += totaltruematches;
	if(totalmatches > 0)
	{
		ans[10] += 1;
		ans[2] += double(totaltruematches)/totalmatches;
	}
	ans[3] += double(sc2) / lent + double(sc) / lens;
	ans[4] += double(fsc2) / lent + double(fsc) / lens;
	ans[5] += double(mc2) / lent + double(mc) / lens;
	ans[6] += double(island2) / lent + double(island1) / lens;
}

int main(int argc, const char * argv[])
{    
    if(argc != 8)
    {
		printf("usage: genSubseq2Seeds.out readFile n k d subsample prek randTableFile\n");
		return 1;
    }

    n = atoi(argv[2]);
    k = atoi(argv[3]);
    d = atoi(argv[4]);
    int subsample = atoi(argv[5]);
    prek = atoi(argv[6]);
    dim1 = (n+1) * (k+1) * d;

    subseq2strobeseeding sub2(n, k, d, subsample, prek);
    sub2.init(argv[7]);

	int chunk_size = sub2.getChunkSize();

	dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
	revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
	h = (int*) malloc(sizeof *h * dim1);
	revh = (int*) malloc(sizeof *revh * dim1);

	ifstream fin(argv[1]);

	string seq, seq2, tmp;
	int x;
	int num = 0;

	while(fin>>seq)
	{
		fin>>seq2;				

		int len = seq.length();
		int lent = seq2.length();

    	vector<int> align;
		for(int i = 0; i < len; i++)
		{
			fin>>x;
			align.push_back(x);
		}

		pseudo_match(seq, seq2, align, sub2);
		num++;
	}	
	printf("%d/%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", n, k, d, subsample,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (2*num), ans[4] / (2*num), ans[5] / (2*num), 
		ans[6] / (2*num), ans[7] / (2*num), ans[8], ans[9]);

    // const char* table_filename = argv[5];

    return 0;
}
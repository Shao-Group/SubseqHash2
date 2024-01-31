#include "util.h"
#include <bits/stdc++.h>

using namespace std;

int n, k, d, subsample;
ssh_index* refindex;

double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 

vector<seedmatch> matches;

vector<vector<bool>> bcover;
double matchdata[5100][2] = {0};


void pseudo_match(int len, int readnum, vector<seed> &seedt, vector<int> &align)
{
	matches.clear();

	int totalmatches = 0;
	int truematches = 0;
	uint64_t x1;

	index_get(refindex, seedt, matches);

	totalmatches += matches.size();

	for(seedmatch m: matches)
	{
		int tp = 0;

		if(m.s1->str != m.s2->str && m.s1->str != m.s2->str_rc)
		{
			totalmatches--;
			continue;
		}

		for(int i = 0; i < n; i++)
			if((m.s2->index>>i) & 1)
			{
				x1 = align[m.s2->st + i]; 
				if(x1 >= m.s1->st && x1 <= m.s1->ed)
					tp += (m.s1->index>>(x1-m.s1->st)) & 1;
			}

		if(2 * tp >= k)
		{
			truematches++;	
			bcover[readnum][m.s2->st/200] = 1;
		}
	}

	matchdata[readnum][0] += totalmatches;
	matchdata[readnum][1] += truematches;
}

void stats(int readnum, int lent)
{
	ans[0] += matchdata[readnum][0];
	ans[1] += matchdata[readnum][1];
	if(matchdata[readnum][0] > 0)
	{
		ans[10] += 1;
		ans[2] += double(matchdata[readnum][1])/matchdata[readnum][0];
	}

	int nblock = lent/200;
	int bc = 0;
	for(int i = 0; i < nblock; i++)
		bc += bcover[readnum][i];

	ans[3] += nblock;
	ans[4] += bc;
}

vector<int> reads;
vector<vector<int>> aligned;
vector<seed> seeds;

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    subsample = atoi(argv[4]);

	ifstream readin(argv[5]);
	int x, num = 0;
	string refseq;

	while(readin>>refseq)
	{
		vector<int> aligninput;

		int len = refseq.length();
		for(int j = 0; j < len; j++)
		{
			readin>>x;
			aligninput.push_back(x);
		}

		reads.push_back(len);
		aligned.push_back(aligninput);

		bcover.push_back(vector<bool>(len/200 + 1, 0));

		num++;
	}	

    uint64_t st, index;
    int64_t hashval;
	kmer kk;

	for(int i = 0; i < subsample; i++)
	{
		string refpath = "./refseed/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
		seeds.clear();
		loadSeedsStr(refpath.c_str(), seeds, k);
		refindex = index_build(seeds);
	
		string qpath = "./readseed/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
		string qpath_rc = "./readseed/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i^1);
	    
	    FILE* fin = fopen(qpath.c_str(), "rb"), *fin_rc;
	    if((i^1) < k)
		    fin_rc = fopen(qpath_rc.c_str(), "rb");

		for(int j = 0; j < num; j++)
		{
			vector<seed> seedt;

			size_t ret = 1;

		    while(ret == 1)
		    {
				ret = fread(&hashval, sizeof(int64_t), 1, fin);
				(void)!fread(&kk, sizeof(kmer), 1, fin);
				(void)!fread(&st, sizeof(uint64_t), 1, fin);
				(void)!fread(&index, sizeof(uint64_t), 1, fin);

				if(hashval == 0 && st == 0 && index == 0)
					break;
				seed tmp;
				tmp.hashval = hashval;
				tmp.st = st;
				tmp.index = index;		
				tmp.str = kk;
				tmp.str_rc = revComp(kk, k);

				uint64_t high1 = ((uint64_t)1)<<63;
				tmp.ed = st + 63;
				while(high1)
					if(index & high1)
						break;
					else
					{
						high1 >>= 1;
						tmp.ed--;
					}

				seedt.push_back(tmp);
			}	

	    	if((i^1) < k)
			    while(ret == 1)
			    {
					ret = fread(&hashval, sizeof(int64_t), 1, fin_rc);
					(void)!fread(&kk, sizeof(kmer), 1, fin_rc);
					(void)!fread(&st, sizeof(uint64_t), 1, fin_rc);
					(void)!fread(&index, sizeof(uint64_t), 1, fin_rc);

					if(hashval == 0 && st == 0 && index == 0)
						break;
					seed tmp;
					tmp.hashval = hashval;
					tmp.st = st;
					tmp.index = index;		
					tmp.str = kk;
					tmp.str_rc = revComp(kk, k);

					uint64_t high1 = ((uint64_t)1)<<63;
					tmp.ed = st + 63;
					while(high1)
						if(index & high1)
							break;
						else
						{
							high1 >>= 1;
							tmp.ed--;
						}
					
					seedt.push_back(tmp);
				}	

			pseudo_match(reads[j], j, seedt, aligned[j]);
		}
	}

	for(int i = 0; i < num; i++)
		stats(i, reads[i]);

	printf("%d/%d/%d/%d, %d, %d, %.4lf, %.4lf, %d, %d, %.4lf\n", n, k, d, subsample, int(ans[0]), int(ans[1]), ans[1]/ans[0], ans[2] / ans[10], int(ans[3]), int(ans[4]), ans[4]/ans[3]);

    return 0;
}
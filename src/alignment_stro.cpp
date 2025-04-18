#include "strobemerseeding.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, w, w_min, w_max, seednum;

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

void pseudo_match(string s, string t, vector<int> &align, strobemerseeding seeding)
{
	int lens = s.length();
	int lent = t.length();

	vector<vector<seed>> seeds(seednum, vector<seed>(0));
	vector<vector<seed>> seedt(seednum, vector<seed>(0));

    fill(scover.begin(), scover.end(), false);
    fill(scover2.begin(), scover2.end(), false);
    fill(mcover.begin(), mcover.end(), false);
    fill(mcover2.begin(), mcover2.end(), false);
    fill(fscover.begin(), fscover.end(), false);
    fill(fscover2.begin(), fscover2.end(), false);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

    clock_t start,end;

    start = clock();
	seeding.get_strobemers(s, seeds);
	seeding.get_strobemers(t, seedt);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    
	for(int j = 0; j < seednum; j++)
	{		
		start = clock();
		ssh_index* ht = index_build(seeds[j]);
		matches.clear();

		index_get(ht, seedt[j], matches);
		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
			
		int truematches = 0;
		int index1, index2, index3, index4;
		for(seedmatch m: matches)
		{			
			if(m.s1->str != m.s2->str)
			{
				totalmatches--;
				continue;
			}

			int tp = 0;

			for(int i = 0; i < k; i++)
				tp += (align[m.s1->st + i] == m.s2->st + i);
			
			index1 = m.s1->index % (1<<10);
			index2 = m.s2->index % (1<<10);

			for(int i = 0; i < k; i++)
				tp += (align[m.s1->st + index1 + i] == m.s2->st + index2 + i);

			if(w == 3)
			{
				index3 = (m.s1->index>>10);
				index4 = (m.s2->index>>10);

				for(int i = 0; i < k; i++)
					tp += (align[m.s1->st + index3 + i] == m.s2->st + index4 + i);
			}

			if(tp >= k)
			{
				truematches++;	

				for(int i = m.s1->st; i <= m.s1->ed; i++)
					mcover[i] = 1;
				for(int i = m.s2->st; i <= m.s2->ed; i++)
					mcover2[i] = 1;
				for(int i = 0; i < k; i++)
				{
					scover[m.s1->st + i] = 1;
					scover[m.s1->st + index1 + i] = 1;
					scover2[m.s2->st + i] = 1;
					scover2[m.s2->st + index2 + i] = 1;
					if(w == 3)
					{
						scover[m.s1->st + index3 + i] = 1;
						scover2[m.s2->st + index4 + i] = 1;
					}
				}
			}
			else
			{	
				for(int i = 0; i < k; i++)
				{
					fscover[m.s1->st + i] = 1;
					fscover[m.s1->st + index1 + i] = 1;
					fscover2[m.s2->st + i] = 1;
					fscover2[m.s2->st + index2 + i] = 1;
					if(w == 3)
					{
						fscover[m.s1->st + index3 + i] = 1;
						fscover2[m.s2->st + index4 + i] = 1;
					}
				}
			}
		}

		ans[7] += double(seeds.size())/s.length() + double(seedt.size())/t.length();
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
    k = atoi(argv[2]);
    w = atoi(argv[3]);
    w_min = atoi(argv[4]);
    w_max = atoi(argv[5]);
    seednum = atoi(argv[6]);

	strobemerseeding seeding(k, w, w_min, w_max, seednum);

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

		pseudo_match(seq, seq2, align, seeding);
		num++;
	}	
	printf("%d/%d/%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", k, w, w_min, w_max, seednum,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (2*num), ans[4] / (2*num), ans[5] / (2*num), 
		ans[6] / (2*num), ans[7] / (2*num), ans[8], ans[9]);

    // const char* table_filename = argv[5];

    return 0;
}
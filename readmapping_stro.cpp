#include "strobemerseeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, w, w_min, w_max, seednum;
vector<string> refname;
vector<vector<ssh_index*>> refindex;

int refnum = 0;
double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 
const int maxlen = 50000;

vector<seedmatch> matches;

vector<bool> scover(maxlen, 0);
vector<bool> scover2(maxlen, 0);
vector<bool> mcover(maxlen, 0);
vector<bool> mcover2(maxlen, 0);
vector<bool> fscover(maxlen, 0);
vector<bool> fscover2(maxlen, 0);

void pseudo_match(string s, vector<int> &align, string name, strobemerseeding seeding, strobemerseeding revseeding)
{
	int lens = s.length();

	vector<vector<seed>> seeds(seednum);

    fill(scover.begin(), scover.end(), false);
    fill(mcover.begin(), mcover.end(), false);
    fill(fscover.begin(), fscover.end(), false);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

    clock_t start,end;

    start = clock();
	seeding.get_strobemers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    int index1, index2;

    for(int j = 0; j < refnum; j++)
    if(refname[j] == name)
    {
    	for(int r = 0; r < ((seednum+1) >> 1); r++)
    	{
			start = clock();

			matches.clear();
			index_get(refindex[j][r], seeds[r], matches);

			end = clock();
			ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
				
			int truematches = 0;
			for(seedmatch m: matches)
			{
				int tp = 0;

				for(int i = 0; i < k; i++)
					tp += (align[m.s2->st + i] == m.s1->st + i);
				
				index1 = m.s1->index % (1<<10);
				index2 = m.s2->index % (1<<10);

				for(int i = 0; i < k; i++)
					tp += (align[m.s2->st + index2 + i] == m.s1->st + index1 + i);

				if(tp >= k)
				{
					truematches++;	

					for(int i = m.s2->st; i <= m.s2->ed; i++)
						mcover[i] = 1;
					
					for(int i = 0; i < k; i++)
					{
						scover[m.s2->st + i] = 1;
						scover[m.s2->st + index2 + i] = 1;
					}
				}
				else
				{	
					for(int i = 0; i < k; i++)
					{
						fscover[m.s2->st + i] = 1;
						fscover[m.s2->st + index2 + i] = 1;
					}
				}
			}

			ans[7] += double(seeds.size())/s.length();
			totalmatches += matches.size();
			totaltruematches += truematches;
		}
	}
	else
	{	
    	for(int r = 0; r < ((seednum+1) >> 1); r++)
		{
			start = clock();
			matches.clear();
			index_get(refindex[j][r], seeds[r], matches);
			end = clock();
			ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

			for(seedmatch m: matches)
			{
				index2 = m.s2->index % (1<<10);

				for(int i = 0; i < k; i++)
				{
					fscover[m.s2->st + i] = 1;
					fscover[m.s2->st + index2 + i] = 1;
				}
			}

			totalmatches += matches.size();
		}
	}

	reverse(s.begin(), s.end());
	for(int i = 0; i < lens; i++)
		s[i] = ALPHABET[3-alphabetIndex(s[i])];

    int fh = (seednum+1) >> 1;
	for(int i = 0; i < fh; i++)
		seeds[i].clear();

    start = clock();
	revseeding.get_strobemers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    
    for(int j = 0; j < refnum; j++)
    if(refname[j] == name)
    {
    	for(int r = 0; r < (seednum >> 1); r++)
    	{
			start = clock();

			matches.clear();
			index_get(refindex[j][fh+r], seeds[r], matches);

			end = clock();
			ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
				
			int truematches = 0;
			for(seedmatch m: matches)
			{
				int tp = 0;
				// cout<<m.s1->st<<" "<<m.s1->ed<<" "<<m.s2->st<<" "<<m.s2->ed<<endl;
				// cout<<lens-m.s2->st-1<<" "<<align[lens-m.s2->st-1]<<endl;
				for(int i = 0; i < k; i++)
					tp += (align[lens - m.s2->st - i - 1] == m.s1->st + i);
				
				index1 = m.s1->index % (1<<10);
				index2 = m.s2->index % (1<<10);

				for(int i = 0; i < k; i++)
					tp += (align[lens - m.s2->st - index2 - i - 1] == m.s1->st + index1 + i);

				if(tp >= k)
				{
					truematches++;	

					for(int i = lens - m.s2->ed - 1; i <= lens - m.s2->st - 1; i++)
						mcover[i] = 1;
					
					for(int i = 0; i < k; i++)
					{
						scover[lens - m.s2->st - i - 1] = 1;
						scover[lens - m.s2->st - index2 - i - 1] = 1;
					}
				}
				else
				{	
					for(int i = 0; i < k; i++)
					{
						fscover[lens - m.s2->st - i - 1] = 1;
						fscover[lens - m.s2->st - index2 - i - 1] = 1;
					}
				}
			}

			ans[7] += double(seeds.size())/s.length();
			totalmatches += matches.size();
			totaltruematches += truematches;
		}
	}
	else
	{	
    	for(int r = 0; r < (seednum >> 1); r++)
		{
			start = clock();
			matches.clear();
			index_get(refindex[j][fh + r], seeds[r], matches);
			end = clock();
			ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

			for(seedmatch m: matches)
			{
				index2 = m.s2->index % (1<<10);

				for(int i = 0; i < k; i++)
				{
					fscover[lens - m.s2->st - i - 1] = 1;
					fscover[lens - m.s2->st - index2 - i - 1] = 1;
				}
			}

			totalmatches += matches.size();
		}
	}
	int sc = 0;
	int fsc = 0;
	int mc = 0;
	int island1 = 0;
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


	ans[0] += totalmatches;
	ans[1] += totaltruematches;
	if(totalmatches > 0)
	{
		ans[10] += 1;
		ans[2] += double(totaltruematches)/totalmatches;
	}
	ans[3] += double(sc) / lens;
	ans[4] += double(fsc) / lens;
	ans[5] += double(mc) / lens;
	ans[6] += double(island1) / lens;
}

vector<vector<vector<seed>>> seedu;

int main(int argc, const char * argv[])
{    
    k = atoi(argv[1]);
    w = atoi(argv[2]);
    w_min = atoi(argv[3]);
    w_max = atoi(argv[4]);
    seednum = atoi(argv[5]);

	strobemerseeding seeding(k, w, w_min, w_max, (seednum+1)>>1);
	strobemerseeding revseeding(k, w, w_min, w_max, seednum>>1);

	string species = argv[6];

	ifstream refin("./ref/" + species);
	string info, refseq;
	refnum = 0;

	while(getline(refin, info))
	{
		getline(refin, refseq);

		vector<vector<seed>> seeds ((seednum+1)>>1);
		seeding.get_strobemers(refseq, seeds);
		seedu.push_back(seeds);

		vector<vector<seed>> seedt (seednum>>1);
		revseeding.get_strobemers(refseq, seedt);

		for(int i = 0; i < (seednum>>1); i++)
			seedu[refnum].push_back(seedt[i]);

		info = info.substr(1, info.find(' ') - 1);

		refname.push_back(info);

		vector<ssh_index*> tmpindex;

		for(int i = 0; i < seednum; i++)
			tmpindex.push_back(index_build(seedu[refnum][i]));
		refindex.push_back(tmpindex);

		refnum++;
	}

	ifstream readin("./reads/" + species);
	int x, num = 0;

	while(readin>>refseq)
	{
		readin>>info;	

		int len = refseq.length();

    	vector<int> align;
		for(int i = 0; i < len; i++)
		{
			readin>>x;
			align.push_back(x);
		}

		pseudo_match(refseq, align, info, seeding, revseeding);
		num++;
	}	
	
	printf("%d/%d/%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", k, w, w_min, w_max, seednum,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (num), ans[4] / (num), ans[5] / (num), 
		ans[6] / (num), ans[7] / (num), ans[8], ans[9]);

    return 0;
}
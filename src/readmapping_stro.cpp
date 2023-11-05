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
double matchdata[1100][3] = {0};
vector<vector<bool>> scover, mcover;
vector<vector<bool>> fscover;

vector<seedmatch> matches;

void pseudo_match(string s, vector<int> &align, string name, strobemerseeding seeding, int readnum)
{
	int lens = s.length();

	vector<vector<seed>> seeds(1);

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
		start = clock();

		matches.clear();
		index_get(refindex[j][0], seeds[0], matches);

		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
			
		int truematches = 0;
		for(seedmatch m: matches)
		{			
			if(m.s1->str != m.s2->str)
			{
				totalmatches--;
				continue;
			}
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
					mcover[readnum][i] = 1;
				
				for(int i = 0; i < k; i++)
				{
					scover[readnum][m.s2->st + i] = 1;
					scover[readnum][m.s2->st + index2 + i] = 1;
				}
			}
			else
			{	
				for(int i = 0; i < k; i++)
				{
					fscover[readnum][m.s2->st + i] = 1;
					fscover[readnum][m.s2->st + index2 + i] = 1;
				}
			}
		}

		matchdata[readnum][2] += double(seeds.size())/s.length();
		totalmatches += matches.size();
		totaltruematches += truematches;
	}
	else
	{	
		start = clock();
		matches.clear();
		index_get(refindex[j][0], seeds[0], matches);
		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

		for(seedmatch m: matches)
		{			
			if(m.s1->str != m.s2->str)
			{
				totalmatches--;
				continue;
			}

			index2 = m.s2->index % (1<<10);

			for(int i = 0; i < k; i++)
			{
				fscover[readnum][m.s2->st + i] = 1;
				fscover[readnum][m.s2->st + index2 + i] = 1;
			}
		}

		totalmatches += matches.size();
	}

	reverse(s.begin(), s.end());
	for(int i = 0; i < lens; i++)
		s[i] = ALPHABET[3-alphabetIndex(s[i])];

	seeds[0].clear();

    start = clock();
	seeding.get_strobemers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    
    for(int j = 0; j < refnum; j++)
    if(refname[j] == name)
    {
		start = clock();

		matches.clear();
		index_get(refindex[j][0], seeds[0], matches);

		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
			
		int truematches = 0;
		for(seedmatch m: matches)
		{
			int tp = 0;

			if(m.s1->str != m.s2->str)
			{
				totalmatches--;
				continue;
			}

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
					mcover[readnum][i] = 1;
				
				for(int i = 0; i < k; i++)
				{
					scover[readnum][lens - m.s2->st - i - 1] = 1;
					scover[readnum][lens - m.s2->st - index2 - i - 1] = 1;
				}
			}
			else
			{	
				for(int i = 0; i < k; i++)
				{
					fscover[readnum][lens - m.s2->st - i - 1] = 1;
					fscover[readnum][lens - m.s2->st - index2 - i - 1] = 1;
				}
			}
		}

		matchdata[readnum][2] += double(seeds.size())/s.length();
		totalmatches += matches.size();
		totaltruematches += truematches;
	}
	else
	{	
		start = clock();
		matches.clear();
		index_get(refindex[j][0], seeds[0], matches);
		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

		for(seedmatch m: matches)
		{			
			if(m.s1->str != m.s2->str)
			{
				totalmatches--;
				continue;
			}
			index2 = m.s2->index % (1<<10);

			for(int i = 0; i < k; i++)
			{
				fscover[readnum][lens - m.s2->st - i - 1] = 1;
				fscover[readnum][lens - m.s2->st - index2 - i - 1] = 1;
			}
		}

		totalmatches += matches.size();
	}

	matchdata[readnum][0] += totalmatches;
	matchdata[readnum][1] += totaltruematches;
}

void stats(int readnum, int lent)
{
	int mc = 0, sc = 0;
	int fsc = 0;
	double island = 0;
	int gap = 0;

	for(int i = 0; i < lent; i++)
	{
		sc += scover[readnum][i];
		fsc += fscover[readnum][i];
		mc += mcover[readnum][i];

		if(!scover[readnum][i])
			gap++;
		else if(gap > 0)
		{
			island += gap * gap;
			gap = 0;
		}
	}

	if(gap > 0)
		island += gap * gap;

	ans[0] += matchdata[readnum][0];
	ans[1] += matchdata[readnum][1];
	if(matchdata[readnum][0] > 0)
	{
		ans[10] += 1;
		ans[2] += double(matchdata[readnum][1])/matchdata[readnum][0];
	}
	ans[3] += double(sc) / lent;
	ans[4] += double(fsc) / lent;
	ans[5] += double(mc) / lent;
	ans[6] += double(island) / lent;
	ans[7] += matchdata[readnum][2];
}

vector<vector<vector<seed>>> seedu;
vector<int> readlen;

int main(int argc, const char * argv[])
{    
    k = atoi(argv[1]);
    w = atoi(argv[2]);
    w_min = atoi(argv[3]);
    w_max = atoi(argv[4]);
    seednum = atoi(argv[5]);

    int num = 0;

    for(int r = 0; r < seednum; r++)
    {
		strobemerseeding seeding(k, w, w_min, w_max, 1);

		for(auto a: refindex)
			for(auto b: a)
				delete b;
			
		refindex.clear();
		seedu.clear();

		string species = argv[6];

		ifstream refin("./ref/" + species);
		string info, refseq;
		refnum = 0;

		while(getline(refin, info))
		{
			getline(refin, refseq);

			vector<vector<seed>> seeds(1);
			seeding.get_strobemers(refseq, seeds);
			seedu.push_back(seeds);

			info = info.substr(1, info.find(' ') - 1);


			vector<ssh_index*> tmpindex;

			tmpindex.push_back(index_build(seedu[refnum][0]));
			refindex.push_back(tmpindex);

			refnum++;

			if(r == 0)
				refname.push_back(info);
		}
		refin.close();

		ifstream readin("./reads/" + species);
		int x, numnum = 0;

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

			if(r == 0)
			{
				readlen.push_back(len);
				scover.push_back(vector<bool>(len, 0));
				mcover.push_back(vector<bool>(len, 0));
				fscover.push_back(vector<bool>(len, 0));
				num++;
			}

			pseudo_match(refseq, align, info, seeding, numnum);
			numnum++;
		}	
		readin.close();
	}
	
	for(int i = 0; i < num; i++)
		stats(i, readlen[i]);

	printf("%d/%d/%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", k, w, w_min, w_max, seednum,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (num), ans[4] / (num), ans[5] / (num), 
		ans[6] / (num), ans[7] / (num), ans[8], ans[9]);

    return 0;
}
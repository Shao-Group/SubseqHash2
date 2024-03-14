#include "strobemerseeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, w, w_min, w_max, seednum;
ssh_index* refindex;

double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 
vector<vector<bool>> bcover;

vector<seedmatch> matches;
double matchdata[5100][2] = {0};

void pseudo_match(string s, vector<int> &align, strobemerseeding seeding, int readnum)
{
	int lens = s.length();

	vector<vector<seed>> seeds(1);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

    clock_t start,end;
	seeding.get_strobemers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;

    int index1, index2;

	start = clock();
	matches.clear();
	index_get(refindex, seeds[0], matches);
	end = clock();

	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
		
	totalmatches += matches.size();

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
			totaltruematches++;	
			bcover[readnum][m.s2->st/200] = 1;
		}
	}

	reverse(s.begin(), s.end());
	for(int i = 0; i < lens; i++)
		s[i] = ALPHABET[3-alphabetIndex(s[i])];

	seeds[0].clear();

	seeding.get_strobemers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;

	start = clock();
	matches.clear();
	index_get(refindex, seeds[0], matches);
	end = clock();

	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
		
	totalmatches += matches.size();

	for(seedmatch m: matches)
	{			
		if(m.s1->str != m.s2->str)
		{
			totalmatches--;
			continue;
		}
		int tp = 0;

		for(int i = 0; i < k; i++)
			tp += (align[lens - m.s2->st - i - 1] == m.s1->st + i);
		
		index1 = m.s1->index % (1<<10);
		index2 = m.s2->index % (1<<10);

		for(int i = 0; i < k; i++)
			tp += (align[lens - m.s2->st - i - 1] == m.s1->st + index1 + i);

		if(tp >= k)
		{
			totaltruematches++;	
			bcover[readnum][(lens - m.s2->ed - 1)/200] = 1;
		}
	}

	matchdata[readnum][0] += totalmatches;
	matchdata[readnum][1] += totaltruematches;
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


vector<vector<seed>> seedu(1);
vector<int> reads;

int main(int argc, const char * argv[])
{    
    k = atoi(argv[1]);
    w = atoi(argv[2]);
    w_min = atoi(argv[3]);
    w_max = atoi(argv[4]);
    seednum = atoi(argv[5]);

    int num = 0;

    string refseq;

    for(int r = 0; r < seednum; r++)
    {
		strobemerseeding seeding(k, w, w_min, w_max, 1);

		delete refindex;
		seedu[0].clear();

		ifstream refin(argv[6]);

		getline(refin, refseq);

		seeding.get_strobemers(refseq, seedu);
		refindex = index_build(seedu[0]);

		refin.close();

		ifstream readin(argv[7]);
		int x, numnum = 0;

		while(readin>>refseq)
		{
			int len = refseq.length();

	    	vector<int> align;
			for(int i = 0; i < len; i++)
			{
				readin>>x;
				align.push_back(x);
			}

			if(r == 0)
			{
				reads.push_back(len);
				bcover.push_back(vector<bool>(len/200 + 1, 0));
				num++;
			}

			pseudo_match(refseq, align, seeding, numnum);
			numnum++;
		}	
		readin.close();
	}

	for(int i = 0; i < num; i++)
		stats(i, reads[i]);
	printf("%d/%d/%d/%d/%d, %d, %d, %.4lf, %.4lf, %d, %d, %.4lf\n", k, w, w_min, w_max, seednum, int(ans[0]), int(ans[1]), ans[1]/ans[0], ans[2] / ans[10], int(ans[3]), int(ans[4]), ans[4]/ans[3]);

    return 0;
}
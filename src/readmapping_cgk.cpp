#include "cgkseeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, w;
ssh_index* refindex;

double ans[20] = {0}; 
//Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 
vector<vector<bool>> bcover;

vector<seedmatch> matches;
double matchdata[5100][2] = {0};

void pseudo_match(string s, vector<int> &align, cgkseeding seeding, int no)
{
	int lens = s.length();

	vector<seed> seeds;


	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

    clock_t start,end;

    start = clock();
	seeding.get_cgk(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;

	start = clock();
	matches.clear();
	index_get(refindex, seeds, matches);
	end = clock();
	
	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
	

	totalmatches += matches.size();

	for(seedmatch m: matches)
	{
		int tp = 0;

		for(int i = 0; i < k; i++)
			tp += (align[m.s2->st + i] == m.s1->st + i);

		if(2 * tp >= k)
		{
			totaltruematches++;	
			bcover[no][m.s2->st/200] = 1;
		}

	}


	reverse(s.begin(), s.end());
	for(int i = 0; i < lens; i++)
		s[i] = ALPHABET[3-alphabetIndex(s[i])];


	seeds.clear();

    start = clock();
	seeding.get_cgk(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;

	start = clock();
	matches.clear();
	index_get(refindex, seeds, matches);
	end = clock();
	
	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

	totalmatches += matches.size();

	for(seedmatch m: matches)
	{
		int tp = 0;

		for(int i = 0; i < k; i++)
			tp += (align[lens - m.s2->ed - 1 + i] == m.s1->ed - i);

		if(2 * tp >= k)
		{
			totaltruematches++;		
			bcover[no][(lens - m.s2->ed - 1)/200] = 1;
		}
	}

	matchdata[no][0] += totalmatches;
	matchdata[no][1] += totaltruematches;
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

vector<seed> seedu;
vector<int> reads;

int main(int argc, const char * argv[])
{    
    k = atoi(argv[1]);
	int repeat = atoi(argv[2]);

	ifstream refin(argv[3]);
	string refseq;
	string line;

	getline(refin, refseq);
	getline(refin, refseq);
	int readnum = 0;

    for(int i = 1; i <= repeat; i++) 
	{
    	cgkseeding seeding(k);
		seedu.clear();
		seeding.get_cgk(refseq, seedu);
		refindex = index_build(seedu);

		ifstream readin(argv[4]);
		int x, num = 0;

		while(readin>>line)
		{
			int len = line.length();
			
			vector<int> align;
			for(int i = 0; i < len; i++)
			{
				readin>>x;
				align.push_back(x);
			}

			if(i == 1)
			{
				reads.push_back(len);
				bcover.push_back(vector<bool>(len/200 + 1, 0));
				readnum++;
			}

			pseudo_match(line, align, seeding, num);

			num++;
		}	
	}

	for(int i = 0; i < readnum; i++)
		stats(i, reads[i]);
	printf("%d/%d, %d, %d, %.4lf, %.4lf, %d, %d, %.4lf\n", k, repeat, int(ans[0]), int(ans[1]), ans[1]/ans[0], ans[2] / ans[10], int(ans[3]), int(ans[4]), ans[4]/ans[3]);
    return 0;
}
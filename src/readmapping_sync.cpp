#include "syncmerseeding.hpp"
#include "util.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, s;
ssh_index* refindex;

double ans[20] = {0}; 
//Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 

vector<seedmatch> matches;

vector<bool> bcover(10000, 0);

void pseudo_match(string s, vector<int> &align, syncmerseeding seeding, int no)
{
	int lens = s.length();

	vector<seed> seeds;

    fill(bcover.begin(), bcover.end(), false);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

    clock_t start,end;

    start = clock();
	seeding.get_syncmers(s, seeds);
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
			bcover[m.s2->st/200] = 1;
		}

	}


	reverse(s.begin(), s.end());
	for(int i = 0; i < lens; i++)
		s[i] = ALPHABET[3-alphabetIndex(s[i])];


	seeds.clear();

    start = clock();
	seeding.get_syncmers(s, seeds);
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
			bcover[(lens - m.s2->ed - 1)/200] = 1;
		}
	}

	ans[0] += totalmatches;
	ans[1] += totaltruematches;
	if(totalmatches > 0)
	{
		ans[10] += 1;
		ans[2] += double(totaltruematches)/totalmatches;
	}

	int nblock = lens/200;
	int bc = 0;
	for(int i = 0; i < nblock; i++)
		bc += bcover[i];

	ans[3] += nblock;
	ans[4] += bc;
}

vector<seed> seedu;

int main(int argc, const char * argv[])
{    
    k = atoi(argv[1]);
    s = atoi(argv[2]);

	syncmerseeding seeding(k, s);

	for(int i = 5; i < argc; i++)
		seeding.add(atoi(argv[i]));


	ifstream refin(argv[3]);
	string refseq;

	getline(refin, refseq);
	getline(refin, refseq);

	seeding.get_syncmers(refseq, seedu);
	refindex = index_build(seedu);

	ifstream readin(argv[4]);
	int x, num = 0;

	while(readin>>refseq)
	{
		int len = refseq.length();

    	vector<int> align;
		for(int i = 0; i < len; i++)
		{
			readin>>x;
			align.push_back(x);
		}

		pseudo_match(refseq, align, seeding, num);

		num++;
	}	
	printf("%d/%d, %d, %d, %.4lf, %.4lf, %d, %d, %.4lf\n", k, s, int(ans[0]), int(ans[1]), ans[1]/ans[0], ans[2] / ans[10], int(ans[3]), int(ans[4]), ans[4]/ans[3]);
    return 0;
}
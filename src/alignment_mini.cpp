#include "minimizerseeding.hpp"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, w;

double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 
const int maxlen = 1<<27;

vector<seedmatch> matches;

vector<bool> scover(maxlen, 0);
vector<bool> scover2(maxlen, 0);
vector<bool> mcover(maxlen, 0);
vector<bool> mcover2(maxlen, 0);
vector<bool> fscover(maxlen, 0);
vector<bool> fscover2(maxlen, 0);

void pseudo_match(string s, string t, vector<int> &align, minimizerseeding seeding)
{
	int lens = s.length();
	int lent = t.length();

	vector<seed> seeds;
	vector<seed> seedt;

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
	seeding.get_minimizers(s, seeds);
	seeding.get_minimizers(t, seedt);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    	
	start = clock();
	ssh_index* ht = index_build(seeds);
	matches.clear();

	index_get(ht, seedt, matches);
	end = clock();
	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
		
	int truematches = 0;
	for(seedmatch m: matches)
	{
		int tp = 0;

		for(int i = 0; i < k; i++)
			tp += (align[m.s1->st + i] == m.s2->st + i);

		if(2 * tp >= k)
		{
			truematches++;	

			for(int i = m.s1->st; i <= m.s1->ed; i++)
				scover[i] = mcover[i] = 1;
			for(int i = m.s2->st; i <= m.s2->ed; i++)
				scover2[i] = mcover2[i] = 1;
		}
		else
		{	
			for(int i = m.s1->st; i <= m.s1->ed; i++)
				fscover[i] = 1;
			for(int i = m.s2->st; i <= m.s2->ed; i++)
				fscover2[i] = 1;
		}
	}

	ans[7] += double(seeds.size())/s.length() + double(seedt.size())/t.length();
	totalmatches += matches.size();
	totaltruematches += truematches;

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
    if(argc != 4)
    {
		printf("usage: alignment_mini.out readFile k w\n");
		return 1;
    }

    k = atoi(argv[2]);
    w = atoi(argv[3]);

	minimizerseeding seeding(k, w);

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
	printf("%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", k+w-1, k, w,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (2*num), ans[4] / (2*num), ans[5] / (2*num), 
		ans[6] / (2*num), ans[7] / (2*num), ans[8], ans[9]);

    // const char* table_filename = argv[5];

    return 0;
}
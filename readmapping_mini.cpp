#include "minimizerseeding.h"
#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int k, w;
vector<string> refname;
vector<ssh_index*> refindex;

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

void pseudo_match(string s, vector<int> &align, string name, minimizerseeding seeding)
{
	int lens = s.length();

	vector<seed> seeds;

    fill(scover.begin(), scover.end(), false);
    fill(mcover.begin(), mcover.end(), false);
    fill(fscover.begin(), fscover.end(), false);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

    clock_t start,end;

    start = clock();
	seeding.get_minimizers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    	
    for(int j = 0; j < refnum; j++)
   	if(refname[j] == name)
    {
		start = clock();

		matches.clear();
		index_get(refindex[j], seeds, matches);

		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
			
		int truematches = 0;
		for(seedmatch m: matches)
		{
			int tp = 0;

			for(int i = 0; i < k; i++)
				tp += (align[m.s2->st + i] == m.s1->st + i);

			if(2 * tp >= k)
			{
				truematches++;	

				for(int i = m.s2->st; i <= m.s2->ed; i++)
					scover[i] = mcover[i] = 1;
			}
			else
			{	
				for(int i = m.s2->st; i <= m.s2->ed; i++)
					fscover[i] = 1;
			}
		}

		ans[7] += double(seeds.size())/s.length();
		totalmatches += matches.size();
		totaltruematches += truematches;
	}
	else
	{		
		start = clock();
		matches.clear();
		index_get(refindex[j], seeds, matches);
		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

		for(seedmatch m: matches)
		{
			for(int i = m.s2->st; i <= m.s2->ed; i++)
				fscover[i] = 1;
		}

		totalmatches += matches.size();
	}

	reverse(s.begin(), s.end());
	for(int i = 0; i < lens; i++)
		s[i] = ALPHABET[3-alphabetIndex(s[i])];

	seeds.clear();

    start = clock();
	seeding.get_minimizers(s, seeds);
    end = clock();
    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;
    	
    for(int j = 0; j < refnum; j++)
    if(refname[j] == name)
    {
		start = clock();

		matches.clear();
		index_get(refindex[j], seeds, matches);

		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;
			
		int truematches = 0;
		for(seedmatch m: matches)
		{
			int tp = 0;

			for(int i = 0; i < k; i++)
				tp += (align[lens - m.s2->ed - 1 + i] == m.s1->ed - i);

			if(2 * tp >= k)
			{
				truematches++;	

				for(int i = lens - m.s2->ed - 1; i <= lens - m.s2->st - 1; i++)
					scover[i] = mcover[i] = 1;
			}
			else
			{	
				for(int i = lens - m.s2->ed - 1; i <= lens - m.s2->st - 1; i++)
					fscover[i] = 1;
			}
		}

		ans[7] += double(seeds.size())/s.length();
		totalmatches += matches.size();
		totaltruematches += truematches;
	}
	else
	{		
		start = clock();
		matches.clear();
		index_get(refindex[j], seeds, matches);
		end = clock();
		ans[9] += (double)(end-start)/CLOCKS_PER_SEC;

		for(seedmatch m: matches)
		{
			for(int i = lens - m.s2->ed - 1; i <= lens - m.s2->st - 1; i++)
				fscover[i] = 1;
		}

		totalmatches += matches.size();
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

	//cout<<totalmatches<<" "<<totaltruematches<<" "<<sc<<" "<<fsc<<endl;
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

vector<vector<seed>> seedu;

int main(int argc, const char * argv[])
{    
    k = atoi(argv[1]);
    w = atoi(argv[2]);

	minimizerseeding seeding(k, w);

	string species = argv[3];

	ifstream refin("./ref/" + species);
	string info, refseq;
	refnum = 0;

	while(getline(refin, info))
	{
		getline(refin, refseq);

		vector<seed> seeds;
		seedu.push_back(seeds);

		seeding.get_minimizers(refseq, seedu[refnum]);

		info = info.substr(1, info.find(' ') - 1);

		refname.push_back(info);
		refindex.push_back(index_build(seedu[refnum]));

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

		pseudo_match(refseq, align, info, seeding);
		num++;
	}	
	printf("%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", k+w-1, k, w,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (num), ans[4] / (num), ans[5] / (num), 
		ans[6] / (num), ans[7] / (num), ans[8], ans[9]);

    return 0;
}
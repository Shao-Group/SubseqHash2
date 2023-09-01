#include "util.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n, k, d, subsample;
vector<string> refname;
ssh_index* refindex[20];

int refnum = 0;
double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 

double matchdata[1100][3] = {0};
vector<seedmatch> matches;

vector<vector<bool>> scover, mcover;
vector<vector<bool>> fscover;

void pseudo_match(int len, int readnum, int refno, vector<seed> &seedt, vector<int> &align)
{
	matches.clear();

	int totalmatches = 0;


    clock_t start,end;


	index_get(refindex[refno], seedt, matches);

	uint64_t x1;

	int truematches = 0;

	totalmatches += matches.size();
	for(seedmatch m: matches)
	{
		int tp = 0;

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
			for(int i = 0; i < n; i++)
				if((m.s2->index>>i) & 1)
					scover[readnum][m.s2->st + i] = 1;

			for(int i = m.s2->st; i <= m.s2->ed; i++)
				mcover[readnum][i] = 1;
		}
		else
		{			
			for(int i = 0; i < n; i++)
				if((m.s2->index>>i) & 1)
					fscover[readnum][m.s2->st + i] = 1;
		}
	}

	matches.clear(); //match to other chromosome
	for(int j = 0; j < refnum; j++)
		if(j != refno)
			index_get(refindex[j], seedt, matches);

	totalmatches += matches.size();
	for(seedmatch m: matches)
	{
		for(int i = 0; i < n; i++)
			if((m.s2->index>>i) & 1)
				fscover[readnum][m.s2->st + i] = 1;
	}


	matchdata[readnum][0] += totalmatches;
	matchdata[readnum][1] += truematches;
	matchdata[readnum][2] += double(seedt.size()) / len;
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

vector<int> reads;
vector<vector<int>> align;
vector<string> trueref;
vector<vector<seed>> seeds(20);

int main(int argc, const char * argv[])
{    
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    d = atoi(argv[3]);
    int subsample = atoi(argv[4]);

	string species = argv[5];

	ifstream refin("./ref/" + species + "list");
	string refseq, info, t;
	refnum = 0;

	while(getline(refin, info))
	{
		refname.push_back(info);
		refnum++;
	}

	ifstream readin("./reads/" + species);
	int x, num = 0;

	while(readin>>refseq)
	{
		vector<int> aligninput;
		readin>>t;

		int len = refseq.length();
		for(int j = 0; j < len; j++)
		{
			readin>>x;
			aligninput.push_back(x);
		}

		reads.push_back(len);
		align.push_back(aligninput);
		trueref.push_back(t);

		scover.push_back(vector<bool>(len, 0));
		mcover.push_back(vector<bool>(len, 0));
		fscover.push_back(vector<bool>(len, 0));

		num++;
	}	

	for(int i = 0; i < subsample; i += 2)
	{
		for(int j = 0; j < refnum; j++)
		{
			string refpath = "./refsub2/" + species + "/" + refname[j] + "_" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);

			seeds[j].clear();
			loadSeeds(refpath.c_str(), seeds[j]);
		}
		if(i + 1 < subsample)
			for(int j = 0; j < refnum; j++)
			{
				//cout<<i<<" "<<j<<endl;
				string refpath = "./refsub2/" + species + "/" + refname[j] + "_" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i+1);

				loadSeeds(refpath.c_str(), seeds[j]);
			}

		for(int j = 0; j < refnum; j++)
		{
			delete refindex[j];
			refindex[j] = index_build(seeds[j]);
		}

		string qpath = "./readsub2/" + species + "/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i);
		string qpath_rc = "./readsub2/" + species + "/" + to_string(n) + "_" + to_string(k) + "_" + to_string(d) + "_" + to_string(i+1);
	    FILE* fin = fopen(qpath.c_str(), "rb"), *fin_rc;

	    if(i + 1 < subsample)
		    fin_rc = fopen(qpath_rc.c_str(), "rb");	

	    uint64_t st, index;
	    int64_t hashval;
		
		for(int j = 0; j < num; j++)
		{
		    vector<seed> seedt;
		   	size_t ret = 1;

		   	//cout<<i<<" "<<j<<endl;
		    while(ret == 1)
		    {
				ret = fread(&hashval, sizeof(int64_t), 1, fin);
				fread(&st, sizeof(uint64_t), 1, fin);
				fread(&index, sizeof(uint64_t), 1, fin);

				if(hashval == 0 && st == 0 && index == 0)
					break;
				seed tmp;
				tmp.hashval = hashval;
				tmp.st = st;
				tmp.index = index;		

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
				
				//cout<<tmp.hashval<<" "<<tmp.st<<" "<<tmp.ed<<" "<<tmp.index<<endl;
				seedt.push_back(tmp);
			}				    
	
	    	if(i + 1 < subsample)
			    while(ret == 1)
			    {
					ret = fread(&hashval, sizeof(int64_t), 1, fin_rc);
					fread(&st, sizeof(uint64_t), 1, fin_rc);
					fread(&index, sizeof(uint64_t), 1, fin_rc);

					if(hashval == 0 && st == 0 && index == 0)
						break;
					seed tmp;
					tmp.hashval = hashval;
					tmp.st = st;
					tmp.index = index;		

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
					
					//cout<<tmp.hashval<<" "<<tmp.st<<" "<<tmp.ed<<" "<<tmp.index<<endl;
					seedt.push_back(tmp);
				}	

			for(int l = 0; l < refnum; l++)
				if(refname[l] == trueref[j])
				{
					pseudo_match(reads[j], j, l, seedt, align[j]);
					break;
				}
		}
	}

	for(int i = 0; i < num; i++)
		stats(i, reads[i]);

	printf("%d/%d/%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.4lf, %.2lf, %.4lf, %.2lf, %.2lf\n", n, k, d, subsample,
		ans[0] / num, ans[1] / num, ans[2] / ans[10], ans[3] / (num), ans[4] / (num), ans[5] / (num), 
		ans[6] / (num), ans[7] / (num), ans[8], ans[9]);

    return 0;
}
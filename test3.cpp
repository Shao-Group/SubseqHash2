#include "subseqhash2seeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

double ans[20] = {0}; //Number of matches, Number of true-matches, precision of seed-matches, sequence cover, false sequence cover, matching cover, island(sc), density
//seeding time, seed-match time 

void pseudo_match(string s, string t, vector<int> &align, subseqhash2seeding & sub2)
{
	int lens = s.length();
	int lent = t.length();

	vector<seedmatch> matches;

	int chunk_size = sub2.getChunkSize();
	int seednum = sub2.getNumPerWindow();

	vector<vector<seed>> seeds(seednum, vector<seed>(0));
	vector<vector<seed>> seedt(seednum, vector<seed>(0));

	DPCell* dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
	DPCell* revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
	int* h = (int*) malloc(sizeof *h * dim1);
	int* revh = (int*) malloc(sizeof *revh * dim1);

    clock_t start,end;

    start = clock();
	sub2.getSubseq2Seeds(s, dp, revdp, h, revh, seeds);
	sub2.getSubseq2Seeds(t, dp, revdp, h, revh, seedt);
    end = clock();

    ans[8] += (double)(end-start)/CLOCKS_PER_SEC;

	vector<bool> scover(lens, 0);
	vector<bool> scover2(lent, 0);
	vector<bool> mcover(lens, 0);
	vector<bool> mcover2(lent, 0);
	vector<bool> fscover(lens, 0);
	vector<bool> fscover2(lent, 0);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

	for(int j = 0; j < seednum; j++)
	{	
    	start = clock();
		ssh_index* ht = index_build(seeds[j]);
		matches.clear();

		index_get(ht, seedt[j], matches);
    	end = clock();
    	ans[9] += (double)(end-start)/CLOCKS_PER_SEC;


		int truematches = 0;
		for(seedmatch m: matches)
		{
			int tp = 0;

			for(int i = 0; i < n; i++)
				if((m.s1->index>>i) & 1)
				{
					x1 = align[m.s1->st + i]; 
					if(x1 >= m.s2->st && x1 <= m.s2->ed)
						tp += (m.s2->index>>(x1-m.s2->st)) & 1;
				}

			if(2 * tp >= k)
			{
				truematches++;	
				for(int i = 0; i < n; i++)
					if((m.s1->index>>i) & 1)
						scover[m.s1->st + i] = 1;
				for(int i = 0; i < n; i++)
					if((m.s2->index>>i) & 1)
						scover2[m.s2->st + i] = 1;

				for(int i = m.s1->st; i <= m.s1->ed; i++)
					mcover[i] = 1;
				for(int i = m.s2->st; i <= m.s2->ed; i++)
					mcover2[i] = 1;
			}
			else
			{			
				for(int i = 0; i < n; i++)
					if((m.s1->index>>i) & 1)
						fscover[m.s1->st + i] = 1;
				for(int i = 0; i < n; i++)
					if((m.s2->index>>i) & 1)
						fscover2[m.s2->st + i] = 1;
			}
		}

		ans[7] += double(seeds[j].size())/s.length() + double(seedt[j].size())/t.length();
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
    if(argc != 7)
    {
		printf("usage: genSubseq2Seeds.out readFile n k d subsample randTableFile\n");
		return 1;
    }

    n = atoi(argv[2]);
    k = atoi(argv[3]);
    d = atoi(argv[4]);
    int subsample = atoi(argv[5]);
    dim1 = (n+1) * (k+1) * d;

    subseqhash2seeding sub2(n, k, d, subsample);
    sub2.init(argv[6]);

	ifstream fin(argv[1]);

	string seq, seq2, tmp;
	int x;
	int num = 0;
	char* str;

	while(fin>>tmp)
	{	
		fin>>tmp;
		fin>>seq;
		fin>>tmp;
		fin>>tmp;

		fin>>seq2;				

		int len = seq.length();
		cout<<seq.length()<<" "<<seq2.length()<<endl;

		int chunk_size = sub2.getChunkSize();
		int seednum = sub2.getNumPerWindow();

		vector<vector<seed>> seeds(seednum, vector<seed>(0));
		vector<vector<seed>> seedt(seednum, vector<seed>(0));

		DPCell* dp = (DPCell*) malloc(sizeof *dp * chunk_size * (dim1));
		DPCell* revdp = (DPCell*) malloc(sizeof *revdp * chunk_size * dim1);
		int* h = (int*) malloc(sizeof *h * dim1);
		int* revh = (int*) malloc(sizeof *revh * dim1);

		sub2.getSubseq2Seeds(seq, dp, revdp, h, revh, seeds);
		sub2.getSubseq2Seeds(seq2, dp, revdp, h, revh, seedt);

		vector<seedmatch> matches;

		for(int j = 0; j < seednum ; j++)
		{	
			ssh_index* ht = index_build(seeds[j]);	
			int sz = seeds[j].size();	
			int sz1 = seedt[j^1].size();		
			matches.clear();

			index_get(ht, seedt[j^1], matches);

			cout<<sz<<" "<<sz1<<" "<<matches.size()<<endl;

			for(int z = 0; z < sz; z++)
				if(seeds[j][z].hashval != seedt[j^1][sz1 -z - 1].hashval)
				{
					for(int y = z-2; y <= z+2; y++)
					{
						cout<<y<<" "<<seeds[j][y].hashval<<" "<<seeds[j][y].st<<" "<<seeds[j][y].ed<<endl;
						cout<<seedt[j^1][sz1 -y - 1].hashval<<" "<<seedt[j^1][sz1 -y - 1].st<<" "<<seedt[j^1][sz1 -y - 1].ed<<endl;
						cout<<seq.substr(seeds[j][y].st, seeds[j][y].ed- seeds[j][y].st + 1)<<endl;
						cout<<seq2.substr(seedt[j^1][sz1 -y - 1].st, seedt[j^1][sz1 -y - 1].ed- seedt[j^1][sz1 -y - 1].st + 1)<<endl;

						printf("%.48s\n", decode(seeds[j][y].str, k, str));
						printf("%.48s\n", decode(seedt[j^1][sz1 -y - 1].str, k, str));
					}
					break;
				}	
		}	
	}

    // const char* table_filename = argv[5];

    return 0;
}
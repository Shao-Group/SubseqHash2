#include "subseqhash2seeding.h"
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int n,k,d, dim1;

double ans[20] = {0};
double ans10[20] = {0};
double subans[50][20] = {0};

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

    ans[4] += (double)(end-start)/CLOCKS_PER_SEC;

	vector<bool> scover(lens, 0);
	vector<bool> scover2(lent, 0);
	vector<bool> fscover(lens, 0);
	vector<bool> fscover2(lent, 0);

	int totalmatches = 0;
	int totaltruematches = 0;

	uint64_t x1;

	for(int j = 0; j < seednum; j++)
	{	
		ssh_index* ht = index_build(seeds[j]);
		matches.clear();

		index_get(ht, seedt[j], matches);

		// cout<<seeds[j].size()<<" "<<seedt[j].size()<<" "<<matches.size()<<endl;

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

		ans[3] += double(seeds[j].size())/s.length() + double(seedt[j].size())/t.length();
		totalmatches += matches.size();
		totaltruematches += truematches;
	}		

	int sc = 0, sc2 = 0;
	int fsc = 0, fsc2 = 0;

	for(int i = 0; i < lens; i++)
	{
		sc += scover[i];
		fsc += fscover[i];
	}

	for(int i = 0; i < lent; i++)
	{
		sc2 += scover2[i];
		fsc2 += fscover[i];
	}

	ans[0] += totalmatches;
	ans[1] += totaltruematches;
	if(totalmatches > 0)
	{
		ans[7] += 1;
		ans[2] += double(totaltruematches)/totalmatches;
	}
	ans[5] += double(sc2) / lent + double(sc) / lens;
	ans[9] += double(fsc2) / lent + double(fsc) / lens;
}

int main(int argc, const char * argv[])
{    
    if(argc != 6)
    {
		printf("usage: genSubseq2Seeds.out readFile n k d randTableFile\n");
		return 1;
    }

    n = atoi(argv[2]);
    k = atoi(argv[3]);
    d = atoi(argv[4]);
    dim1 = (n+1) * (k+1) * d;

    subseqhash2seeding sub2(n, k, d);
    sub2.init(argv[5]);

	ifstream fin(argv[1]);

	string seq, seq2, tmp;
	int x;
	int num = 0;

	while(fin>>seq)
	{
		fin>>seq2;				

		int len = seq.length();

    	vector<int> align;
		for(int i = 0; i < len; i++)
		{
			fin>>x;
			align.push_back(x);
		}

		pseudo_match(seq, seq2, align, sub2);
		num++;
	}	
	printf("%d/%d, %.2lf, %.2lf, %.4lf, %.4lf, %.4lf, %.2lf, %.2lf\n", n, k, ans[0] / num, ans[1] / num, ans[2] / ans[7], ans[5] / (2*num), ans[9] / (2*num), ans[4], ans[3] / (2*num));

    // const char* table_filename = argv[5];

    return 0;
}
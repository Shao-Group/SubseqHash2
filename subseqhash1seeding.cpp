#include "subseqhash1seeding.h"

const double INF=1e15;

//accessing DPCell[n+1][k+1][d]
inline int subseqhash1seeding::dpIndex(int d1, int d2, int d3)
{
    return d1 * dim1 + d2 * dim2 + d3;
}

void subseqhash1seeding::init(const char* table_filename)
{
    FILE *filein = fopen(table_filename, "r");
    
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%lf ", &A[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%d,%d", &B1[i][j][q], &B2[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C[i][j]);
}

void subseqhash1seeding::DP(std::string s, DPCell* dp, int* h, std::vector<seed>& seeds)
{
    int dp_index;

	int len = s.length();
    int del = n - k;

    memset(h, 0, sizeof *h * dim1 * (n+1));
	for(int i = 0; i <= n; i++)
	{
	    dp_index = dpIndex(i, 0, 0); //[i][0][0]

	    dp[dp_index].f_max = 0;
	    dp[dp_index].f_min = 0;
	    
	    dp_index += 1;
	    for(int j = 1; j < d; ++j, ++dp_index)
	    {
			dp[dp_index].f_max = -INF;
			dp[dp_index].f_min = INF;
	    }
	}

    for(int st = 0; st + n <= len; st++)
    {
		for(int i = 0; i <= n; i++)
		{
			dp_index = dpIndex(i, 0, 0);
		    h[dp_index] = st + 1;
		}
	
		int d1 = C[0][alphabetIndex(s[st])];
		dp_index = dpIndex(1, 1, d1);

		dp[dp_index].f_min = dp[dp_index].f_max = B2[0][alphabetIndex(s[st])][d1] * A[0][alphabetIndex(s[st])][d1];
		dp[dp_index].g_min = dp[dp_index].g_max = 0;
		h[dp_index] = st + 1;

		//printf("%d %d %d %.2lf %.2lf %d %d\n", st, d1, dp_index, dp[dp_index].f_max, dp[dp_index].f_min, dp[dp_index].g_max, dp[dp_index].g_min);
		double v1, v2;

		for(int i = 2; i <= n; i++)
		{
		    int minj = std::max(1, i - del);
		    int maxj = std::min(i, k);

		    for(int j = minj; j <= maxj; j++)
		    {
				int now = alphabetIndex(s[st + i - 1]); //Current char
				int v = C[j-1][now];

				for(int q = 0; q < d; q++)
				{
				    int z = (q + v) % d; // q: previous d, z: current d
				    
				    int idx = dpIndex(i, j, z);
				    
				    int laidx1 = dpIndex(i-1, j, z);// = index_cal(i-1,j,z);
				    int laidx2 = dpIndex(i-1, j-1, q);// = index_cal(i-1,j-1,q);

				    if(h[idx] != st + 1)
				    {
						dp[idx].f_min = INF;
						dp[idx].f_max = -INF;
				    }
				    
				    if(h[laidx1] == st + 1)
				    {
						if(dp[laidx1].f_min < dp[idx].f_min)
						{
						    dp[idx].f_min = dp[laidx1].f_min;
						    dp[idx].g_min = dp[laidx1].g_min;
						}
						if(dp[laidx1].f_max > dp[idx].f_max)
						{
						    dp[idx].f_max = dp[laidx1].f_max;
						    dp[idx].g_max = dp[laidx1].g_max;
						}
						h[idx] = st + 1;
				    }

				    if(h[laidx2] == st + 1)
				    {
						if(B1[j-1][now][z] == -1)
						{
						    v1 = -dp[laidx2].f_min;
						    v2 = -dp[laidx2].f_max;

						    if(B2[j-1][now][z] == -1)
						    {
								v1 -= A[j-1][now][z];
								v2 -= A[j-1][now][z];
						    }
						    else
						    {
								v1 += A[j-1][now][z];
								v2 += A[j-1][now][z];
						    }

						    if(v1 > dp[idx].f_max)
						    {
								dp[idx].f_max = v1;
								dp[idx].g_max = i-1;
						    }
						    if(v2 < dp[idx].f_min)
						    {
								dp[idx].f_min = v2;
								dp[idx].g_min = i-1;
						    }
						}

						else
						{
						    v1 = dp[laidx2].f_min;
						    v2 = dp[laidx2].f_max;

						    if(B2[j-1][now][z] == -1)
						    {
								v1 -= A[j-1][now][z];
								v2 -= A[j-1][now][z];
						    }
						    else
						    {
								v1 += A[j-1][now][z];
								v2 += A[j-1][now][z];
						    }

						    if(v1 < dp[idx].f_min)
						    {
								dp[idx].f_min = v1;
								dp[idx].g_min = i-1;
						    }
						    if(v2 > dp[idx].f_max)
						    {
								dp[idx].f_max = v2;
								dp[idx].g_max = i-1;
						    }
						}
						
						h[idx] = st + 1;
				    }			

				    //printf("%d %d %d %d %d %.2lf %.2lf %d %d\n", i,j,z,idx,h[idx],dp[idx].f_max, dp[idx].f_min, dp[idx].g_max, dp[idx].g_min);
				}
			}
		}

		seed tmp;
		int mod;
		dp_index = dpIndex(n, k, 0);

		for(int i = 0; i < d; i++)
			if(h[dp_index + i] == st + 1)
			{
				if(fabs(dp[dp_index + i].f_min) > dp[dp_index + i].f_max)
					tmp.hashval = uint64_t(dp[dp_index + i].f_min * 32768);
				else
					tmp.hashval = uint64_t(dp[dp_index + i].f_max * 32768);

				mod = i;
				break;
			}

		kmer hashval = 0;
	    std::vector<size_t> index;

	    int x = n;
	    int y = k;
	    size_t nextpos;
	    int nextval, nextd;

	    bool z = 1;
	    if(tmp.hashval < 0)
			z = 0;

	    while(y > 0)
	    {
			dp_index = dpIndex(x, y, mod);
			if(z == 1)
			{
			    nextpos = st + dp[dp_index].g_max;
			    nextval = alphabetIndex(s[nextpos]);

			    index.push_back(nextpos);
			    hashval = (hashval<<2) | nextval;

			    if(B1[y-1][nextval][mod] == -1)
					z = 0;

			    x = dp[dp_index].g_max;
			    y--;
			    mod = (mod + d - C[y][nextval]) % d;
			}
			else
			{
			    nextpos = st + dp[dp_index].g_min;
			    nextval = alphabetIndex(s[nextpos]);

			    index.push_back(nextpos);
			    hashval = (hashval<<2) | nextval;

			    if(B1[y-1][nextval][mod] == -1)
					z = 1;

			    x = dp[dp_index].g_min;
			    y--;
			    mod = (mod + d - C[y][nextval]) % d;
			}
		}

	    tmp.st = index.back();
	    tmp.ed = index[0];

	    for(size_t a: index)
	    	tmp.index |= ((uint64_t)1)<<(a - tmp.st);

	    if(seeds.size() > 0 && tmp.hashval == seeds.back().hashval)
	    {
	    	if(tmp.st <= seeds.back().st && tmp.ed >= seeds.back().ed)
				continue;
			else
				if(tmp.st >= seeds.back().st && tmp.ed <= seeds.back().ed)
				{
					seeds.back().st = tmp.st;
					seeds.back().index = tmp.index;
					continue;
				}
	    }

	    tmp.str = hashval;
	    tmp.str_rc = revComp(hashval, k);
	    seeds.push_back(tmp);
	}
} // end of DP



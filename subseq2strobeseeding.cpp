#include "subseq2strobeseeding.h"
const int64_t INF=(int64_t)1<<62;

//accessing DPCell[chunk_size][n+1][k+1][d]
inline int subseq2strobeseeding::dpIndex(int d1, int d2, int d3, int d4)
{
    return d1 * dim1 + d2 * dim2 + d3 * dim3 + d4;
}

//accessing int[n+1][k+1][d]
inline int subseq2strobeseeding::hIndex(int d2, int d3, int d4)
{
    return d2 * dim2 + d3 * dim3 + d4;
}

void subseq2strobeseeding::init(const char* table_filename)
{
    FILE *filein = fopen(table_filename, "r");

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%" SCNd64, &A[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%d,%d", &B1[i][j][q], &B2[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C1[i][j]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%" SCNd64, &revA[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	   		for(int q = 0; q < d; q++)
				fscanf(filein, "%d,%d", &revB1[i][j][q], &revB2[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C2[i][j]);
	
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d,%d", &combine1[i][j], &combine2[i][j]);
	

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%" SCNd64, &A3[i][j]);
	
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	   		fscanf(filein, "%d ", &combine3[i][j]);
	
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C3[i][j]);


    int x;
    int subsample = num_valid;

    if((num_valid&1) && (k&1))
    {
    	subsample--;
    	valid[k/2] = subsample+1;
    }

    for(int i = 0; i < k-2; i++)
    {
		fscanf(filein, "%d", &x);
		if(i < subsample)
			valid[x+1] = i+1;
    }
}

void subseq2strobeseeding::DP(std::string s, size_t start, size_t end, DPCell* dp, int* h)
{
    int dp_index, h_index;

    int len = end - start + 1;
    int del = n - k;

    memset(h, 0, sizeof *h * dim1);

    for(int st = 0; st < len; st++)
    {
		int st_base_index = st * dim1;
		size_t real_index_st = st + start;
		int en = std::min(st + n - 1, len - 1);
		int interval_len = en - st + 1;

		for(int i = 0; i <= interval_len; i++)
		{
		    h_index = hIndex(i, 0, 0); //[i][0][0]
		    dp_index = st_base_index + h_index; //[st][i][0][0]

		    dp[dp_index].f_max = 0;
		    dp[dp_index].f_min = 0;
		    h[h_index] = st + 1;
		    
		    dp_index += 1;
		    for(int j = 1; j < d; ++j, ++dp_index)
		    {
				dp[dp_index].f_max = -INF;
				dp[dp_index].f_min = INF;
		    }
		}

	
		int d1 = C1[0][alphabetIndex(s[real_index_st])];
		h_index = hIndex(1, 1, 0);
		dp_index = st_base_index + h_index; //[st][1][1][0]
		int idx1 = dp_index + d1;

		for(int i = 0; i < d; ++i, ++dp_index)
		{
			dp[dp_index].f_max = -INF;
			dp[dp_index].f_min = INF;
		}


		dp[idx1].f_min = dp[idx1].f_max = B2[0][alphabetIndex(s[real_index_st])][d1] * A[0][alphabetIndex(s[real_index_st])][d1];
		dp[idx1].g_min = dp[idx1].g_max = 0;
		h_index += d1; //[1][1][d1]
		h[h_index] = st + 1;

		int64_t v1, v2;

		for(int i = 2; i <= interval_len; i++)
		{
		    int minj = std::max(1, i - del);
		    int maxj = std::min(i, k - 2);

		    for(int j = minj; j <= maxj; j++)
		    {
				int now = alphabetIndex(s[real_index_st + i - 1]); //Current char
				int v = C1[j-1][now];

				for(int q = 0; q < d; q++)
				{
				    int z = (q + v) % d; // q: previous d, z: current d
				    
				    h_index = hIndex(i, j, z);
				    dp_index = st_base_index + h_index;
				    int idx = dp_index;
				    
				    int laidx1;// = index_cal(st,i-1,j,z);
				    int laidx2;// = index_cal(st,i-1,j-1,q);

				    if(h[h_index] != st + 1)
				    {
						dp[idx].f_min = INF;
						dp[idx].f_max = -INF;
				    }

				    int prev_h_index = h_index;
				    h_index = hIndex(i-1, j, z);
				    laidx1 = st_base_index + h_index; 
				    
				    if(h[h_index] == st + 1)
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
						h[prev_h_index] = st + 1;
				    }

				    h_index = hIndex(i-1, j-1, q);
				    laidx2 = st_base_index + h_index;

				    if(h[h_index] == st + 1)
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
						
						h[prev_h_index] = st + 1;
				    }					
				}
			}
		}
	}
} // end of DP


void subseq2strobeseeding::revDP(std::string s, size_t start, size_t end, DPCell* revdp, int* revh)
{
    int dp_index, h_index;
    
    int len = end - start + 1;
    int del = n - k;

    memset(revh, 0, sizeof *revh * dim1);

    for(int st = len - 1; st >= 0; st--)
    {
		int st_base_index = st * dim1;
		size_t real_index_st = st + start;	
		int en = std::max(st - n + 1, 0);
		int interval_len = st - en + 1;


		for(int i = 0; i <= interval_len; i++)
		{
		    h_index = hIndex(i, 0, 0);
		    dp_index = st_base_index + h_index;
		    
		    revdp[dp_index].f_max = 0;
		    revdp[dp_index].f_min = 0;
		    revh[h_index] = st + 1;

		    dp_index += 1;
		    for(int j = 1; j < d; ++j, ++dp_index)
		    {
				revdp[dp_index].f_max = -INF;
				revdp[dp_index].f_min = INF;
		    }
		}

		int d1 = C2[0][alphabetIndex(s[real_index_st])];
		h_index = hIndex(1, 1, 0);
		dp_index = st_base_index + h_index;
		int idx1 = dp_index + d1;

		for(int i = 0; i < d; ++i, ++ dp_index)
		{
		    revdp[dp_index].f_max = -INF;
		    revdp[dp_index].f_min = INF;
		}

		revdp[idx1].f_min = revdp[idx1].f_max = revB2[0][alphabetIndex(s[real_index_st])][d1] * revA[0][alphabetIndex(s[real_index_st])][d1];
		revdp[idx1].g_min = revdp[idx1].g_max = 0;
		h_index += d1;
		revh[h_index] = st + 1;

		int64_t v1, v2;

		for(int i = 2; i <= interval_len; i++)
		{
		    int minj = std::max(1, i - del);
		    int maxj = std::min(i, k - 2);

		    for(int j = minj; j <= maxj; j++)
		    {
				int now = alphabetIndex(s[real_index_st - i + 1]);
				int v = C2[j-1][now];

				for(int q = 0; q < d; q++)
				{
				    int z = (q + v) % d; // q: previous d, z: current d
				    
				    h_index = hIndex(i, j, z);
				    dp_index = st_base_index + h_index;
				    int idx = dp_index;
				    
				    int laidx1; // = index_cal(st,i-1,j,z);
				    int laidx2; // = index_cal(st,i-1,j-1,q);

				    if(revh[h_index] != st + 1)
				    {
						revdp[idx].f_min = INF;
						revdp[idx].f_max = -INF;
				    }

				    int prev_h_index = h_index;
				    h_index = hIndex(i-1, j, z);
				    laidx1 = st_base_index + h_index;

				    if(revh[h_index] == st + 1)
				    {
						if(revdp[laidx1].f_min < revdp[idx].f_min)
						{
						    revdp[idx].f_min = revdp[laidx1].f_min;
						    revdp[idx].g_min = revdp[laidx1].g_min;	
						}

						if(revdp[laidx1].f_max > revdp[idx].f_max)
						{
						    revdp[idx].f_max = revdp[laidx1].f_max;
						    revdp[idx].g_max = revdp[laidx1].g_max;
						}
						revh[prev_h_index] = st + 1;
				    }

				    h_index = hIndex(i-1, j-1, q);
				    laidx2 = st_base_index + h_index;

				    if(revh[h_index] == st + 1)
				    {
						if(revB1[j-1][now][z] == -1)
						{
						    v1 = -revdp[laidx2].f_min;
						    v2 = -revdp[laidx2].f_max;

						    if(revB2[j-1][now][z] == -1)
						    {
								v1 -= revA[j-1][now][z];
								v2 -= revA[j-1][now][z];
						    }
						    else
						    {
								v1 += revA[j-1][now][z];
								v2 += revA[j-1][now][z];
						    }

						    if(v1 > revdp[idx].f_max)
						    {
								revdp[idx].f_max = v1;
								revdp[idx].g_max = i-1;
						    }
						    if(v2 < revdp[idx].f_min)
						    {
								revdp[idx].f_min = v2;
								revdp[idx].g_min = i-1;
						    }
						}

						else
						{
						    v1 = revdp[laidx2].f_min;
						    v2 = revdp[laidx2].f_max;

						    if(revB2[j-1][now][z] == -1)
						    {
								v1 -= revA[j-1][now][z];
								v2 -= revA[j-1][now][z];
						    }
						    else
						    {
								v1 += revA[j-1][now][z];
								v2 += revA[j-1][now][z];
						    }

						    if(v1 < revdp[idx].f_min)
						    {
								revdp[idx].f_min = v1;
								revdp[idx].g_min = i-1;
						    }
						    if(v2 > revdp[idx].f_max)
						    {
								revdp[idx].f_max = v2;
								revdp[idx].g_max = i-1;
						    }
						}
						revh[prev_h_index] = st + 1;
				    }
				}
			}
		}
    }
} // end of revDP

void subseq2strobeseeding::combine(std::string s, size_t start, size_t end, DPCell* dp, DPCell* revdp, std::vector<std::vector<seed>>& seeds)
{
    int dp_index, dp_index2;
    int len = end - start + 1;
    int del = n - k;
    int num;

    int64_t ans1[MAXK];
    int ans2[MAXK];
    int ans3[MAXK];
    std::vector<std::vector<seed>> seedtmp(num_valid, std::vector<seed>(0));

    for(int st = 0; st + n - 1 < len; st++)
    {
		for(int j = 1; j <= k - 2; j++)
		{
		    ans1[j] = -INF;
		    ans2[j] = -1;
		}

		for(int i = 1; i < n-1; i++)
		{
		    int minj = std::max(1, i - del);
		    int maxj = std::min(i, k - 2);

		    int nt = alphabetIndex(s[start + i + st]); // mid nucleotide

		    for(int j = minj; j <= maxj; j++)
			if(valid[j])
			{
			    for(int q = 0; q < d; q++)
			    {
					int d3 = C3[j][nt];
					dp_index = dpIndex(st+i+1, n-i-1, k-j-1, q);
					dp_index2 = dpIndex(st+i-1, i, j, (d-d3-q+d)%d);

					if(dp[dp_index].f_min >= INF || revdp[dp_index2].f_min >= INF)
						continue;

					int64_t v = combine3[j][nt] * A3[j][nt];

					if(combine1[j][nt] == 1)
					    v += dp[dp_index].f_max;
					else
					    v -= dp[dp_index].f_min;

					if(combine2[j][nt] == 1)
					    v += revdp[dp_index2].f_max;
					else
					    v -= revdp[dp_index2].f_min;

					if(v > ans1[j])
					{
					    ans1[j] = v;
					    ans2[j] = i;
					    ans3[j] = q;
					}
			    }
			}
		}	

		for(int j = 1; j <= k - 2; j++)
		{
		    seed tmp;
		    if(!valid[j])
				continue;
			if(ans2[j] == -1) //make sure each position has a seed
			{
				tmp.hashval = 0;
				seedtmp[valid[j] - 1].push_back(tmp);
				continue;
			}

		    kmer hashval = 0;
		    std::vector<size_t> index, tmp1;

		    int x = ans2[j];
		    int y = j;
		    int nt = alphabetIndex(s[start + x + st]);

		    int d1 = ans3[j];
		    int d3 = C3[j][nt];
		    int d2 = (d + d - d3 - d1) % d;
		    size_t nextpos;
		    int nextval, nextd;

		    bool z = 1;
		    if(combine2[j][nt] < 0)
				z = 0;

		    while(y > 0)
		    {
				dp_index = dpIndex(st+ans2[j]-1, x, y, d2);
				if(z == 1)
				{
				    nextpos = start + st + ans2[j] - revdp[dp_index].g_max - 1;
				    nextval = alphabetIndex(s[nextpos]);

				    index.push_back(nextpos);
				    hashval = (hashval<<2) | nextval;

				    if(revB1[y-1][nextval][d2] == -1)
						z = 0;

				    nextd = (d2 + d - C2[y-1][nextval]) % d;
				    x = revdp[dp_index].g_max;
				    y--;
				    d2 = nextd;
				}
				else
				{
				    nextpos = start + st + ans2[j] - revdp[dp_index].g_min - 1;
				    nextval = alphabetIndex(s[nextpos]);

				    index.push_back(nextpos);
				    hashval = (hashval<<2) | nextval;

				    if(revB1[y-1][nextval][d2] == -1)
						z = 1;
							
				    nextd = (d2 + d - C2[y-1][nextval]) % d;
				    x = revdp[dp_index].g_min;
				    y--;
				    d2 = nextd;
				}
		    }

		    index.push_back(start + st + ans2[j]);
		    hashval = (hashval<<2) | alphabetIndex(s[start + st + ans2[j]]);

		    x = n - ans2[j] - 1;
		    y = k - j - 1;

		    z = 1;
		    if(combine1[j][nt] < 0)
				z = 0;

		    //cout<<j<<" "<<ans1[j]<<" "<<ans2[j]<<endl;
		    while(y > 0)
		    {
				dp_index = dpIndex(st+ans2[j]+1, x, y, d1);
				if(z == 1)
				{	
				    nextpos = start + st + ans2[j] + dp[dp_index].g_max + 1;
				    nextval = alphabetIndex(s[nextpos]);

				    tmp1.push_back(nextpos);

				    if(B1[y-1][nextval][d1] == -1)
						z = 0;
							
				    nextd = (d1 + d - C1[y-1][nextval]) % d;
				    x = dp[dp_index].g_max;
				    y--;
				    d1 = nextd;
				}
				else
				{								
				    nextpos = start + st + ans2[j] + dp[dp_index].g_min + 1;
				    nextval = alphabetIndex(s[nextpos]);

				    tmp1.push_back(nextpos);

				    if(B1[y-1][nextval][d1] == -1)
						z = 1;
							
				    nextd = (d1 + d - C1[y-1][nextval]) % d;
				    x = dp[dp_index].g_min;
				    y--;
				    d1 = nextd;
				}
		    }


		    tmp.st = start + st;
		    tmp.ed = tmp1[0];
		    tmp.index = 0;

		    for(size_t a: index)
		    	tmp.index |= 1<<(a - tmp.st);

		    for(int a = k - j - 2; a >= 0; a--)
		    {
		    	tmp.index |= 1<<(tmp1[a] - tmp.st);			
				hashval = (hashval<<2) | alphabetIndex(s[tmp1[a]]);
		    }

		    num = valid[j] - 1;

		    tmp.hashval = ans1[j];
		    tmp.str = hashval;
		    //tmp.str_rc = revComp(hashval, k);
		    seedtmp[num].push_back(tmp);
		}
    }


    for(int st = 0; st + w * n + prek <= len; st++)
    {
    	for(int j = 1; j <= k-2; j++)
    		if(valid[j])
	    	{
	    		int num = valid[j] - 1;

	    		if(seedtmp[num][st + prek].hashval == 0)
	    			break;
	    		if(w == 2 && seedtmp[num][st + prek + n].hashval == 0)
	    			break;

	    		seed tmp;
	    		tmp.st = st + start;

	    		tmp.ed = seedtmp[num][st + prek + (w-1) * n].ed;

	    		tmp.index = ((1<<prek) - 1);
	    		for(int window = 0; window < w; window++)
	    			tmp.index |= (seedtmp[num][st + prek + window * n].index << (prek + window * n));

	    		kmer enc = 0lu;
				for(int i = 0; i < prek; i+=1)
					enc = (enc << 2) | alphabetIndex(s[st + start + i]);

				tmp.hashval = enc;

				for(int window = 0; window < w; window++)
					enc = (enc << (2*k)) + seedtmp[num][st + prek + window * n].str;

				tmp.str = enc;

				for(int window = 0; window < w; window++)
					tmp.hashval += seedtmp[num][st + prek + window * n].hashval;

				seeds[num].push_back(tmp);
	    	}
    }

} // end of combine

void subseq2strobeseeding::getSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, int* h, int* revh,
		     std::vector<std::vector<seed>>& seeds)
{	
    int len = s.length();

    int st = 0, en = 0;

    while(en < len - 1)
    {
		en = st + chunk_size - 1;
		if(en >= len)
		    en = len - 1;

		DP(s, st, en, dp, h);
		revDP(s, st, en, revdp, revh);
		combine(s, st, en, dp, revdp, seeds);

		st = en - 2 * n  - k + 2;
    }
}

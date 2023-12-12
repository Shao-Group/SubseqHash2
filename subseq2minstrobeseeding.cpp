#include "subseq2minstrobeseeding.h"
const int64_t INF=(int64_t)1<<62;

//accessing DPCell[chunk_size][n+1][k+1][d]
inline int subseq2minstrobeseeding::dpIndex(int d1, int d2, int d3, int d4)
{
    return d1 * dim1 + d2 * dim2 + d3 * dim3 + d4;
}

//accessing int[n+1][k+1][d]
inline int subseq2minstrobeseeding::hIndex(int d2, int d3, int d4)
{
    return d2 * dim2 + d3 * dim3 + d4;
}

void subseq2minstrobeseeding::init(const char* table_filename)
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

    for(int i = 0; i < k; i++)
    {
	fscanf(filein, "%d", &x);
	if(i < subsample)
	    valid[x] = i+1;
    }
}

void subseq2minstrobeseeding::DP(std::string s, size_t start, size_t end, DPCell* dp, int* h)
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
	    int maxj = std::min(i, k - 1);

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


void subseq2minstrobeseeding::revDP(std::string s, size_t start, size_t end, DPCell* revdp, int* revh)
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
	    int maxj = std::min(i, k - 1);

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

void subseq2minstrobeseeding::combine(std::string s, size_t start, size_t end, DPCell* dp, DPCell* revdp, std::vector<std::vector<seed>>& seeds)
{
    int dp_index, dp_index2;
    int len = end - start + 1;
    int del = n - k;
    int num;

    int64_t ans1[MAXK];
    int ans2[MAXK];
    int ans3[MAXK];
    int ans4[MAXK];

    std::vector<std::vector<seed>> seedtmp(num_valid, std::vector<seed>(0));

    for(int st = 0; st + n - 1 < len; st++)
    {
	for(int j = 0; j < k; j++)
	{
	    ans1[j] = -INF;
	    ans2[j] = -1;
	    ans4[j] = d;
	}

	for(int i = 0; i < n; i++)
	{
	    int nt = alphabetIndex(s[start + i + st]); // mid nucleotide

	    if(n - i >= k && valid[0])
	    {
		int d3 = C3[0][nt];
		int d1 = (d-d3+d)%d;
		dp_index = dpIndex(st+i+1, n-i-1, k-1, d1);

		while(dp[dp_index].f_min >= INF)
		{	
		    d1 = (d1+1) % d;
		    dp_index = dpIndex(st+i+1, n-i-1, k-1, d1);
		}

		int64_t v = combine3[0][nt] * A3[0][nt];

		if(combine1[0][nt] == 1)
		    v += dp[dp_index].f_max;
		else
		    v -= dp[dp_index].f_min;

		if((v > ans1[0] && (d1+d3) % d == ans4[0]) || (d1+d3) % d < ans4[0])
		{
		    ans1[0] = v;
		    ans2[0] = i;
		    //ans3[0] = (d-d3+d)%d;
		    ans3[0] = d1;
		    ans4[0] = (d1+d3) % d; //omega, position, j, psi
		}
	    }

	    if(i >= k && valid[k-1])
	    {
		int d3 = C3[k-1][nt];
		int d1 = (d-d3+d)%d;
		dp_index2 = dpIndex(st+i-1, i, k-1, d1);

		while(revdp[dp_index2].f_min >= INF)
		{	
		    d1 = (d1+1) % d;
		    dp_index2 = dpIndex(st+i-1, i, k-1, d1);
		}			

		int64_t v = combine3[k-1][nt] * A3[k-1][nt];

		if(combine2[k-1][nt] == 1)
		    v += revdp[dp_index2].f_max;
		else
		    v -= revdp[dp_index2].f_min;

		if((v > ans1[k-1] && (d1+d3) % d == ans4[k-1]) || (d1+d3) % d < ans4[k-1])
		{
		    ans1[k-1] = v;
		    ans2[k-1] = i;
		    ans3[k-1] = 0;
		    ans4[k-1] = (d1+d3) % d;
		}
	    }

	    int minj = std::max(1, i - del);
	    int maxj = std::min(i, k - 2);

	    for(int j = minj; j <= maxj; j++)
		if(valid[j])
		{				
		    dp_index = dpIndex(st+i+1, n-i-1, k-j-1, 0);
		    dp_index2 = dpIndex(st+i-1, i, j, 0);
		    int d2num = 0, totald = 0;

		    for(int q = 0; q < d; q++)
			if(revdp[dp_index2 + q].f_min < INF)
			    d2num += (1<<q);

		    for(int q = 0; q < d; q++)
			if(dp[dp_index + q].f_min < INF)
			    totald |= (d2num<<q);

		    int targetd = 0;
		    int d3 = C3[j][nt], d1 = (d - d3) % d;

		    while(targetd < d)
		    {
			if(((totald>>d1)&1) || ((totald>>(d1+d)&1)))
			    break;
			d1 = (d1+1) % d;
			targetd++;
		    }

		    for(int q = 0; q < d; q++)
		    {
			dp_index = dpIndex(st+i+1, n-i-1, k-j-1, q);
			dp_index2 = dpIndex(st+i-1, i, j, (targetd-d3-q+2*d)%d);

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

			if((v > ans1[j] && targetd == ans4[j]) || targetd < ans4[j]) 
			{
			    ans1[j] = v;
			    ans2[j] = i;
			    ans3[j] = q;
			    ans4[j] = targetd;
			}
		    }
		}
	}	

	for(int j = 0; j < k; j++)
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
	    int d2 = (ans4[j] + 2 * d - d1 - d3) % d;
	    size_t nextpos;
	    int nextval, nextd;

	    bool z = 1;
	    if(combine2[j][nt] < 0)
		z = 0;

	    if(j > 0)
	    {
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
	    }

	    index.push_back(start + st + ans2[j]);
	    hashval = (hashval<<2) | alphabetIndex(s[start + st + ans2[j]]);

	    if(j < k-1)
	    {
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
	    }

	    tmp.st = start + st;
	    if(j < k-1)
		tmp.ed = tmp1[0];
	    else
		tmp.ed = index[k-1];

	    tmp.index = 0;

	    for(size_t a: index)
		tmp.index |= ((uint64_t)1)<<(a - tmp.st);

	    for(int a = k - j - 2; a >= 0; a--)
	    {
		tmp.index |= ((uint64_t)1)<<(tmp1[a] - tmp.st);			
		hashval = (hashval<<2) | alphabetIndex(s[tmp1[a]]);
	    }

	    num = valid[j] - 1;

	    tmp.hashval = ans1[j];
	    tmp.str = hashval;
	    tmp.psi = ans4[j];
	    //tmp.str_rc = revComp(hashval, k);
	    seedtmp[num].push_back(tmp);
	}
    }


    for(int st = 0; st + w * n + prek <= len; st++)
    {
    	for(int j = 0; j < k; j++)
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

		tmp.index = (((uint64_t)1)<<prek) - 1;
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

		tmp.psi = seedtmp[num][st + prek].psi;
		seeds[num].push_back(tmp);
	    }
    }

} // end of combine

/*
  Minimizer implementation from minimap2.
  !!This hash function can only handle kmers of length <= 32!!
*/
static inline uint64_t hash64(kmer s, uint64_t mask)
{
    uint64_t key = (uint64_t) s;
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

struct MMEntry{
    size_t st;
    uint64_t hashval;
    MMEntry(int st, uint64_t hashval): st(st), hashval(hashval){}
};

/*
  Subsampling all_seeds by prek minimizers in mm_w.
*/
static void minimizerSketchSeeds(const std::string& s,
				 const int k, const int w, const int num_valid,
				 const std::vector<std::vector<seed>>& all_seeds,
				 std::vector<std::vector<seed>>& seeds){
    //get minimizer positions
    std::vector<size_t> mm_pos;
    kmer mask = (1ULL<<(k<<1)) - 1;
    size_t len = s.length();
    std::list<MMEntry> candidates;
    size_t j = len - k + 1;
    j = (j < w ? j : w);

    char cur[k];
    s.copy(cur, k, 0);
    kmer now = encode(cur, k);
    uint64_t hashval = hash64(now, mask);

    candidates.push_back(MMEntry(0, hashval));

    for(size_t i=1; i<j; ++i){
	now = (now<<2) & mask;
	now |= alphabetIndex(s[i+k-1]);
	hashval = hash64(now, mask);

	while(!candidates.empty() && hashval<=candidates.back().hashval){
	    candidates.pop_back();
	}
	candidates.push_back(MMEntry(i, hashval));
    }

    for(size_t i=w; i<= len - k; ++i){
	while(!candidates.empty() && candidates.front().st < i-w){
	    candidates.pop_front();
	}

	if(mm_pos.empty() || candidates.front().st != mm_pos.back()){
	    mm_pos.push_back(candidates.front().st);
	}

	now = (now<<2) & mask;
	now |= alphabetIndex(s[i+k-1]);
	hashval = hash64(now, mask);

	while(!candidates.empty() && hashval<=candidates.back().hashval){
	    candidates.pop_back();
	}

	candidates.push_back(MMEntry(i, hashval));
    }

    while(!candidates.empty() && candidates.front().st <= len - k - w){
	candidates.pop_front();
    }

    if(mm_pos.empty() || candidates.front().st != mm_pos.back()){
	mm_pos.push_back(candidates.front().st);
    }

    //subsample seeds
    std::vector<size_t> it(num_valid, 0);
    for(size_t st : mm_pos){
	for(int i=0; i<num_valid; ++i){
	    size_t& c = it[i];
	    while(c < all_seeds[i].size() && all_seeds[i][c].st < st) ++c;
	    if(c < all_seeds[i].size() && all_seeds[i][c].st == st){
		seeds[i].push_back(all_seeds[i][c]);
	    }
	}
    }
}

void subseq2minstrobeseeding::getSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, int* h, int* revh,
					   std::vector<std::vector<seed>>& seeds)
{	
    int len = s.length();

    int st = 0, en = 0;

    std::vector<std::vector<seed>> all_seeds(num_valid, std::vector<seed>(0));

    while(en < len - 1)
    {
	en = st + chunk_size - 1;
	if(en >= len)
	    en = len - 1;

	DP(s, st, en, dp, h);
	revDP(s, st, en, revdp, revh);
	combine(s, st, en, dp, revdp, all_seeds);

	st = en - 2 * n  - k + 2;
    }
    minimizerSketchSeeds(s, prek, mm_w, num_valid, all_seeds, seeds);
}

double subseq2minstrobeseeding::getSeeds(std::string& s, const size_t s_idx,
				      const char* output_dir, const int dir_len,
				      DPCell* dp, DPCell* revdp,
				      int* h, int* revh){
    std::vector<std::vector<seed>> seeds(num_valid, std::vector<seed>(0));
    getSubseq2Seeds(s, dp, revdp, h, revh, seeds);
    int seed_len = prek + k;

    double density = 0.0;
    char output_filename[500];
    for(int i=0; i<num_valid; ++i){
	sprintf(output_filename, "%.*s/%d-%zu.ss2msseed",
		dir_len, output_dir, i<<1, s_idx);
	//saveSeeds(output_filename, seed_len, seeds[i]);
	saveSeedsWithScore(output_filename, seed_len, seeds[i]);
	density += seeds[i].size();
	seeds[i].clear();
    }

    getSubseq2Seeds(revComp(s), dp, revdp, h, revh, seeds);
    for(int i=0; i<num_valid; ++i){
	sprintf(output_filename, "%.*s/%d-%zu.ss2msseed",
		dir_len, output_dir, (i<<1)+1, s_idx);
	//saveSeeds(output_filename, seed_len, seeds[i]);
	saveSeedsWithScore(output_filename, seed_len, seeds[i]);
	density += seeds[i].size();
    }
    
    
    return density/(s.length()*num_valid*2);
}

#pragma GCC target("avx512vl")
// #pragma GCC optimize("O3")

#include "subseq2simdseeding.hpp"


const short INF= 32767;

//accessing DPCell[chunk_size][n+1][k+1][d]
inline int subseqhash2SIMDseeding::dpIndex(int d1, int d2, int d3)
{
    return d1 * dim1 + d2 * dim2 + d3;
}

//accessing int[n+1][k+1][d]
inline int subseqhash2SIMDseeding::hIndex(int d2, int d3, int d4)
{
    return MAXD * (d2 * dim2 + d3) + d4;
}

void subseqhash2SIMDseeding::init(const char* table_filename)
{
    FILE *filein = fopen(table_filename, "r");

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				(void)!fscanf(filein, "%hd", &A[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				(void)!fscanf(filein, "%hd,%hd", &B1[i][j][q], &B2[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	(void)!fscanf(filein, "%hd", &C1[i][j]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				(void)!fscanf(filein, "%hd", &revA[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	   		for(int q = 0; q < d; q++)
				(void)!fscanf(filein, "%hd,%hd", &revB1[i][j][q], &revB2[i][j][q]);

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	(void)!fscanf(filein, "%hd", &C2[i][j]);
	
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	(void)!fscanf(filein, "%hd,%hd", &combine1[i][j], &combine2[i][j]);
	

    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	(void)!fscanf(filein, "%hd", &A3[i][j]);
	
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	   		(void)!fscanf(filein, "%hd ", &combine3[i][j]);
	
    for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
	    	(void)!fscanf(filein, "%hd", &C3[i][j]);


    int x;
    int subsample = num_valid;

    // if((num_valid&1) && (k&1))
    // {
    // 	subsample--;
    // 	valid[k/2] = subsample+1;
    // }

    for(int i = 0; i < k; i++)
    {
		(void)!fscanf(filein, "%d", &x);
		if(i < subsample)
			valid[x] = i+1;
    }
}

void subseqhash2SIMDseeding::DP(std::string s, size_t start, size_t end, DPCell* dp, short* h)
{
    int dp_index, h_index;
	unsigned short cho;

    int len = end - start + 1;

    for(int i = 0; i < dim1 * MAXD; i++)
    	h[i] = 0;

	__m512i z32 = _mm512_setzero_si512();
	__m512i maxv = _mm512_set1_epi16(INF);
	__m512i minv = _mm512_set1_epi16(-INF);

    for(short st = 0; st < len; st++)
    {
		size_t real_index_st = st + start;

		__m512i stmk = _mm512_set1_epi16(st+1);

		if(s[real_index_st] == 'N')
			continue;

		short en = std::min(st + n - 1, len - 1);
		short interval_len = en - st + 1;

		for(short i = 0; i <= interval_len; i++)
		{
		    h_index = hIndex(i, 0, 0); //[i][0][0]
		    dp_index = dpIndex(st, i, 0); //[st][i][0][0]

		    dp[dp_index].f_max[0] = 0;
		    dp[dp_index].f_min[0] = 0;
		    h[h_index] = st + 1;
		    
		    for(int j = 1; j < d; ++j)
		    {
				dp[dp_index].f_max[j] = -INF;
				dp[dp_index].f_min[j] = INF;
		    }
		}

		int d1 = C1[0][alphabetIndex(s[real_index_st])];
		
		h_index = hIndex(1, 1, d1); //[1][1][d1]
		h[h_index] = st + 1;

		dp_index = dpIndex(st, 1, 1); //[st][1][1]

		for(int i = 0; i < d; ++i)
		{
			dp[dp_index].f_max[i] = -INF;
			dp[dp_index].f_min[i] = INF;
		}

		dp[dp_index].f_min[d1] = dp[dp_index].f_max[d1] = B2[0][alphabetIndex(s[real_index_st])][d1] * A[0][alphabetIndex(s[real_index_st])][d1];
		dp[dp_index].g_min[d1] = dp[dp_index].g_max[d1] = 0;

		int64_t v1, v2;

		for(int i = 2; i <= interval_len; i++)
		{
		    if(s[real_index_st + i - 1] == 'N')
		    	break;

		    int minj = std::max(1, i - n + k);
		    int maxj = std::min(i, k - 1);

		    for(int j = minj; j <= maxj; j++)
		    {
				int now = alphabetIndex(s[real_index_st + i - 1]); //Current char
				int v = C1[j-1][now];

	    		__m512i index = _mm512_load_si512(&shift[d-v]);

				h_index = hIndex(i-1, j-1, 0);

				__m512i nh = _mm512_load_si512(&h[h_index]);
				nh = _mm512_permutexvar_epi16(index, nh);
	    		__mmask32 mk2 = _mm512_cmpeq_epi16_mask(nh, stmk);

				dp_index = dpIndex(st, i-1, j-1);
			
				__m512i v1 = _mm512_load_si512(&dp[dp_index].f_min[0]);
				v1 = _mm512_permutexvar_epi16(index, v1);

				__m512i v2 = _mm512_load_si512(&dp[dp_index].f_max[0]);
				v2 = _mm512_permutexvar_epi16(index, v2);
				

	    		__m512i a = _mm512_load_si512(&A[j-1][now][0]);
	    		__m512i b = _mm512_load_si512(&B2[j-1][now][0]);

	    		// __m512i odd = _mm512_mul_epi16(a, b);
	    		// __m512i even = _mm512_mul_epi16(_mm512_srli_epi64(a,32), _mm512_srli_epi64(b,32));
	    		// a = _mm512_mask_shuffle_epi16(odd, 0xaaaa, even, _MM_SHUFFLE(2,2,0,0));	

	    		a = _mm512_mullo_epi16(a, b);

	    		b = _mm512_load_si512(&B1[j-1][now][0]);
	    		// odd = _mm512_mullo_epi16(v1, b);
	    		// even = _mm512_mul_epi16(_mm512_srli_epi64(v1,32), _mm512_srli_epi64(b,32));
	    		// v1 = _mm512_mask_shuffle_epi16(odd, 0xaaaa, even, _MM_SHUFFLE(2,2,0,0));
	    		v1 = _mm512_mullo_epi16(v1, b);

				v1 = _mm512_add_epi16(v1, a);

				// odd = _mm512_mul_epi16(v2, b);
	    		// even = _mm512_mul_epi16(_mm512_srli_epi64(v2,32), _mm512_srli_epi64(b,32));
	    		// v2 = _mm512_mask_shuffle_epi16(odd, 0xaaaa, even, _MM_SHUFFLE(2,2,0,0));
	    		v2 = _mm512_mullo_epi16(v2, b);

				v2 = _mm512_add_epi16(v2, a);
				
				__mmask32 sign = _mm512_cmp_epi16_mask(b, z32, 5);

				__m512i vmin = _mm512_mask_mov_epi16(v2, sign, v1);
				__m512i vmax = _mm512_mask_mov_epi16(v1, sign, v2);

				vmin = _mm512_mask_mov_epi16 (maxv, mk2, vmin);
				vmax = _mm512_mask_mov_epi16 (minv, mk2, vmax);

			    h_index = hIndex(i-1, j, 0);
			    dp_index = dpIndex(st, i-1, j);
				int dest_index = dpIndex(st, i, j);

	    		__m512i ch = _mm512_load_si512(&h[h_index]);
	    		sign = _mm512_cmpeq_epi16_mask(ch, stmk);

				__m512i f_prev = _mm512_load_si512(&dp[dp_index].f_min[0]);
				f_prev = _mm512_mask_mov_epi16 (maxv, sign, f_prev);				unsigned int valid1, valid2;
				__mmask32 choice = _mm512_cmp_epi16_mask(vmin, f_prev, 1);
	    		
				v1 = _mm512_mask_mov_epi16(f_prev, choice, vmin);
	    		_mm512_stream_si512((__m512i*)dp[dest_index].f_min, v1);

    			f_prev = _mm512_load_si512(&dp[dp_index].f_max[0]);
				f_prev = _mm512_mask_mov_epi16 (minv, sign, f_prev);
				__mmask32 choice2 = _mm512_cmp_epi16_mask(vmax, f_prev, 6);

				v2 = _mm512_mask_mov_epi16(f_prev, choice2, vmax);
	    		_mm512_stream_si512((__m512i*)dp[dest_index].f_max, v2);

				__m512i g = _mm512_set1_epi16(i-1);
				__m512i g_prev = _mm512_load_si512(&dp[dp_index].g_min[0]);
				v1 = _mm512_mask_mov_epi16(g_prev, choice, g);
    			g_prev = _mm512_load_si512(&dp[dp_index].g_max[0]);
				v2 = _mm512_mask_mov_epi16(g_prev, choice2, g);

	    		_mm512_stream_si512((__m512i*)dp[dest_index].g_min, v1);
	    		_mm512_stream_si512((__m512i*)dp[dest_index].g_max, v2);

				h_index = hIndex(i, j, 0);
				ch = _mm512_mask_mov_epi16 (ch, mk2, stmk);
	    		_mm512_store_si512(h + h_index, ch);

	    		// printf("%d %d\n", i, j);
	    		// for(int z = 0; z < d; z++)
	    		// 	printf("%d %d\n", dp[dest_index].f_min[z], dp[dest_index].f_max[z]);
			}
		}
	}
} // end of DP


void subseqhash2SIMDseeding::revDP(std::string s, size_t start, size_t end, DPCell* dp, short* h)
{
    int dp_index, h_index;

    int len = end - start + 1;

    for(int i = 0; i < dim1 * MAXD; i++)
    	h[i] = 0;

	__m512i z32 = _mm512_setzero_si512();
	__m512i maxv = _mm512_set1_epi16(INF);
	__m512i minv = _mm512_set1_epi16(-INF);

    for(short st = len - 1; st >= 0; st--)
    {
		size_t real_index_st = st + start;

		__m512i stmk = _mm512_set1_epi16(st+1);

		if(s[real_index_st] == 'N')
			continue;

		short en = std::max(st - n + 1, 0);
		short interval_len = st - en + 1;

		for(short i = 0; i <= interval_len; i++)
		{
		    h_index = hIndex(i, 0, 0); //[i][0][0]
		    dp_index = dpIndex(st, i, 0); //[st][i][0][0]

		    dp[dp_index].f_max[0] = 0;
		    dp[dp_index].f_min[0] = 0;
		    h[h_index] = st + 1;
		    
		    for(int j = 1; j < d; ++j)
		    {
				dp[dp_index].f_max[j] = -INF;
				dp[dp_index].f_min[j] = INF;
		    }
		}

		int d1 = C2[0][alphabetIndex(s[real_index_st])];

		h_index = hIndex(1, 1, d1); //[1][1][d1]
		h[h_index] = st + 1;

		dp_index = dpIndex(st, 1, 1); //[st][1][1]
		int idx1 = dp_index + d1;

		for(int i = 0; i < d; ++i)
		{
			dp[dp_index].f_max[i] = -INF;
			dp[dp_index].f_min[i] = INF;
		}

		dp[dp_index].f_min[d1] = dp[dp_index].f_max[d1] = revB2[0][alphabetIndex(s[real_index_st])][d1] * revA[0][alphabetIndex(s[real_index_st])][d1];
		dp[dp_index].g_min[d1] = dp[dp_index].g_max[d1] = 0;
		// for(int z = 0; z < d; z++)
		// 	printf("%d %d %d\n", z, dp[dp_index].f_min[z], dp[dp_index].f_max[z]);

		int64_t v1, v2;

		for(int i = 2; i <= interval_len; i++)
		{
		    if(s[real_index_st + i - 1] == 'N')
		    	break;

		    int minj = std::max(1, i - n + k);
		    int maxj = std::min(i, k - 1);

		    for(int j = minj; j <= maxj; j++)
		    {
				int now = alphabetIndex(s[real_index_st - i + 1]); //Current char
				int v = C2[j-1][now];

	    		__m512i index = _mm512_load_si512(&shift[d-v]);

				h_index = hIndex(i-1, j-1, 0);

				__m512i nh = _mm512_load_si512(&h[h_index]);
				nh = _mm512_permutexvar_epi16(index, nh);
	    		__mmask32 mk2 = _mm512_cmpeq_epi16_mask(nh, stmk);

				dp_index = dpIndex(st, i-1, j-1);
			
				__m512i v1 = _mm512_load_si512(&dp[dp_index].f_min[0]);
				v1 = _mm512_permutexvar_epi16(index, v1);

				__m512i v2 = _mm512_load_si512(&dp[dp_index].f_max[0]);
				v2 = _mm512_permutexvar_epi16(index, v2);
				

	    		__m512i a = _mm512_load_si512(&revA[j-1][now][0]);
	    		__m512i b = _mm512_load_si512(&revB2[j-1][now][0]);

	    		// __m512i odd = _mm512_mul_epi16(a, b);
	    		// __m512i even = _mm512_mul_epi16(_mm512_srli_epi64(a,32), _mm512_srli_epi64(b,32));
	    		// a = _mm512_mask_shuffle_epi16(odd, 0xaaaa, even, _MM_SHUFFLE(2,2,0,0));	
	    		a = _mm512_mullo_epi16(a, b);


	    		b = _mm512_load_si512(&revB1[j-1][now][0]);
	    		// odd = _mm512_mul_epi16(v1, b);
	    		// even = _mm512_mul_epi16(_mm512_srli_epi64(v1,32), _mm512_srli_epi64(b,32));
	    		// v1 = _mm512_mask_shuffle_epi16(odd, 0xaaaa, even, _MM_SHUFFLE(2,2,0,0));
	    		v1 = _mm512_mullo_epi16(v1, b);
				v1 = _mm512_add_epi16(v1, a);

				// odd = _mm512_mul_epi16(v2, b);
	    		// even = _mm512_mul_epi16(_mm512_srli_epi64(v2,32), _mm512_srli_epi64(b,32));
	    		// v2 = _mm512_mask_shuffle_epi16(odd, 0xaaaa, even, _MM_SHUFFLE(2,2,0,0));
	    		v2 = _mm512_mullo_epi16(v2, b);
				v2 = _mm512_add_epi16(v2, a);
				
				__mmask32 sign = _mm512_cmp_epi16_mask(b, z32, 5);

				__m512i vmin = _mm512_mask_mov_epi16(v2, sign, v1);
				__m512i vmax = _mm512_mask_mov_epi16(v1, sign, v2);

				vmin = _mm512_mask_mov_epi16 (maxv, mk2, vmin);
				vmax = _mm512_mask_mov_epi16 (minv, mk2, vmax);

			    h_index = hIndex(i-1, j, 0);
			    dp_index = dpIndex(st, i-1, j);
				int dest_index = dpIndex(st, i, j);

	    		__m512i ch = _mm512_load_si512(&h[h_index]);
	    		sign = _mm512_cmpeq_epi16_mask(ch, stmk);

				__m512i f_prev = _mm512_load_si512(&dp[dp_index].f_min[0]);
				f_prev = _mm512_mask_mov_epi16 (maxv, sign, f_prev);
				__mmask32 choice = _mm512_cmp_epi16_mask(vmin, f_prev, 1);
	    		
				v1 = _mm512_mask_mov_epi16(f_prev, choice, vmin);
	    		_mm512_stream_si512((__m512i*)dp[dest_index].f_min, v1);

    			f_prev = _mm512_load_si512(&dp[dp_index].f_max[0]);
				f_prev = _mm512_mask_mov_epi16 (minv, sign, f_prev);
				__mmask32 choice2 = _mm512_cmp_epi16_mask(vmax, f_prev, 6);

				v2 = _mm512_mask_mov_epi16(f_prev, choice2, vmax);
	    		_mm512_stream_si512((__m512i*)dp[dest_index].f_max, v2);

				__m512i g = _mm512_set1_epi16(i-1);
				__m512i g_prev = _mm512_load_si512(&dp[dp_index].g_min[0]);
				v1 = _mm512_mask_mov_epi16(g_prev, choice, g);
    			g_prev = _mm512_load_si512(&dp[dp_index].g_max[0]);
				v2 = _mm512_mask_mov_epi16(g_prev, choice2, g);

	    		_mm512_stream_si512((__m512i*)dp[dest_index].g_min, v1);
	    		_mm512_stream_si512((__m512i*)dp[dest_index].g_max, v2);


				h_index = hIndex(i, j, 0);
				ch = _mm512_mask_mov_epi16 (ch, mk2, stmk);
	    		_mm512_store_si512(h + h_index, ch);
			}
		}
	}
} // end of revDP

void subseqhash2SIMDseeding::combine(std::string s, size_t start, size_t end, DPCell* dp, 
		DPCell* revdp, std::vector<std::vector<seed>>& seeds)
{
    int dp_index, dp_index2;
    int len = end - start + 1;
    int del = n - k;
    int num;

    int64_t ans1[MAXK];
    int ans2[MAXK];
    int ans3[MAXK];
    int ans4[MAXK];
	__m512i maxv = _mm512_set1_epi16(INF);

    for(int st = 0; st + n - 1 < len; st++)
    {
    	bool skip = 0;
    	int nst;

		for(int i = st; i < st + n; i++)
			if(s[start + i] == 'N')
			{
				nst = i;
				skip = 1;
			}	

		if(skip)
		{
			st = nst;
			continue;
		}

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
				dp_index = dpIndex(st+i+1, n-i-1, k-1);

				while(dp[dp_index].f_min[d1] >= INF)
					d1 = (d1+1) % d;

				int v = combine3[0][nt] * A3[0][nt];

				if(combine1[0][nt] == 1)
				    v += dp[dp_index].f_max[d1];
				else
				    v -= dp[dp_index].f_min[d1];

				if((v > ans1[0] && (d1+d3) % d == ans4[0]) || (d1+d3) % d < ans4[0])
				{
				    ans1[0] = v;
				    ans2[0] = i;
				   //ans3[0] = (d-d3+d)%d;
					ans3[0] = d1;
					ans4[0] = (d1+d3) % d; //omega, position, j, psi
				}
		    }

		    if(i >= k-1 && valid[k-1])
		    {
	    		int d3 = C3[k-1][nt];
	    		int d1 = (d-d3+d)%d;
				dp_index2 = dpIndex(st+i-1, i, k-1);

				while(revdp[dp_index2].f_min[d1] >= INF)
					d1 = (d1+1) % d;

				int v = combine3[k-1][nt] * A3[k-1][nt];

				if(combine2[k-1][nt] == 1)
				    v += revdp[dp_index2].f_max[d1];
				else
				    v -= revdp[dp_index2].f_min[d1];

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
				dp_index = dpIndex(st+i+1, n-i-1, k-j-1);
				dp_index2 = dpIndex(st+i-1, i, j);
				int d2num = 0, totald = 0;

				for(int q = 0; q < d; q++)
					if(revdp[dp_index2].f_min[q] < INF)
						d2num += (1<<q);

				for(int q = 0; q < d; q++)
					if(dp[dp_index].f_min[q] < INF)
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


				short midv = combine3[j][nt] * A3[j][nt];
				__m512i v = _mm512_set1_epi16(midv);

				dp_index = dpIndex(st+i+1, n-i-1, k-j-1);

				__m512i fmax_f = _mm512_loadu_si512(&dp[dp_index].f_max[0]);
				__m512i fmin_f = _mm512_loadu_si512(&dp[dp_index].f_min[0]);

				if(combine1[j][nt] == 1)
					v = _mm512_add_epi16(v, fmax_f);
				else
					v = _mm512_sub_epi16(v, fmin_f);

				dp_index2 = dpIndex(st+i-1, i, j);

				__m512i index = _mm512_load_si512(&comp[(targetd-d3+2*d)%d][0]);
				
				__m512i fmax_r = _mm512_loadu_si512(&revdp[dp_index2].f_max[0]);
				fmax_r = _mm512_permutexvar_epi16(index, fmax_r);
				__m512i fmin_r = _mm512_loadu_si512(&revdp[dp_index2].f_min[0]);
				fmin_r = _mm512_permutexvar_epi16(index, fmin_r);

				if(combine2[j][nt] == 1)
					v = _mm512_add_epi16(v, fmax_r);
				else
					v = _mm512_sub_epi16(v, fmin_r);
 				
 				__mmask32 va1 = _mm512_cmp_epi16_mask(fmin_f, maxv, 1);
				__mmask32 va2 = _mm512_cmp_epi16_mask(fmin_r, maxv, 1);

				unsigned int valid1, valid2;

				_store_mask32((__mmask32*)&valid1, va1);
				_store_mask32((__mmask32*)&valid2, va2);

				valid1 &= valid2;

				short result[32];

	    		_mm512_storeu_si512(&result[0], v);

				// if(st == 0 && j == 1)
				// {
				// 	printf("%d ", i);
				// 	printf("%u\n", valid1);
			    // 	for(int q = 0; q < d; q++)
			    // 		printf("%d ", result[q]);
			    // 	printf("\n");

				// }

			    for(int q = 0; q < d; q++)
					if((valid1>>q)&1)
						if((result[q] > ans1[j] && targetd == ans4[j]) || targetd < ans4[j]) 
						{
						    ans1[j] = result[q];
						    ans2[j] = i;
						    ans3[j] = q;
						    ans4[j] = targetd;
						}

				// if(dp[dp_index].f_min >= INF || revdp[dp_index2].f_min >= INF)
				// 	continue;
			}
		}	
		for(int j = 0; j < k; j++)
		{

		    if(!valid[j] || ans2[j] == -1 || ans1[j] < threshold)
				continue;

		    seed tmp;
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
					dp_index = dpIndex(st+ans2[j]-1, x, y);
					if(z == 1)
					{
					    nextpos = start + st + ans2[j] - revdp[dp_index].g_max[d2] - 1;
					    nextval = alphabetIndex(s[nextpos]);

					    index.push_back(nextpos);
					    hashval = (hashval<<2) | nextval;

					    if(revB1[y-1][nextval][d2] == -1)
							z = 0;

					    nextd = (d2 + d - C2[y-1][nextval]) % d;
					    x = revdp[dp_index].g_max[d2];
					    y--;
					    d2 = nextd;
					}
					else
					{
					    nextpos = start + st + ans2[j] - revdp[dp_index].g_min[d2] - 1;
					    nextval = alphabetIndex(s[nextpos]);

					    index.push_back(nextpos);
					    hashval = (hashval<<2) | nextval;

					    if(revB1[y-1][nextval][d2] == -1)
							z = 1;
								
					    nextd = (d2 + d - C2[y-1][nextval]) % d;
					    x = revdp[dp_index].g_min[d2];
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
					dp_index = dpIndex(st+ans2[j]+1, x, y);
					// printf("%d %d %d\n", st, ans2[j], d1);
					// printf("%d %d %d %d %d\n", x, y, dp_index, g_min[dp_index], g_max[dp_index]);

					if(z == 1)
					{	
					    nextpos = start + st + ans2[j] + dp[dp_index].g_max[d1] + 1;
					    nextval = alphabetIndex(s[nextpos]);

					    tmp1.push_back(nextpos);

					    if(B1[y-1][nextval][d1] == -1)
							z = 0;
								
					    nextd = (d1 + d - C1[y-1][nextval]) % d;
					    x = dp[dp_index].g_max[d1];
					    y--;
					    d1 = nextd;
					}
					else
					{								
					    nextpos = start + st + ans2[j] + dp[dp_index].g_min[d1] + 1;
					    nextval = alphabetIndex(s[nextpos]);

					    tmp1.push_back(nextpos);

					    if(B1[y-1][nextval][d1] == -1)
							z = 1;
								
					    nextd = (d1 + d - C1[y-1][nextval]) % d;
					    x = dp[dp_index].g_min[d1];
					    y--;
					    d1 = nextd;
					}
			    }
			}

		    tmp.st = index[0];
		    tmp.index = 0;
		    if(j < k-1)
		    	tmp.ed = tmp1[0];
			else
				tmp.ed = index[k-1];
		    
		    
		    for(size_t a: index)
		    {
		    	tmp.index |= (((uint64_t)1)<<(a - tmp.st));
		    }

		    for(int a = k - j - 2; a >= 0; a--)
		    {
		    	tmp.index |= (((uint64_t)1)<<(tmp1[a] - tmp.st));
				hashval = (hashval<<2) | alphabetIndex(s[tmp1[a]]);
		    }

		    num = valid[j] - 1;

		    // if(seeds[num].size() > 0 && ans1[j] == seeds[num].back().hashval && tmp.st == seeds[num].back().st 
		    // 	&& seeds[num].back().index == tmp.index)
			//     continue;
		    if(seeds[num].size() > 0 && hashval == seeds[num].back().str && ans4[j] == seeds[num].back().psi)
		    {
		    	if(tmp.st <= seeds[num].back().st && tmp.ed >= seeds[num].back().ed)
					continue;
				else
					if(tmp.st >= seeds[num].back().st && tmp.ed <= seeds[num].back().ed)
					{
						seeds[num].back().st = tmp.st;
						seeds[num].back().index = tmp.index;
						continue;
					}
		    }

		    // if(tmp.index > ((uint64_t)1<<n))
		    // 	printf("%d %d %d %d %llu\n", j, num, tmp.st, tmp.ed, tmp.index);

			//tmp.hashval = ans1[j];
		    
		    tmp.str = hashval;
		    tmp.str_rc = revComp(hashval, k);
		    tmp.hashval = kmerHash{}(std::min(hashval, tmp.str_rc));
		    tmp.psi = ans4[j];
		    // printf("%d %lld %d %d %llu\n", j, tmp.hashval, tmp.st, tmp.ed, tmp.index);

		    seeds[num].push_back(tmp);
		}
    }
} // end of combine

void subseqhash2SIMDseeding::getSubseq2Seeds(std::string s, DPCell* dp, short* h, DPCell* revdp, 
	short* revh, std::vector<std::vector<seed>>& seeds)
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

		st = en - n + 2;
    }
}

void subseqhash2SIMDseeding::writeSubseq2Seeds(std::string s, DPCell* dp, DPCell* revdp, short* h, short* revh,
		     std::vector<std::vector<seed>>& seeds, std::vector<FILE*> fout)
{	
    int len = s.length();

    int st = 0, en = 0;
    uint64_t pos[2];//st, index

    while(en < len - 1)
    {
		en = st + chunk_size - 1;
		if(en >= len)
		    en = len - 1;

		DP(s, st, en, dp, h);
		revDP(s, st, en, revdp, revh);
		combine(s, st, en, dp, revdp, seeds);

		for(int i = 0; i < num_valid; i++)
		{
		    for(auto s : seeds[i])
		    {
				fwrite(&(s.hashval), sizeof(int64_t), 1, fout[i]);
				fwrite(&(s.str), sizeof(kmer), 1, fout[i]);
				pos[0] = s.st;
				pos[1] = s.index;
				fwrite(pos, sizeof(uint64_t), 2, fout[i]);
		    }

		    seeds[i].clear();
		}

		st = en - n + 2;
    }


	for(int i = 0; i < num_valid; i++)
    	fclose(fout[i]);
}

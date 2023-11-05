#include <iostream>
#include <vector>
#include <set>
#include <chrono>
#include <climits>
#include <random>
#include <fstream>
#include <map>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <inttypes.h>

using namespace std;

random_device rd;
const int64_t INF = (int64_t)1<<62;

map<char, int> dict;
int b, d, num_valid;

int64_t f_max[100][100][100][50];
int64_t f_min[100][100][100][50];
bool h[100][100][50];

int64_t revf_max[100][100][100][50];
int64_t revf_min[100][100][100][50];
bool revh[100][100][50];

int64_t A[100][4][50];
int B1[100][4][50];
int B2[100][4][50];

int64_t revA[100][4][50];
int revB1[100][4][50];
int revB2[100][4][50];

int64_t A3[100][4];

int combine1[100][4];
int combine2[100][4];
int combine3[100][4];

int C1[100][4];
int C2[100][4];
int C3[100][4];

int valid[50] = {0};

void init(string path)
{
	dict['A'] = 0;
	dict['C'] = 1;
	dict['G'] = 2;
	dict['T'] = 3;


	FILE *filein = fopen(path.c_str(), "r");

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%" SCNd64, &A[i][j][q]);

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%d,%d", &B1[i][j][q], &B2[i][j][q]);

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C1[i][j]);

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	for(int q = 0; q < d; q++)
				fscanf(filein, "%" SCNd64, &revA[i][j][q]);

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	   		for(int q = 0; q < d; q++)
				fscanf(filein, "%d,%d", &revB1[i][j][q], &revB2[i][j][q]);

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C2[i][j]);
	
    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d,%d", &combine1[i][j], &combine2[i][j]);
	

    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%" SCNd64, &A3[i][j]);
	
    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	   		fscanf(filein, "%d ", &combine3[i][j]);
	
    for(int i = 0; i < b; i++)
		for(int j = 0; j < 4; j++)
	    	fscanf(filein, "%d", &C3[i][j]);


    int x;
    int subsample = num_valid;

    for(int i = 0; i < b; i++)
    {
		fscanf(filein, "%d", &x);
		if(i < subsample)
			valid[x] = i+1;
    }


	fclose(filein);
}

void DP(string s)
{
	int slen = s.length();
	int len = slen;
	int del = len - b;

	for(int st = 0; st < slen; st++)
	{
		memset(h, 0, sizeof(h));

		for(int i = 0; i <= len; i++)
		{
			f_max[st][i][0][0] = 0;
			f_min[st][i][0][0] = 0;
			h[i][0][0] = 1;
			for(int j = 1; j < d; j++)
			{
				f_max[st][i][0][j] = -INF;
				f_min[st][i][0][j] = INF;
			}
		}

		int d1 = C1[0][dict[s[st]]];
		
		for(int i = 0; i < d; i++)
		{
			f_max[st][1][1][i] = -INF;
			f_min[st][1][1][i] = INF;
		}
		f_min[st][1][1][d1] = f_max[st][1][1][d1] = B2[0][dict[s[st]]][d1] * A[0][dict[s[st]]][d1];
		h[1][1][d1] = 1;

		int64_t v1, v2, v3, v4;
		for(int i = 2; i <= len; i++)
		{
			int minj = max(1, i - del);
			int maxj = min(i, b);

			for(int j = minj; j <= maxj; j++)
			{
				int now = dict[s[st + i - 1]];
				int v = C1[j-1][now];

				for(int q = 0; q < d; q++)
				{
					int z = (q + v) % d; // q: previous d, z: current d
					
					if(h[i][j][z] == 0)
					{
						f_min[st][i][j][z] = INF;
						f_max[st][i][j][z] = -INF;
					}

					if(h[i-1][j][z] == 1)
					{
						f_min[st][i][j][z] = min(f_min[st][i][j][z], f_min[st][i-1][j][z]);
						f_max[st][i][j][z] = max(f_max[st][i][j][z], f_max[st][i-1][j][z]);
						h[i][j][z] = 1;
					}

					if(h[i-1][j-1][q])
					{
						if(B1[j-1][now][z] == -1)
						{
							v1 = -f_min[st][i-1][j-1][q];
							v2 = -f_max[st][i-1][j-1][q];

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

							if(v1 > f_max[st][i][j][z])
								f_max[st][i][j][z] = v1;
							if(v2 < f_min[st][i][j][z])
								f_min[st][i][j][z] = v2;
						}

						else
						{
							v1 = f_min[st][i-1][j-1][q];
							v2 = f_max[st][i-1][j-1][q];

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

							if(v1 < f_min[st][i][j][z])
								f_min[st][i][j][z] = v1;
							if(v2 > f_max[st][i][j][z])
								f_max[st][i][j][z] = v2;
						}

						h[i][j][z] = 1;
					}
				}
			}
		}
		len--;
		del += 1;
	}
}

void revDP(string s)
{
	int slen = s.length();
	int len = slen;
	int del = len - b;

	for(int st = slen-1; st >= 0; st--)
	{
		memset(revh, 0, sizeof(revh));

		for(int i = 0; i <= len; i++)
		{
			revf_max[st][i][0][0] = 0;
			revf_min[st][i][0][0] = 0;
			revh[i][0][0] = 1;

			for(int j = 1; j < d; j++)
			{
				revf_min[st][i][0][j] = INF;
				revf_max[st][i][0][j] = -INF;
			}
		}

		int d1 = C2[0][dict[s[st]]];


		for(int i = 0; i < d; i++)
		{
			revf_min[st][1][1][i] = INF;
			revf_max[st][1][1][i] = -INF;
		}

		revf_min[st][1][1][d1] = revf_max[st][1][1][d1] = revB2[0][dict[s[st]]][d1] * revA[0][dict[s[st]]][d1];
		revh[1][1][d1] = 1;

		int64_t v1, v2, v3, v4;
		for(int i = 2; i <= len; i++)
		{
			int minj = max(1, i - del);
			int maxj = min(i, b);

			for(int j = minj; j <= maxj; j++)
			{
				int now = dict[s[st - i + 1]];
				int v = C2[j-1][now];

				for(int q = 0; q < d; q++)
				{
					int z = (q + v) % d; // q: previous d, z: current d
					
					if(revh[i][j][z] == 0)
					{
						revf_min[st][i][j][z] = INF;
						revf_max[st][i][j][z] = -INF;
					}

					if(revh[i-1][j][z] == 1)
					{
						revf_min[st][i][j][z] = min(revf_min[st][i][j][z], revf_min[st][i-1][j][z]);
						revf_max[st][i][j][z] = max(revf_max[st][i][j][z], revf_max[st][i-1][j][z]);
						revh[i][j][z] = 1;
					}

					if(revh[i-1][j-1][q])
					{
						if(revB1[j-1][now][z] == -1)
						{
							v1 = -revf_min[st][i-1][j-1][q];
							v2 = -revf_max[st][i-1][j-1][q];

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

							if(v1 > revf_max[st][i][j][z])
								revf_max[st][i][j][z] = v1;
							if(v2 < revf_min[st][i][j][z])
								revf_min[st][i][j][z] = v2;
						}

						else
						{
							v1 = revf_min[st][i-1][j-1][q];
							v2 = revf_max[st][i-1][j-1][q];

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

							if(v1 < revf_min[st][i][j][z])
								revf_min[st][i][j][z] = v1;
							if(v2 > revf_max[st][i][j][z])
								revf_max[st][i][j][z] = v2;
						}

						revh[i][j][z] = 1;
					}
				}
			}
		}
		len--;
		del += 1;
	}
}

int tp = 0;
int skip = 0;
int cnt[100] = {0};
int globalcnt = 0;

void match(string s, string t)
{
	if(b > t.length())
	{
		skip++;
		return;
	}
	
	DP(s);
	revDP(s);

	int len = s.length();

	int64_t ans1[50];
	int64_t ans2[50];
	int ansd1[50];
	int ansd2[50];

	for(int j = 0; j < b; j++)
	{
		ansd1[j] = ansd2[j] = d;
		ans1[j] = ans2[j] = -INF;
	}

	for(int i = 0; i < len; i++)
	{
		int nt = dict[s[i]];

	    if(len - i >= b && valid[0])
	    {
			int d3 = C3[0][nt];
			int d1 = (d-d3)%d;

			//if(f_min[i+1][len-i-1][b-1][d1] >= INF)
			//	continue;
			while(f_min[i+1][len-i-1][b-1][d1] >= INF)
				d1 = (d1 + 1) % d;

			int64_t v = combine3[0][nt] * A3[0][nt];

			if(combine1[0][nt] == 1)
			    v += f_max[i+1][len-i-1][b-1][d1];
			else
			    v -= f_min[i+1][len-i-1][b-1][d1];

			if((d1+d3) % d < ansd1[0] || (((d1+d3) % d == ansd1[0]) && v > ans1[0]))
			{
			    ans1[0] = v;
			    ansd1[0] = (d1+d3) % d;
			}
	    }

	    if(i >= b-1 && valid[b-1])
	    {
    		int d3 = C3[b-1][nt];
			int d1 = (d-d3)%d;

			//if(revf_min[i-1][i][b-1][(d - d3 + d) % d] >= INF)
			//	continue;

			while(revf_min[i-1][i][b-1][d1] >= INF)
				d1 = (d1 + 1) % d;

			int64_t v = combine3[b-1][nt] * A3[b-1][nt];

			if(combine2[b-1][nt] == 1)
			    v += revf_max[i-1][i][b-1][d1];
			else
			    v -= revf_min[i-1][i][b-1][d1];

			if((d1+d3) % d < ansd1[b-1] || (((d1+d3) % d == ansd1[b-1]) && v > ans1[b-1]))
			{
			    ans1[b-1] = v;
			    ansd1[b-1] = (d1+d3) % d;
			}
		}

		int minj = max(1, b + i - len);
		int maxj = min(i, b-2);

		for(int j = minj; j <= maxj; j++)
		if(valid[j])
		{
			int64_t v = combine3[j][nt] * A3[j][nt];
			int d3 = C3[j][nt];
							
			int d2num = 0, totald = 0;

			for(int q = 0; q < d; q++)
				if(revf_min[i-1][i][j][q] < INF)
					d2num += (1<<q);

			for(int q = 0; q < d; q++)
				if(f_min[i+1][len-i-1][b-j-1][q] < INF)
					totald |= (d2num<<q);

			int targetd = 0;
			int d1 = (d - d3) % d;

			while(targetd < d)
			{
				if(((totald>>d1)&1) || ((totald>>(d1+d)&1)))
					break;
				d1 = (d1+1) % d;
				targetd++;
			}

			for(int q = 0; q < d; q++)
				if(f_min[i+1][len-i-1][b-j-1][q] < INF && revf_min[i-1][i][j][(targetd - d3 - q + 2*d) % d] < INF)
				{
					int64_t v1 = v;
					if(combine1[j][nt] == 1)
						v1 += f_max[i+1][len-i-1][b-j-1][q];
					else
						v1 -= f_min[i+1][len-i-1][b-j-1][q];

					if(combine2[j][nt] == 1)
						v1 += revf_max[i-1][i][j][(targetd - d3 - q + 2*d) % d];
					else
						v1 -= revf_min[i-1][i][j][(targetd - d3 - q + 2*d) % d];

					//cout<<i<<" "<<targetd<<" "<<q<<" "<<(targetd - d3 - q + 2*d) % d<<" "<<v1<<endl;
					if(targetd < ansd1[j] || ((targetd == ansd1[j]) && v1 > ans1[j]))
					{				
						ans1[j] = v1;
						ansd1[j] = targetd;
					}
				}
		}
	}

	DP(t);
	revDP(t);

	len = t.length();
	
	for(int i = 0; i < len; i++)
	{
		int nt = dict[t[i]];

	    if(len - i >= b && valid[0])
	    {
			int d3 = C3[0][nt];
			int d1 = (d-d3)%d;

			//if(f_min[i+1][len-i-1][b-1][d1] >= INF)
			//	continue;
			while(f_min[i+1][len-i-1][b-1][d1] >= INF)
				d1 = (d1 + 1) % d;

			int64_t v = combine3[0][nt] * A3[0][nt];

			if(combine1[0][nt] == 1)
			    v += f_max[i+1][len-i-1][b-1][d1];
			else
			    v -= f_min[i+1][len-i-1][b-1][d1];

			if((d1+d3) % d < ansd2[0] || (((d1+d3) % d == ansd2[0]) && v > ans2[0]))
			{
			    ans2[0] = v;
			    ansd2[0] = (d1+d3) % d;
			}
	    }

	    if(i >= b && valid[b-1])
	    {
    		int d3 = C3[b-1][nt];
			int d1 = (d-d3)%d;

			//if(revf_min[i-1][i][b-1][(d - d3 + d) % d] >= INF)
			//	continue;

			while(revf_min[i-1][i][b-1][d1] >= INF)
				d1 = (d1 + 1) % d;

			int64_t v = combine3[b-1][nt] * A3[b-1][nt];

			if(combine2[b-1][nt] == 1)
			    v += revf_max[i-1][i][b-1][d1];
			else
			    v -= revf_min[i-1][i][b-1][d1];

			if((d1+d3) % d < ansd2[b-1] || (((d1+d3) % d == ansd2[b-1]) && v > ans2[b-1]))
			{
			    ans2[b-1] = v;
			    ansd2[b-1] = (d1+d3) % d;
			}
		}

		int minj = max(1, b + i - len);
		int maxj = min(i, b-2);

		for(int j = minj; j <= maxj; j++)
		if(valid[j])
		{
			int64_t v = combine3[j][nt] * A3[j][nt];
			int d3 = C3[j][nt];
							
			int d2num = 0, totald = 0;

			for(int q = 0; q < d; q++)
				if(revf_min[i-1][i][j][q] < INF)
					d2num += (1<<q);

			for(int q = 0; q < d; q++)
				if(f_min[i+1][len-i-1][b-j-1][q] < INF)
					totald |= (d2num<<q);

			int targetd = 0;
			int d1 = (d - d3) % d;

			while(targetd < d)
			{
				if(((totald>>d1)&1) || ((totald>>(d1+d)&1)))
					break;
				d1 = (d1+1) % d;
				targetd++;
			}

			for(int q = 0; q < d; q++)
				if(f_min[i+1][len-i-1][b-j-1][q] < INF && revf_min[i-1][i][j][(targetd - d3 - q + 2*d) % d] < INF)
				{
					int64_t v1 = v;
					if(combine1[j][nt] == 1)
						v1 += f_max[i+1][len-i-1][b-j-1][q];
					else
						v1 -= f_min[i+1][len-i-1][b-j-1][q];

					if(combine2[j][nt] == 1)
						v1 += revf_max[i-1][i][j][(targetd - d3 - q + 2*d) % d];
					else
						v1 -= revf_min[i-1][i][j][(targetd - d3 - q + 2*d) % d];

					//cout<<i<<" "<<targetd<<" "<<q<<" "<<(targetd - d3 - q + 2*d) % d<<" "<<v1<<endl;
					if(targetd < ansd2[j] || ((targetd == ansd2[j]) && v1 > ans2[j]))
					{				
						ans2[j] = v1;
						ansd2[j] = targetd;
					}
				}
		}
	}
	globalcnt++;

	for(int j = 0; j < b; j++)
	{
		// if(valid[j])
		// 	cout<<globalcnt<<" "<<ans1[j]<<" "<<ansd1[j]<<" "<<ans2[j]<<" "<<ansd2[j]<<endl;

		if(valid[j] && ans1[j] > -INF && (ans1[j] == ans2[j]) && (ansd1[j] == ansd2[j]))
		{
			//cout<<globalcnt<<" "<<ans1[j]<<" "<<ansd1[j]<<endl;
			tp++;
			break;
		}		
	}
	for(int j = 0; j < b; j++)
		if(valid[j] && ans1[j] > -INF && (ans1[j] == ans2[j]) && (ansd1[j] == ansd2[j]))
		{
			cnt[j]++;
		}	
}


int main(int argc, const char * argv[])
{
	b = stoi(argv[1]);
	d = stoi(argv[2]);
	num_valid = stoi(argv[3]);
	string save = argv[4];

	init(save);

	string s, t;
	int n = 0;


	while(cin>>s)
	{
		cin>>t;
		match(s, t);
		n++;
	}
	printf("%.2lf\n", tp * 100.0 / (n-skip));

	// for(int i = 0; i < b; i++)
	// 	cout<<cnt[i]<<" ";
	// cout<<endl;
	return 0;
}
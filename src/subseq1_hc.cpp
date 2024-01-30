#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <map>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <inttypes.h>


using namespace std;

double f_max[100][100][101];
double f_min[100][100][101];
bool h[100][100][101];

map<char, int> dict;


double A[50][100][4][100];
int B1[50][100][4][100];
int B2[50][100][4][100];
int C[50][100][4];

int b, p, r;

void init(string path)
{
	dict['A'] = 0;
	dict['C'] = 1;
	dict['G'] = 2;
	dict['T'] = 3;

	for(int ri = 0; ri < r; ri++)
	{
		string table_filename = path + "/" + to_string(ri+1);
    	FILE *filein = fopen(table_filename.c_str(), "r");
    
	    for(int i = 0; i < b; i++)
			for(int j = 0; j < 4; j++)
		    	for(int q = 0; q < p; q++)
					(void)!fscanf(filein, "%lf ", &A[ri][i][j][q]);

	    for(int i = 0; i < b; i++)
			for(int j = 0; j < 4; j++)
		    	for(int q = 0; q < p; q++)
					(void)!fscanf(filein, "%d,%d", &B1[ri][i][j][q], &B2[ri][i][j][q]);

	    for(int i = 0; i < b; i++)
			for(int j = 0; j < 4; j++)
		    	(void)!fscanf(filein, "%d", &C[ri][i][j]);
	}
}

pair<int, double> DP(string s, int ri)
{
	int len = s.length();
	int del = len - b;

	memset(h, 0, sizeof(h));

	f_max[0][0][0] = f_min[0][0][0] = 0;
	h[0][0][0] = 1;
	for(int i = 1; i <= len; i++)
	{
		f_max[i][0][0] = 0;
		f_min[i][0][0] = 0;	    
		for(int j = 1; j < p; ++j)
	    {
			f_max[i][0][j] = -1e15;
			f_max[i][0][j] = 1e15;
	    }
		h[i][0][0] = 1;
	}

	int d1 = C[ri][0][dict[s[0]]];
	f_min[1][1][d1] = f_max[1][1][d1] = B2[ri][0][dict[s[0]]][d1] * A[ri][0][dict[s[0]]][d1];
	h[1][1][d1] = 1;

	double v1, v2, v3, v4;
	for(int i = 2; i <= len; i++)
	{
		int minj = max(1, i - del);
		int maxj = min(i, b);

		for(int j = minj; j <= maxj; j++)
		{
			int now = dict[s[i - 1]];
			int v = C[ri][j-1][now];

			for(int k = 0; k < p; k++)
			{
				int z = (k + v) % p;
				if(h[i][j][z] == 0)
				{
					f_min[i][j][z] = 1e15;
					f_max[i][j][z] = -1e15;
				}

				if(h[i-1][j][z] == 1)
				{
					f_min[i][j][z] = min(f_min[i][j][z], f_min[i-1][j][z]);
					f_max[i][j][z] = max(f_max[i][j][z], f_max[i-1][j][z]);
					h[i][j][z] = 1;
				}

				if(h[i-1][j-1][k])
				{
					if(B1[ri][j-1][now][z] == -1)
					{
						v1 = -f_min[i-1][j-1][k];
						v2 = -f_max[i-1][j-1][k];

						if(B2[ri][j-1][now][z] == -1)
						{
							v1 -= A[ri][j-1][now][z];
							v2 -= A[ri][j-1][now][z];
						}
						else
						{
							v1 += A[ri][j-1][now][z];
							v2 += A[ri][j-1][now][z];
						}

						if(v1 > f_max[i][j][z])
							f_max[i][j][z] = v1;
						if(v2 < f_min[i][j][z])
							f_min[i][j][z] = v2;
					}

					else
					{
						v1 = f_min[i-1][j-1][k];
						v2 = f_max[i-1][j-1][k];

						if(B2[ri][j-1][now][z] == -1)
						{
							v1 -= A[ri][j-1][now][z];
							v2 -= A[ri][j-1][now][z];
						}
						else
						{
							v1 += A[ri][j-1][now][z];
							v2 += A[ri][j-1][now][z];
						}

						if(v1 < f_min[i][j][z])
							f_min[i][j][z] = v1;
						if(v2 > f_max[i][j][z])
							f_max[i][j][z] = v2;
					}

					// cout<<i<<" "<<j<<" "<<k<<" "<<v<<" "<<z<<endl;
					// cout<<h[i-1][j][z]<<" "<<h[i-1][j-1][k]<<endl;
					// cout<<f_min[i-1][j][z]<<" "<<f_max[i-1][j][z]<<" "<<f_min[i-1][j-1][k]<<" "<<f_max[i-1][j-1][k]<<" "<<f_min[i][j][z]<<" "<<f_max[i][j][z]<<endl;
					h[i][j][z] = 1;
				}
			}
		}
	}

	double ans = 0;

	pair<int, double> ret;

	for(int i = 0; i < p; i++)
		if(h[len][b][i])
		{
			ans = max(ans, abs(f_min[len][b][i]));
			ans = max(ans, f_max[len][b][i]);

			if(ans > 0)
			{
				ret.first = i;
				ret.second = ans;
				return ret;
			}
		}

	return ret;
}


int tp = 0;
int skip = 0;

void match(string s, string t)
{
	if(b > t.length())
	{
		skip++;
		return;
	}
	

	for(int i = 0; i < r; i++)
	{
		pair<int, double> v1 = DP(s, i);
		pair<int, double> v2 = DP(t, i);

		if(v1 == v2)
		{
			tp ++;
			break;
		}
	}
}


int main(int argc, const char * argv[])
{
	b = stoi(argv[1]);
	p = stoi(argv[2]);
	r = stoi(argv[3]);
	string path = argv[4];

	string s, t;
	init(path);

	int n = 0;

	while(cin>>s)
	{
		cin>>t;
		match(s, t);
		n++;
	}
	printf("%.2lf\n", tp * 100.0 / (n-skip));


	return 0;
}
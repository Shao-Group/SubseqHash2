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

map<char, int> dict;


int A[100][4][50];
int B1[100][4][50];
int B2[100][4][50];

int A3[100][4];

int combine1[100][4];
int combine2[100][4];
int combine3[100][4];

int C1[100][4];
int C3[100][4];

void init(int k, int d, string path)
{
	dict['A'] = 0;
	dict['C'] = 1;
	dict['G'] = 2;
	dict['T'] = 3;

	vector<int> possign;
	vector<int> sign;

	for(int i = 0; i < 4; i++)
		possign.push_back(i);
	for(int i = 0; i < d; i++)
		sign.push_back(i);

	if(d < 4)
		for(int i = d; i < 4; i++)
			sign.push_back(i%d);
		
	unsigned seed;
	mt19937 generator(rd());
	uniform_int_distribution<int> distribution(1<<5, 1<<10);	

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				A[i][j][q] = distribution(generator);

			A3[i][j] = distribution(generator);
		}


		for(int q = 0; q < d; q++)
		{
			seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
			shuffle(possign.begin(), possign.end(), default_random_engine(seed));

			for(int j = 0; j < 4; j++)
			{
				B1[i][j][q] = (possign[j] % 2) ? 1: -1;
				B2[i][j][q] = (possign[j] / 2) ? 1: -1;
			}
		}

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(sign.begin(), sign.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
			C1[i][j] = sign[j];		

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(sign.begin(), sign.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
			C3[i][j] = sign[j];

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(possign.begin(), possign.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
		{
			combine1[i][j] = (possign[j] % 2) ? 1: -1;
			combine2[i][j] = (possign[j] / 2) ? 1: -1;
		}

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(possign.begin(), possign.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
		{
			combine3[i][j] = (possign[j] % 2) ? 1: -1;
		}
	}	

	// vector<int> samples, sampledpos;
	// if(k <= sample_size + 2)
	// {
	//     //not enough, select everything (except the middle position if k is odd)
	//     sample_size = (k >> 1) - 1;
	//     int i;
	//     for(i = 0; i < sample_size; i++)
	// 		sampledpos.push_back(i);
	    
	//     sample_size <<= 1;
	//     if(sample_size < k-2) 
	//     	++i;

	//     for(; i<=sample_size; i++)
	// 		sampledpos.push_back(i);
	// }
	// else
	// {
	//     int i = (k>>1) - 1;
	//     //candidates
	//     for(int j=0; j<i; ++j)
	// 		samples.push_back(j);
	    
	//     seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
	//     shuffle(samples.begin(), samples.end(), default_random_engine(seed));

	//     for(int j=0; j<(sample_size>>1); ++j)
	//     {
	// 		sampledpos.push_back(samples[j]);
	// 		sampledpos.push_back(k-3-samples[j]);
	// 	}
	// }
	
	vector<int> samples, sampledpos; //priority of positions in subsampling

    int sample_size = k;
    int half = sample_size>>1;
    for(int j=0; j<half; ++j)
		samples.push_back(j);
    
    seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
    shuffle(samples.begin(), samples.end(), default_random_engine(seed));

    for(int j=0; j<half; ++j)
    {
		sampledpos.push_back(samples[j]);
		sampledpos.push_back(k-1-samples[j]);
    }
	
	if(k&1)
		sampledpos.push_back(half);
	
	FILE *file = fopen(path.c_str(), "w");



	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				fprintf(file, "%d ", A[i][j][q]);
			fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				fprintf(file, "%d,%d\t", B1[i][j][q], B2[i][j][q]);
			fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", C1[i][j]);
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				fprintf(file, "%d ", A[i][3-j][q]);
			fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				fprintf(file, "%d,%d\t", B1[i][3-j][q], B2[i][3-j][q]);
			fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", C1[i][3-j]);
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	fprintf(file, "\n");

	for(int i = 0; i < k/2; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d,%d\t", combine1[i][j], combine2[i][j]);
		fprintf(file, "\n");
	}	

	if(k&1)
	{
		fprintf(file, "%d,%d\t", combine1[k/2][0], combine2[k/2][0]);
		fprintf(file, "%d,%d\t", combine1[k/2][1], combine2[k/2][1]);
		fprintf(file, "%d,%d\t", combine2[k/2][1], combine1[k/2][1]);
		fprintf(file, "%d,%d\t", combine2[k/2][0], combine1[k/2][0]);
		fprintf(file, "\n");		
	}

	for(int i = (k+1)/2; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d,%d\t", combine2[k-1-i][3-j], combine1[k-1-i][3-j]);
		fprintf(file, "\n");
	}	
	fprintf(file, "\n");
	fprintf(file, "\n");

	for(int i = 0; i < k/2; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", A3[i][j]);
		fprintf(file, "\n");
	}

	if(k&1)
	{
		fprintf(file, "%d ", A3[k/2][0]);
		fprintf(file, "%d ", A3[k/2][1]);
		fprintf(file, "%d ", A3[k/2][1]);
		fprintf(file, "%d ", A3[k/2][0]);
		fprintf(file, "\n");		
	}

	for(int i = (k+1)/2; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", A3[k-1-i][3-j]);
		fprintf(file, "\n");
	}
	fprintf(file, "\n");

	for(int i = 0; i < k/2; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d\t", combine3[i][j]);
		fprintf(file, "\n");
	}		
	if(k&1)
	{
		fprintf(file, "%d\t", combine3[k/2][0]);
		fprintf(file, "%d\t", combine3[k/2][1]);
		fprintf(file, "%d\t", combine3[k/2][1]);
		fprintf(file, "%d\t", combine3[k/2][0]);
		fprintf(file, "\n");		
	}
	for(int i = (k+1)/2; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d\t", combine3[k-1-i][3-j]);
		fprintf(file, "\n");
	}	
	fprintf(file, "\n");

	for(int i = 0; i < k/2; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", C3[i][j]);
		fprintf(file, "\n");
	}
	if(k&1)
	{
		fprintf(file, "%d ", C3[k/2][0]);
		fprintf(file, "%d ", C3[k/2][1]);
		fprintf(file, "%d ", C3[k/2][1]);
		fprintf(file, "%d ", C3[k/2][0]);
		fprintf(file, "\n");		
	}
	for(int i = (k+1)/2; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", C3[k-1-i][3-j]);
		fprintf(file, "\n");
	}
	fprintf(file, "\n");

	// fprintf(file, "%d\t", sample_size);
	for(int i = 0; i < sample_size; i++)
		fprintf(file, "%d\t", sampledpos[i]);
	fprintf(file, "\n");
}

int main(int argc, const char * argv[])
{
	int k = stoi(argv[1]);
	int p = stoi(argv[2]);
	string save = argv[3];

	init(k, p, save);

	
	return 0;
}

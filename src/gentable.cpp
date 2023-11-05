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

int sample_size;

random_device rd;

map<char, int> dict;


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
	uniform_int_distribution<int64_t> distribution((int64_t)1<<45, (int64_t)1<<55);	

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
			{
				A[i][j][q] = distribution(generator);
				revA[i][j][q] = distribution(generator);
			}

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
			seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
			shuffle(possign.begin(), possign.end(), default_random_engine(seed));

			for(int j = 0; j < 4; j++)
			{
				revB1[i][j][q] = (possign[j] % 2) ? 1: -1;
				revB2[i][j][q] = (possign[j] / 2) ? 1: -1;
			}
		}

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(sign.begin(), sign.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
			C1[i][j] = sign[j];		

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(sign.begin(), sign.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
			C2[i][j] = sign[j];

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
				fprintf(file, "%" PRId64 " ", A[i][j][q]);
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
				fprintf(file, "%" PRId64 " ", revA[i][j][q]);
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
				fprintf(file, "%d,%d\t", revB1[i][j][q], revB2[i][j][q]);
			fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", C2[i][j]);
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d,%d\t", combine1[i][j], combine2[i][j]);
		fprintf(file, "\n");
	}	
	fprintf(file, "\n");
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%" PRId64 " ", A3[i][j]);
		fprintf(file, "\n");
	}
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d\t", combine3[i][j]);
		fprintf(file, "\n");
	}	
	fprintf(file, "\n");

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			fprintf(file, "%d ", C3[i][j]);
		fprintf(file, "\n");
	}
	fprintf(file, "\n");

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

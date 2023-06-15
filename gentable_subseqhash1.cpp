#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <climits>
#include <map>
#include <cstring>
#include <unordered_set>

using namespace std;
random_device rd;
map<char, int> dict;

double A[100][4][100];
int B1[100][4][100];
int B2[100][4][100];
int C[100][4];

ofstream tout;
void init(int k, int d, string path)
{
	dict['A'] = 0;
	dict['C'] = 1;
	dict['G'] = 2;
	dict['T'] = 3;

	vector<int> pos;
	vector<int> possign;

	for(int i = 0; i < 4; i++)
		possign.push_back(i);

	for(int i = 0; i < d; i++)
		pos.push_back(i);
	if(d < 4)
		for(int i = d; i < 4; i++)
			pos.push_back(i%d);

	unsigned seed;
  	default_random_engine generator(rd());
	uniform_real_distribution<double> distribution((int64_t)1<<30, (int64_t)1<<31);


	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			for(int q = 0; q < d; q++)
				A[i][j][q] = distribution(generator);

		seed = unsigned(chrono::system_clock::now().time_since_epoch().count());
		shuffle(pos.begin(), pos.end(), default_random_engine(seed));

		for(int j = 0; j < 4; j++)
			C[i][j] = pos[j];

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
	}

	ofstream tout(path.c_str());

	for(int i = 0; i < k; i++)
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				tout << setprecision(4) << fixed << A[i][j][q] << " ";
			tout<<endl;
		}
	tout<<endl;

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			for(int q = 0; q < d; q++)
				tout<<B1[i][j][q]<<","<<B2[i][j][q]<<" ";
			tout<<endl;
		}
		tout<<endl;
	}
	tout<<endl;

	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			tout<<C[i][j]<<" ";
		tout<<endl;
	}
	tout<<endl;
	tout<<endl;
}

int main(int argc, const char * argv[])
{
	int k, d, num;
	string path;

	k = stoi(argv[1]);
	d = stoi(argv[2]);
	num = stoi(argv[3]);
	path = argv[4];
	

	for(int i = 1; i <= num; i++)
		init(k, d, path + '/' + to_string(i));
	return 0;
}
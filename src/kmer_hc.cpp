#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <climits>
#include <map>

using namespace std;


map<char,int> dict;

uint64_t get_hash(string s, int pos, int len)
{
	uint64_t ret = 0;
	for(int i = 0; i < len; i++)
		ret = ret * 4 + dict[s[pos + i]];

	return ret;
}

vector<uint64_t> get_minimizers(string s, int k)
{
	int len = s.length();

	vector<uint64_t> mins;

	for(int i = 0; i <= len - k; i++)
	{
		uint64_t v = get_hash(s, i, k);
		mins.push_back(v);
	}

	return mins;
}

int tp = 0;
void minimizer_match(string s, string t, int k)
{
	vector<uint64_t> v1 = get_minimizers(s, k);
	vector<uint64_t> v2 = get_minimizers(t, k);

	for(auto i: v1)
		for(auto j: v2)
			if(i == j)
			{
				tp++;
				return;
			}
}


int main(int argc, const char * argv[])
{
	int k = stoi(argv[1]);

	string s, t;
	int n = 0;

    dict['A'] = 0;
    dict['C'] = 1;
    dict['G'] = 2;
    dict['T'] = 3;

	while(cin>>s)
	{
		cin>>t;
		minimizer_match(s, t, k);
		n++;
	}
	printf("%.2lf\n", tp * 100.0 / n);

	return 0;
}
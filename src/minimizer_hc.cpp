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

uint64_t getblock64_mini ( const uint64_t * p, int i )
{
    return p[i];
}

inline uint64_t ROTL64_mini ( uint64_t x, int8_t r )
{
    return (x << r) | (x >> (64 - r));
}

uint64_t fmix64_mini ( uint64_t k )
{
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;

    return k;
}

uint64_t MurmurHash3_x64_128_mini (const void * key, const uint64_t seed)
{
    const uint8_t * data = (const uint8_t*)key;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = 0x87c37b91114253d5;
    const uint64_t c2 = 0x4cf5ad432745937f;


    const uint8_t * tail = (const uint8_t*)(data);

    uint64_t k1 = 0;
    uint64_t k2 = 0;

    k1 ^= ((uint64_t)tail[ 7]) << 56;
    k1 ^= ((uint64_t)tail[ 6]) << 48;
    k1 ^= ((uint64_t)tail[ 5]) << 40;
    k1 ^= ((uint64_t)tail[ 4]) << 32;
    k1 ^= ((uint64_t)tail[ 3]) << 24;
    k1 ^= ((uint64_t)tail[ 2]) << 16;
    k1 ^= ((uint64_t)tail[ 1]) << 8;
    k1 ^= ((uint64_t)tail[ 0]) << 0;
    k1 *= c1; 
    k1 = ROTL64_mini(k1,31); 
    k1 *= c2; 
    h1 ^= k1;

    //----------
    // finalization

    h1 ^= 8; h2 ^= 8;

    h1 += h2;
    h2 += h1;

    h1 = fmix64_mini(h1);
    h2 = fmix64_mini(h2);

    h1 += h2;
    h2 += h1;

    return h1;
}

uint64_t get_value(string s, int pos, int len)
{
	uint64_t ret = 0;
	for(int i = 0; i < len; i++)
		ret = ret * 4 + dict[s[pos + i]];

	return ret;
}

uint64_t get_hash(string s, int pos, int len, uint64_t seed)
{
	uint64_t ret = 0;
	for(int i = 0; i < len; i++)
		ret = ret * 4 + dict[s[pos + i]];

	// return ret; //lexicographic
    // uint64_t m = 0xFFFFFFFFFFFFFFFF;

    // ret = (~ret + (ret << 21)) & m;
    // ret = ret ^ ret >> 24;
    // ret = (ret + (ret<<3) + (ret<<8)) & m;
    // ret = ret ^ ret >> 14;
    // ret = (ret + (ret<<2) + (ret<<4)) & m;
    // ret = ret ^ ret >> 28;
    // ret = (ret + (ret<<31)) & m;

    // //minimap2-style order
	// return ret;

	return MurmurHash3_x64_128_mini(&ret, seed);
}

uint64_t get_minimizers(string s, int k, uint64_t seed)
{
	int len = s.length();
	vector<pair<int, uint64_t>> ret;
	vector<uint64_t> mins;

	for(int i = 0; i <= len - k; i++)
	{
		uint64_t v = get_hash(s, i, k, seed);
		mins.push_back(v);
	}

	uint64_t minimum = LLONG_MAX;
	int pos = 0;
	for(int i = 0; i <= len - k; i++)
	{
		//cout<<i<<" "<<mins[i]<<endl;
		if(mins[i] < minimum)
		{	
			minimum = mins[i];
			pos = i;
		}
	}
	//cout<<minimum<<" "<<pos<<endl;
	//cout<<s.substr(pos, k)<<" ";
	return get_value(s, pos, k);
}

int tp = 0;
int skip = 0;
int repeat;
vector<uint64_t> randseeds;

void minimizer_match(string s, string t, int k)
{
	if(k > t.length())
	{
		skip++;
		return;
	}

	for(int i = 0; i < repeat; i++)
	{
		uint64_t v1 = get_minimizers(s, k, randseeds[i]);
		uint64_t v2 = get_minimizers(t, k, randseeds[i]);

		if(v1 == v2)
		{
			tp++;
			break;
		}
	}
}


int main(int argc, const char * argv[])
{
	int k = stoi(argv[1]);
	repeat = stoi(argv[2]);

	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_int_distribution<int64_t> distribution(0, (uint64_t)-1);

	for(int i = 0; i < repeat; i++)
		randseeds.push_back(distribution(generator));

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
	printf("%.2lf\n", tp * 100.0 / (n - skip));

	return 0;
}
#include "util.h"

kmer encode(const char* s, const int k)
{
    kmer enc = 0lu;

    for(int i=0; i<k; i+=1)
		enc = (enc << 2) | alphabetIndex(s[i]);
    
    return enc;
}

char* decode(const kmer enc, const int k, char* str)
{
    if(str == NULL)
		str = (char*)malloc(sizeof *str *k);

    kmer enc_copy = enc;

    for(int i=k-1; i>=0; i-=1)
    {
		str[i] = ALPHABET[enc_copy & 3];
		enc_copy >>= 2;
    }

    return str;
}

kmer revComp(kmer s, int k)
{
    kmer r = 0;
    for(int i=0; i<k; ++i)
    {
		r <<= 2;
		r |= (s & 3) ^ 3;
		s >>= 2;
    }
    return r;
}

ssh_index* index_build(std::vector<seed>& seeds)
{
	ssh_index* ht = new ssh_index;
	seedinfo tmp;

	int sz = seeds.size();
	for(int i = 0; i < sz; i++)
	{
		tmp.s = &seeds[i];

		if(ht->find(seeds[i].hashval) == ht->end())
			ht->emplace(seeds[i].hashval, std::vector<seedinfo>(1, tmp));
		else
			(*ht)[seeds[i].hashval].push_back(tmp);
	}

	return ht;
}


void index_get(ssh_index* ht, std::vector<seed>& seeds, std::vector<seedmatch>& matches)
{
	int sz = seeds.size();
	for(int i = 0; i < sz; i++)
	{
		ssh_index::const_iterator got = ht->find(seeds[i].hashval);
		if(got != ht->end())
		{
			for(seedinfo y: got->second)
				matches.push_back((seedmatch){y.s, &seeds[i]});
		}
	}	
}

void saveSeeds(const char* filename, int k, const std::vector<seed>& seeds)
{
    FILE* fout = fopen(filename, "wb");
    fwrite(&k, sizeof(k), 1, fout);
    int pos[2];//st, ed
    for(auto s : seeds)
    {
		fwrite(&(s.str), sizeof(kmer), 1, fout);
		pos[0] = s.st;
		pos[1] = s.ed;
		fwrite(pos, sizeof(int), 2, fout);
    }

    fclose(fout);
}

void saveSeedsPosition(const char* filename, const std::vector<seed>& seeds)
{
    FILE* fout = fopen(filename, "wb");
    uint64_t pos[2];//st, ed
    for(auto s : seeds)
    {
		fwrite(&(s.hashval), sizeof(int64_t), 1, fout);
		pos[0] = s.st;
		pos[1] = s.index;
		fwrite(pos, sizeof(uint64_t), 2, fout);
    }

    fclose(fout);
}

void saveSeedsStrPosition(const char* filename, const std::vector<seed>& seeds)
{
    FILE* fout = fopen(filename, "wb");
    uint64_t pos[2];//st, ed
    for(auto s : seeds)
    {
		fwrite(&(s.hashval), sizeof(int64_t), 1, fout);
		fwrite(&(s.str), sizeof(kmer), 1, fout);
		pos[0] = s.st;
		pos[1] = s.index;
		fwrite(pos, sizeof(uint64_t), 2, fout);
    }

    fclose(fout);
}

void loadSeeds(const char* filename, std::vector<seed> &seeds)
{
    FILE* fin = fopen(filename, "rb");
    size_t ret = 1;
    uint64_t st, index;
	int64_t hashval;

    while(ret == 1)
    {
		ret = fread(&hashval, sizeof(uint64_t), 1, fin);
		fread(&st, sizeof(uint64_t), 1, fin);
		fread(&index, sizeof(uint64_t), 1, fin);

		if(!ret)
			break;
		seed tmp;
		tmp.hashval = hashval;
		tmp.st = st;
		tmp.index = index;
		
		uint64_t j = ((uint64_t)1)<<63;
		tmp.ed = st + 63;
		while(j)
			if(index & j)
				break;
			else
			{
				j >>= 1;
				tmp.ed--;
			}

		seeds.push_back(tmp);
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}

void loadSeedsStr(const char* filename, std::vector<seed> &seeds)
{
    FILE* fin = fopen(filename, "rb");
    size_t ret = 1;
    uint64_t st, index;
	int64_t hashval;
	kmer kk;

    while(ret == 1)
    {
		ret = fread(&hashval, sizeof(uint64_t), 1, fin);
		fread(&kk, sizeof(kmer), 1, fin);
		fread(&st, sizeof(uint64_t), 1, fin);
		fread(&index, sizeof(uint64_t), 1, fin);

		if(!ret)
			break;
		seed tmp;

		tmp.str = kk;
		tmp.hashval = hashval;
		tmp.st = st;
		tmp.index = index;
		
		uint64_t j = ((uint64_t)1)<<63;
		tmp.ed = st + 63;
		while(j)
			if(index & j)
				break;
			else
			{
				j >>= 1;
				tmp.ed--;
			}

		seeds.push_back(tmp);
    }
    if(ferror(fin)){
	fprintf(stderr, "Error reading %s\n", filename);
    }
    fclose(fin);
}
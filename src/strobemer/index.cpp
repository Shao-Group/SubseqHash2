//
//  index.cpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#include "index.hpp"
#include <iostream>
#include <math.h>       /* pow */
#include <bitset>



/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from http://www.cse.yorku.ca/~oz/hash.html:

uint64_t hash(std::string kmer)
{
    unsigned long hash = 5381;
    int c;
    for (std::string::size_type i=0; i< kmer.length(); i++) {
        c = kmer[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    return hash;
}

/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from minimap2:sketch.c :
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}//hash64

uint64_t getblock64 ( const uint64_t * p, int i )
{
    return p[i];
}

inline uint64_t ROTL64 ( uint64_t x, int8_t r )
{
    return (x << r) | (x >> (64 - r));
}

uint64_t fmix64 ( uint64_t k )
{
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;

    return k;
}

uint64_t MurmurHash3_x64_128 (const void * key, const uint64_t seed)
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
    k1 = ROTL64(k1,31); 
    k1 *= c2; 
    h1 ^= k1;

    //----------
    // finalization

    h1 ^= 8; h2 ^= 8;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    return h1;
}


static unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table

static inline uint64_t kmer_to_uint64(std::string &kmer, uint64_t kmask)
{
    uint64_t bkmer = 0;

    for (char i : kmer) {
        int c = seq_nt4_table[(uint8_t)i];
        bkmer = (bkmer << 2 | c) & kmask;

    }
    return bkmer;
}


mers_vector construct_flat_vector_three_pos(three_pos_index &tmp_index, uint64_t &unique_elements){
    mers_vector flat_vector;
    for (auto &it : tmp_index)  {
        for (auto &t : it.second) // it.second is the vector of k-mers, t is a tuple
        {
            flat_vector.push_back(t);
        }
    }
    //    flat_array sort
    std::sort(flat_vector.begin(), flat_vector.end());

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;
    unique_elements = 1;
    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k != prev_k){
            unique_elements ++;
        }
        prev_k = curr_k;
    }

    return flat_vector;
}


void index_vector_three_pos(mers_vector  &flat_vector, kmer_lookup &mers_index){
//    robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index;
    uint64_t offset = 0;
    uint64_t prev_offset = 0;
    unsigned int count = 0;

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int, unsigned short int, unsigned short int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;

    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k == prev_k){
            count ++;
        }
        else {
            std::tuple<uint64_t, unsigned int> s(prev_offset, count);
            mers_index[prev_k] = s;
            count = 1;
            prev_k = curr_k;
            prev_offset = offset;
        }
        offset ++;
    }

    // last k-mer
    std::tuple<uint64_t, unsigned int> s(prev_offset, count);
    mers_index[curr_k] = s;

//    return mers_index;
}




// initialize queue and current minimum and position
static inline void initialize_window(std::vector<uint64_t> &string_hashes, std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, int w_min, int w_max, int k){
    robin_hood::hash<std::string> robin_hash;
    for (int i = w_min; i < w_max; i++) {
//        auto strobe = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe, kmask);
//        q.push_back(bstrobe);
        uint64_t hv = string_hashes[i];
        q.push_back(hv);
        if (hv < q_min_val) {
            q_min_val = hv;
            q_min_pos = i;
        }
    }
}

// update queue and current minimum and position
static inline void update_window(std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, uint64_t new_strobe_hashval, int w_min, int w_max, int i, int seq_length ){
    uint64_t popped_val;
    popped_val = q.front();
    q.pop_front();
    q.push_back(new_strobe_hashval);
//    std::cout << "LOL " << popped_val << " " << q_min_val  << std::endl;
    if (popped_val == q_min_val){ // we popped the minimum value, find new brute force
//        std::cout << "OK "  << std::endl;
        q_min_val = UINT64_MAX;
//        q_min_pos = -1;
//        q_min_pos = i+w_min;
        q_min_pos = ((i+w_min + seq_length) - abs(i+w_min - seq_length))/2;
//        std::cout << "OK " << q_min_val  << std::endl;
//        std::cout << "OK " << q_min_pos  << std::endl;
        for (int j = 0; j <= q.size()-1; j++) {
//            std::cout << "(" << q[j] << " " << j << " " << i + w_min << "), ";
            if (q[j] < q_min_val) {
                q_min_val = q[j];
                q_min_pos = i + w_min + j + 1;
            }
        }
//        std::cout << "Changed: " << q_min_pos << " " << q_min_val << std::endl;
    }
    else if ( new_strobe_hashval < q_min_val ) { // the new value added to queue is the new minimum
        q_min_val = new_strobe_hashval;
        q_min_pos = i + w_max;
    }
}



static inline void make_string_to_hashvalues2(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k, uint64_t kmask, uint64_t seed) {
    robin_hood::hash<uint64_t> robin_hash;
//    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
    unsigned int hash_count = 0;
    int l;
    uint64_t x = 0;
    for (int i = l = 0; i < seq.length(); i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
//                uint64_t hash_k = x;
//                uint64_t hash_k = robin_hash(x);
//                uint64_t hash_k = hash64(x, kmask);
                uint64_t hash_k = MurmurHash3_x64_128(&x, seed);
                string_hashes.push_back(hash_k);
                pos_to_seq_choord.push_back( i - k + 1);
//                pos_to_seq_choord[hash_count] = i - k + 1;
                hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
}



static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
//    unsigned int min_pos;
//    min_pos = -1;
    for (auto i = w_start; i <= w_end; i++) {
//        uint64_t res = (strobe_hashval + string_hashes[i]) & q ;
        uint64_t res = (strobe_hashval ^ string_hashes[i]) ;
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}

mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index)
{
    mers_vector kmers;
    int l;
    int i;
    uint64_t mask=(1ULL<<2*k) - 1;
    uint64_t x = 0;
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    int cnt = 0;
    for (int i = l = 0; i <= seq.length()-1; i++) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = x;
                std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_k, ref_index, i-k+1, i-k+1, i-k+1);
                kmers.push_back(s);
                cnt ++;
//                if ((cnt % 10000000) == 0 ){
//                    std::cout << cnt << " kmers created." << std::endl;
//                }
            }
        }
        else {
            l = 0, x = 0; // if there is an "N", restart
        }

    }
    return  kmers;
}

mers_vector seq_to_randstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, uint64_t seed)
{
    mers_vector randstrobes2;

    if (seq.length() < w_max) {
        return randstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, seed);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


       // auto strobe1 = seq.substr(i, k);
       // std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}


mers_vector seq_to_randstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, uint64_t seed)
{
    mers_vector randstrobes3;

    if (seq.length() < 2*w_max) {
        return randstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
//    std::vector<uint64_t> pos_to_seq_choord;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, seed);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return randstrobes3;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " randstrobes created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        unsigned int strobe_pos_next1;
        uint64_t strobe_hashval_next1;
        unsigned int strobe_pos_next2;
        uint64_t strobe_hashval_next2;

        if (i + 2*w_max <= seq_length - 1){
            unsigned int w1_start = i+w_min;
            unsigned int w1_end = i+w_max;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

            unsigned int w2_start = i+w_max + w_min;
            unsigned int w2_end = i+2*w_max;
            uint64_t conditional_next = (strobe_hash ^ strobe_hashval_next1);
            get_next_strobe(string_hashes, conditional_next, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
        }
        else if ((i + 2*w_min + 1 < seq_length) && (seq_length <= i + 2*w_max) ){
//            unsigned int w_start = i+w_min;
//            unsigned int w_end = seq_length -1;
//            uint64_t strobe_hash;
//            strobe_hash = string_hashes[i];
//            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w_start, w_end, q);

            int overshot;
            overshot = i + 2*w_max - seq_length;
            unsigned int w1_start = i+w_min;
            unsigned int w1_end = i+w_max - overshot/2;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

            unsigned int w2_start = i+w_max - overshot/2 + w_min;
            unsigned int w2_end = i+2*w_max - overshot;
            uint64_t conditional_next = (strobe_hash ^ strobe_hashval_next1);
            get_next_strobe(string_hashes, conditional_next, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
        }
        else{
            return randstrobes3;
        }

//        uint64_t hash_randstrobe3 = (((strobe_hash << k) ^ strobe_hashval_next1) << k) ^ strobe_hashval_next2;
        uint64_t hash_randstrobe3 = (strobe_hash/3) + (strobe_hashval_next1/4) + (strobe_hashval_next2/5);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 =  pos_to_seq_choord[strobe_pos_next1]; //seq_pos_strobe1 + (strobe_pos_next1 - i); //
        unsigned int seq_pos_strobe3 =  pos_to_seq_choord[strobe_pos_next2]; //seq_pos_strobe1 + (strobe_pos_next2 - i); //
//        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << pos_to_seq_choord.size() << std::endl;

        // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
        if (strobe_pos_next2 ==  seq_length){
//            std::cout << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
            seq_pos_strobe3 = seq_length-1;
        }
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
        randstrobes3.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next1 - (i+k), ' ') << std::string(k, 'X') << std::string(strobe_pos_next2 - strobe_pos_next1 - k, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_length << std::endl;
    }
    return randstrobes3;
}



mers_vector seq_to_minstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, uint64_t seed)
{
    mers_vector minstrobes2;
    if (seq.length() < w_max) {
        return minstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, seed);
    unsigned int seq_length = string_hashes.size();
//    unsigned int seq_length = seq.length();

    // initialize the deque
    std::deque <uint64_t> q;
    uint64_t q_min_val = UINT64_MAX;
    int q_min_pos = -1;
    initialize_window(string_hashes, q, q_min_val, q_min_pos, w_min, w_max, k);

//    std::cout << seq << std::endl;

    // create the minstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);

        uint64_t bstrobe = string_hashes[i];
//        uint64_t hash_minstrobe2 = (bstrobe << k) ^ q_min_val;
        uint64_t hash_minstrobe2 = (bstrobe/2) + (q_min_val/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[q_min_pos];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_minstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        minstrobes2.push_back(s);

        // update queue and current minimum and position
        if (i + w_max <= seq_length){
//            auto new_strobe = seq.substr(i + w_max, k);
//            uint64_t new_bstrobe = kmer_to_uint64(new_strobe, kmask);
//            uint64_t new_strobe_hashval = hash64(new_bstrobe, kmask);
            uint64_t new_strobe = string_hashes[i + w_max];

            update_window(q, q_min_val, q_min_pos, new_strobe, w_min, w_max, i, seq_length );
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length < i + w_max) ){
//            uint64_t new_strobe_hashval =  UINT64_MAX;
            uint64_t new_strobe = UINT64_MAX;
            update_window(q, q_min_val, q_min_pos, new_strobe, w_min, w_max, i, seq_length );
        }
        else{
            return minstrobes2;
        }

//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(q_min_pos - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return minstrobes2;
}



mers_vector seq_to_hybridstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, uint64_t seed)
{
    mers_vector hybridstrobes2;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < w_max) {
        return hybridstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, seed);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;

    // initialize deque1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k);

    // initialize deque2
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k);

    // initialize deque3
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i];

        int strobe2_pos;
        uint64_t strobe2_val;
        unsigned int r =  bstrobe % 3;
        if (r == 0){
            strobe2_val = q1_min_val;
            strobe2_pos = q1_min_pos;
        }
        else if (r == 1){
            strobe2_val = q2_min_val;
            strobe2_pos = q2_min_pos;
        }
        else{
            strobe2_val = q3_min_val;
            strobe2_pos = q3_min_pos;
        }

        uint64_t hash_hybridstrobe2;
//        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
        hash_hybridstrobe2 = (bstrobe/2) + (strobe2_val/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe2_pos];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_hybridstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        hybridstrobes2.push_back(s);


        // update queue and current minimum and position
        if (i + w_max < seq_length){
//            auto new_strobe1 = seq.substr(i + w_min+x, k);
//            uint64_t new_bstrobe1 = kmer_to_uint64(new_strobe1, kmask);
//            uint64_t new_strobe_hashval1 = hash64(new_bstrobe1, kmask);
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);

//            auto new_strobe2 = seq.substr(i + w_min+2*x, k);
//            uint64_t new_bstrobe2 = kmer_to_uint64(new_strobe2, kmask);
//            uint64_t new_strobe_hashval2 = hash64(new_bstrobe2, kmask);
            uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);

//            auto new_strobe3 = seq.substr(i + w_max, k);
//            uint64_t new_bstrobe3 = kmer_to_uint64(new_strobe3, kmask);
//            uint64_t new_strobe_hashval3 = hash64(new_bstrobe3, kmask);
            uint64_t new_strobe_hashval3 = string_hashes[i + w_max];
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval, w_min, w_min+x, i, seq_length);
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i, seq_length);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min+x < seq_length) && (seq_length <= i + w_min+2*x)  ){
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i, seq_length);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min+2*x < seq_length) && (seq_length <= i + w_max) ){
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);
            uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else{
            return hybridstrobes2;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return hybridstrobes2;
}




mers_vector seq_to_hybridstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, uint64_t seed)
{
    mers_vector hybridstrobes3;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < 2*w_max) {
        return hybridstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, seed);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;

    // initialize deque1 window1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k);

    // initialize deque2 window1
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k);

    // initialize deque3 window1
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k);

    // initialize deque1 window2
    std::deque <uint64_t> q4;
    uint64_t q4_min_val = UINT64_MAX;
    int q4_min_pos = -1;
    initialize_window(string_hashes, q4, q4_min_val, q4_min_pos, w_max+ w_min, w_max + w_min+ x, k);

    // initialize deque2 window2
    std::deque <uint64_t> q5;
    uint64_t q5_min_val = UINT64_MAX;
    int q5_min_pos = -1;
    initialize_window(string_hashes, q5, q5_min_val, q5_min_pos, w_max+ w_min+x, w_max+ w_min+2*x, k);

    // initialize deque3 window2
    std::deque <uint64_t> q6;
    uint64_t q6_min_val = UINT64_MAX;
    int q6_min_pos = -1;
    initialize_window(string_hashes, q6, q6_min_val, q6_min_pos, w_max+ w_min+2*x, 2*w_max, k);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i];

        int strobe2_pos;
        uint64_t strobe2_val;
        unsigned int r =  bstrobe % 3;
        if (r == 0){
            strobe2_val = q1_min_val;
            strobe2_pos = q1_min_pos;
        }
        else if (r == 1){
            strobe2_val = q2_min_val;
            strobe2_pos = q2_min_pos;
        }
        else{
            strobe2_val = q3_min_val;
            strobe2_pos = q3_min_pos;
        }

        uint64_t hash_hybridstrobe2;
//        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
        hash_hybridstrobe2 = (bstrobe/3) + (strobe2_val/4);

        int strobe3_pos;
        uint64_t strobe3_val;
        unsigned int r2 =  hash_hybridstrobe2 % 3;
        if (r2 == 0){
            strobe3_val = q4_min_val;
            strobe3_pos = q4_min_pos;
        }
        else if (r2 == 1){
            strobe3_val = q5_min_val;
            strobe3_pos = q5_min_pos;
        }
        else{
            strobe3_val = q6_min_val;
            strobe3_pos = q6_min_pos;
        }

        uint64_t hash_hybridstrobe3;
//        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
        hash_hybridstrobe3 = hash_hybridstrobe2 + (strobe3_val/5);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe2_pos];
        unsigned int seq_pos_strobe3 = pos_to_seq_choord[strobe3_pos];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_hybridstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
        hybridstrobes3.push_back(s);


        // update queue and current minimum and position
        uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
        update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);

        uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
        update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);

        uint64_t new_strobe_hashval3 = string_hashes[i + w_max];
        update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i, seq_length);

        if (i + 2*w_max < seq_length){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);

            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval5, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval6 = string_hashes[i + 2*w_max];
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval6, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min + 1 < seq_length) && (seq_length <= i + w_max + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval, w_max + w_min, w_max + w_min+x, i, seq_length);
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+x < seq_length) && (seq_length <= i + w_max + w_min+2*x)  ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+2*x < seq_length) && (seq_length <= i + 2*w_max) ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval5, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else{
            return hybridstrobes3;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k) - 1, ' ') << std::string(k, 'X') << std::string(strobe3_pos - strobe2_pos - k, ' ') << std::string(k, 'X') << std::endl;

//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return hybridstrobes3;
}









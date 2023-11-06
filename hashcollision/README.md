
### Introduction
Codes and instructions for comparing the Probability of Hash Collision between SubseqHash2, SubseqHash, and minimizer. 

### Usage
- To simulate pairs of strings:
  ```
  g++ -o simulate simulate.cpp -std=c++11
  ./simulate k n e data_directory
  ```
    - Run `./simulate k n e data_directory` to simulate pairs of strings. `k` is the length of random generated string in each pair. `n` is the number of pairs. `e` is maximal edit distance between two strings in each pair. `data_directory` is the path of the directory that save `e` different files named from 1 to e. In file `i`, there are 2 * n lines, every two lines contain two strings, the first one is of lenth `k` and the edit distance of the two strings is `i`.
    - 10 sample data files are in folder `./data`, each file has 100 pairs of strings with k = 20. The 10 files are named from 1 to 10 which is the edit distance between the string pairs in that file.

- To generate random tables for SubseqHash/SubseqHash2
  ```
  ../src/gentable.out k d ./table_saving_path
  ../src/gentable_subseqhash1.out k d t ./table_saving_directory
  ```
    - ` ../src/gentable.out k d ./table_saving_path` generates random tables for SubseqHash2, where `k` is the length of subsequence, `d` is a parameter, and `./table_saving_path` is the path of the file saving the tables. A sample file is at `./table_subseq2/15_11` with k=15 and d = 11.
    - `../src/gentable_subseqhash1.out k d t ./table_saving_directory`generates random tables for SubseqHash, where `k` is the length of subsequence, `d` is a parameter, `t` is the number of repeating, and `./table_saving_directory` is the directory saving the files. There will be `t` files in that directory, and each one has one set of random tables for SubseqHash.  A sample folder is at `./table_subseq1/15_11` with k = 15, d = 11 and 10 files with different tables in it.
	

- To get the probability of hash collision:
  - Run `../src/subseq2_hc.out k d t table_file < data_file` to get the probability of hash collision of SubseqHash2 on the pairs of strings in the input file `data_file`. `k` and `d` are the parameters for SubseqHash2 tables. The table file must be generated with the same parameters. `t` is the time of repeats. The output format is single number of `probability of hash collision` * 100.  An example usage: `../src/subseq2_hc.out 15 11 1 ./table_subseq2/15_11 <./data/1` 
  -	 Run `../src/subseq1_hc.out k d t table_directory < data_file` to get the probability of hash collision of SubseqHash, an example usage:	`../src/subseq1_hc.out 15 11 1 ./table_subseq1/15_11 < ./data/1`
  -	 Run `../src/minimizer_hc.out k t < data_file` to get the probability of hash collision of minimizer repeating `t` times with length `k`, an example usage:	`../src/minimizer_hc.out 15 1 < ./data/1`
  -	 Run `../src/kmer_hc.out k < data_file` to get the probability of hash collision of using all kmers with length `k`, an example usage:	`../src/kmer_hc.out 15 < ./data/1`
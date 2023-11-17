
### Introduction
Codes and instructions for using SubseqHash2 for seeding. 

### Usage

- To generate random tables for SubseqHash2
  ```
  ../src/gentable.out k d ./table_saving_path
  ```
    - ` ../src/gentable.out k d ./table_saving_path` generates random tables for SubseqHash2, where `k` is the length of subsequence, `d` is a parameter, and `./table_saving_path` is the path of the file saving the tables. A sample file is at `./table_subseq2/15_11` with k=15 and d = 11.
	

- To get the probability of hash collision:
  - Run `../src/sample.out n k d t table_file data_file output_file` to do seeding for sequences in data file using SubseqHash2. `k` and `d` are the parameters for SubseqHash2 tables. The table file must be generated with the same parameters. `t` is the time of repeats. The output format is single number of `probability of hash collision` * 100.  An example usage: `../src/subseq2_hc.out 15 11 1 ./table_subseq2/15_11 <./data/1` 
  - The output includes the seeds for all the sequences from the input file. For each sequence, seeds are categoried into each order from 1 to t. For each order, one line contains one seed including the score and optimal subsequence.
  - `../src/sample.out 20 15 11 5 ./table_subseq2/15_11 sample_sequences output`

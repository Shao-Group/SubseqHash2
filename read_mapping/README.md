### Introduction
Codes and instructions for comparing SubseqHash2 (SIMD version), minimizer, syncmer on mapping reads to reference genome. 

### Usage
    
- To generate random tables for SubseqHash2
  ```
  ../src/gentable_rc.out k d ./table_saving_path
  ```
	 - ` ../src/gentable_rc.out k d ./table_saving_path` generates random tables for SubseqHash2, where `k` is the length of subsequence, `d` is a parameter, and `./table_saving_path` is the path of the file saving the tables. A sample file is at `./table_subseq2/25_11` with k=25 and d = 11. Slightly from the two previous expeirments, here we use the variant SubseqHash2r which requires to generate random tables that can produce the same seed for a string and its reverse complement.

- Data format:
  There are two data files. The first one is the reference genome of fasta format. It only has one sequence. There is an example at `./ref/sample`. 
The second file is the read dataset with ground-truth. Each read takes two lines, the first line is the sequence. 
The second line is the ground-truth alignment, the i-th number `$x_i$` in the second line means the i-th char of the read should be aligned to the position `$x_i$`. 
There is also an example `./reads/sample`.

- Get results of SubseqHash2:
  - Seeding for reference genome `../src/refseeding_sub2.out n k d t table_file ref_file`
    `n`, `k`, `d` are the parameter for SubseqHash2. ``t` is the number of repeating times. The table file must be generated with the same parameters. 
`ref_file` is the name of the reference genome. The output seeds are store at ./refseed/.

  - Seeding for reads `../src/readseeding_sub2.out n k d t table_file read_file`.
`read_file` is the name of the reference genome. The output seeds are store at ./readseed/.

  - Get seed matching results `../src/readmapping_sub2.out n k d t read_file`. It will take the seeds from `./refseed` and `./readseed`. The output format is parameters, 
number of seed-matches, number of true seed-matches, ratio of true seed-matches, the average ratio of true seed-matches per read, 
number of 200bp segments, number of 200bp segments with true seed-matches, sensitivity of segments with true seed-matches.

  - An example on sample data:
    ```
    ../src/refseeding_sub2.out 30 25 11 2 ./table_subseq2/25_11 ./ref/sample
    ../src/readseeding_sub2.out 30 25 11 2 ./table_subseq2/25_11 ./reads/sample
    ../src/readmapping_sub2.out 30 25 11 2 ./reads/sample
    ```

- Get results of Other methods:
  - An example of minimizer on sample data. The arguments are k w ref_file read_file.
      ```
      ../src/readmapping_mini.out 25 5 ./ref/sample ./reads/sample
      ```
  - An example of syncmer on sample data. The arguments are k s ref_file read_file pos_1 pos_2 .. pos_s.
      ```
      ../src/readmapping_sync.out 20 10 ./ref/sample ./reads/sample 0 10
      ```

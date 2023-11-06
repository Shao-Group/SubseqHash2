### Introduction
Codes and instructions for comparing SubseqHash2 (SubseqHash2r variant), SubseqHash2w, SubseqHash, minimizer, syncmer and strobemer on mapping reads to reference genome. 

### Usage
    
- To generate random tables for SubseqHash/SubseqHash2
  ```
  ../src/gentable_rc.out k d ./table_saving_path
  ../src/gentable_subseqhash1.out k d t ./table_saving_directory
  ```

	 - ` ../src/gentable_rc.out k d ./table_saving_path` generates random tables for SubseqHash2, where `k` is the length of subsequence, `d` is a parameter, and `./table_saving_path` is the path of the file saving the tables. A sample file is at `./table_subseq2/25_11` with k=25 and d = 11. Slightly from the two previous expeirments, here we use the variant SubseqHash2r which requires to generate random tables that can produce the same seed for a string and its reverse complement.
	- `../src/gentable_subseqhash1.out k d t ./table_saving_directory`generates random tables for SubseqHash, where `k` is the length of subsequence, `d` is a parameter, `t` is the number of repeating, and `./table_saving_directory` is the directory saving the files. There will be `t` files in that directory, and each one has one set of random tables for SubseqHash.  A sample folder is at `./table_subseq1/25_11` with k = 25, d = 11 and 10 files with different tables in it.

- Data format:
  There are two three files. The first one is the reference genome of fasta format. It can have multiple sequences. There is an example at `./ref/sample`. The second file only has the names of all the genome sequences, one name per line which is the same as the first world in each description line of fasta file, and an example file is at `./ref/samplelist`. The last file is the read dataset with ground-truth. Each read takes three lines, the first line is the sequence. The scecond line is the name of reference genome sequence that read should be mapped to (the name must be the same as the one in reference genome file). The third line is the ground-truth alignment, the i-th number `$x_i$` in the second line means the i-th char of the read should be aligned to the position `$x_i$`. There is also an example `./reads/sample`. If you want to use other datasets, suppose the name of dataset is X, the three files must be at `./ref/X`, `./ref/Xlist`, `./reads/X`.

- Get results of SubseqHash2:
  - Seeding for reference genome `../src/refseeding_sub2.out n k d t table_file X`
    `n`, `k`, `d` are the parameter for SubseqHash2. ``t` is the number of repeating times. The table file must be generated with the same parameters. `X` is the name of the dataset. The output seeds are store at ./refsub2/X.

  - Seeding for reads `../src/readseeding_sub2.out n k d t table_file X`. The output seeds are store at ./readsub2/X.

  - Get seed matching results `../src/readmapping_sub2.out n k d t X`. It will take the seeds from `./refsub2/X` and `./readsub2/X`. The output format is parameters, number of seed-matches, number of true seed-matches, the average ratio of true seed-matches, coverage of true seed-matches, coverage of false seed-matches, matching coverage of true seed-matches, the island size, the density of seeds.

  - An example on sample data:
    ```
    ../src/refseeding_sub2.out 30 25 11 2 ./table_subseq2/25_11 sample
    ../src/readseeding_sub2.out 30 25 11 2 ./table_subseq2/25_11 sample
    ../src/readmapping_sub2.out 30 25 11 2 sample
    ```

- Get results of Other methods:
  - An example of SubseqHash2w on sample data. The arguments are n k d t k0 table_file X.
    ```
    ../src/refseeding_substro.out 30 25 11 2 5 ./table_subseq2/25_11 sample
    ../src/readseeding_substro.out 30 25 11 2 5 ./table_subseq2/25_11 sample
    ../src/readmapping_substro.out 30 25 11 2 5 sample
    ```
  - An example of SubseqHash on sample data. The arguments are n k d t table_file X.
    ```
    ../src/refseeding_sub1.out 30 25 11 0 ./table_subseq1/25_11/0 sample
    ../src/refseeding_sub1.out 30 25 11 1 ./table_subseq1/25_11/1 sample
    ../src/readseeding_sub1.out 30 25 11 0 ./table_subseq1/25_11/0 sample
    ../src/readseeding_sub1.out 30 25 11 1 ./table_subseq1/25_11/1 sample
    ../src/readmapping_sub1.out 30 25 11 2 sample
    ```
  - An example of minimizer on sample data. The arguments are k w X.
      ```
      ../src/readmapping_mini.out 25 5 sample
      ```
  - An example of syncmer on sample data. The arguments are k s X pos_1 pos_2 .. pos_s.
      ```
      ../src/readmapping_sync.out 20 10 sample 0 10
      ```
  - An example of strobemer on sample data. The arguments are k, n, w_min, w_max, t, X.
      ```
      ../src/readmapping_stro.out 10 2 10 20 sample
      ```

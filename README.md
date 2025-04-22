### Introduction
SubseqHash2 employs a seeding approach that utilizes minimized subsequences as seeds. 
Comparing to SubseqHash, SubseqHash2 not only produces high-quality seeds for input genome sequences but also significantly improves running time by algorithmic innovation and applying SIMD instructions.
This repository includes the source code for four distinct experiments, each housed in its respective folder:
folder `hash_collision` contains code for computing the probability of hash collision of different methods.
Folder `alignment` has code for comparing different methods in pairwise sequence alignment. Folder `read_mapping` has code for comparing different methods in mapping reads to reference genome.
Folder `overlap_detection` provides code for comparing SubseqHash with other methods on the seeding step of overlap detection.

### Usage

The repository can be installed with following command. The detailed instruction of each experiment can be found in the corresponding folder.

```
git clone https://github.com/Shao-Group/subseqhash2.git
cd src
make product
``` 

To use SubseqHash2 for seed generation, please refer to the example in `./src/sample.cpp`. You can run the following commands:

- ` ../src/gentable.out k d ./table_saving_path` generates a random table for SubseqHash2, where `k` is the length of subsequence, `d` is a parameter, and `./table_saving_path` is the path of the file saving the tables. A sample file is at `./table/6_5` with k = 6 and d = 5.

- Run `../src/sample.out n k d t table_file data.fasta` to get the sequence alignment results using SubseqHash2. `n`, `k` and `d` are the parameters, `t` is the number of repeating times, `table_file` is the path to the previously generated table file (with matching parameters). The data file must be in FASTA format. An example is `../src/sample.out 10 6 5 6 ./table/6_5 ./sample.fa`. The output gives the starting position, hash value, and the actual subsequence for each generated seed.

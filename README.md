### Introduction
SubseqHash2 employs a seeding approach that utilizes minimized subsequences as seeds. 
Comparing to SubseqHash, SubseqHash2 not only produces high-quality seeds for input genome sequences but also demonstrates significantly improved running time.
This repository includes the source code for four distinct experiments, each housed in its respective folder:
folder `hashcollision` contains code for computing the probability of hash collision of different methods.
`alignment` has code for comparing different methods in pairwise sequence alignment. `readmaping` has code for comparing different methods in mapping reads to reference genome.
`overlap-detection` provides code for comparing SubseqHash with other methods on the seeding step of overlap detection.

### Usage

The repository can be installed with following command. The detailed instruction of each experiment can be found in the corresponding folder.

```
git clone https://github.com/Shao-Group/subseqhash2.git
cd src
make product
``` 
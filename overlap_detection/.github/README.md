### Introduction
Code for comparing SubseqHash2 with minimizer, syncmer, and strobemer on the seeding step of overlap detection. For this purpose, sequences that share at least one common seed are considered overlapping.

All the seeding programs take ordinary fasta files as input.
A sample file with 3 reads, [SRX533603-sample3.efa](../sample-reads/SRX533603-sample3.efa), is provided.

### Installation
```
git clone https://github.com/Shao-Group/SubseqHash2.git
cd SubseqHash2/overlap_detection
make product
```
### Usage
- To generate random tables used by SubseqHash2 and SubseqHash2w: Run `../src/gentable_rc.out k d tables/table_name` where `k` and `d` are parameters for SubseqHash2. Two sample tables are provided in [tables/](tables).

- All seeding programs by default use 15 threads to generate seeds in parallel (one thread per read), the number of threads used can be changed by modifying `NUMTHREADS` in [src/seedfactory.h](../src/seedfactory.h). 

- To generate seeds with SubseqHash2: Run `./genSubseq2Seeds.out readFile n k d randTableFile`. Example: `./genSubseq2Seeds.out sample-reads/SRX533603-sample3.efa 60 48 11 tables/48_11` which corresponds to the second left most point of the `SH2r 60 11 k` curve in Fig.5 of our paper. This call creates a directory `sample-reads/SRX533603-sample3.efa-locSeeds-n60-k48-d11/` containing all the seeds. SubseqHash2 generates `k` sets of seeds for each read.

- To generate seeds with minimizer: Run `./genMinimizerSeeds.out readFile w k`. Note that `w` here is the number of `k`-mers in a window. Example: `./genMinimizerSeeds.out sample-reads/SRX533603-sample3.efa 4 22` which corresponds to the second left most point of the `minimizer 25` curve. This call creates a directory `sample-reads/RX533603-sample3.efa-minimizerSeeds-w4-k22/` containing all the seeds. Two seed files per read are stored, one for the forward direction, the other for the reverse complement. 

- To generate seeds with syncmer: Run `./genSyncmerSeeds.out readFile k s` to select `k`-mers according to the location of their smallest `s`-mers. Example: `./genSyncmerSeeds.out sample-reads/SRX533603-sample3.efa 14 11` which corresponds to the second right most point of the `syncmer 14` curve. This call creates a directory `sample-reads/SRX533603-sample3.efa-syncmerSeeds-k14-s11/` containing all the seeds. Similar as minimizer, two seed files per read are stored.

- To generate seeds with strobemer: Run `./genStrobemerSeeds.out readFile k w w_min w_max seednum` to generate seeds with `w` strobes of length `k` each, the second strobe is picked from the window `[w_min, w_max]` relative to the beginning of the seed; the `seednum` parameter controls the number of repeats applied, namely, `2*seednum` sets of seeds with different hash functions are generated for each read. Example: `./genStrobemerSeeds.out sample-reads/SRX533603-sample3.efa 22 2 22 38 12` which corresponds to the left most point of the `strobemer 60 k` curve. This call creates a directory `sample-reads/SRX533603-sample3.efa-strobemerSeeds-k22-w2-w_min22-w_max38/`.

- To generate seeds with SubseqHash2w: Run `./genSubseq2strobeSeeds.out readFile n k d randTableFile w pre_k` to produce seeds with `w+1` parts -- a leading `pre_k`-mer followed by `w` subsequneces of length `k` picked from a window of size `n`. Example: `./genSubseq2strobeSeeds.out sample-reads/SRX533603-sample3,efa 44 32 11 tables/32_11 1 16` which corresponds to the second left most point of the `SH2w 60 11 k` curve. This call creates a directory `sample-reads/SRX533603-sample3.efa-subseq2strobeSeeds-n44-k32-d11-w1-pre_k16/`. SubseqHash2w generates `2k` sets of seeds for each read.

- To get overlapping pairs with the seeds generated using any of the above methods: Run `overlapBySeeds.out seedsDir k repeatId numReads fileExt` where `seedsDir` is the directory containing the seed files; `k` is the length of each seed; `repeatId` specifies seeds from which repeat is to be used; only the top `numReads` reads are considered for overlapping; and `fileExt` is the filename extension of the seed files, all possible choices are listed in the following table.
| Seeding Method | fileExt |
| --- | --- |
| SubseqHash2 | subseqseed2 |
| minimizer | mmseed |
| syncmer | syncmerseed |
| strobemer | strobemerseed |
| SubseqHash2w | ss2sseed |
For seeding methods without repeat (minimizer, syncmer), the `repeatId` should be 0; for others, it starts from 0. Each call creates a file `overlap-n{numReads}-r{repeatId}.all-pair` within the seed directory. In this file, each line represents an overlapping pair in the format `read_id1 read_id2 number_of_shared_seeds`.

- The `overlapFileOp.js` script can be used to take union/intersection of overlap detection results. It requires the [k8](https://github.com/attractivechaos/k8) library.
  - Example for obtaining the union of two result files (used for repeating SubseqHash2): `./overlapFileOp.js union sample-reads/SRX533603-sample3.efa-locSeeds-n60-k48-d11/overlap-n3-r{0,1}.all-pair > union2.all-pair`.
  - Example for obtaining the intersection of two files (used for comparing reported overlapping pairs versus the ground truth): `./overlapFileOp.js intersection sample-reads/SRX533603-sample3.truepairs sample-reads/SRX533603-sample3.efa-minimizerSeeds-w4-k22/overlap-n3-r0.all-pair` where the truepairs file is in the same format as result files except that the third column is the number of overlapping base pairs.


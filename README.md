# eMCDP

A tool built on top of mc-overlaps [1] (mcdp2 [2] predecessor), allowing for colocalization analysis on genomic windows. Written in C++.

## How to use this?

### Basic data set

```bash
make clean && make run_sample # to run the basic algorithm on basic data set
```

#### Running a large (|R|=2500, |Q|=5000) dataset

```bash
make clean && make run_simple_pvalue REF_PATH=data/02-synth-data/g24_8.ref.tsv QUERY_PATH=data/02-synth-data/g24_8.query.tsv CHR_SIZES_PATH=data/02-synth-data/g24_sizes.tsv OUTPUT_PATH=data/output/02-synth-data-g24_8.txt
```

Logs will be in the `data/output` directory.

### TODO

- figure out a way to implement merge two consecutive P_DP tables, implement it into a function (possibly can be done with `joint_logprobs`, but not sure, is next talking point)

> [1] Askar Gafurov, Broňa Brejová, Paul Medvedev,
> Markov chains improve the significance computation of overlapping genome annotations,
> Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i203–i211, https://doi.org/10.1093/bioinformatics/btac255

> [2] Gafurov, A., Vinař, T., Medvedev, P., Brejová, B. (2024). Efficient Analysis of Annotation Colocalization Accounting for Genomic Contexts. In: Ma, J. (eds) Research in Computational Molecular Biology. RECOMB 2024. Lecture Notes in Computer Science, vol 14758. Springer, Cham. https://doi.org/10.1007/978-1-0716-3989-4_3

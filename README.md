# GALA

A tool built on top of mc-overlaps [1] (mcdp2 [2] predecessor), allowing for colocalization analysis on genomic windows. Written in C++.

#### How to use this?

```
bash run.sh
```

### TODO

- finish implementing `eval_probs_single_chr_direct` from [eval_probs_single_direct_lm](https://github.com/fmfi-compbio/mc-overlaps/blob/bdc19acc5e0a8da200870442a8391cc0633392ef/src/simple_model.py#L78)
- implement `calculate_joint_pvalue` function from [joint_pvalue](https://github.com/fmfi-compbio/mc-overlaps/blob/bdc19acc5e0a8da200870442a8391cc0633392ef/src/simple_model.py#L78)
- figure out a way to implement merge two consecutive P_DP tables, implement it into a function

> [1] Askar Gafurov, Broňa Brejová, Paul Medvedev,
> Markov chains improve the significance computation of overlapping genome annotations,
> Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i203–i211, https://doi.org/10.1093/bioinformatics/btac255

> [2] Gafurov, A., Vinař, T., Medvedev, P., Brejová, B. (2024). Efficient Analysis of Annotation Colocalization Accounting for Genomic Contexts. In: Ma, J. (eds) Research in Computational Molecular Biology. RECOMB 2024. Lecture Notes in Computer Science, vol 14758. Springer, Cham. https://doi.org/10.1007/978-1-0716-3989-4_3

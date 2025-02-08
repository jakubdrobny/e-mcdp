# eMCDP - e(xtended) [MCDP2](https://github.com/fmfi-compbio/mcdp2)

A tool built on top of MCDP2 [1], allowing for colocalization analysis on genomic windows. Written in C++.

#### How to use this?

```
bash run.sh
```

### TODO

- check what is the difference between [eval_probs_single_chr_direct](https://github.com/fmfi-compbio/mc-overlaps/blob/master/src/simple_model.py#L78) and [eval_probs_single_chr_eigen](https://github.com/fmfi-compbio/mc-overlaps/blob/master/src/simple_model.py#L137)
- implement both if they are not the same, otherwise just the direct one
- implement `calculate_joint_pvalue` function from [joint_pvalue](https://github.com/fmfi-compbio/mc-overlaps/blob/master/src/helpers.py#L119)

> [1] Askar Gafurov, Broňa Brejová, Paul Medvedev,
> Markov chains improve the significance computation of overlapping genome annotations,
> Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i203–i211, https://doi.org/10.1093/bioinformatics/btac255

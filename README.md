# eMCDP

A tool built on top of mc-overlaps [1] (mcdp2 [2] predecessor), allowing for colocalization analysis on genomic windows. Written in C++.

## How to use this?

Make sure to have [make](https://www.gnu.org/software/make/manual/make.html) and [cmake](https://cmake.org/) installed. Then first install the dependencies
```bash
sudo apt install libgtest-dev googletest libomp-dev
```

And the you can install the program with following commands:
```bash
mkdir build && cd build
cmake ..
make
sudo make install
```

Now you have the `emcdp` executable installed and it can be run from the command line.

The program provides a set of flags to operate it:

- `--r <path-to-your-ref-intervals-file>` - REQUIRED, tells the program where to find the file with reference annotation intervals
- `--q <path-to-your-query-intervals-file>` - REQUIRED, tells the program where to find the file with query annotation intervals
- `--chs <path-to-your-chromosome-sizes-file>` - REQUIRED, tells the program where to find the file with chromosome sizes
- `--log <name-of-log-file>` - logs will be written into this file, which will be located in the `data/logs/` directory, *if not set the logs will be written to stdout*
- `--o <name-of-results-file>` - results of the program will be written into this file, which will be located in the `data/output/` directory, *if not set the output will be written to stdout*
- `--windows.source <file|basic|dense>` - specifies the windows source of windows set
- `--windows.path <path-to-your-windows-file>` - required with the `--windows.source file` flag, tells the program the location of the window set file
- `--windows.size <windows-size>` - required with the `--windows.source <basic|dense>` flags, tells the program the size of windows to generate
- `--windows.step <windows-step>` - required with the `--windows.source dense` flag, tells the program the shift when generating overlapping set of windows
- `--algorithm <naive|slow_bad|slow|fast_bad|fast>` - defaults to naive, is used to choose algorithm when evaluating windows
- `--significance <enrichment|depletion|combined>` - defaults to enrichment, is used to choose whether to measure enrichment or depletion, combined measures enrichment if observed overlap is larger than mean and depletion otherwise
- `--test` - if this flag is specified, all other flags (except `--help`) are ignored and all the tests in the `src/Tests` are ran and then the program quits
- `--help` - if this flag is specified, all other flags are ignored and a help text will be shown

## References

> [1] Askar Gafurov, Bronislava Brejová, Paul Medvedev.
> Markov chains improve the significance computation of overlapping genome annotations,
> Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i203–i211, https://doi.org/10.1093/bioinformatics/btac255

> [2] Askar Gafurov, Tomáš Vinař, Paul Medvedev, Bronislava Brejová. Efficient Analysis of Annotation Colocalization Accounting for Genomic Contexts. In: Ma, J. (eds) Research in Computational Molecular Biology. RECOMB 2024. Lecture Notes in Computer Science, vol 14758. Springer, Cham. https://doi.org/10.1007/978-1-0716-3989-4_3

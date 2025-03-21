# eMCDP

A tool built on top of mc-overlaps [1] (mcdp2 [2] predecessor), allowing for colocalization analysis on genomic windows. Written in C++.

## How to use this?

Make sure to have [make](https://www.gnu.org/software/make/manual/make.html) installed, then first build the binary with:

```bash
make clean && make build
```

Now there is an `e-mcdp` executable placed in the `bin/` directory in the project root.

The program provides a set of flags to operate it:

- `--r <path-to-your-ref-intervals-file>` - REQUIRED, tells the program where to find the file with reference annotation intervals
- `--q <path-to-your-query-intervals-file>` - REQUIRED, tells the program where to find the file with query annotation intervals
- `--chs <path-to-your-chromosome-sizes-file>` - REQUIRED, tells the program where to find the file with chromosome sizes
- `--log <name-of-log-file>` - logs will be written into this file, which will be located in the `data/logs/` directory
- `--o <name-of-results-file>` - results of the program will be written into this file, which will be located in the `data/output/` directory
- `--windows.source <file|basic|dense>` - specifies the windows source of windows set
- `--windows.path <path-to-your-windows-file>` - required with the `--windows.source file` flag, tells the program the location of the window set file
- `--windows.size` - required with the `--windows.source <basic|dense>` flags, tells the program the size of windows to generate
- `--windows.step` - required with the `--windows.source dense` flag, tells the program the shift when generating overlapping set of windows
- `--algorithm <naive|fast>` - defaults to naive, is used to choose algorithm when evaluating windows

From here you can either run the executable with flags manually or by using any of the `Makefile` directives:

- Samples:
  - `make run_sample` - to run the basic algorithm on basic data set
  - `make run_sample_basic_windows WINDOWS_SIZE=<your-windows-size>` - to run the algorithm with a set of consecutive non-overlapping windows on a basic data set
  - `make run_sample_dense_windows WINDOWS_SIZE=<your-windows-size> WINDOWS_STEP=<your-windows-step>` - to run the algorithm with a set of overlapping windows on a basic data set
- Custom data:
  - arguments `REF_PATH=<ref-path>`, `QUERY_PATH=<query-path>`, `CHR_SIZES_PATH=<chr-sizes-path>`, `OUTPUT_PATH=<output-path>`, `LOG_PATH=<log-path>`, `WINDOWS_PATH=<windows-path>`, `ALGORITHM=<algorithm>` are supposed to replace `--r <ref-path>`, `--q <query-path>`, `--chs <chr-sizes-path>`, `--o <output-path>`, `--log <log-path>`, `--windows.path <windows-path>`, `--algorithm <algorithm>` flags respecitvely
  - `make run_simple_pvalue` - run the basic algorithm on the whole genome
  - `make run_basic_windows` - run the basic algorithm on a basic set of windows
  - `make run_dense_windows` - run the basic algorithm on a dense set of windows
  - `make run_windows_from_file` - run he basic algorithm on a window set specified from file

If you want to redirect the output and logs to stdout, just append `_console` to any of the directives names.

Check the Makefile for more details.

#### Running a large (|R|=2500, |Q|=5000) dataset, p-value for whole genome

```bash
make run_simple_pvalue REF_PATH=data/02-synth-data/g24_8.ref.tsv QUERY_PATH=data/02-synth-data/g24_8.query.tsv CHR_SIZES_PATH=data/02-synth-data/g24_sizes.tsv OUTPUT_PATH=data/output/02-synth-data-g24_8.txt
```

### TODO

- make `CMakeLists.txt` file so everyone can use it easily
- add --h/--help flag and print help info for all flags
- implement a flow to run fast algorithm for overlapping sets of windows:
  - specify --algorithm naive/fast
  - split overlapping windows into sections
  - calculate all sections (4 sets of probs, starting/ending in 0/1)
  - write a function to merge 4 sets of probs for 2 sections
  - for each window merge necessary sections

> [1] Askar Gafurov, Broňa Brejová, Paul Medvedev,
> Markov chains improve the significance computation of overlapping genome annotations,
> Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i203–i211, https://doi.org/10.1093/bioinformatics/btac255

> [2] Gafurov, A., Vinař, T., Medvedev, P., Brejová, B. (2024). Efficient Analysis of Annotation Colocalization Accounting for Genomic Contexts. In: Ma, J. (eds) Research in Computational Molecular Biology. RECOMB 2024. Lecture Notes in Computer Science, vol 14758. Springer, Cham. https://doi.org/10.1007/978-1-0716-3989-4_3

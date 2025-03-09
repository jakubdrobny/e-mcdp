CXX := g++
CXXFLAGS := -Wall -O2 -std=c++20 -fopenmp

SOURCES := $(wildcard src/*.cpp) $(wildcard src/**/*.cpp)
BIN := bin/e-mcdp

# enable sequential prerequisites execution
.NOTPARALLEL:

all: $(BIN)

$(BIN): $(SOURCES)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $^ -o $@ -lgtest -lpthread

install: $(BIN)
	@cp $(BIN) /usr/local/bin
	@echo "Installed e-mcdp"

clean:
	@rm -rf bin
	@echo "Cleaned bin directory."

run_simple_pvalue: clean $(BIN) 
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --o $(OUTPUT_PATH) --log $(LOG_PATH)

run_simple_pvalue_console: clean $(BIN) 
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH)

run_sample_basic_windows_1000000_console: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --windows.source basic --windows.size 1000000

run_sample_console: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv

run_sample: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --o data/output/01-sample-data-hirt.txt --log data/logs/01-sample-data-hirt.txt

.PHONY: all clean install test run_simple_pvalue_console run_simple_pvalue run_sample_console run_sample

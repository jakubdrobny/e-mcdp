CXX := g++
CXXFLAGS := -Wall -O2 -std=c++20 -fopenmp

SOURCES := $(wildcard src/*.cpp) $(wildcard src/**/*.cpp)
BIN := bin/e-mcdp

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

run_simple_pvalue: $(BIN) 
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --o $(OUTPUT_PATH)

run_simple_pvalue_console: $(BIN) 
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH)

run_sample: $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv

.PHONY: all clean install test

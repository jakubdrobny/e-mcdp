CXX := g++
CXXFLAGS := -Wall -O3 -std=c++20 -fopenmp

SOURCES := $(wildcard src/*.cpp) $(wildcard src/**/*.cpp)
BIN := bin/e-mcdp

# enable sequential prerequisites execution
.NOTPARALLEL:

build: $(BIN)

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

run_basic_windows: clean $(BIN)
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --o $(OUTPUT_PATH) --log $(LOG_PATH) --windows.source basic --windows.size $(WINDOWS_SIZE) --algorithm $(ALGORITHM) 

run_basic_windows_console: clean $(BIN)
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --windows.source basic --windows.size $(WINDOWS_SIZE)

run_dense_windows: clean $(BIN)
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --o $(OUTPUT_PATH) --log $(LOG_PATH) --windows.source dense --windows.size $(WINDOWS_SIZE) --windows.step $(WINDOWS_STEP)

run_dense_windows_console: clean $(BIN)
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --windows.source dense --windows.size $(WINDOWS_SIZE) --windows.step $(WINDOWS_STEP)

run_windows_from_file: clean $(BIN)
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --o $(OUTPUT_PATH) --log $(LOG_PATH) --windows.source file --windows.path $(WINDOWS_PATH)

run_windows_from_file_console: clean $(BIN)
	@./bin/e-mcdp --r $(REF_PATH) --q $(QUERY_PATH) --chs $(CHR_SIZES_PATH) --windows.source file --windows.path $(WINDOWS_PATH)

run_sample: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --o data/output/01-sample-data-hirt.txt --log data/logs/01-sample-data-hirt.txt

run_sample_console: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv

run_sample_basic_windows: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --windows.source basic --windows.size $(WINDOWS_SIZE) --log data/logs/01-sample-data-hirt-windows.txt --o data/output/01-sample-data-hirt-windows.tsv

run_sample_basic_windows_console: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --windows.source basic --windows.size $(WINDOWS_SIZE)

run_sample_dense_windows: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --windows.source dense --windows.size $(WINDOWS_SIZE) --windows.step $(WINDOWS_STEP) --log data/logs/01-sample-data-hirt-windows.txt --o data/output/01-sample-data-hirt-windows.tsv

run_sample_dense_windows_console: clean $(BIN)
	@./bin/e-mcdp --r data/01-sample-data/tcga-ref-intervals.tsv --q data/01-sample-data/hirt-query-intervals.tsv --chs data/01-sample-data/chr-sizes.tsv --windows.source dense --windows.size $(WINDOWS_SIZE) --windows.step $(WINDOWS_STEP)

.PHONY: all clean install test run_sample run_sample_console run_simple_pvalue run_simple_pvalue_console run_basic_windows run_basic_windows_console run_dense_windows run_dense_windows_console run_windows_from_file run_windows_from_file_console run_sample_basic_windows run_sample_basic_windows_console run_sample_dense_windows run_sample_dense_windows_console

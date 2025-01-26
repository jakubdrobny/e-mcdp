CXX := g++
CXXFLAGS := -Wall -O2

SOURCES := $(wildcard src/*.cpp) $(wildcard src/**/*.cpp)
BIN := bin/e-mcdp

all: $(BIN)

$(BIN): $(SOURCES)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $^ -o $@

install: $(BIN)
	@cp $(BIN) /usr/local/bin
	@echo "Installed e-mcdp"

clean:
	@rm -rf bin
	@echo "Cleaned bin directory."

.PHONY: all clean install

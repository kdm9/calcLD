
.PHONY: all
all:  bin/calcld bin/binld

bin/calcld: main.cc vcf_reader.cc ld_calculator.cc
	@mkdir -p bin
	g++ -O3 -g -march=native -std=c++17 -o $@ $^ -lhts

bin/binld: bin.cc
	@mkdir -p bin
	g++ -O3 -g -march=native -std=c++17 -o $@ $^

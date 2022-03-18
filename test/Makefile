CWD = $(shell pwd)
INCDIR = $(CWD)/..
SDSLDIR = $(CWD)/../external/sdsl-lite/include/
SIMDEDIR = $(CWD)/../simde/
CXX ?= g++
CPPFLAGS = -std=c++14 -m64 -g -O3 -I$(INCDIR) -I$(SDSLDIR) -I$(SIMDEDIR)


all: test benchmark

.PHONY: clean

test: test.cpp $(INCDIR)/wfa_lm.hpp $(INCDIR)/wfa_lm_st.hpp
	$(CXX) $(CPPFLAGS) test.cpp -o test

benchmark: benchmark.cpp $(INCDIR)/wfa_lm.hpp $(INCDIR)/wfa_lm_st.hpp
	$(CXX) $(CPPFLAGS) benchmark.cpp -o benchmark

clean:
	rm test
	rm benchmark

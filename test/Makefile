CWD = $(shell pwd)
INCDIR = $(CWD)/..
SDSLDIR = $(CWD)/../external/sdsl-lite/include/
CXX ?= g++
CPPFLAGS = -std=c++14 -m64 -g -O3 -I$(INCDIR) -I$(SDSLDIR)


all: test benchmark

.PHONY: clean $(INCDIR)/wfalm.hpp $(INCDIR)/wfalm_st.hpp

test: test.cpp $(INCDIR)/wfalm.hpp $(INCDIR)/wfalm_st.hpp
	$(CXX) $(CPPFLAGS) test.cpp -o test

benchmark: benchmark.cpp
	$(CXX) $(CPPFLAGS) benchmark.cpp -o benchmark

clean:
	rm test
	rm benchmark

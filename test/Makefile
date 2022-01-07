CWD = $(shell pwd)
INCDIR = $(CWD)/..
SDSLDIR = $(CWD)/../external/sdsl-lite/include/
CXX ?= g++
CPPFLAGS = -std=c++14 -m64 -g -O3 -I$(INCDIR) -I$(SDSLDIR)


all: test benchmark

.PHONY: clean $(INCDIR)/wfalm.hpp

test: test.cpp 
	$(CXX) $(CPPFLAGS) test.cpp -o test

benchmark: benchmark.cpp
	$(CXX) $(CPPFLAGS) benchmark.cpp -o benchmark

clean:
	rm test
	rm benchmark
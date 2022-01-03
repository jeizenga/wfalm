CWD = $(shell pwd)
INCDIR = $(CWD)/..
CXX ?= g++
CPPFLAGS = -std=c++11 -m64 -g -O3 -I$(INCDIR)


all: 
	make test

.PHONY: clean $(INCDIR)/wfalm.hpp

test: test.cpp 
	$(CXX) $(CPPFLAGS) test.cpp -o test

clean:
	rm test
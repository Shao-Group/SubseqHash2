CPP=g++
CFLAGS+= -m64 -g -Wall -std=c++11

.PHONY: all
all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))

.PHONY: product
product: CFLAGS = -O3 -std=c++14
product: $(patsubst %.cpp, %.out, $(wildcard *.cpp))

%.out: %.cpp ../seedfactory.h $(wildcard ../*seeding.o) ../strobemer/index.o ../util.o makefile
	g++ $(CFLAGS) -o $@ $(filter-out makefile ../seedfactory.h, $^) -pthread
#g++ -Wall -g -std=c++14 -o $@ $(filter-out makefile,$^) -pthread

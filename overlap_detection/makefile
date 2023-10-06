CPP=g++
CFLAGS+= -m64 -g -Wall -std=c++11

.PHONY: dependency-%
dependency-%:
	$(MAKE) -C .. $(subst dependency-,,$@)

.PHONY: all
all: dependency-all $(patsubst %.cpp, %.out, $(wildcard *.cpp))

.PHONY: product
product: CFLAGS = -O3 -std=c++14
product: dependency-product $(patsubst %.cpp, %.out, $(wildcard *.cpp))

DEPENDENCIES = ../seedfactory.h $(wildcard ../*seeding.o) ../strobemer/index.o ../util.o

%.out: %.cpp $(DEPENDENCIES) makefile 
	g++ $(CFLAGS) -o $@ $(filter-out makefile ../seedfactory.h, $^) -pthread
#g++ -Wall -g -std=c++14 -o $@ $(filter-out makefile,$^) -pthread

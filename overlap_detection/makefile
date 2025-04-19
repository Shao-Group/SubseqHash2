CPP=g++
CFLAGS+= -m64 -g -Wall -std=c++11

.PHONY: dependency-%
dependency-%:
	$(MAKE) -C ../src $(subst dependency-,,$@)

.PHONY: all
all: dependency-all $(patsubst %.cpp, %.out, $(wildcard *.cpp))

.PHONY: product
product: CFLAGS = -O3 -std=c++14
product: dependency-product $(patsubst %.cpp, %.out, $(wildcard *.cpp))

DEPENDENCIES = ../src/seedfactory.hpp $(patsubst %.cpp, %.o, $(wildcard ../src/*seeding.cpp)) ../src/strobemer/index.o ../src/util.o
.SECONDARY: $(DEPENDENCIES)

%.out: %.cpp $(DEPENDENCIES) makefile 
	g++ $(CFLAGS) -o $@ $(filter-out makefile ../seedfactory.hpp, $^) -pthread
#g++ -Wall -g -std=c++14 -o $@ $(filter-out makefile,$^) -pthread

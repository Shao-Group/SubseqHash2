CC=gcc
CPP=g++
CFLAGS+= -m64 -g -Wall -std=c++11 -march=native
LIBS=
ALLDEP:= $(patsubst %.hpp,%.o,$(filter-out seedfactory.hpp, $(wildcard *.hpp))) seedfactory.hpp ./strobemer/index.o

.PHONY: all
all: strobemer $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.hpp,%.cpp,$(wildcard *.hpp)), $(wildcard *.cpp)))

.PHONY: product-c++11
product-c++11: CFLAGS = -O3 -std=c++11 -march=native
product-c++11: strobemer $(ALLDEP) $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.o,%.cpp,$(ALLDEP)) $(ALLILP), $(wildcard *.cpp)))

.PHONY: product
product: CFLAGS = -O3 -std=c++14 -march=native
product: strobemer $(ALLDEP) $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.o,%.cpp,$(ALLDEP)) $(ALLILP), $(wildcard *.cpp)))

./strobemer/index.o: ./strobemer/index.cpp
	$(CPP) $(CFLAGS) -MMD -c ./strobemer/index.cpp -o ./strobemer/index.o

%.out: %.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $(filter-out %.hpp, $^) $(LIBS)
%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.cpp %.hpp makefile
	$(CPP) $(CFLAGS) -MMD -c $< -o $@
%.o: %.c %.hpp makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d
	rm ./strobemer/*.o ./strobemer/*.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     

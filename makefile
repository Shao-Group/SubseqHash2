CC=gcc
CPP=g++
CFLAGS+= -m64 -g -Wall -std=c++14
LDFLAGS= -L$$GUROBI_HOME/lib -lgurobi91
LIBS=./strobemer/index.o
INC= $$GUROBI_HOME/include/
ALLDEP:= $(patsubst %.h,%.o,$(wildcard *.h))
ALLILP:= $(wildcard *_ILP.c)

.PHONY: all
all: strobemer $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.h,%.cpp,$(wildcard *.h)), $(wildcard *.cpp)))

.PHONY: product-c++11
product-c++11: CFLAGS = -O3 -std=c++11
product-c++11: strobemer $(ALLDEP) $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.o,%.cpp,$(ALLDEP)) $(ALLILP), $(wildcard *.cpp)))

.PHONY: product
product: CFLAGS = -O3 -std=c++14
product: strobemer $(ALLDEP) $(patsubst %.cpp,%.out,$(filter-out $(patsubst %.o,%.cpp,$(ALLDEP)) $(ALLILP), $(wildcard *.cpp)))

test%.out: test%.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $^ $(LIBS) -pthread
test2%.out: test%.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $^ $(LIBS) -pthread

strobemer: ./strobemer/index.cpp
	$(CPP) $(CFLAGS) -MMD -c ./strobemer/index.cpp -o ./strobemer/index.o

%.out: %.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $^ $(LIBS)
%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.cpp %.h makefile
	$(CPP) $(CFLAGS) -MMD -c $< -o $@
%.o: %.c %.h makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

#gurobi make
%_ILP: %_ILP.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ -I$(INC) $(LDFLAGS) -lm

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     

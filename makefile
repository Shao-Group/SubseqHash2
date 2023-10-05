CC=gcc
CPP=g++
CFLAGS+= -m64 -g -Wall -std=c++11
LIBS=
ALLDEP:= $(patsubst %.h,%.o,$(filter-out seedfactory.h, $(wildcard *.h))) seedfactory.h ./strobemer/index.o

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

./strobemer/index.o: ./strobemer/index.cpp
	$(CPP) $(CFLAGS) -MMD -c ./strobemer/index.cpp -o ./strobemer/index.o

%.out: %.cpp $(ALLDEP)
	$(CPP) $(CFLAGS) -o $@ $(filter-out %.h, $^) $(LIBS)
%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.cpp %.h makefile
	$(CPP) $(CFLAGS) -MMD -c $< -o $@
%.o: %.c %.h makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     

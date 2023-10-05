%.out: %.cpp ../seedfactory.h ../seeding.o ../subseqhash2seeding.o ../util.o makefile
	g++ -O3 -std=c++14 -o $@ $(filter-out makefile ../seedfactory.h, $^) -pthread
#g++ -Wall -g -std=c++14 -o $@ $(filter-out makefile,$^) -pthread

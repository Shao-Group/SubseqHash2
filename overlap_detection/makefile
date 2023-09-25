%.out: %.cpp ../seedfactory.o ../seeding.o ../subseqhash2seeding.o ../util.o makefile
	g++ -O3 -std=c++14 -o $@ $(filter-out makefile,$^) -pthread
#g++ -Wall -g -std=c++14 -o $@ $(filter-out makefile,$^) -pthread

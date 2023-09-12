%.out: %.cpp ../seedfactory.o ../seeding.o ../subseqhash2seeding.o ../util.o
	g++ -O3 -std=c++14 -o $@ $^ -pthread

CXXFLAGS = g++ --std=c++11 -O2
CXXLDFLAGS = -lpthread

all: pi.o
	$(CXXFLAGS) -o pi.out pi.o $(CXXLDFLAGS)  


pi.o: pi.cpp
	$(CXXFLAGS) -c pi.cpp

run:
	./pi

check:
	./pi > result.txt

clean:
	rm -rf *.out *.o result.txt
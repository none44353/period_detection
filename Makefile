CFLAGS=-g --std=c++11 -O2

all:test

test:testDataset.o
	g++ ${CFLAGS} -o $@ $+

.cpp.o:
	g++ -c ${CFLAGS} -o $@ $+

clean:
	rm -f *.o test

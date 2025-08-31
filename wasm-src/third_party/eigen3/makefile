all: qrep

qrep: qrep.o
	g++ -o qrep qrep.o

mst-backbone.o: qrep.cpp
	g++ -c qrep.cpp -o qrep.o -O3 -I. -std=c++14

clean:
	rm qrep qrep.o



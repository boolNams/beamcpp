prog: main.o func.o
	g++ -o prog main.o func.o

main.o: main.cpp header.h
	g++ -c main.cpp

func.o: main.cpp header.h
	g++ -c func.cpp
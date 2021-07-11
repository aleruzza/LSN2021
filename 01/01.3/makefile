CC = g++
CFLAGS = -Wall -O3 --std=c++11

main: build/main.o build/random.o
	$(CC) build/random.o build/main.o -o main
build/main.o: src/main.cpp
	$(CC) -c src/main.cpp -o build/main.o $(CFLAGS)
build/random.o: src/random.cpp src/random.h
	$(CC) -c src/random.cpp -o build/random.o $(CFLAGS)
clean:
	rm build/*.o main seed.out

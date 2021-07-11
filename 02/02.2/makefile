CC = g++
CFLAGS = -Wall -O3 --std=c++11

main: build/main.o build/random.o build/utils.o src/Walker.h build/Walker.o build/WalkerLattice.o build/WalkerContinuum.o
	$(CC) build/random.o build/main.o build/utils.o build/Walker.o build/WalkerLattice.o build/WalkerContinuum.o -o main
build/main.o: src/main.cpp
	$(CC) -c src/main.cpp -o build/main.o $(CFLAGS)
build/random.o: src/random.cpp src/random.h
	$(CC) -c src/random.cpp -o build/random.o $(CFLAGS)
build/utils.o: src/utils.cpp src/utils.h
	$(CC) -c src/utils.cpp -o build/utils.o $(CFLAGS)
build/Walker.o: src/Walker.cpp src/Walker.h
	$(CC) -c src/Walker.cpp -o build/Walker.o $(CFLAGS)
build/WalkerLattice.o: src/WalkerLattice.cpp src/WalkerLattice.h
	$(CC) -c src/WalkerLattice.cpp -o build/WalkerLattice.o $(CFLAGS)
build/WalkerContinuum.o: src/WalkerContinuum.cpp src/WalkerContinuum.h
	$(CC) -c src/WalkerContinuum.cpp -o build/WalkerContinuum.o $(CFLAGS)
clean:
	rm build/*.o main seed.out

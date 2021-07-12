CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main: build/main.o build/random.o build/utils.o build/IncubatorGA.o build/TSPIndividual.o build/TSProblem.o
	$(CC) build/random.o build/main.o build/utils.o build/IncubatorGA.o build/TSPIndividual.o build/TSProblem.o -o main
build/main.o: src/main.cpp
	$(CC) -c src/main.cpp -o build/main.o $(CFLAGS)
build/TSProblem.o: src/TSProblem.cpp
	$(CC) -c src/TSProblem.cpp -o build/TSProblem.o $(CFLAGS)
build/IncubatorGA.o: src/IncubatorGA.cpp
	$(CC) -c src/IncubatorGA.cpp -o build/IncubatorGA.o
build/TSPIndividual.o: src/TSPIndividual.cpp
	$(CC) -c src/TSPIndividual.cpp -o build/TSPIndividual.o $(CFLAGS)
build/random.o: src/random.cpp src/random.h
	$(CC) -c src/random.cpp -o build/random.o $(CFLAGS)
build/utils.o: src/utils.cpp src/utils.h
	$(CC) -c src/utils.cpp -o build/utils.o $(CFLAGS)
clean:
	rm build/*.o main seed.out

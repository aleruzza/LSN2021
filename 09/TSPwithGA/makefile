CC = g++
CFLAGS = -Wall -O3 --std=c++11

main: build/main.o build/random.o build/utils.o build/IncubatorGA.o build/TSPIndividual.o build/TSProblem.o src/TSPsquare.h src/TSPcircumference.h
	$(CC) build/random.o build/main.o build/utils.o build/IncubatorGA.o build/TSPIndividual.o build/TSProblem.o -o main $(CFLAGS)
build/main.o: src/main.cpp src/random.h src/TSPsquare.h src/TSPcircumference.h
	$(CC) -c src/main.cpp -o build/main.o $(CFLAGS)
build/TSProblem.o: src/TSProblem.cpp src/TSProblem.h
	$(CC) -c src/TSProblem.cpp -o build/TSProblem.o $(CFLAGS)
build/IncubatorGA.o: src/IncubatorGA.cpp src/IncubatorGA.h
	$(CC) -c src/IncubatorGA.cpp -o build/IncubatorGA.o $(CFLAGS)
build/TSPIndividual.o: src/TSPIndividual.cpp src/TSPIndividual.h
	$(CC) -c src/TSPIndividual.cpp -o build/TSPIndividual.o $(CFLAGS)
build/random.o: src/random.cpp src/random.h
	$(CC) -c src/random.cpp -o build/random.o $(CFLAGS)
build/utils.o: src/utils.cpp src/utils.h
	$(CC) -c src/utils.cpp -o build/utils.o $(CFLAGS)
clean:
	rm build/*.o main seed.out

#ifndef UTILS
#define UTILS

#include "random.h"
#include <cmath>
struct point 
{
    double x;
    double y;
    double z;
};

//funzione che inizializza il generatore di numeri casuali utilizzando file di input e output
// Primes seed.in seed.out
void InitRandomGenerator(Random& rnd);

double getR(point p);

point getPoint(double x, double y, double z);

#endif

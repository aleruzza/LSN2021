#ifndef UTILS
#define UTILS

#include "random.h"
#include <cmath>

//punto in 1D
struct point 
{
    double x;
   // double y;
   // double z;
};

struct measure
{
    double val;
    double err;
    double acceptance_rate;
};

//funzione che inizializza il generatore di numeri casuali utilizzando file di input e output
// Primes seed.in seed.out
void InitRandomGenerator(Random& rnd);

double getR(point p);

point getPoint(double x, double y, double z);

#endif

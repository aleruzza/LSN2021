#include "utils.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void InitRandomGenerator(Random& rnd)
{
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    rnd.SaveSeed();
}

double getR(point p)
{
    return sqrt(p.x*p.x+ p.y*p.y + p.z*p.z);
}

point getPoint(double x, double y, double z)
{
    point p;
    p.x = x; p.y = y, p.z=z;
    return p;
}

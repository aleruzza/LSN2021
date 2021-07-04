#ifndef PROBDIST
#define PROBDIST

#include "utils.h"
#include <cmath>

class ProbabilityDist
{
    public:
        virtual double eval(point p) = 0;
};

class HydrogenAtom1S : public ProbabilityDist
{
    public:
        double eval(point p)
        {
            double r = getR(p);
            return exp(-2*r)/M_PI;
        }
};

class HydrogenAtom2P : public ProbabilityDist
{
    public:
        double eval(point p)
        {
            double r = getR(p);
            return pow(exp(-r/2)*sqrt(2/M_PI)*r*cos_theta/8, 2);
        }
};

#endif

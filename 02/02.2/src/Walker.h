#ifndef WALKERH
#define WALKERH

#include <cmath>
#include "random.h"

class Walker
{
    public:
        Walker(Random* rnd);
        
        Walker(int dim, Random* rnd);
        
        int getSteps();
        
        virtual void walk() = 0;
        
        double getR2();
        
        double getR();
        
        void restart();
        
    protected:
        double* m_pos;
        
        int m_dim;
        int steps;
        Random *m_rnd;
};

#endif

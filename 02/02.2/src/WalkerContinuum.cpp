#include "WalkerContinuum.h"
#include <cmath>
#include <iostream>

#define X 0
#define Y 1
#define Z 2

WalkerContinuum::WalkerContinuum(Random* rnd, double step_length) : Walker(DIM, rnd), m_step_l(step_length){}

WalkerContinuum::WalkerContinuum(Random* rnd) : Walker(DIM, rnd), m_step_l(DEF_STEP_L){}

void WalkerContinuum::walk()
{
    double theta = randomTheta();
    double phi = randomPhi();
    
    m_pos[X] += m_step_l*sin(theta)*cos(phi);
    m_pos[Y] += m_step_l*sin(theta)*sin(phi);
    m_pos[Z] += m_step_l*cos(theta);
    steps++;
}

double WalkerContinuum::randomTheta()
{
    double y = acos(1-2*m_rnd->Rannyu(0,1));
    //double y = m_rnd->Rannyu(0, M_PI);
    return y;
}

double WalkerContinuum::randomPhi()
{
    double ph = m_rnd->Rannyu(0, 2*M_PI);
    return ph;
}

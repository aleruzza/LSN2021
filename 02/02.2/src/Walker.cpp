#include "Walker.h"
#include "random.h"

Walker::Walker(Random* rnd) : Walker(3, rnd)
{}

Walker::Walker(int dim, Random* rnd)
{
    m_rnd = rnd;
    m_dim = dim;
    m_pos = new double[m_dim];
    restart();
}

void Walker::restart()
{
    for(int i=0;i<m_dim;i++)
        m_pos[i] = 0;
    steps=0;
}

double Walker::getR2()
{
    double r2=0;
    for(int i=0;i<m_dim;i++)
        r2+=pow(m_pos[i], 2);
    return r2;
}

double Walker::getR()
{
    return sqrt(getR2());
}

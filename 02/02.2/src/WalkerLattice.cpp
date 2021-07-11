#include "WalkerLattice.h"

//costruttori
WalkerLattice::WalkerLattice(Random* rnd) : Walker(3, rnd), m_latt_a(DEF_LATT_A)
{}

WalkerLattice::WalkerLattice(Random* rnd, double latt_a) : Walker(3, rnd),  m_latt_a{latt_a}
{}

void WalkerLattice::walk()
{
    int dir = randomDir();
    int sign = randomSign();
    m_pos[dir] += sign*m_latt_a;
    steps++;
}

int WalkerLattice::randomDir()
{
    int dir = m_dim*(m_rnd->Rannyu());
    return dir;
}

int WalkerLattice::randomSign()
{
    if((m_rnd->Rannyu())<0.5)
        return -1;
    else
        return 1;
}


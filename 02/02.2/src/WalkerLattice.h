#ifndef WALKERLATT
#define WALKERLATT
#include "Walker.h"
#define DEF_LATT_A 1
class WalkerLattice: public Walker 
{
    public:
        
        WalkerLattice(Random* rnd);
        
        //constructor used to set the lattice step
        WalkerLattice(Random* rnd, double latt_a);
        
        //Walker virtual methods that need to be implemented
        virtual void walk();
        
    private:
        double m_latt_a;
        
        //functions used internally
        int randomDir();
        int randomSign();
};
#endif

#ifndef WALKERCONT
#define WALKERCONT
#include "Walker.h"

#define DEF_STEP_L 1
#define DIM 3
class WalkerContinuum: public Walker 
{
    public:
        
        WalkerContinuum(Random* rnd);
        
        //constructor used to set the step length
        WalkerContinuum(Random* rnd, double step_length);
        
        //metodi virtuali di Walker da implementare
        virtual void walk();
        
    private:
        
        double m_step_l;
        
        //funzioni utili internamente
        double randomTheta();
        double randomPhi();
        void RandomPos(double& phi, double& theta);
};
#endif

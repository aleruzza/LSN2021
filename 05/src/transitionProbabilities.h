#ifndef TRANSPROB
#define TRANSPROB

#include "utils.h"

class TransitionProbability
{
    public:
        virtual point  nextPoint(point pi) = 0;
};

class UniformTP : public TransitionProbability
{
    public:
        UniformTP(double r, Random& rnd): m_r(r), m_rnd(rnd){}
        
        point nextPoint(point pi)
        {
            point pf;
    
            pf.x = pi.x + m_r*(m_rnd.Rannyu()-0.5);
            pf.y = pi.y + m_r*(m_rnd.Rannyu()-0.5);
            pf.z = pi.z + m_r*(m_rnd.Rannyu()-0.5);
            
            return pf;
        }
        
    private:
        double m_r;
        Random m_rnd;
};

class GaussTP : public TransitionProbability
{
    public:
        GaussTP(double sigma, Random& rnd): m_sigma(sigma), m_rnd(rnd){}
        
        point nextPoint(point pi)
        {
            point pf;
            
            pf.x = m_rnd.Gauss(pi.x, m_sigma);
            pf.y = m_rnd.Gauss(pi.y, m_sigma);
            pf.z = m_rnd.Gauss(pi.z, m_sigma);
            
            return pf;
        }
        
    private:
        double m_sigma;
        Random m_rnd;
};
#endif //TRANSPROB

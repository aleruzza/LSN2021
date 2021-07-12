#ifndef TRANSPROB
#define TRANSPROB

#include "utils.h"
#include <iostream>
#include <vector>
using namespace std;

class TransitionProbability
{
    public:
        virtual point nextPoint(point pi) = 0;
};

class UniformTP : public TransitionProbability
{
    public:
        UniformTP(double r, Random& rnd): m_r(r), m_rnd(rnd){
            for(int i=0;i<100;i++)
            {
                random_numbers.push_back(m_rnd.Rannyu());
            }

            random_index=m_rnd.Rannyu();
        }
        
        point nextPoint(point pi)
        {
            point pf;
            
            int int_random_index = random_index*100;
           
            //pf.x = m_rnd.Rannyu(pi.x-m_r, pi.x+m_r);
            
            pf.x = pi.x + 2.*m_r*(random_numbers[int_random_index]-0.5);
            random_index = random_numbers[int_random_index];
            random_numbers[int_random_index] = m_rnd.Rannyu();
            // for(int i=0;i<m_rnd.Rannyu()*5;i++)
            //m_rnd.Rannyu();
            // pf.y = pi.y + m_r*(m_rnd.Rannyu()-0.5);
            // pf.z = pi.z + m_r*(m_rnd.Rannyu()-0.5);
            
            return pf;
        }
        
    private:
        double m_r;
        Random m_rnd;
        vector<double> random_numbers;
        double random_index;
};

class GaussTP : public TransitionProbability
{
    public:
        GaussTP(double sigma, Random& rnd): m_sigma(sigma), m_rnd(rnd){}
        
        point nextPoint(point pi)
        {
            point pf;
            
            pf.x = m_rnd.Gauss(pi.x, m_sigma);
            //pf.y = m_rnd.Gauss(pi.y, m_sigma);
            //pf.z = m_rnd.Gauss(pi.z, m_sigma);
            
            return pf;
        }
        
    private:
        double m_sigma;
        Random m_rnd;
};
#endif //TRANSPROB

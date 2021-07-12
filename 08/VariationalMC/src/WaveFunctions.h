#ifndef WAVEFUNCT
#define WAVEFUNCT

#include "utils.h"
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

//kept the use of the struct point even if now is just a double
//to allow easier adaptation of the code to a 3D case

//virtual class of a generic WaveFunction
class WaveFunction
{
    public:
        ~WaveFunction(){}
        virtual double eval(point p) = 0;
        virtual double evalAbsSquared(point p) = 0;
        virtual double evalFirstDerivative(point p) = 0;
        virtual double evalSecondDerivative(point p) = 0;
        virtual std::string getParametersInString() = 0;
};

//Parametric 1D Wave Function used for the exercise
//NN stands for Not Normalized, not needed for M(RT)^2
class NNDoubleGaussWF : public WaveFunction
{
    public:
        NNDoubleGaussWF(double mu, double sigma): m_sigma(sigma), m_mu(mu){}
        
        double eval(point p)
        {
            double psi = exp(-0.5*pow((p.x-m_mu)/m_sigma, 2)) + exp(-0.5*pow((p.x+m_mu)/m_sigma, 2));
            return psi;
        }
        
        double evalAbsSquared(point p)
        {
            //this wf is real
            //double psi2 = exp(-pow((p.x-m_mu)/m_sigma, 2)) + exp(-pow((p.x+m_mu)/m_sigma, 2)) + 2*exp(-0.5*pow((p.x-m_mu)/m_sigma, 2)-0.5*pow((p.x+m_mu)/m_sigma, 2));
            return pow(eval(p), 2);

            /*
            if(p.x>-1 && p.x<1)
                return 1;
            else
                return 0;*/
           // return exp(-0.5*p.x*p.x);
            //return psi2;
        }

        double evalFirstDerivative(point p)
        {
            double dpsi = -(exp(-pow((p.x-m_mu)/m_sigma, 2)/2)*(p.x-m_mu) 
                                + exp(-pow((p.x+m_mu)/m_sigma, 2)/2)*(p.x+m_mu))/(m_sigma*m_sigma);
            return dpsi;
        }

        double evalSecondDerivative(point p)
        {
            double ddpsi = (exp(-pow((p.x-m_mu)/m_sigma, 2)/2)*(pow((p.x-m_mu)/m_sigma, 2)-1)
                                + exp(-pow((p.x+m_mu)/m_sigma, 2)/2)*(pow((p.x+m_mu)/m_sigma, 2)-1))/(m_sigma*m_sigma);
            return ddpsi;
        }

        std::string getParametersInString()
        {
            std::ostringstream strs;
            strs << std::setw(12) << m_mu << std::setw(12) << m_sigma;
            std::string str = strs.str();
            return str;
        }
        
    private:
        double m_sigma, m_mu;
};


#endif

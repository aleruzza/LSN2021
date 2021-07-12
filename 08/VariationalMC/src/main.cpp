#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <string>
#include <algorithm>
#include "random.h"
#include "utils.h"
#include "TransitionProbabilities.h"
#include "WaveFunctions.h"

using namespace std;
 
#define UNIFORM 1
#define GAUSS 2

bool accept(point p1, point p2, WaveFunction* p, Random& rnd);
double HOverPsi(WaveFunction* psi, point p);
measure computeEnergy(WaveFunction* psi, int nblock, int nstep, Random& rnd, bool saveMeans, bool savePoints, int npoints);
void Input();

//steps
int n_block, n_step, n_step_eq;

//ranges of WF parameters
double mu_min, mu_max, mu_n;
double sigma_min, sigma_max, sigma_n;

//flags
bool save_means, save_points;

//number of sampled points to save
int npoints_tosave;

int main (int argc, char *argv[]){
    
    //Initialization of the random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    //Get the parameters from io/input.dat
    Input();

    /*
        Open the output file for final measures
        refering to simulations done with different values of 
        sigma and mu for the wave function.
        Every line will have the following data
            mu  sigma   energy  energy_error   acceptance_rate
    */
    ofstream fout("io/final_measures.output");
    

    //the (mu, sigma) space is explored with a grid sampling
    double mu_current, mu_step=(mu_max-mu_min)/(double) mu_n;
    double sigma_current, sigma_step=(sigma_max-sigma_min)/(double) sigma_n;
    
    WaveFunction *psi;
    measure meas;
    int totsim = (mu_n)*(sigma_n);
    int current_sim = 0;
    int i_mu, i_sigma;

    for(i_mu = 0;i_mu < mu_n; i_mu++)
    {
        mu_current = mu_min + mu_step*i_mu;
        for(i_sigma = 0; i_sigma < sigma_n; i_sigma++)
        {
            sigma_current = sigma_min + sigma_step*i_sigma;
            current_sim += 1;

            //Wave function for current simulation
            psi = new NNDoubleGaussWF(mu_current, sigma_current);

            //Compute energy with a MC simulation using data blocking. 
            //save_means == true => save the updated mean and error after every block
            //save_points == true => will save some (npoints_tosave) of the sampled points
            meas = computeEnergy(psi, n_block, n_step, rnd, save_means, save_points, npoints_tosave);

            //output of the final measure on file
            fout << psi->getParametersInString() << setw(15) << meas.val << setw(15) << meas.err << setw(15) << meas.acceptance_rate << endl;

            //terminal output to show progress
            //if(!current_sim%10)
            { 
                cout << setw(12) << current_sim << "/" << totsim << setw(12) << "( " << setw(6) << mu_current << ", " << setw(6) << sigma_current << " ): ";
                cout << setw(12) << meas.val << setw(8) << "+- " << setw(10) << meas.err << setw(20) << "(a.r.: " << meas.acceptance_rate << ")" << endl;
            }
            delete psi;
        }
    }
    
    fout.close();
    
    return 0;
}

void Input()
{
    ifstream input("io/input.dat");

    input >> n_block;
    input >> n_step;
    input >> n_step_eq;
    input >> mu_min;
    input >> mu_max;
    input >> mu_n;
    input >> sigma_min;
    input >> sigma_max;
    input >> sigma_n;
    input >> save_means;
    input >> save_points;
    input >> npoints_tosave;

}

measure computeEnergy(WaveFunction* psi, int nblock, int nstep, Random& rnd, bool saveMeans, bool savePoints, int npoints)
{

    //variables used in the saving of sampled points
    int dist_points = 0;
    if(savePoints && npoints!=0)
        dist_points = nstep/npoints;

    //output files, open them only if needed
    ofstream f_blkMeans;
    ofstream f_sampledPoints;
    if(saveMeans)
        f_blkMeans.open("io/blkmeans.output", ios::app);
    if(savePoints)
        f_sampledPoints.open("io/points.output", ios::app);

    //Transition function used in M(RT)^2 algorithm
    //UniformTP, here in 1D, generate x_new 
    //with uniform probability between x_current-r and x_current+r
    //r is chosen by trial and error to obtain appr. 50% acceptance rate 
    double r = 2;
    TransitionProbability *trP = new UniformTP(r, rnd);

    //points used during the simulation, in this 1D code only p.x is used
    point p1, p2;
    
    //starting point
    p2.x = 0;
    
    //equilibration phase
    for(int i=0;i<n_step_eq;i++)
    {
        p1 = p2;
        p2 = trP->nextPoint(p1);
        if(!accept(p1, p2, psi, rnd))
            p2 = p1;
    }
    
    //measurement phase
    int nstep_per_block = nstep/nblock;
    double blockAcc, meansAcc = 0.0, meansSquaredAcc = 0.0;
    double val=0.0, unc=0.0;
    long int accepted = 0, rejected = 0;
    
    for(int iblk=0;iblk<nblock;iblk++)
    {
        blockAcc = 0.0;
        for(int istep=0; istep<nstep_per_block; istep++)
        {
            p1 = p2;
            p2 = trP->nextPoint(p1);
            if(!accept(p1, p2, psi, rnd))
            {
                p2 = p1;
                rejected++;
            }
            else
                accepted++;
            
            //update block accumulator
            blockAcc += HOverPsi(psi, p2);
            
            //save the coordinate of sampled points on file if requested
            if(savePoints)
            {
                if(!( ((iblk*nstep_per_block)+istep)%dist_points ))
                    f_sampledPoints << psi->getParametersInString() << setw(20) << p2.x << endl;
            }
        }
        
        //update block accumulators
        blockAcc /= (double) nstep_per_block;
        meansAcc += blockAcc;
        meansSquaredAcc += blockAcc*blockAcc;

        //compute global mean at current block with its uncertainty
        val = meansAcc/(iblk+1);
        if(iblk>0)
            unc = sqrt((meansSquaredAcc/(iblk+1) - val*val)/iblk);

        //save block mean, last global mean and uncertainty (if requested)
        if(saveMeans)
            f_blkMeans << psi->getParametersInString() << setw(12) << blockAcc << setw(12) << val << setw(12) << unc << endl;
    }
    
    //close files used, only if opened in the first time
    if(saveMeans)
        f_blkMeans.close();
    if(savePoints)
        f_sampledPoints.close();

    //prepare results of simulation to return 
    measure m;
    m.val = val;
    m.err = unc;
    m.acceptance_rate = accepted/(double) (accepted+rejected);

    return m;
}

double HOverPsi(WaveFunction* psi, point p)
{
    //compute H(psi)/psi in p (p.x)
    //this function implements a specific potential
    //code units: m=1, hbar = 1
    double psi_in_p = psi->eval(p);
    double t=0, v=0;
    if(psi_in_p!=0)
    {
        t = -psi->evalSecondDerivative(p)/(2.*psi_in_p);
        v = (pow(p.x, 4) - 5*p.x*p.x/2);
    }
    return t+v;
}

//function used to decide whether to accept the moves proposed 
//in the M(RT)^2 algorithm, p1 -> p2.
//Implemented for symmetric transition probabilities
bool accept(point p1, point p2, WaveFunction* p, Random& rnd)
{
    if(p->evalAbsSquared(p1)==0)
        return true;
    if(p->evalAbsSquared(p2)==0)
        return false;

    double alpha = min(1., p->evalAbsSquared(p2)/p->evalAbsSquared(p1));

    //accept with probability alpha
    return (rnd.Rannyu()<alpha);
}


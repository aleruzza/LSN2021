
/* Compute an integral using Monte Carlo methods.
 * The same integral can be computed using both uniform sampling 
 * or importance sampling. The pdf used for importance 
 * sampling has been chosen looking at the function to integrate.
 * 
 * Usage:
 * -------
 * ./main <mode> <file_output>
 * mode: 1-uniform sampling
 *       2-importance sampling
 * 
 * Output:
 * ---------
 * In each output file the computed results are stored in this format:
 * each row refers to a different block and contains the best estimate 
 * of the integral at that block with its  statistical uncertainty
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "utils.h"

//Monte Carlo steps
#define M 1000000
//number of blocks
#define K 100

using namespace std;
 
double func(double x);
double rand_lin(Random& rnd);
double d(double x);

int main (int argc, char *argv[]){
    
    //open files for the output 
    if(argc!=3)
    {
        cout << "Usage: ./main <mode> <file_output>" << endl;
        cout << "mode: 1- uniform sampling " << endl << "2- importance sampling" << endl;
        exit(1);
    }
    
    ofstream fout(argv[2]);
    
    //initialize random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    double block_i;
    double sum_i=0, sumi_2=0;
    double x;
    double last_i_est;
    double unc;
        
    const int N = M/K;
    if(atoi(argv[1])==1)
    {
        //compute integral with uniform sampling
        //cycle over blocks
        for(int i=0;i<K;i++)
        {
            //estimate integral with M/K=N throws
            block_i=0;
            for(int j=0;j<N;j++)
            {
                //estraggo la x tra 0 e 1 
                x =rnd.Rannyu();
                block_i += func(x)/N;
            }
            
            //compute last estimate of the integral
            sum_i += block_i;
            last_i_est = sum_i/(i+1);
            
            //estimate statistical uncertainty
            sumi_2 += pow(block_i, 2);
            if(i>1)
                unc = sqrt((sumi_2/(i+1) - pow(last_i_est, 2))/i);
            
            fout << last_i_est << " " << unc << endl;
            cout << "blocco " << i << " " << last_i_est << " +- " << unc << endl;
        }
    }
    else
        if(atoi(argv[1])==2)
        {
            for(int i=0;i<K;i++)
            {
                //importance sampling with pdf f(x) = 2(1-x)
                //estimate integral with M/K=N throws
                block_i=0;
                for(int j=0;j<N;j++)
                {
                    //extract x between 0 and 1 with exponential probability distribution
                    x = rand_lin(rnd);
                    block_i += func(x)/(N*lin_pdf(x));
                }
                
                //estimate the integral
                sum_i += block_i;
                last_i_est = sum_i/(i+1);
                
                //estimate statistical uncertainty 
                sumi_2 += pow(block_i, 2);
                if(i>1)
                    unc = sqrt((sumi_2/(i+1) - pow(last_i_est, 2))/i);
                
                fout << last_i_est << " " << unc << endl;
                cout << "blocco " << i << " " << last_i_est << " +- " << unc << endl;
            }
        }
    
    return 0;
}

//function to integrate
double func(double x)
{
    double f = 0.5*M_PI*cos(M_PI*x*0.5);
    return f;
}

//pdf used for importance sampling
double lin_pdf(double x)
{
    double d;
    if(x>=0 && x<1)
        d = 2*(1-x);
    else
        d = 0;
    return d;
}

//generate random number between 0 and 1
//according to the pdf f(x) = 2(1-x) 
double rand_lin(Random& rnd)
{
    double n;
    n = (1-sqrt(1-rnd.Rannyu()));
    return n;
}


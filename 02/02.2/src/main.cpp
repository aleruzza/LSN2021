/* Simulates a random walk and computes the mean distance 
 * from the origin travelled as a function of the number of steps. 
 * This code implements both the walk on a lattice
 * and the walk in the continuum
 * 
 * Usage:   
 * ---------
 *  ./main <mode> <output>
 *          mode: l - random walk on a lattice
 *                c - random walk in the continuum
 * 
 * Output
 * ---------
 * In each output file the computed results are stored in this format:
 * each row refers to a different moment during the walk characterized 
 * by the number of steps made. The first column reports the number of steps,
 * the second reports the squared distance of the walker from the origin
 * averaged over the M simulations of random walk divided in K blocks. 
 * The last columns reports the statistical uncertainty of the mean
 * distance value computed.
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "utils.h"
#include "WalkerLattice.h"
#include "WalkerContinuum.h"
#include "Walker.h"

using namespace std;

#define M 1000000  //simulations of random walk
#define K 100    // number of blocks
#define STEPS 100 // walker's steps
#define DIM 3 //dimensions of space

void resetToZero(double* vett, int dim);

int main (int argc, char *argv[]){
    
    //open output files
    if(argc!=3)
    {
        cout << "Usage: ./main <mode> <output>" << endl;
        cout << "mode: l - random walk on a lattice" << endl;
        cout << "      c - random walk in the continuum" << endl;
        exit(1);
    }
    ofstream fout(argv[2]);
    
    //initializing the random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    //instantiating the walker
    Walker* walker;
    if(argv[1][0]=='l')
        walker = new WalkerLattice(&rnd);
    else
        if(argv[1][0]=='c')
            walker = new WalkerContinuum(&rnd);
        else
            exit(2);
    
    //number of random walks for each block
    const int N = M/K;
    
    //accumulators
    double block_r2[STEPS];
    double acc_r2[STEPS];
    double acc_r4[STEPS];
    double unc=0;
    
    resetToZero(acc_r2, STEPS);
    resetToZero(acc_r4, STEPS);
    //cycle over blocks
    for(int i=0;i<K;i++)
    {
        if(i%10==0)
            cout << "sono al blocco " << i << endl;
        //reset block accumulator
        resetToZero(block_r2, STEPS);
        
        //for each block running N random walks
        for(int j=0;j<N;j++)
        {
            //each random walk consists in STEPS steps,
            //after each step the value of r^2 is accumulated in 
            //block_r2[step]
            walker->restart();
            for(int step=0; step<STEPS; step++)
            {
                walker->walk();
                //dividing by N to compute the mean value
                block_r2[step] += walker->getR2()/N; 
            }
        }
        
        //when all the N random walks have been performed, the array block_r2
        //contains the block averages for each step
        //updating global accumulators    
        for(int k=0; k<STEPS;k++)
        {
            acc_r2[k] += block_r2[k];
            acc_r4[k] += pow(block_r2[k], 2);
        }
    }
    
    //when all the blocks have been computed the <r^2> values can be found in acc_r2
    //for each step computing statistical uncertainty and printing values to the output file
    for(int t=0; t<STEPS;t++)
    {
        if(t>0)
            unc = sqrt((acc_r4[t]/K-pow(acc_r2[t]/K, 2))/(K-1));
        fout << t+1 << " " << acc_r2[t]/K << " " << unc << endl;
    }
    return 0;
}



void resetToZero(double* vett, int dim)
{
    for(int k=0;k<dim;k++)
            vett[k] = 0;
}

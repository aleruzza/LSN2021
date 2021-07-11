/**************************************************************************************
   
    Compute price for plain-vanilla put/call options.
    The asset price is modelled with a geometric brownian motion and sampled.
    Options price is computed from asset price estimating holder's profit and 
    discounting the risk-free interest.

    Usage:
    -------
    ./main <output> <mode> <type>
    mode: 0 - campionamento diretto, 1 - campionamento discretizzato
    type: 0 - call option, 1 - put option

    Output:
    --------
    Each row of the output files has three numerical values. In order:
    - the number i of the block considered
    - the best estimate up to the block i of the option price that is being computed
    - the statistical uncertainty of the value above

 * ************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include "random.h"
#include "utils.h"
#include "params.h"

using namespace std;

//utility constants
#define DIR_SAMPL 0
#define DISCR_SAMPL 1
#define CALL 0
#define PUT 1

void print_help();

int main (int argc, char *argv[]){
    
    //open files for output
    if(argc!=4)
        print_help();
    
    ofstream fout(argv[1]);
    
    int mode = atoi(argv[2]);
    int type = atoi(argv[3]);
    
    if(mode>1 || type>1)
        print_help();
    
    //initialize random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    const int L = M/N; //iterations per block
    int i, l, w;
    
    //direct sampling call option
    double block_acc = 0;
    double p_acc = 0;
    double p2_acc = 0;
    
    double z, s, p, unc;
    
    //cycle over blocks
    for(i=0;i<N;i++)
    {
        block_acc = 0;
        //compute l times the price
        for(l=0;l<L;l++)
        {

            //I need to sample the asset price at T
            //to estimate the call(put)-option price
            //This can be done directly or with 
            //discretized sampling
            if(mode==DIR_SAMPL)
            {
                //direct sampling
                z = rnd.Gauss(0, 1);
                s = S_0*exp((R-pow(VOL, 2)/2.)*T+VOL*z*sqrt(T));
            }
            else
                if(mode==DISCR_SAMPL)
                {
                    //discretized sampling
                    s = S_0;
                    for(w=0;w<W;w++)
                    {
                        z = rnd.Gauss(0, 1);
                        s = s*exp((R-pow(VOL, 2)/2.)*(T/W)+VOL*z*sqrt(T/W));
                    }
                }
                
            //compute the option price
            //This is the only part that differs 
            //between call and put options
            if(type==CALL)
                p = exp(-R*T)*max(0., s-K);
            else
                if(type==PUT)
                    p = exp(-R*T)*max(0., K-s);
            
            //update block accumulator
            //dividing here for numerical stability
            block_acc += p/L;
        }
        
        //block_acc contains the block mean
        
        //updating global accumulators
        p_acc += block_acc;
        p2_acc += pow(block_acc, 2);
        
        //computing uncertainty
        if(i>1)
            unc = sqrt((p2_acc/(i+1) - pow(p_acc/(i+1), 2))/i);
        
        //output to file
        fout << i+1 << " " <<p_acc/(i+1) << " " << unc << endl;
    }
    
    return 0;
}

void print_help()
{
    cout << "Usage: ./main <output> <mode> <type>" << endl;
        cout << "mode: 0 - campionamento diretto, 1 - campionamento discretizzato" << endl;
        cout << "type: 0 - call option, 1 - put option" << endl;
        exit(1);
}

/*********************************************************
Sample probability distributions from wave functions
using the Metropolis algorithm. Two different functions
can be used as transition probability distributions

Usage: 
---------
./main <transition_probability_type> <conf> <fout> [<output> <n_punti>]
 - transition_probability_type\t: u -> uniforme, g -> gaussiana
 - conf\t\t\t\t: 1 -> 1s, 2 -> 2p
 - gli argomenti tra [] sono opzionali, 
    permettono di salvare <n_punti> tra i punti campionati su <output>

Output:
---------
The output files generated from the block averaging contain 
best estimates up to the block related to the line, and its 
statistical uncertainty.
The output files used to store the sampled points contain
in each row the three coordinates of a different point

**********************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <string>
#include <algorithm>
#include "random.h"
#include "utils.h"
#include "transitionProbabilities.h"
#include "ProbabilityDist.h"

using namespace std;

//number of M(RT)^2 steps used in the simulation
#define N 1e7
//number of blocks
#define K 100
//number of M(RT)^2 steps used for equilibration
#define NEQ 100
 
#define UNIFORM 1
#define GAUSS 2

bool accept(point p1, point p2, ProbabilityDist* p, Random& rnd);
void printHelp();
void valutaRaggiDistribuzione(Random& rnd, int type, string fileName);

int main (int argc, char *argv[]){
    
    if(argc<4)
        printHelp();
    
    //initialize random number generator
    Random rnd;
    InitRandomGenerator(rnd);
    
    //the program can be run as
    //./main r <type> <file_di_output>
    //in this case the acceptance rate is studied
    //varying the parameter of the transition probability
    // specified with the argument
    //type: 1-uniform, 2-gaussian
    if(argv[1][0]=='r')
    {
        int type = atoi(argv[2]);
        if(type==UNIFORM || type==GAUSS)
        {
            string fn(argv[3]);
            valutaRaggiDistribuzione(rnd, type, fn);
        }
    }
    
    //open output file
    ofstream fout(argv[3]);

    //number of iterations per block
    const int M = N/K;
    
    //uniform distribution steps
    double r = 0;
    if(atoi(argv[2])==1) //1s config
        r = 2.5;
    else
        if(atoi(argv[2])==2) //2p config
            r = 6;
        else
            printHelp();
    

    //sigma of the gaussian distribution
    double s = 0;
    if(atoi(argv[2])==1) //1s config
        s = 0.8;
    else
        if(atoi(argv[2])==2) //2p config
            s = 1.9;
        else
            printHelp();
    
    //instantiate the transition function
    TransitionProbability* trP;
    switch(argv[1][0])
    {
        case 'u':
            trP = new UniformTP(r, rnd);
            break;
        case 'g':
            trP = new GaussTP(s, rnd);
            break;
        default:
            printHelp();
    }
    
    //istantiate the probability distribution
    ProbabilityDist* p;
    if(atoi(argv[2])==1)
        p = new HydrogenAtom1S();
    else
        if(atoi(argv[2])==2)
            p = new HydrogenAtom2P();
        else
            printHelp();
        
    //if 5 arguments are passed to the program 
    //then some sampled points must be saved
    ofstream foutPoints;
    int pointsInt = 0; //every this much points a point is saved
    bool printPoints = false;
    if(argc>5)
    {
        printPoints = true;
        foutPoints.open(argv[4]);
        pointsInt = N/atoi(argv[5]);
    }
    
    //points used during the simulation
    point p1, p2;
    
    //set the starting point
    p2 = getPoint(0,0,0);
    
    //equilibration phase, NEQ steps
    for(int i=0;i<NEQ;i++)
    {
        p1 = p2;
        p2 = trP->nextPoint(p1);
        if(!accept(p1, p2, p, rnd))
            p2 = p1;
    }
    
    //measure phase
    double blockAcc, meansAcc = 0, meansSquaredAcc = 0;
    double val, unc=0;
    long int accepted = 0, rejected = 0;
    
    for(int k=0;k<K;k++)
    {
        blockAcc = 0;
        for(int m=0; m<M; m++)
        {
            p1 = p2;
            p2 = trP->nextPoint(p1);
            if(!accept(p1, p2, p, rnd))
            {
                p2 = p1;
                rejected++;
            }
            else
                accepted++;
            
            blockAcc += getR(p2)/M;
            
            //print sampled points if requested
            if(printPoints)
            {
                if(!((k*M+m)%pointsInt))
                    foutPoints << p2.x << " " << p2.y << " " << p2.z << endl;
            }
        }
        
        meansAcc += blockAcc;
        meansSquaredAcc += pow(blockAcc, 2);
        val = meansAcc/(k+1);
        if(k>0)
            unc = sqrt((meansSquaredAcc/(k+1) - val*val)/k);
        fout << val << " " << unc << endl;
    }
    
    cout << "points accepted: " << accepted << endl;
    cout << "points rejected: " << rejected << endl;
    cout << "acceptance rate: " << accepted/(double) (accepted + rejected) << endl;
    cout << "results: r = (" << val << " +- " << unc << ")*a_0" << endl; 
    
    fout.close();
    if(printPoints)
        foutPoints.close();
    
    return 0;
}

//function used to determine wheter or not
//a point must be accepted according to the
//metropolis algorithm. Implemented for symmetrical 
//transition functions.
bool accept(point p1, point p2, ProbabilityDist* p, Random& rnd)
{
    double alpha = min(1., p->eval(p2)/p->eval(p1));
    return (rnd.Rannyu()<alpha);
}

void printHelp()
{
        cout << endl << "**********************************************************************************" << endl << endl;
        cout << "Usage: ./main <transition_probability_type> <conf> <fout> [<output> <n_punti>]" << endl << endl;
        cout << "\t - transition_probability_type\t: u -> uniforme, g -> gaussiana" << endl;
        cout << "\t - conf\t\t\t\t: 1 -> 1s, 2 -> 2p" << endl;
        cout << "\t - gli argomenti tra [] sono opzionali, \n \t  permettono di salvare <n_punti> tra i punti campionati su <output>" << endl << endl;
        cout << "**********************************************************************************" << endl << endl;
        exit(1);
}

//function used to evaluate the acceptance rate
//varying the parameter of the transition function.
//Each acceptance rate is computed using 1e6 steps
//of the Metropolis algorithm.
void valutaRaggiDistribuzione(Random& rnd, int type, string fileName)
{
    ofstream fout(fileName.c_str());
    point p1, p2, q1, q2;
    TransitionProbability* trP;
    ProbabilityDist *pr1, *pr2;
    
    pr1 = new HydrogenAtom1S();
    pr2 = new HydrogenAtom2P();
    
    p2 = getPoint(0,0,0);
    q2 = getPoint(0,0,0);
    
    
    double start = 0.1;
    double stop = 7;
    double step = 0.1;
    double r = start;
    
    cout << "Valuto i raggi per la distribuzione scelta" << endl;
    cout << "da: " << start << endl;
    cout << "a: " << stop << endl;
    cout << "step: " << step << endl;
    cout << "scrivo i risultati sul file " << fileName << endl;
    
    long int accepted1S=0, accepted2P=0;
    
    //cycle over the range of parameters 
    //that needs to be studied
    while(r < stop)
    {
        accepted1S=0, accepted2P=0;
        
        if(type==UNIFORM)
            trP = new UniformTP(r, rnd);
        else
            if(type==GAUSS)
                trP = new GaussTP(r, rnd);
        for(int i=0; i<1e6; i++)
        {
            p1 = p2;
            p2 = trP->nextPoint(p1);
            if(!accept(p1, p2, pr1, rnd))
                p2 = p1;
            else
                accepted1S++;
            
            q1 = q2;
            q2 = trP->nextPoint(q1);
            if(!accept(q1, q2, pr2, rnd))
                q2 = q1;
            else
                accepted2P++;
        }
        
        fout << r << " " << accepted1S/(double) (1e6) << " " <<  accepted2P/(double) (1e6) << endl;
        delete trP;
        r += step;
    }
    
    fout.close();
    
 exit(0);   
}

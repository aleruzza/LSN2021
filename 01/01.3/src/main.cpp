/* Exercise 01.3: Estimate Pi simulating Buffon's experiment
   Use block averaging to estimate statistical uncertainty

   Usage:
   -------
   ./main <output_file>

   Output:
   --------
   Each row of the file refers to a different block and reports the 
   best estimate up to that block with its statistical uncertainty
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

#define N 100 //number of blocks
#define M 1000000 //number of throws
#define D 0.1 //lines spacing
#define L 0.09 //bar length (must be less than D)

using namespace std;

double RandomCosine(Random& rnd);
int main (int argc, char *argv[]){
    
    const int K=M/N; //throws per block
    
    //open output files
    if(argc!=2)
    {
        cout << "Usage: ./main <output_file>" << endl;
        exit(1);
    }
    ofstream fout(argv[1]);
    
    //initialize random numbers generator
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    //Simulating Buffon's experiment
    int i,j, t, count=0;
    double center, cos, a, b, sum=0, sum_squares=0, unc=0, block_estime, last_estime;
    //cycle over blocks
    for(i=0;i<N;i++)
    {
        count=0;
        //estimate pi throwing the bar K times
        for(j=0;j<K;j++)
        {
            //simulate the throw randomly extracting the center of the bar and the angle
            center = rnd.Rannyu();
            cos  = RandomCosine(rnd);
            
            //compute the extremities 
            a = center-L*abs(cos)/2;
            b = center+L*abs(cos)/2;
            
            //check if a line has been hit
            t = b/D; //index of first line at the left of b
            if(a <= t*D)
                count++;
        }
        
        block_estime = 2*L*K/(count*D);
        sum += block_estime;
        last_estime = sum/(i+1);
        sum_squares += pow(block_estime,2);
        if(i>1)
            unc = sqrt( (sum_squares/(i+1)-pow(sum/(i+1),2))/i );
        fout << last_estime << " " << unc << endl;
        cout << "blocco " << i << " : " << last_estime << " pm " << unc << endl;
    }
    
    fout.close();

    rnd.SaveSeed();
    return 0;
}

//function used for generating a random angle 
//and thus cosine, without using pi
double RandomCosine(Random& rnd)
{
    double x, y, d;
    do
    {
        x = rnd.Rannyu(-1,1);
        y = rnd.Rannyu(-1,1);
        d = pow(x,2)+pow(y,2);
    }while(d>1);
    
    return x/sqrt(d);
}

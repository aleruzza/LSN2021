/**************************************************************************************
   Exercise 01.2 - Verify the central limit theorem generating SETS sets of K means each 
    computed extracting increasing numbers of values. The same thing is done drawing numbers 
    with different probability distributions.
 
    Usage:
    ----------
    ./main <fout_std> <fout_exp> <fout_lor>
    the three required arguments specify file names for the output of the results 
    obtained from each distribution
    
    Output:
    ----------
    For each file created every column refers to a different set and contains the 
    K means computed
 ************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "utils.h"

//number of means to compute for each set
#define K 10000

//parameters of the pdf used
#define LAMBDA 1
#define M 0
#define GAMMA 1

#define SETS 4

using namespace std;
 
int main (int argc, char *argv[]){
    
    //open files for output
    if(argc!=4)
    {
        cout << "Usage: ./main <output_std_dice> <output_exp_dice> <output_lor_dice>" << endl;
        exit(1);
    }
    ofstream fout_std(argv[1]);
    ofstream fout_exp(argv[2]);
    ofstream fout_lor(argv[3]);
    
    //initialize the random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    //N[] contains the number of throws used to compute the means in each set
    int N[SETS] = {1, 2, 10, 100};
    int k, i, j;
    double sum_std, sum_exp, sum_lor;
    
    //throwning the 3 dices simultaneously 
    //i'm using this particular structure of cycles 
    //because makes output on files, in the specified format, simpler
    for(k=0;k<K;k++)
    {
        for(i=0;i<SETS;i++)
        {
            sum_std = 0;
            sum_exp = 0;
            sum_lor = 0;
            for(j=0;j<N[i];j++)
            {
                sum_std += rnd.Rannyu()/N[i];
                sum_exp += rnd.Exp(LAMBDA)/N[i];
                sum_lor += rnd.Cauchy(M, GAMMA)/N[i];
            }
            fout_std << sum_std << " ";
            fout_exp << sum_exp << " ";
            fout_lor << sum_lor << " ";
        }
        fout_std << endl;
        fout_exp << endl;
        fout_lor << endl;
    }
    
    fout_std.close();
    fout_exp.close();
    fout_lor.close();
    
    return 0;
}

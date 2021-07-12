/**************************************************************************************
   
    Solves Travel Salesman Problem with simulated annealing.
    The algorithm is implements in the class IncubatorGA.
    In the current version the program can be used to solve TSP with
    - a randomly generated map of cities on the unitary circumference
    - a randomly generated map of cities inside a square

    Usage:
    -------
    ./main <best_path.output> [<loss_foreach_beta.out> <every_loss.out>]
    All the input parameters of the simulation are read
    from input.dat. See Init() for the order and format in
    which they have to be stored.

    Output:
    ---------
    - <best_path.output>: each row is related to a precise
        city and contains its coordinates. The order of the rows 
        indicates the path in which the cities must be visited
        to minimize the loss function, according to the solution found.
    - <loss_foreach_beta.out>: created only if specified as an argument.
        Each row refers to a different beta of the simulated annealing 
        and contains:
            * beta
            * number n of steps simulated at this beta
            * minimal loss obtained in these n steps
            * loss of the last individual after the n steps
            * acceptance rate
    - <every_loss.out>: created only if specified as an argument,
        contains the loss function value for every individual
        ever created during the simulation

 **************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "utils.h"
#include "IncubatorSA.h"
#include "TSProblem.h"
#include "TSPIndividual.h"
#include "TSPcircumference.h"
#include "TSPsquare.h"

#define DEBUG
#define CIRC 0
#define SQUARE 1
using namespace std;

int n_cities;
int map_type;
double beta_min;
double beta_max;
int spacing_type;
int n_betasteps;
int mcsteps_for_beta;
string city_map_filename;
string best_path_filename;
bool save_bestloss_foreach_beta;
string loss_foreach_beta_filename;
bool save_every_loss;
string every_loss_filename;
 
void Init();
void printHelp();

int main (int argc, char *argv[]){
    
    Random rnd;
    InitRandomGenerator(rnd);

    Init();

    if(argc<2)
        printHelp();

    ofstream of_pathEnd(argv[1]);

    cout << "Generating problem: ";
    TSProblem* problem;
    cout << city_map_filename << endl;
    if(map_type==CIRC)
    {
        cout << n_cities << " cities placed randomly on a circunference of radius 1" << endl;
        problem = new TSPcircumference(n_cities, &rnd);
    }
    else
        if(map_type==SQUARE)
        {
            cout << n_cities << " cities placed randomly inside [0;1) x [0;1)" << endl;
            problem = new TSPsquare(n_cities, &rnd);
        }

    cout << "Instanzio l'incubatore " << endl;
    IncubatorSA incubator(problem, beta_min, beta_max, spacing_type, mcsteps_for_beta, n_betasteps, rnd);
    
    if(argc>2)
        incubator.print_bestloss_foreach_beta(argv[2]);
    if(argc>3)
        incubator.print_every_loss(argv[3]);

    cout << "Evolvo" << endl;
    TSPIndividual best = incubator.evolve();

    of_pathEnd << "loss: " << best.loss() << endl;
    problem->printPathCoordinates(of_pathEnd, best.getSolution());

    of_pathEnd.close();
    return 0;
}

void Init()
{
    ifstream fin("input.dat");
    fin >> n_cities;
    fin >> map_type;
    fin >> beta_min;
    fin >> beta_max;
    fin >> spacing_type;
    fin >> n_betasteps;
    fin >> mcsteps_for_beta;
    fin.close();
}

void printHelp()
{
    cout << "Usage: ./main <best_path.output> [<loss_foreach_beta.out> <every_loss.out>]" << endl;
    exit(1);
}
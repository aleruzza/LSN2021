/**************************************************************************************
    Solves Travel Salesman Problem with a genetic algorithm.
    The algorithm is implements in the class IncubatorGA.
    In the current version the program can be used to solve TSP with
    - a randomly generated map of cities on the unitary circumference
    - a randomly generated map of cities inside a square

    Usage:
    -------
    ./main <best_path.output> [<loss_foreach_gen.out>]
    All the input parameters of the simulation are read
    from input.dat. See Init() for the order and format in
    which they have to be stored.

    Output:
    ---------
    - <best_path.output>: each row is related to a precise
        city and contains its coordinates. The order of the rows 
        indicates the path in which the cities must be visited
        to minimize the loss function, according to the solution found.
    - <loss_foreach_gen.out>: created only if specified as an argument.
        Each row refers to a different generation and contains the values
        of the loss function evaluated on every individual of that generation.

 * ************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "utils.h"
#include "IncubatorGA.h"
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
double ind_for_generation;
double q_exp;
double p_crossover;
double p_mut;
unsigned int n_generations;


void Init();
void printHelp();

int main (int argc, char *argv[]){
    
    //initialize random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    Init();

    if(argc<2)
        printHelp();

    ofstream of_pathEnd(argv[1]);

    cout << endl << "Generating problem: ";
    TSProblem* problem = NULL;
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

    cout << endl <<"Solving with a genetic algorithm with hyperparameters:" << endl;
    cout << "number of individuals for generation: " << ind_for_generation << endl;
    cout << "crossover probability: " << p_crossover << endl;
    cout << "mutation probability: " << p_mut << endl;
    cout << "index of selection power law: " << q_exp << endl;
    IncubatorGA incubator(problem, ind_for_generation, q_exp, p_crossover, p_mut, rnd);
    
    if(argc>2)
        incubator.print_loss_during_evolution(argv[2]);

    cout << endl << "Evolving for " << n_generations << " generations..." << endl;
    TSPIndividual best = incubator.evolve(n_generations);

    //print info of best individual
    cout << endl << "Best solution after " << n_generations << " generations has loss: " << best.loss() << endl;
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
    fin >> ind_for_generation;
    fin >> q_exp;
    fin >> p_crossover;
    fin >> p_mut;
    fin >> n_generations;
    fin.close();
}

void printHelp()
{
    cout << "Usage: ./main <best_path.output> [<loss_foreach_gen.out>]" << endl;
    exit(1);
}
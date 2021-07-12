/**************************************************************************************
   
    Solves Travel Salesman Problem using a genetic algorithm parallelized using 
    MPI libraries.
    Each process represents a fictitious different continent. They can exchange
    their best individuals every n_gen_exch_interval generations.

    Usage:
    -------
    mpiexec -np <n> ./main
    The mpiexec utility must be used. No arguments required.
    All the input and output files are stored in the io directory.
    The input parameters are read from the file input.dat, see Init() for 
    the order and format in which they have to be stored.

    Output:
    ---------
    - N_Rloss.output: contains the losses evalueated on every individual generated
        in the process N. Each row refers to a different generation.
    
    - best_path.output: each row is related to a precise
        city and contains its coordinates. The order of the rows 
        indicates the path in which the cities must be visited
        to minimize the loss function, according to the solution found.

 * ************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "mpi.h"
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
#define ROOT_PROC 0
using namespace std;

int n_cities;
int n_coordinates;
int map_type;
int ind_for_generation;
double q_exp;
double p_crossover;
double p_mut;
unsigned int n_generations;
unsigned int ngen_exch_interval;

void Init();
void printHelp();

int main (int argc, char *argv[]){
    
    //init MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    //getting number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //getting current process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //status 
    MPI_Status stat;

    //Initialize random number generator
    Random rnd;
    InitRandomGenerator(rnd);

    TSProblem* problem;

    ofstream of_pathEnd;

    //if root get parameters, initialize problem and broadcast everything 
    if(rank==ROOT_PROC)
    {
        //get input parameters
        Init();
    
        //generate TSP problem, only ROOT_PROC
        cout << endl << "Generating problem in root process: " << endl;
        
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

        n_coordinates = problem->get_n_coordinates();

        //displaying some informations
        cout << endl <<"Solving with a genetic algorithm with hyperparameters:" << endl;
        cout << "number of individuals for generation: " << ind_for_generation << endl;
        cout << "crossover probability: " << p_crossover << endl;
        cout << "mutation probability: " << p_mut << endl;
        cout << "index of selection power law: " << q_exp << endl;
    }
   
    //Broadcast parameters and problem
    MPI_Bcast(&n_coordinates, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&n_cities, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&ind_for_generation, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&n_generations, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&ngen_exch_interval, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&q_exp, 1, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&p_crossover, 1, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&p_mut, 1, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

    //n_generations must be divided between processes
    n_generations /= size;

    //Bcast cities coordinates
    vector<double> cities_coordinates(n_cities*n_coordinates, 0);
    if(rank==ROOT_PROC)
        cities_coordinates = problem->getCityCoordinatesFORMPI();
    MPI_Bcast(&cities_coordinates.front(), cities_coordinates.size(), MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

    //generate problem for secondary processes
    if(rank!=ROOT_PROC)
    {
        problem = new TSProblem(n_cities,n_coordinates, cities_coordinates, rnd);
    }

    //Init Incubator
    IncubatorGA incubator(problem, ind_for_generation, q_exp, p_crossover, p_mut, rnd);
    
    //change file name foreach rank
    string loss_filename = "_Rloss.output";
    string prefix = to_string(rank);
    incubator.print_loss_during_evolution(("io/" + prefix+loss_filename).c_str());

    if(rank == ROOT_PROC)
        cout << endl << "Starting evolution..." << endl;

    int iter = n_generations/ngen_exch_interval;
    TSPIndividual best;
    for(int i_iter = 0; i_iter < iter; i_iter++)
    {
        best = incubator.evolve(ngen_exch_interval);

        //now exchange best individual randomly
        //root_process determines number of exchanges and broadcast it
        int n_exch = 0;
        if(rank==ROOT_PROC)
            n_exch = rnd.Rannyu(0, 2*size);
        MPI_Bcast(&n_exch, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);

        if(rank==ROOT_PROC)
            cout << "Continent migration " << i_iter << endl;
        //generate exchanges schedule 
        int* ex_schedule_from = new int[n_exch];
        int* ex_schedule_to = new int[n_exch];
        if(rank==ROOT_PROC)
            for(int i = 0 ;i<n_exch;i++)
            {
                ex_schedule_from[i] = rnd.Rannyu(0, size);
                do
                    ex_schedule_to[i] = rnd.Rannyu(0, size);
                while(ex_schedule_to[i]==ex_schedule_from[i]);
            }
            
        MPI_Bcast(ex_schedule_from, n_exch, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);
        MPI_Bcast(ex_schedule_to, n_exch, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD);

        //now can do exchanges
        for(int i_exch=0;i_exch<n_exch;i_exch++)
        {
            if(rank==ex_schedule_from[i_exch])
            {
                vector<int> ex_data = incubator.selectIndividual().getGene();
                MPI_Send(&ex_data.front(), ex_data.size(), MPI_INTEGER, ex_schedule_to[i_exch], 1, MPI_COMM_WORLD);
            }
            else
            if(rank==ex_schedule_to[i_exch])
            {
                vector<int> ex_data(n_cities-1, 0);
                MPI_Recv(&ex_data.front(), ex_data.size(), MPI_INTEGER, ex_schedule_from[i_exch], 1, MPI_COMM_WORLD, &stat);
                incubator.injectIndividual(TSPIndividual(ex_data, problem, &rnd));
            }

        }
    }

    //decide which process has the best individual
    double* losses = new double[size];
    double loss_cp = best.loss();
    MPI_Gather(&loss_cp, 1, MPI_DOUBLE, losses, 1, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
    double i_min_loss = 0;
    if(rank==ROOT_PROC)
        for(int i=0;i<size;i++)
            if(losses[i]<loss_cp)
            {
                i_min_loss = i;
                loss_cp = losses[i];
            }
    
    int best_proc = 0;
    if(rank==ROOT_PROC)
        best_proc = i_min_loss;

    //comunicate to everyone who has the best individual
    MPI_Bcast(&best_proc, 1, MPI_INTEGER, ROOT_PROC, MPI_COMM_WORLD );

    if(rank==best_proc)
    {
        of_pathEnd.open("io/best_path.output");
        //print info of best individual
        cout << endl << "Best solution after " << n_generations << " generations (per process) has loss: " << best.loss() << endl;
        cout << "It was found in process " << best_proc << "/" << size << endl;
        of_pathEnd << "loss: " << best.loss() << endl;
        problem->printPathCoordinates(of_pathEnd, best.getSolution());

        of_pathEnd.close();
    }
    
    //closing MPI
    MPI_Finalize();
    return 0;
} 

void Init()
{
    ifstream fin("io/input.dat");
    fin >> n_cities;
    fin >> map_type;
    fin >> ind_for_generation;
    fin >> q_exp;
    fin >> p_crossover;
    fin >> p_mut;
    fin >> n_generations;
    fin >> ngen_exch_interval;
    fin.close();
}

void printHelp()
{
    cout << "Usage: ./main <city_map.output> <best_path.output> [<loss_foreach_gen.out>]" << endl;
    exit(1);
}
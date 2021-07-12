#ifndef INCUBATORGACLASS
#define INCUBATORGACLASS

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "TSProblem.h"
#include "TSPIndividual.h"
#include "random.h"

using namespace std;
class IncubatorGA{

    public:
        //costruttore: inizializzare generazione 
        IncubatorGA(TSProblem *p, int ind_for_generation, double q_exp, double p_crossover, double p_mut, Random& rnd);
        TSPIndividual evolve(unsigned int n_generations);
        void apply_selection();
        void printGeneration(int gen);
        void mutate_generation();
        void print_loss_during_evolution(string filename);
        TSPIndividual selectIndividual();
        void injectIndividual(TSPIndividual individual);

    private:
        vector<TSPIndividual> current_generation;
        vector<TSPIndividual> new_generation;
        TSProblem* problem;
        int individuals_for_generation;

        //hyperparameters
        double q = 3; //hyperparameter used in selecting individuals for breeding
        double p_c = 0.5;
        double p_m = 0.3;
        
        Random m_rnd;

        ofstream fout_loss;
        bool saveLoss = false;
        int wd = 15;
};

#endif //INCUBATORGACLASS
#ifndef INCUBATOR
#define INCUBATOR

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "TSProblem.h"
#include "TSPIndividual.h"
#include "random.h"


#define LOG_SPACING 0
#define LIN_SPACING 1

using namespace std;
class IncubatorSA{

    public:
        //costruttore: inizializzare generazione 
        IncubatorSA(TSProblem *p, double beta_start, double beta_end,
                         int spacing_type, int mcstep_for_beta, int n_betasteps, Random& rnd);
        
        TSPIndividual evolve();
        bool accept(double beta);
        
        void print_bestloss_foreach_beta(string filename);
        void print_every_loss(string filename);
       

    private:
        TSPIndividual current_ind;
        TSPIndividual current_step_best_ind;
        TSPIndividual new_ind_proposal;

        vector<double> beta;
        vector<int>  n_moves;

        TSProblem* m_problem;
        
        Random m_rnd;

        ofstream fout_loss;
        ofstream fout_bestloss;
        bool saveLoss = false;
        bool savebestLoss = false;
        int accepted;

        double current_best_loss;

        int wd = 15;
};

#endif //INCUBATOR
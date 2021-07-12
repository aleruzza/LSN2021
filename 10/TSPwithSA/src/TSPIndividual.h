#ifndef TSPINDIVIDUALCLASS
#define TSPINDIVIDUALCLASS

#include <vector>
#include "random.h"
#include <cmath>
#include <fstream>
#include "TSProblem.h"

/*
    Implementation of Individual for the Travelling salesman problem
    The cities coordinates are saved in the instance of TSPproblem referenced
    m_gene contains the proposed solution without the first (and last) city 
    which is always '0' to avoid degenaracies. getSolution() adds the '0' and gives the complete solution. 
*/

using namespace std;

class TSPIndividual
{
    public:
        TSPIndividual();
        TSPIndividual(vector<int> gene, TSProblem* problem, Random * rnd);
        TSPIndividual(TSProblem* prob, Random *rnd);
        //getter and setter
        void setGene(vector<int> gene);
        vector<int> getSolution() const;
        vector<int> getGene();

        //loss function
        double loss() const;
        bool isValid();

        //operator < needed to use algorithm::sort
        bool operator < (const TSPIndividual& ind2) const;
        
        //mutation operators
        void pairPermutationMutation();
        void shiftMutation();
        void mPermutationMutation();
        void inversionMutation();

        void randomMutation();
        
        //crossover operator
        void Crossoover(TSPIndividual &parent2);


        void printInfo();
        void printInfo(vector<int> gene);
        void printInfo(ofstream &str);

    private:
        vector<int> m_gene; 
        TSProblem* m_problem;
        Random* m_rnd;
};
#endif
#ifndef TSPROBLEMCLASS
#define TSPROBLEMCLASS
#include <vector>
#include <fstream>
#include <algorithm>
#include "random.h"
#include <iostream>
#include <iomanip>

using namespace std;

class TSProblem
{
     public:
        TSProblem(int ncities, int ncoordinates, vector<vector<double>> cities_coordinates, Random &rnd );
        TSProblem(int ncities, int ncoordinatess, Random &rnd );
        TSProblem(int ncities, int ncoordinates, vector<double> cities_coordinates, Random &rnd );
        //getter
        double getCityCoordinate(int city, int coordinate);
        int get_n_cities();
        int get_n_coordinates();

        void printPathCoordinates(ofstream &fout);
        void printPathCoordinates(ofstream &fout, vector<int> indexes);
        vector<double> getCityCoordinatesFORMPI();

     protected:
        Random m_rnd;
        int n_cities;
        int n_coordinates;
        vector<vector<double>> m_cities_coordinates;
 };

#endif
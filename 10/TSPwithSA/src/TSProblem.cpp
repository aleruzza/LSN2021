#include "TSProblem.h"
#include <iostream>

TSProblem::TSProblem(int ncities, int ncoordinates, vector<vector<double>> cities_coordinates, Random &rnd ):
        n_cities(ncities),
        n_coordinates(ncoordinates),
        m_cities_coordinates(cities_coordinates),
        m_rnd(rnd){}

TSProblem::TSProblem(int ncities, int ncoordinates, Random &rnd ):
        n_cities(ncities),
        n_coordinates(ncoordinates),
        m_rnd(rnd){}

double TSProblem::getCityCoordinate(int city, int coordinate)
{
    return m_cities_coordinates[city][coordinate];
}

int TSProblem::get_n_cities(){ return n_cities;}

int TSProblem::get_n_coordinates(){return n_coordinates;}

void TSProblem::printPathCoordinates(ofstream &fout, vector<int> indexes)
{
    for(int i=0;i<indexes.size();i++)
    {
            
        for(int i_coord=0;i_coord<n_coordinates;i_coord++)
            fout << setw(15) << getCityCoordinate(indexes[i], i_coord);
        fout << endl;
    }
}

void TSProblem::printPathCoordinates(ofstream &fout)
{
    for(int i_city=0;i_city<n_cities;i_city++)
    {
        for(int i_coord=0;i_coord<n_coordinates;i_coord++)
            fout << setw(15) << getCityCoordinate(i_city, i_coord);
        fout << endl;
    }
}
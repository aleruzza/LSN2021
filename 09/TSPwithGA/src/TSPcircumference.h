#ifndef TSPCIRC
#define TSPCIRC

#include "TSProblem.h"
#include <cmath>
#include <iostream>
class TSPcircumference : public TSProblem
{
    public:
        TSPcircumference(int n_cities, Random* rnd): TSProblem(n_cities, 2, *rnd)
        {
            vector<double> sing_coord;
            double theta;
            m_cities_coordinates.clear();
            for(int i=0;i<n_cities;i++)
            {
                sing_coord.clear();
                theta = rnd->Rannyu(0,2*M_PI);
                sing_coord.push_back(cos(theta));
                sing_coord.push_back(sin(theta));
                m_cities_coordinates.push_back(sing_coord);
            }
        }
};
#endif
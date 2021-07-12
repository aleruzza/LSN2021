#ifndef TSPSQUARE
#define TSPSQUARE

#include "TSProblem.h"

class TSPsquare : public TSProblem
{
    public:
        TSPsquare(int n_cities, Random* rnd): TSProblem(n_cities, 2, *rnd)
        {
            vector<double> sing_coord;
            m_cities_coordinates.clear();
            for(int i=0;i<n_cities;i++)
            {
                sing_coord.clear();
                sing_coord.push_back(rnd->Rannyu());
                sing_coord.push_back(rnd->Rannyu());
                m_cities_coordinates.push_back(sing_coord);
            }
        }
};
#endif
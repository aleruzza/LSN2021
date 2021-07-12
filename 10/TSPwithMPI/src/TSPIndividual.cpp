#include "TSPIndividual.h"
#include <iostream>

TSPIndividual::TSPIndividual(){}

TSPIndividual::TSPIndividual(TSProblem* prob, Random* rnd):
m_problem(prob), m_rnd(rnd)
{
    m_gene.clear();
    
    for(int i=1;i<prob->get_n_cities();i++)
    {
        m_gene.push_back(i);
    }
    
    //shuffle vector generating 20 random couples and swapping
    for(int k=0;k<200;k++)
        pairPermutationMutation();
    
}

//constructor from given data
TSPIndividual::TSPIndividual(vector<int> gene, TSProblem* problem, Random *rnd):
m_gene(gene), m_problem(problem), m_rnd(rnd)
{
    //TODO: maybe check if data in gene is valid
}

vector<int> TSPIndividual::getSolution() const
{
    //add 0 at the beginning and the end 
    vector<int> gene = m_gene;
    gene.insert(gene.begin(),0);
    gene.push_back(0);
    return gene;
}

//square sum
double TSPIndividual::loss() const
{
    vector<int> route = getSolution();
    int n_cities = m_problem->get_n_cities();
    int n_coordinates = m_problem->get_n_coordinates();


    double loss=0.0;
    double r2;
    for(int i=0;i<n_cities;i++)
    {
        r2=0.0;
        for(int j=0;j<n_coordinates;j++)
            r2 += pow(m_problem->getCityCoordinate(route[i],j)-m_problem->getCityCoordinate(route[i+1],j), 2);
        loss += r2;
    }
    return loss;
}

vector<int> TSPIndividual::getGene()
{
    return m_gene;
}

bool TSPIndividual::isValid()
{
    //to be valid the vector route obtained from gene.getm_gene() must
    // -have lenght n_cities+1
    // -start and end with '0'
    // -contain only numbers between 0 and n_cities-1 
    // -contain only one instance for each number (except for 0)
    bool result = true;
    vector<int> route = getSolution();
    int n_cities = m_problem->get_n_cities();

    //check lenght
    if(route.size()!=(n_cities+1))
        result = false;
    
    //check if starts and ends with 0
    if(route[0]!=0 || route[n_cities]!=0)
        result = false;

    //check if some number is > n_cities-1, if every city is visited 
    //and if every city is visited just once
    vector<int> counter;
    for(int i=0;i<n_cities;i++)
    {
        counter[i] = 0;
    }
    for(int i=0;i<n_cities;i++)
    {
        counter[route[i]]++;
    }
    for(int i=0;i<n_cities;i++)
    {
        if(counter[i]!=1)
            result = false;
    }
    
    return result;
}

bool TSPIndividual::operator < (const TSPIndividual& ind2) const
{
    return (loss() < ind2.loss());
}

void TSPIndividual::Crossoover(TSPIndividual &parent2)
{
    int position_of_cut = m_rnd->Rannyu(1, m_gene.size());
    vector<int> old_gene_two = parent2.getGene();
    vector<int> old_gene_one_copy = m_gene;
    vector<int> old_gene_two_copy = old_gene_two;

    vector<int> new_gene_one;
    vector<int> new_gene_two;

    #ifdef VERBOSE
        cout << "position of cut " << position_of_cut << endl;
    #endif

    unsigned int j;
    for(int i=0;i<position_of_cut;i++)
    {
        new_gene_one.push_back(m_gene.at(i));
        new_gene_two.push_back(old_gene_two.at(i));

        j=0;
        while(j<old_gene_two_copy.size() && new_gene_one[i]!=old_gene_two_copy[j])
            j++;
        if(j<old_gene_two_copy.size())
            old_gene_two_copy.erase(old_gene_two_copy.begin()+j);
        
        j=0;
        while(j<old_gene_one_copy.size() && new_gene_two[i]!=old_gene_one_copy[j])
            j++;
        if(j<old_gene_one_copy.size())
            old_gene_one_copy.erase(old_gene_one_copy.begin()+j);
    }

    for(unsigned int i=0;i<old_gene_one_copy.size();i++)
    {
        new_gene_one.push_back(old_gene_two_copy[i]);
        new_gene_two.push_back(old_gene_one_copy[i]);
    }

    #ifdef VERBOSE
        cout << "new gene one" << endl;
        printInfo(new_gene_one);
    #endif

    m_gene = new_gene_one;
    parent2.setGene(new_gene_two);
}

void TSPIndividual::setGene(vector<int> gene)
{
    m_gene = gene;
}

void TSPIndividual::printInfo()
{
    cout << "0 - ";
    for(unsigned int i=0;i<m_gene.size();i++)
        cout << m_gene.at(i) << " - ";
    cout << " 0" << endl;
}

void TSPIndividual::printInfo(ofstream &str)
{
    str << "0 - ";
    for(unsigned int i=0;i<m_gene.size();i++)
        str << m_gene.at(i) << " - ";
    str << " 0" << endl;
}

void TSPIndividual::printInfo(vector<int> gene)
{
    cout << "0 - ";
    for(unsigned int i=0;i<gene.size();i++)
        cout << gene.at(i) << " - ";
    cout << " 0" << endl;
}

void TSPIndividual::pairPermutationMutation()
{
    int a, b, t;
    a = m_rnd->Rannyu()*m_gene.size();
    b = m_rnd->Rannyu()*m_gene.size();

    //swap
    t = m_gene[a];
    m_gene[a] = m_gene[b];
    m_gene[b] = t;
}

void TSPIndividual::shiftMutation()
{
    int i_start = m_rnd->Rannyu(0, m_gene.size()-2);
    int dim = 2;
    int n_shift = m_rnd->Rannyu(0, m_gene.size()-i_start-dim);

    rotate(m_gene.begin()+i_start, m_gene.begin()+i_start+2, m_gene.begin()+i_start+2+n_shift);

    #ifdef DEBUG
    if(!isValid())
        cout << "error: generated not valid individual in shiftMutation()" << endl;
    #endif
}

void TSPIndividual::mPermutationMutation()
{
    int m = m_rnd->Rannyu(0, m_gene.size()/2);
    int i1 = m_rnd->Rannyu(0, m_gene.size()/2-m);
    int i2 = m_rnd->Rannyu(m_gene.size()/2, m_gene.size()-m);

    int t, j;
    for(j=0;j<m;j++)
    {
        t = m_gene[i1+j];
        m_gene[i1+j] = m_gene[i2+j];
        m_gene[i2+j] = t;
    }

    #ifdef DEBUG
    if(!isValid())
        cout << "error: generated not valid individual in mPermutationMutation()" << endl;
    #endif
}

void TSPIndividual::inversionMutation()
{
    reverse(m_gene.begin(), m_gene.end());
}

void TSPIndividual::randomMutation()
{
    int mut_type = m_rnd->Rannyu(1,5);
    switch(mut_type)
    {
        case 1:
            pairPermutationMutation();
            break;
        case 2:
            shiftMutation();
            break;
        case 3:
            mPermutationMutation();
            break;
        case 4:
            inversionMutation();
            break;
        default:
            break;
    }
}



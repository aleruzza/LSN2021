#include "IncubatorGA.h"


IncubatorGA::IncubatorGA(TSProblem *p, int ind_for_generation, double q_exp, double p_crossover, double p_mut, Random& rnd)
: m_rnd(rnd),
problem(p),
individuals_for_generation(ind_for_generation),
q(q_exp),
p_c(p_crossover),
p_m(p_mut)
{
    //initialize generation
    for(unsigned int i=0;i<individuals_for_generation;i++)
        current_generation.push_back(TSPIndividual(p, &rnd));
    
    //order original generation
    apply_selection();

    #ifdef VERBOSE
        for(int i=0;i<current_generation.size();i++)
        {
            current_generation.at(i).printInfo();
        }
    #endif
}


TSPIndividual IncubatorGA::evolve(unsigned int n_generations)
{
    int i_firstInd, i_secondInd;

    for(unsigned int i_gen=0;i_gen<n_generations;i_gen++)
    {
        #ifdef VERBOSE
            printGeneration(i_gen);
        #endif

        //clean new_generation vector
        new_generation.clear();
        while(new_generation.size() < individuals_for_generation)
        {
            
            i_firstInd = current_generation.size()*pow(m_rnd.Rannyu(), q);
            i_secondInd = current_generation.size()*pow(m_rnd.Rannyu(), q);
            
            if(m_rnd.Rannyu()<p_c)               
                current_generation.at(i_firstInd).Crossoover(current_generation.at(i_secondInd));

            new_generation.push_back(current_generation.at(i_firstInd));
            new_generation.push_back(current_generation.at(i_secondInd));
        }
        
        //moving new_generation in current_generation
        current_generation = new_generation;

        //apply mutations and then selection
        mutate_generation();

        apply_selection();

        //print the loss() value of every individual if requested
        if(saveLoss)
        {
            for(unsigned int i_ind=0;i_ind<current_generation.size(); i_ind++)
            {
                fout_loss << setw(wd) << current_generation[i_ind].loss();
            }
            fout_loss << endl;
        }

        if(i_gen%(n_generations/10)==0)
            cout  << "generation: " << setw(10) << i_gen << "     loss: " << setw(10) << current_generation.at(0).loss() << endl;
    }

    //the first GAIndividual in current_generation is the best
    return  current_generation.at(0);
}

void IncubatorGA::apply_selection()
{
    //must order current_generation from the best to the worst 
    //valueing each GAIndividual based on the value of problem.loss()
    //which has to be minimized
    std::sort(current_generation.begin(), current_generation.end());
}

void IncubatorGA::printGeneration(int gen)
{
    ofstream fileo;
    fileo.open("generations.out", ios::app);
    fileo << "generazione " << gen << endl;
    for(unsigned int i=0;i<current_generation.size();i++)
        current_generation.at(i).printInfo(fileo);
    fileo << endl<<endl;
}

void IncubatorGA::mutate_generation()
{
    for(unsigned int i_ind=0;i_ind<current_generation.size();i_ind++)
        {
            if(m_rnd.Rannyu()<p_m)
            {
                current_generation[i_ind].randomMutation();
            }
        }
}

void IncubatorGA::print_loss_during_evolution(string filename)
{
    saveLoss = true;
    fout_loss = ofstream(filename.c_str());
    fout_loss << "individuals for generation: " << individuals_for_generation;
    fout_loss << ",  q: " << q << ",  crossover probability: " << p_c;
    fout_loss << ",  mutation probability: " << p_m << endl;
}
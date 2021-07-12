#include "IncubatorSA.h"

IncubatorSA::IncubatorSA(TSProblem *p, double beta_start, double beta_end, int spacing_type, int mcstep_for_beta, int n_betasteps, Random& rnd)
: current_ind(p, &rnd), m_rnd(rnd), m_problem(p)
{
    //generate annealing schedule
    double beta_spacing;
    if(spacing_type==LIN_SPACING)
            beta_spacing = (beta_end-beta_start)/(double) (n_betasteps-1);
        else
            if(spacing_type==LOG_SPACING)
                beta_spacing = (log10(beta_end)-log10(beta_start))/(double) (n_betasteps-1);
    
    for(int i=0;i<n_betasteps;i++)
    {
        n_moves.push_back(mcstep_for_beta);

        if(spacing_type==LIN_SPACING)
            beta.push_back(beta_start+i*beta_spacing);
        else
            if(spacing_type==LOG_SPACING)
                beta.push_back(pow(10, log10(beta_start)+i*beta_spacing));
    }
}


TSPIndividual IncubatorSA::evolve()
{

    cout << "Inizio l'evoluzione" << endl;
    for(unsigned int i_beta=0;i_beta< beta.size();i_beta++)
    {
        accepted=0;
        current_step_best_ind = current_ind;
        current_best_loss = current_step_best_ind.loss();
        cout << "beta " << beta[i_beta] << ", mcsteps: " << n_moves[i_beta] << endl;
        for(int i_mc=0;i_mc < n_moves[i_beta]; i_mc++)
        {
            
            //propose new individual from current
            new_ind_proposal = current_ind;
            new_ind_proposal.randomMutation();

            if(accept(beta[i_beta]))
            {
                current_ind = new_ind_proposal;
                accepted++;
            }
            if(saveLoss)
                fout_loss << setw(wd) << beta[i_beta] << setw(wd) << current_ind.loss() << endl;
        }

        //print the loss() value of every individual if requested
        if(savebestLoss)
            fout_bestloss << setw(wd) << beta[i_beta] << setw(wd) << n_moves[i_beta] << setw(wd) << current_best_loss << setw(wd) << current_ind.loss() <<  setw(wd) << accepted/(double) n_moves[i_beta] << endl;
    }

    //the first GAIndividual in current_generation is the best
    return  current_step_best_ind;
}

bool IncubatorSA::accept(double beta)
{
    double new_ind_loss = new_ind_proposal.loss();
    double current_ind_loss = current_ind.loss();
    if(new_ind_loss < current_ind_loss)
    {
        if(new_ind_loss<current_best_loss)
        {
            current_best_loss = new_ind_loss;
            current_step_best_ind = new_ind_proposal;
        }
        return true;
    }
    else
    {
        double alpha = exp(-beta*(new_ind_loss-current_ind_loss));
        return m_rnd.Rannyu()<alpha;
    }
}

void IncubatorSA::print_bestloss_foreach_beta(string filename)
{
    savebestLoss = true;
    fout_bestloss = ofstream(filename.c_str());
    fout_bestloss << endl;
}

void IncubatorSA::print_every_loss(string filename)
{
    saveLoss = true;
    fout_loss = ofstream(filename.c_str());
    fout_loss << endl;
}


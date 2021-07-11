/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/***************************************************************
Run a simulation of the 1D Ising model. Both Metropolis and Gibbs 
algorithm are used to sample configurations of the system.

Usage:
--------
No arguments required. Input parameters are read from io/input.dat
If correctly specified in that file it is possible to restart a 
previous simulation from the configuration stored in config.final

Output:
---------
The following files are created
  * io/config.final - contains the final configuration
  * io/output.[feature].0 - contains all the intermediate values 
                          of block averaging computed for 
                          measuring [feature]
  * io/output.allFinal.0 - contains the final values of 
                           every feature measured
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
    Input(); //load config
    for(int isim = 0; isim < n_sim; isim++)
    {
        //set temp of the simulation
        if(n_sim!=1)
          temp = temp_i + isim*(temp_f-temp_i)/(double) (n_sim-1);
        else
          temp = temp_i;
        beta = 1./temp;
        
        cout << "Running simulation number " << isim+1 << " at temperature " << temp << endl;
        Init(); //initialization and equilibration
        if(metro)
          cout << "format: block number - acceptance rate" << endl;
          
        for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; ++istep)
            {
                Move(metro);
                Measure();
                Accumulate(); //Update block averages
            }
            Averages(iblk);   //Print results for current block
        }
        ConfFinal(); //Write final configuration
    }
    
  return 0;
}


void Input(void)
{
    ifstream ReadInput;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    
    //Read input informations
    ReadInput.open("io/input.dat");
    
    ReadInput >> temp_i;
    temp=temp_i;
    beta = 1.0/temp;

    ReadInput >> temp_f;
    
    ReadInput >> n_sim;
    
    cout << "Running " << n_sim << " simulations with temperature from " << temp_i << "K to " << temp_f << "K" << endl;
    
    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
        
    ReadInput >> metro; // if=1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;

    ReadInput >> restart;

    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();

    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input;
    if(restart)
      input.open("seed.out");
    else
      input.open("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();

    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility
    
    n_props = 4; //Number of observables 
}

//generate initial configuration, compute initial energy
//and equilibrate with eq_steps steps.
void Init()
{
    //initial configuration
    if(restart)
    {
      ifstream conf("config.final");
      for (int i=0; i<nspin; ++i)
      {
        conf >> s[i];
      }
    }
    else
    {
      for (int i=0; i<nspin; ++i)
      {
          if(rnd.Rannyu() >= 0.5) s[i] = 1;
          else s[i] = -1;
      }
    }
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
    
    //equilibration
    for(int i=0;i<eq_steps;i++)
    {
        Move(metro);
    }
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
        //evaluate acceptance probability for flipping the spin s[o]
        energy_old = Boltzmann(s[o], o);
        energy_new = Boltzmann(flip(s[o]), o);
        p = min(exp(-beta*(energy_new-energy_old)),1.);
        if(rnd.Rannyu()<p)
        {
            s[o] = flip(s[o]);
            accepted++;
        }
        attempted++;
    }
    else //Gibbs sampling
    {
        energy_up = Boltzmann(+1, o);
        energy_down = Boltzmann(-1, o);

        //p is the probability for setting s[o] to +1
        p = 1./(1+exp(beta*(energy_up-energy_down)));

        if(rnd.Rannyu()<p)
            s[o] = +1;
        else
            s[o] = -1;
    }
  }
}


double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


void Measure()
{
  double u = 0.0, m = 0.0;

  //cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    //measure of energy
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     
     //measure of magnetization
     m+= s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[im] = m;
  walker[ix] = pow(m,2);
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
  
  if(metro)
  {
    cout << setw(5) << iblk << " - " << setw(10) << accepted/attempted;
    if(!(iblk%3))
        cout << endl;
  }
  else
  {
    if(!(iblk%(nblk/20)))
    {
      cout << '#';
      cout.flush();
    }
  }

  Ene.open("io/output.ene.0",ios::app);
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);
  Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  Ene.close();
  
  Heat.open("io/output.heat.0",ios::app);
  stima_c = pow(beta,2)*(blk_av[ic]/blk_norm - pow(nspin*stima_u, 2))/(double)nspin; //Heat capacity
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  Heat.close();
  
  Mag.open("io/output.mag.0",ios::app);
  stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);
  Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  Mag.close();

  Chi.open("io/output.chi.0",ios::app);
  stima_x = beta*(blk_av[ix]/blk_norm/(double)nspin); //Magnetic susceptibility
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  Chi.close();
  
}


void ConfFinal(void)
{
  ofstream WriteConf;
  ofstream resultsMultSim;

  cout << endl << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("io/config.final");
  
  resultsMultSim.open("io/output.allFinal.0",ios::app);
  resultsMultSim << setw(wd) << temp << " ";
  resultsMultSim << setw(wd) << glob_av[iu]/(double)nblk << " " << setw(wd) << err_u << " ";
  resultsMultSim << setw(wd) << glob_av[ic]/(double)nblk << " " << setw(wd) << err_c << " ";
  resultsMultSim << setw(wd) << glob_av[im]/(double)nblk << " " << setw(wd) << err_m << " ";
  resultsMultSim << setw(wd) << glob_av[ix]/(double)nblk << " " << setw(wd) << err_x << endl;
  
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double flip(double spin)
{
    if(spin==-1)
        return 1;
    else
        return -1;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

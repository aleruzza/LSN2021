/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/****************************************************************
  Molucular dynamics simulation using Verlet algorithm
  In this implementation the following physical features are measured:
  - kinetic energy
  - potential energy
  - total energy
  - temperature
  - pressure
  - g(r)

  Final measurements are provided with uncertainties, computed
  using block averaging.
  
  Usage:
  -------
  The program requires no arguments. The simulation parameters 
  are read from io/input.dat. The initial configuration is read from 
  io/config.0, if the simulation is set to restart from old configuration
  then the configuration at t-dt is read from old.0.


  Output:
  --------
  All the output files are created in the io/ directory. 
  Can be created the folloing output files:
  * output_<property>.dat: contains every measure made of <property>
  * output_<property>_bm.dat: contains block averages and uncertainties
                              in the folowing order:
                              - number i of block
                              - block mean
                              - mean of block means from 0 to i
                              - uncertainty 
  * config.final: configuration at final t
  * old.final: configuration at final t-dt
*****************************************************************/

#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <iomanip>    
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

#define N 10 //number of blocks
int main(){ 

  //Initialization
  Input();             

  int nconf = 1;       //counter for measured configurations

  //motion integration steps 
  //block averaging implemented
  for(int iblk = 0; iblk<nblock;iblk++)
  {
    Reset();
    for(int istep=1; istep <= nstep/(double) nblock; ++istep)
    {
      Move();           //Move particles with Verlet algorithm
      if((istep+iblk*nstep/nblock)%iprint == 0) cout << "Number of time-steps: " << (istep+iblk*nstep/nblock) << endl;
      if(istep%10 == 0)
      {
        Measure();     //Properties measurement
        //ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk+1);
  }  

  ConfFinal();         //Write final configuration to restart

  return 0;
}

/*
 * Read input parameters, configure random number generator,
 *   print some information about the simulation and
 *   prepare the starting configuration
 */
void Input(void)
{ 
  //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("io/input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  //rho is the numerical density
  ReadInput >> rho;
  cout << "Numerical density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut; //null potential for r > rcut
  ReadInput >> delta; //time step for ODE integration
  ReadInput >> nstep; //total number of steps simulated
  ReadInput >> nblock; //number of blocks
  ReadInput >> iprint; //steps between output to stout
  ReadInput >> nbins; //number of bins for g(r) computation
  ReadInput >> restart; //if true restart from old config
  
  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of bins for g(r) = " << nbins << endl;
  if(restart)
      cout << "restarting from old config files" << endl;
  cout << endl;
  
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  ip = 4;
  
  n_props = 5; //Number of observables


  //Prepare everything to measure g(r)
  walker_gr = new double[nbins];
  blk_av_gr = new double[nbins];
  acc_means_gr = new double[nbins];
  acc_means2_gr = new double[nbins];
  bin_size = (box/2.0)/(double)nbins;


//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("io/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  if(restart)
  {
    //read r(t-dt)
    cout << "Read initial configuration at t-dt from file old.0 " << endl << endl;
    ReadConf.open("io/old.0");
    for (int i=0; i<npart; ++i){
        ReadConf >> xold[i] >> yold[i] >> zold[i];
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;
    }
    ReadConf.close();
    
    //compute r(t+dt) with Verlet algorithm, then compute velocities
    //in x+dt is stored in x, x is stored in xold
    Move();
    
    //rescale velocities to match desired temperature
    double sumv2=0;
    for (int i=0; i<npart; ++i) sumv2 += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    sumv2 /=(double) npart;
    double fs = sqrt(3. * temp / sumv2);  
    for (int i=0; i<npart; ++i)
    {
      
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;
    
        xold[i] = Pbc(x[i] - vx[i] * delta);
        yold[i] = Pbc(y[i] - vy[i] * delta);
        zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  else
  {
        //if the code is not restarting from a previous configuration
        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rand()/double(RAND_MAX) - 0.5;
            vy[i] = rand()/double(RAND_MAX) - 0.5;
            vz[i] = rand()/double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3. * temp / sumv2);   // fs = velocity scale factor 
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }


    //reset the accumulators
    for(int i=0;i<n_props;i++)
    {
      accumMeansSquared[i] = 0.0;
      accumMeans[i] = 0.0;
    }
    for(int i=0;i<nbins;i++)
    {
      acc_means_gr[i] = 0.0;
      acc_means2_gr[i] = 0.0;
    }

   return;
}

//Move particles with Verlet algorithm
void Move(void)
{ 
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i)
  { 
    //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i)
  { 

    //Verlet integration scheme
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

//Compute forces as -Grad_ip V(r)
double Force(int ip, int idir)
{ 
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
       // distance ip-i in pbc
      dvec[0] = Pbc( x[ip] - x[i] ); 
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      //l'if è per il cutoff del potenziale
      if(dr < rcut)
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
    }
  }
  
  return f;
}

//Properties measurement
void Measure()
{

  double t, vij, wij;
  double v=0.0, w=0.0;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;
  int bin;

  Epot.open("io/output_epot.dat",ios::app);
  Ekin.open("io/output_ekin.dat",ios::app);
  Temp.open("io/output_temp.dat",ios::app);
  Etot.open("io/output_etot.dat",ios::app);
  Pres.open("io/output_pres.dat",ios::app);
  
  v = 0.0; //reset observables
  t = 0.0;

  //reset the hystogram of g(r)
  for (int k=0; k<nbins; ++k) walker_gr[k]=0.0;

  //cycle over pairs of particles -m  only one time for pair
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //for dr > rcut particles do not interact
     if(dr < rcut)
     {
      vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
      wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
      //Potential energy
      v += vij;
      w += wij;
     }

     //update the histogram of g(r)
      if(dr < box/2)
      {
        //look for the histogram to fill
        bin = dr/bin_size;
        walker_gr[bin] += 2;
      }
    }          
  }

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  
  //update block accumulators

  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  stima_p = 48.0 * w / (3.0*(double)npart);

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  Pres << stima_p << endl;

  accumBlk[iv] += stima_pot;
  accumBlk[ik] += stima_kin;
  accumBlk[ie] += stima_etot;
  accumBlk[it] += stima_temp;
  accumBlk[ip] += stima_p;
  
  for(int i=0;i<nbins;i++)
  {
    blk_av_gr[i] += walker_gr[i];
  }

  block_norm += 1.0;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();
    return;
}

void Averages(int iblk)
{
  ofstream fileBM[5];

  //compute block averages and uncertanties

  fileBM[iv].open("io/output_epot_bm.dat",ios::app);
  fileBM[ik].open("io/output_ekin_bm.dat",ios::app);
  fileBM[it].open("io/output_temp_bm.dat",ios::app);
  fileBM[ie].open("io/output_etot_bm.dat",ios::app);
  fileBM[ip].open("io/output_pres_bm.dat",ios::app);

  double block_mean, last_mean, unc;

  for(int i=0;i<n_props; i++)
  {
    unc=0.0;
    block_mean = accumBlk[i]/block_norm;
    
    //pressure needs further operations
    if(i==4)
      block_mean = rho * temp + (block_mean) / vol;

    accumMeans[i] += block_mean;
    accumMeansSquared[i] += pow(block_mean,2);
    last_mean = accumMeans[i]/(double) iblk;

    fileBM[i] << setw(12) << iblk << setw(12) << block_mean << setw(12) << last_mean;
    
    if(iblk>1)
      unc = sqrt(((accumMeansSquared[i]/(double) iblk) - pow(last_mean,2))/(double) (iblk-1));
    
    fileBM[i] << setw(12) << unc << endl;

    fileBM[i].close();
  }
  
  //compute g(r), save only final result
  double r, gdir, err_gdir, volgdir;
  ofstream fileG;
  if(iblk==nblock)
    fileG.open("io/output_gdir_final.dat", ios::app);
  
  for(int i=0; i< nbins; i++)
  {
    r = (i)*bin_size+bin_size/2;
    volgdir = 4*M_PI*(pow(r+bin_size/2, 3)-pow(r-bin_size/2,3))/3;
    gdir = blk_av_gr[i]/(npart*rho*volgdir*block_norm);
    acc_means_gr[i] += gdir;
    acc_means2_gr[i] += gdir*gdir;
    
    if(iblk==nblock)
    {
      err_gdir = Error(acc_means_gr[i], acc_means2_gr[i], iblk);
      fileG << setw(12) << r << setw(12) << acc_means_gr[i]/(double)iblk << setw(12) << err_gdir << endl;
    }
      
  }

  fileG.close();

}

void Reset()
{
  for(int i=0;i<n_props;i++)
    accumBlk[i] = 0.0;
  for(int i=0;i<nbins;i++)
    blk_av_gr[i] = 0.0;
  block_norm = 0.0;
}

//Write final configuration
void ConfFinal(void)
{ 
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("io/config.final");

  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  
  cout << "Print final configuration of t-dt to file old.final " << endl;
  WriteConf.open("io/old.final");

  for (int i=0; i<npart; ++i)
  {
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

//Write configuration in .xyz format
void ConfXYZ(int nconf)
{ 
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//Algorithm for periodic boundary conditions with side L=box
double Pbc(double r)
{ 
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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

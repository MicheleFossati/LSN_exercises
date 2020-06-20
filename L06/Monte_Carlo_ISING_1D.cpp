/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  vector<double> temperatures;
  for (double t=0.5; t<=2; t=t+0.1)
  {
    temperatures.push_back(t);
  }
  
  for(auto it : temperatures)
  {
    temp = it;
    beta = 1./(it);
    PrintInput();
    for(int iblk=1; iblk <= nblk; ++iblk) //cicle over blocks
    {
      Reset(iblk);   //Reset block averages
      Equilibrate(100);
      for(int istep=1; istep <= nstep; ++istep) //cicle in a block
      {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //updates data for statistics
  
    }
    PrintGlobStat(nblk);
  //  ConfFinal(); //Write final configuration
  }
  rnd.SaveSeed();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;
 /* cout << "///////////////////////////////////" << endl;
  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  cout << "///////////////////////////////////" << endl << endl;
*/
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");
  
  ReadInput >> temp;
  beta = 1.0/temp;
  ReadInput >> nspin;
  ReadInput >> J;
  ReadInput >> h;
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;
  
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

  }

void PrintInput()
{
  cout << endl << "INPUT PARAMETERS" << endl;
  cout << "Temperature = " << temp << endl;
  cout << "Number of spins = " << nspin << endl;
  cout << "Exchange interaction = " << J << endl;
  cout << "External field = " << h << endl;
  if(metro==1) cout << "Metropolis moves" << endl;
  else cout << "Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Energy = " << walker[iu]/(double)nspin << endl;


}
void Equilibrate( int eq_steps)
{
    double u, m;
    const int wd = 12;
    ofstream of("equilibration.dat");
    for(int i=0; i<eq_steps; i++)
    {
        u=0;
        m=0;
        Move(metro);
        //cycle over spins
        for (int i=0; i<nspin; ++i)
        {
            u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
            m += s[i];
        }
        of << setw(wd) << u/m_spin << setw(wd) << m/m_spin << endl;
        
    }
}
void Move(int metro)
{
  int o;
 // double p, energy_old, energy_new, sm;
 // double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
        double pr_switch = exp(-beta* EnergyCost(o, s[o]));
        if ( pr_switch > 1)
        {
            s[o] = -s[o];
            accepted++;
        }
        else
            if (rnd.Rannyu() < pr_switch)
            {
                s[o] = -s[o];
                accepted++;
            }
    }
      
    else //Gibbs sampling
    {
      double p_up = 1./( 1 + exp(-beta*EnergyCost(o, 1))); // calculating p_up, so spin "goes" from -1 to 1.
      if (rnd.Rannyu() < p_up) s[o] = 1;
      else s[o] = -1;
    }
    attempted++;
  }

}

double EnergyCost(int site, double spin_val)
{
    return 2 * spin_val * ( J * (s[Pbc(site-1)] + s[Pbc(site+1)]) + h);
}
                         
double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure() // measures observables at the current step, saves data in walker
{
 // int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[im] = m;
  walker[ic] = u*u;
  walker[ix] = m*m;

// INCLUDE YOUR CODE HERE
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
   blk_norm = blk_norm + 1.0; //blk_norm contains size of block
}


void Averages(int iblk) //Print results for current block
{
    cout.flush() << iblk << '\r';
    if (metro==1) cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy in this block
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);

    stima_m = blk_av[im]/blk_norm/(double)nspin;
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);
    
    stima_c = beta*beta* (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2) ) /(double)nspin;
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);
  
    //stima_x = beta* (blk_av[ix]/blk_norm - pow(blk_av[im]/blk_norm,2) ) /(double)nspin;
    stima_x = beta* (blk_av[ix]/blk_norm  ) /(double)nspin;
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);
}

void PrintGlobStat(int iblk) //call it after Averages!
{
  ofstream of;
  const int wd=12;
  
  if (metro == 1) //metropolis --> printi in file.m
  {
    of.open("output.ene.m",ios::app);
    of << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    of.close();
    
    of.open("output.mag.m", ios::app);
    of  << setw(wd) << temp << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    of.close();
    
    of.open("output.heat.m", ios::app);
    of  << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    of.close();
    
    of.open("output.chi.m", ios::app);
    of  << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    of.close();
  }
  
  else // gibbs --> print in file.g
  {
    of.open("output.ene.g",ios::app);
    of << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    of.close();
    
    of.open("output.mag.g", ios::app);
    of  << setw(wd) << temp << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    of.close();
    
    of.open("output.heat.g", ios::app);
    of  << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    of.close();
    
    of.open("output.chi.g", ios::app);
    of  << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    of.close();
  }
  
}
void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    double err = sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2)) /(double)iblk );
    if (isnan(err) || err<0)
    {
        return 0;
    }
    else
    {
        return err;
    }
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

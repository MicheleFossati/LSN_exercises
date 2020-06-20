#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "walker.h"
#include <vector>
#include <cmath>


using namespace std;
int main(){

   //random number generator setup
	Random random_gen;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
	  Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

  ifstream input("seed.out");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			random_gen.SetRandom(seed,p1,p2);
		}

	input.close();
	} else cerr << "PROBLEM: Unable to open seed.out" << endl;

    Random* rnd_pointer;
    rnd_pointer = &random_gen;
    int n_walkers = 10000;
    int n_steps = 500;
    vector<Walker> walkers; // we will store walkers in vector
    vector<double> ave_r2(n_steps);
    fill(ave_r2.begin(), ave_r2.end(), 0.);
       
    
    for(int i=0; i<n_walkers; i++) // initializing walkers
    {
        Walker sample_walker(rnd_pointer,1, 0, 0, 0);
        walkers.push_back(sample_walker);
    }
    
    for(int i=0; i<n_steps; i++)
    {
        double step_ave2 = 0;
        for(auto it = walkers.begin(); it != walkers.end(); it++) //evolvin all walkers
        {
            (*it).lattice_evol();
            step_ave2 += (*it).R2();
        }
        ave_r2[i] = step_ave2/n_walkers; //saving average square radius
        
    }
    
    // printing on file
    ofstream output;
    output.open("random_walk.dat");
    for(auto it = ave_r2.begin(); it != ave_r2.end(); it++)
    {
        output << *it << endl;
    }
    output.close();
    
    //resetting walkers
    for(auto it = walkers.begin(); it != walkers.end(); it++)
    {
        (*it).set_pos(0., 0., 0.);
    }
    
    fill(ave_r2.begin(), ave_r2.end(), 0.); // resetting average
    for(int i=0; i<n_steps; i++) // evolving walkers in spherical-symmetric way
       {
           double step_ave2 = 0;
           for(auto it = walkers.begin(); it != walkers.end(); it++)
           {
               (*it).cont_evol();
               step_ave2 += (*it).R2();
           }
           ave_r2[i] = step_ave2/n_walkers;
           
       }
    
    output.open("random_walk_cont.dat"); //printing new data
    for(auto it = ave_r2.begin(); it != ave_r2.end(); it++)
    {
        output << *it << endl;
    }
    output.close();
    return 0;
}

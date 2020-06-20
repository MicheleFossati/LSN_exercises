#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>


using namespace std;
int main(){

   // random number generator setup
	Random rnd;
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
			rnd.SetRandom(seed,p1,p2);
		}

	input.close();
	} else cerr << "PROBLEM: Unable to open seed.out" << endl;
 // ////////////////
	
	const int n_blk=1000, n_th=100;
	ofstream of("integral.dat");
				
	for (int i_blk=0; i_blk<n_blk; i_blk++)
	{
		double integral_unif=0, integral_straight=0;
		
		for(int i_th=0; i_th<n_th; i_th++)
		{
			integral_unif += M_PI*0.5*cos(M_PI*rnd.Rannyu()*0.5);
			
			double x = rnd.Straight();
			integral_straight += M_PI*0.5*cos(M_PI*x*0.5) / (2*(1-x));
		}
		
		integral_unif = integral_unif/n_th;
		integral_straight = integral_straight/n_th;
		
		of << setw(13) << integral_unif << setw(13) << integral_straight << endl;
	}

	of.close();
	rnd.SaveSeed();
	return 0;
}

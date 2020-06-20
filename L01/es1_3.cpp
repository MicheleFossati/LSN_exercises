
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>

using namespace std;
int main(){
    //random number generator setup
    Random rnd;
     int seed[4];
     int p1, p2;
     ifstream Primes("Primes");
     if (Primes.is_open()){
        Primes >> p1 >> p2 ;
     } else cerr << "PROBLEM: Unable to open Primes" << endl;
     Primes.close();

    ifstream input("seed.in");
     string property;
     if (input.is_open())
	 {
		  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		  rnd.SetRandom(seed,p1,p2);
		  input.close();
	 }
		
	else cerr << "PROBLEM: Unable to open seed.in" << endl;

   
	const int n_blk=10000, n_th=1000, len=1, dist=2;
	int n_hit=0;
	ofstream of("buffon_pi.dat");
	
	for(int i_blk=0; i_blk<n_blk; i_blk++)
	{
		for(int i_th=0; i_th<n_th; i_th++)
		{
			// throw needle
			double x = dist*rnd.Rannyu();
			double theta = rnd.Unif0toPi();
			
			if(x + len*0.5*sin(theta) > dist || x - len*0.5*sin(theta) <= 0)
				n_hit++;
		}
	
		of << 2*len*double(n_th)/(double(n_hit)*dist) << endl;
		
		n_hit=0;
	}

	of.close();
  

	
	rnd.SaveSeed();
    return 0;
}


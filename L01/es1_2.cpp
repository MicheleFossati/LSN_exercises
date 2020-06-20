
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
// //////////
	
    int copies = 10000;
	int N;
	
    ofstream myfile;
    myfile.open ("data_2_S1.txt");
    for(int i=0; i<copies; i++){
        myfile << rnd.Rannyu() << " " << rnd.Exp(1.) << "   " << rnd.Lorentz(0.,1.) << endl;
        }
    myfile.close();
   
	
   myfile.open ("data_2_S2.txt");
   N = 2;
   for(int i=0; i<copies; i++){
	   double unif = 0, exp = 0, lor = 0;
	   for(int j=0; j<N; j++)
	   {
		   unif += rnd.Rannyu();
		   exp += rnd.Exp(1.);
		   lor += rnd.Lorentz(0., 1.);
	   }
	   unif = unif/N;
	   exp = exp/N;
	   lor = lor/N;
	   
       myfile << unif  << " " << exp << "   " << lor << endl;
	}
    myfile.close();
    
	
    myfile.open ("data_2_S10.txt");
	N = 10;
    for(int i=0; i<copies; i++){
        double unif = 0, exp = 0, lor = 0;
		for(int j=0; j<N; j++)
		{
		  unif += rnd.Rannyu();
		  exp += rnd.Exp(1.);
		  lor += rnd.Lorentz(0., 1.);
		}
		unif = unif/N;
		exp = exp/N;
		lor = lor/N;
		
		myfile << unif  << " " << exp << "   " << lor << endl;
	}
    myfile.close();
    
	
    myfile.open ("data_2_S100.txt");
	N = 100;
    for(int i=0; i<copies; i++){
		double unif = 0, exp = 0, lor = 0;
		for(int j=0; j<N; j++)
		{
		  unif += rnd.Rannyu();
		  exp += rnd.Exp(1.);
		  lor += rnd.Lorentz(0., 1.);
		}
		unif = unif/N;
		exp = exp/N;
		lor = lor/N;
		
		myfile << unif  << " " << exp << "   " << lor << endl;
	}
	
    myfile.close();
    
	rnd.SaveSeed();
    return 0;
}

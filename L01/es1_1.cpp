#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>

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
	
	// /////////////////

    int M = 100000;
    int M1 = 1000000;
    
    
    ofstream myfile;
    myfile.open ("data.txt");
    for(int i=0; i<M; i++){
        myfile << rnd.Rannyu() << endl;
        }
    myfile.close();
    
   
    myfile.open ("data1.txt");
    for(int i=0; i<M1; i++){
        myfile << rnd.Rannyu() << endl;
        }
    myfile.close();
    rnd.SaveSeed();
    return 0;
}

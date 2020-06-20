
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>


using namespace std;
int main()
{
    //random number generator setup
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

    // ///////////////////
    double S_0 = 100;
    double T = 1;
    double mu = 0.1;
    double sigma = 0.25;
    int throws = 100000;
    int n_times = 100;
    double d_t = T/n_times;
    ofstream output("price_dir.dat");
	
	//sampling GRW at time T directly, by the "right formula"
    for (int i= 0; i< throws; i++)
    {
        output << S_0*exp((mu - 0.5*pow(sigma,2))*T + sigma* rnd.Gauss(0,T)) << endl;
    }
    output.close();
    
    
    vector<double> times(n_times);
    
    output.open("price_walk.dat");
    for(int i=0; i< n_times; i++)
    {
        times[i] = static_cast<double>((i+1))*d_t;
    }

	//sampling GRW at time T "step by step"
    for(int thr =0; thr<throws; thr++)
    {
        double S_t = S_0;
        for(int t=0; t<n_times; t++) //evolve
        {
            S_t = S_t*exp( (mu - 0.5*pow(sigma,2))*d_t + sigma*rnd.Gauss(0,1)*sqrt(d_t) );
        }
        output << S_t << endl;
    }
    output.close();
    return 0;
}

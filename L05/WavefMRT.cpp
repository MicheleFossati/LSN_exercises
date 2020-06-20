#include <iostream>
#include <fstream>
#include "WavefMRT.h"

using namespace std;

int main()
{
  if(!Init()){
    cout << "Quitting..." << endl;
    return 1;
  }
  mode = "unif";
  level = "1s";
  unif_side = 2.4;
  
  Params();
  
  of.open("pos1sUnif.dat", ios::app);
  for(int istep=1; istep<=nsteps; ++istep){
    MRT();
    MeasInst();
  }
  of.close();
  
  cout << " acceptance " << double(accepted)/double(attempted) << endl << endl;
  
  rnd.SaveSeed();

 // ///////////
  
  if(!Init()){
    cout << "Quitting..." << endl;
    return 1;
  }
  mode = "gauss";
  level = "1s";
  sigma = 0.75;
  
  Params();
  
  of.open("pos1sGauss.dat", ios::app);
  for(int istep=1; istep<=nsteps; ++istep){
    MRT();
    MeasInst();
  }
  of.close();
  
  cout << " acceptance " << double(accepted)/double(attempted) << endl << endl;
  
  rnd.SaveSeed();

  // /////////////////
  if(!Init()){
    cout << "Quitting..." << endl;
    return 1;
  }
  mode = "unif";
  level = "2p";
  unif_side = 6;
  
  Params();
  
  of.open("pos2pUnif.dat", ios::app);
  for(int istep=1; istep<=nsteps; ++istep){
    MRT();
    MeasInst();
  }
  of.close();
  
  cout << " acceptance " << double(accepted)/double(attempted) << endl << endl;
  
  rnd.SaveSeed();

  // /////////////////
  if(!Init()){
    cout << "Quitting..." << endl;
    return 1;
  }
  mode = "gauss";
  level = "2p";
  sigma = 1.9;
  
  Params();
  
  of.open("pos2pGauss.dat", ios::app);
  for(int istep=1; istep<=nsteps; ++istep){
    MRT();
    MeasInst();
  }
  of.close();
  
  cout << " acceptance " << double(accepted)/double(attempted) << endl << endl;
  
  rnd.SaveSeed();
  return 0;
}

// //////////////////////////////
//     FUNCTIONS
// //////////////////////////////

bool
Init()
{
  // loading random gen
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
    }
   else {
     cerr << "PROBLEM: Unable to open seed.in" << endl;
     return false;
   }
  
  // loading params from file
  input.open("input.dat");
   
   if(!input){
     cout << "Cannot open input.dat" << endl;
     input.close();
     return false;
   }
   
   else{
     input >> x;
     input >> y;
     input >> z;
     input >> nsteps;
     
     accepted = 0;
     attempted = 0;
   }
  input.close();
   return true;
}

void
Params()
{
  cout << "PARAMETERS" << '\n'
       << "initial posizion = (" << x << ", " << y << ", " << z << ")" << '\n'
       << mode << " mode" << '\n';
  
  if (mode == "unif")
    cout  << "unif side = " << unif_side << '\n';
  if (mode == "gauss")
    cout  << "sigma = " << sigma << '\n';
    
   cout << "sampling " << level << " level" << '\n'
        << "n. of steps = " << nsteps << '\n' ;
}

void
Run()
{
  
}
  
void
MRT(){
  double x_step, y_step, z_step;

  if( mode == "unif"){
    x_step = x + unif_side*(rnd.Rannyu()-0.5) ;
    y_step = y + unif_side*(rnd.Rannyu()-0.5) ;
    z_step = z + unif_side*(rnd.Rannyu()-0.5) ;
  }
  
  else if( mode == "gauss"){
    x_step = x + rnd.Gauss(0,sigma) ;
    y_step = y + rnd.Gauss(0,sigma) ;
    z_step = z + rnd.Gauss(0,sigma) ;
  }
  
  else {
    cout << "insert valid mode: unif or gauss" << endl;
    return;
  }
  
  double r = EuclNorm(x,y,z);
  double r_s = EuclNorm(x_step,y_step,z_step);
  double q;

  if (level == "1s") q = exp(-2*(EuclNorm(x_step, y_step, z_step) - EuclNorm(x,y,z)));
  
  else if (level == "2p") q = exp(r-r_s) * z_step/z * z_step/z;
  else{
    cout << "insert valid level: 1s or 2p"  << endl;
    return;
  }
  
  if ( q > 1.) {
    x = x_step;
    y = y_step;
    z = z_step;
    accepted++;
  }
  else if (rnd.Rannyu() < q){
    x = x_step;
    y = y_step;
    z = z_step;
    accepted++;
  }
  attempted++;
}


void
MeasInst()
{
  r_inst = EuclNorm(x,y,z);
  r_block += r_inst;
  
  of << setw(wd) << x << setw(wd) << y << setw(wd) << z << setw(wd) << r_inst << endl;
  
}

void
MeasBlock()
{
  ofstream of("r_block.dat", ios::app);
  of << r_block/nsteps << endl;
  
  //reset
  attempted = 0;
  accepted = 0;
  r_block = 0;
}

double
EuclNorm(double mx, double my, double mz)
{
  return sqrt(mx*mx +my*my + mz*mz);
}

// //////////////////////////
// Variational Monte Carlo //
//    Michele Fossati      //
// //////////////////////////
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"

using namespace std;


//random n generator
Random rnd;

// position
double pos = 0;

//simulation params
int nsteps=10000, nblocks;
double len, mu, sigma;

// other params
const int wd = 13;

//statistics
int attempted=0, accepted=0;
double E_block=0;

//functions for simulation
bool InitRnd();
void Params();
void Equilibrate();
void MRT();
void MeasInst();
void PrintResult(ostream &);

//useful functions
double Psi(const double &);
double Pot(const double &);
double Acceptance();
  

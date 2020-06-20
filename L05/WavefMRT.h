#include <string>
#include "random.h"
using namespace std;

//random n generator
Random rnd;

//ofstream
ofstream of;

// position
double x, y, z;

//simulation params
int nsteps;
string mode, level;
double unif_side, sigma;

// other params
const int wd = 13;

//statistics
int attempted, accepted;
double r_inst, r_block;

//functions for simulation
bool Init();
void Params();
void MRT();
void MeasInst();
void Run();

//useful functions
double EuclNorm(double, double, double);
  



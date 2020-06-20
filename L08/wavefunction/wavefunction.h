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

class point
{
public:
	point(){};
	point(double mx, double my, double mz){ x=mx; y=my; z=mz;}
	~point(){};
	
	double x;
	double y;
	double z;
	
	double r() const {return sqrt(x*x + y*y + z*z);}
};


//random n generator
Random rnd;



//ofstream
//ofstream of_block("block.dat");

// position
double pos = 0;

//simulation params
int nsteps, nblocks;
double len, mu, sigma;

// other params
const int wd = 13;

//statistics
int attempted=0, accepted=0;
double pos_block, E_block;

//functions for simulation
bool InitRnd();
void Params();
void Equilibrate();
void MRT();
void MeasInst(ostream &);
void PrintResult(ostream &);

//useful functions
double Psi(const double &);
double Pot(const double &);
double Acceptance();
  

point operator-(const point &);
point operator+(const point &, const point &);
point operator-(const point &, const point &);
point operator+(const point &, const double &);
point operator-(const point &, const double &);
ostream &operator<< (ostream &, const point &);
double operator*(const point &, const point &);




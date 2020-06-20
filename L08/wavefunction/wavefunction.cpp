#include "wavefunction.h"
using namespace std;

int main(int argc, char * argv[])
{
  if(!InitRnd()){
    cout << "Bye bye..." << endl;
    return 1;
  }
	
	// loading mu-s and sigma-s
	if (argc != 4)
	{
		cout << "Error. Insert valid arguments" << endl;
		cout << "Syntax: mu, sigma, nstep" << endl;
		return 1;
	}
	
	mu = atof(argv[1]);
  sigma = atof(argv[2]);
	nsteps = atof(argv[3]);
	len = 2;

	// open output stream
	ofstream outObs("wavef.dat");

	Params();
	Equilibrate();
	for(int istep=1; istep<=nsteps; ++istep)
	{
		MRT();
		MeasInst(outObs);
	}
	
	cout << "acceptance " << Acceptance() << endl << endl;
	//PrintResult(outPar);
			
	


	// closing routine
	outObs.close();
	//of_block.close();
  rnd.SaveSeed();

  return 0;
}



// /////////////////////////////////////////
//              FUNCTIONS
// /////////////////////////////////////////

bool
InitRnd()
{
  // loading random gen
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
       Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
     while ( !input.eof() ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
     }
       
     input.close();
		 return true;
    }
   else {
     cerr << "PROBLEM: Unable to open seed.in" << endl;
     return false;
   }
}

void
Params()
{
  cout << "PARAMETERS" << '\n'
	<< "initial posizion = " << pos << '\n'
	<< "mu = " << mu << '\n'
	<< "sigma = " << sigma << '\n'
	<< "nsteps = " << nsteps << '\n' ;
}

void
Equilibrate()
{
	bool loop = true;
	while (loop)
	{
		for(int i=0; i<500; i++) MRT();
		
		if (Acceptance() < 0.45) len = 0.9*len;
		else if (Acceptance() > 0.55)	len = 1.1*len;
		else loop = false;
		
		cout << "Equilibrating... acceptance = " << Acceptance() << endl;
		accepted=0;
		attempted=0;
		pos=0;
	}
}

void
MRT(){
	double pos_step = pos + len*(rnd.Rannyu()-0.5);

	double psir = Psi(pos);
	double psir_step = Psi(pos_step);

  double q = (psir_step * psir_step) / (psir * psir);
  
  if ( q > 1.) {
		pos = pos_step;
    accepted++;
  }
  else if (rnd.Rannyu() < q){
		pos = pos_step;
    accepted++;
  }
  attempted++;
}


void
MeasInst(ostream & os)
{
	double x_plus = (pos + mu) / sigma;
	double x_minus = (pos - mu) / sigma;
	double T_psi = 0.5/(sigma*sigma) * (
										 Psi(pos)
									 - pow(x_minus,2)* exp(-0.5 * pow(x_minus,2))
									 - pow(x_plus,2)* exp(-0.5 * pow(x_plus,2))
										 );
	
	double E_inst = T_psi/Psi(pos) + Pot(pos);
  
	// print inst
  os << pos << setw(wd) << E_inst << endl;
	
  pos_block += pos;
	E_block += E_inst;
  
}

void
PrintResult(ostream & os)
{
	//print mu, sigma and Energy
  os << setw(wd) << mu << setw(wd) << sigma << setw(wd) << E_block/nsteps << endl;
  
  //reset
  attempted = 0;
  accepted = 0;
	
	pos = 0;
  pos_block = 0;
	E_block = 0;
}

double
Psi(const double & x)
{
	return exp(-(x - mu)*(x - mu) / (2*sigma*sigma)) +
				 exp(-(x + mu)*(x + mu) / (2*sigma*sigma));
}

double Pot(const double & x)
{
	return pow(x,4) - 2.5*pow(x,2);
}

double Acceptance()
{
	return double(accepted)/double(attempted);
}

// wrong! 3D was not required!
/*
// /////////////////////////////////////////////////
//         basic operators
// /////////////////////////////////////////////////

// opposite
point operator-(const point & a)
{
	point p(- a.x, - a.y, - a.z);
	return p;
}

// sum between points
point operator+(const point & a, const point & b)
{
	point p(a.x + b.x, a.y + b.y, a.z + b.z);
	return p;
}

point operator-(const point & a, const point & b)
{
	point p(a.x - b.x, a.y - b.y, a.z - b.z);
	return p;
}

point operator+(const point & a, const double & scalar)
{
	point p(a.x + scalar, a.y + scalar, a.z + scalar);
	return p;
}

point operator-(const point & a, const double & scalar)
{
	point p(a.x - scalar, a.y - scalar, a.z - scalar);
	return p;
}

// scalar product
double operator*(const point & a, const point & b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

ostream &operator<< (ostream &os, const point & p)
{
	os  << setw(wd) << p.x << setw(wd) << p.y << setw(wd) << p.z << setw(wd) << p.r() ;
	return os;
}

*/

#include "optimize.h"
using namespace std;

int main(int argc, char * argv[])
{
  if(!InitRnd()){
    cout << "Bye bye..." << endl;
    return 1;
  }
	
	// loading mu-s and sigma-s
	if (argc != 6)
	{
		cout << "Error. Insert valid arguments" << endl;
		cout << "Syntax: mu_min, mu_max, sigma_min, sigma_max, param_step" << endl;
		return 1;
	}
	
	double mu_min = atof(argv[1]);
	double mu_max = atof(argv[2]);
	double sigma_min = atof(argv[3]);
	double sigma_max = atof(argv[4]);
	double param_step = atof(argv[5]);

	
	vector<double> mus, sigmas;
	for(int i=0; i< param_step; i++)
	{
		mus.push_back(mu_min + i*(mu_max-mu_min)/ param_step);
		sigmas.push_back(sigma_min + i*(sigma_max-sigma_min)/param_step);
	}
	
	// open output stream
	ofstream outPar("mu_sigma_E.dat");
	
	for(int i_mu=0; i_mu < param_step; i_mu++)
	{
		for(int i_sigma=0; i_sigma < param_step; i_sigma++)
		{
			mu = mus[i_mu];
			sigma = sigmas[i_sigma];
			len = 1;
			Params();
			Equilibrate();
			for(int istep=1; istep<=nsteps; ++istep)
			{
				MRT();
				MeasInst();
			}
			
			cout << "acceptance " << Acceptance() << endl << endl;
			PrintResult(outPar);
			
		}
	}

	// closing routine
	outPar.close();
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

   ifstream input("seed.out");
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
MeasInst()
{
	double x_plus = (pos + mu) / sigma;
	double x_minus = (pos - mu) / sigma;
	double T_psi = 0.5/(sigma*sigma) * (
										 Psi(pos)
									 - pow(x_minus,2)* exp(-0.5 * pow(x_minus,2))
									 - pow(x_plus,2)* exp(-0.5 * pow(x_plus,2))
										 );
	
	double E_inst = T_psi/Psi(pos) + Pot(pos);
  
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

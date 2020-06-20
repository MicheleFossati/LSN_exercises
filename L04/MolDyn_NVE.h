/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000, nbins=100, print_wd =13;
int n_props;
int iv,ik,it,ie, igofr;
double gofr_norm[nbins];

// averages
double acc,att;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_gofr[nbins];
double block_pot=0, block_kin=0, block_etot=0, block_temp=0, block_gofr[nbins];
std::ofstream blockprint("block_obs.dat"),
         instprint("inst_obs.dat"),
         gofrprint("block_gofr.dat");
//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut, bin_size;

// simulation
int nstep, nblock, iprint, seed, intervalmeas;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void PrintInputPar();
void CorrVel();
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(bool);
void PrintBlock();
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  for(int iblock=1; iblock <= nblock; ++iblock){
    for(int istep=1; istep <= nstep; ++istep){
    Move();           //Move particles with Verlet algorithm
      if(istep%iprint == 0) cout.flush() << "Block n." << iblock << " step n." << istep << " of " << nstep << '\r';
      if(istep%intervalmeas == 0){
        Measure(true);     //Properties measurement
       // ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
      }
    }
    PrintBlock();
  }
  
  ConfFinal();         //Write final configuration to restart
  blockprint.close();
  instprint.close();
  gofrprint.close();
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;
  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblock;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> intervalmeas;
  
  ReadInput.close();

  PrintInputPar();
  
//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  //reset measure value
  for(int i=0;i<nbins;i++){
    block_gofr[i] =0 ;
    stima_gofr[i] = 0;
  }
  
//evaluate normalization factors for g
  bin_size = (box/2.0)/(double)nbins;
  for(int i=0; i<nbins; i++){
    gofr_norm[i] = rho*npart*4*pi/3*pow(bin_size,3) * (pow(i+1,3) - pow(i,3));
  }
  
  // here we print the values of r for g(r)
  for(int i=0; i<nbins; i++){
    gofrprint << setw(print_wd) << bin_size*(i+0.5);
  }
  gofrprint << endl;
  
  
//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  ReadConf.open("old.0");
  if(!ReadConf){
    cout << "Cannot open old.0, preparing random velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
  else{
    cout << "Read old configurations from old.0" << endl;
    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();
    cout << "Correcting velocities..." << endl;
    CorrVel();
  }
  
  return;
}

void PrintInputPar()
{
 /* cout << "INFO" << endl;
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
  cout << "The program uses Lennard-Jones units " << endl; */
  cout << "PARAMS:" << endl;
  cout << "Temperature = " << temp << endl;
  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "Volume of the simulation box = " << vol << endl;
  //cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of blocks = " << nblock << endl;
  cout << "Number of steps for block = " << nstep << endl;
  cout << "Measure every " << intervalmeas << " steps" << endl << endl;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void CorrVel()
{
  double alpha; //correction factor
  
  Move();
  Measure(false);
  cout << "T after one step: " << stima_temp << endl;
  alpha = temp/stima_temp;
  for (int i=0; i< npart; i++)
  {
    vx[i] = vx[i]*sqrt(alpha); //correct velocities
    vy[i] = vy[i]*sqrt(alpha);
    vz[i] = vz[i]*sqrt(alpha);
    
    xold[i] = Pbc(x[i] - vx[i]*delta); //correct positions
    yold[i] = Pbc(y[i] - vy[i]*delta);
    zold[i] = Pbc(z[i] - vz[i]*delta);
  }
  Measure(false); //measure but not print on file!
  cout << "T corrected: " << stima_temp << endl;
}

void Measure(bool print){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  t = 0.0;
  //reset the hystogram of g(r)
  for (int i=0; i<nbins; ++i) stima_gofr[i]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

    //update of the histogram of g(r), maybe you should restrict r to r_cut
     if(dr<box/2)
       stima_gofr[static_cast<int>(dr/bin_size)] +=2;
 
    //Potential energy
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij;
     }
    }
  }
  
  //normalize histogram
  for(int i=0; i<nbins;i++) stima_gofr[i] = stima_gofr[i]/gofr_norm[i];

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
 
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle


  if(print==true){
    block_pot += stima_pot;
    block_kin += stima_kin;
    block_etot += stima_etot;
    block_temp += stima_temp;
	  
    for(int i=0;i<nbins;i++) block_gofr[i] += stima_gofr[i];
 
    instprint << setw(13) << stima_pot << setw(13) << stima_kin << setw(13) << stima_etot << setw(13) << stima_temp << endl;
  }
  return;
}

void PrintBlock()
{
	double pot = block_pot/(double)nstep*(double)intervalmeas;
	double kin = block_kin/(double)nstep*(double)intervalmeas;
	double etot = block_etot/(double)nstep*(double)intervalmeas;
	double temp = block_temp/(double)nstep*(double)intervalmeas;
  
	blockprint << setw(13) << pot << setw(13) << kin << setw(13) << etot << setw(13) << temp << endl;
  
  for(int i=0;i<nbins;i++){
    block_gofr[i] = block_gofr[i]/(double)nstep*intervalmeas;
    gofrprint << setw(print_wd) << block_gofr[i];
  }
  gofrprint << endl;
  
  // reset block ave
  block_pot=0;
  block_kin=0;
  block_etot=0;
  block_temp=0;
  for(int i=0;i<nbins;i++) block_gofr[i] =0;
}
void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

 // cout << "Print final configuration to file config.final and old.final " << endl << endl;
  
  WriteConf.open("config.final");
  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  
  WriteConf.open("old.final");
  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "random.h"
#include "mpi.h"

using namespace std;


class city
{
public:
	city() {x=0; y=0;};
	city (double i_x, double i_y): x(i_x), y(i_y) {};
	~city(){};
	
	double x;
	double y;
};

using path = vector<city>;

bool operator==(const city &, const city &);
city operator-(const city &, const city &);
double operator*(const city &, const city &); //scalar product
double Length(const path &);
bool operator < (const path &, const path &);

std::ostream &operator<< (ostream &, const path &);
std::ostream &operator<< (ostream &, const vector<path> &);





int main(int argc, char ** argv)
{
        MPI::Init(argc,argv);
        int size = MPI::COMM_WORLD.Get_size();
        int rank = MPI::COMM_WORLD.Get_rank();

	int ncities, npaths, ngen, nmigr;
	int rankA, rankB;
	string mode;
	double length, prob_swap, prob_shift, prob_reverse, prob_sex, sup_ave_len=0;
	
	path cities, best_path;
	vector<path> popul;
		
	// loading random gen
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open())
	{
		for(int i=0; i<rank; i++) Primes.ignore(20, '\n');
		Primes >> p1 >> p2 ;
	}
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
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
		return 1;
	}

	//cout << "rank" << rank << "primes: " << p1 << " " << p2 << endl;

	// loading params
	input.open("input.dat");
	input >> ncities;
	input >> npaths;
	input >> ngen;
	input >> nmigr;
	input >> mode;
	input >> length;
	input >> prob_swap;
	input >> prob_shift;
	input >> prob_reverse;
	input >> prob_sex;
	input.close();
	
	// generating cities. they must be the same for all the cores
	vector<double> cities_x;
	vector<double> cities_y;
	if (rank == 0)
	{
	if (mode == "CIRCLE"){
		for(int i=0; i<ncities; i++) {
			double theta = 2*M_PI*rnd.Rannyu();
			city m_city( length*cos(theta), length*sin(theta) );
			cities.push_back(m_city);
		}
	}

	else if (mode == "SQUARE"){
	for (int i=0; i<ncities; i++) {
	city m_city(length*rnd.Rannyu(), length*rnd.Rannyu());
	cities.push_back(m_city);
	}
	}

	else {
	cout << "Insert valid mode CIRCLE or SQUARE" << endl;
	return 1;
	}
	


	for(int i=0; i<cities.size(); i++)
	{
		cities_x.push_back(cities[i].x);
		cities_y.push_back(cities[i].y);
	}

	}

	if(rank !=0)
	{
		cities_x.resize(ncities);
		cities_y.resize(ncities);
	}

	MPI_Bcast(cities_x.data(), cities_x.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(cities_y.data(), cities_y.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank !=0)
	{
		for(int i=0; i<ncities; i++)
		{
			city m_city(cities_x[i], cities_y[i]);
			cities.push_back(m_city);
		}
	}	


	ofstream of_cities(to_string(rank) + "cities.dat"); 
	of_cities << cities  << endl;
	of_cities.close();
	
	// generating paths
	path sample_path = cities;

	popul.push_back(sample_path);
	for(int i=1; i< npaths; i++)
	{
		random_shuffle(next(sample_path.begin()), sample_path.end(), [&rnd](int n) {
			return static_cast<int>(rnd.Rannyu(0, n));
			});
		
		popul.push_back(sample_path);
	}
		
	
	// DEBUG
	//cout << "POPULATION" << endl;
	//cout << popul;
	
	// statistics
	ofstream of(mode + to_string(rank) + "_len_best_ave.dat");

	best_path = *min_element(popul.begin(), popul.end());
	sort(popul.begin(), popul.end()); //using overloaded  <
	for(int i=0; i<npaths/2; i++)
	{
		sup_ave_len += Length(popul[i]);
	}
	sup_ave_len = sup_ave_len/ (npaths/2);
	of <<std::setw(13) << Length(best_path) << setw(13) << sup_ave_len << endl;

	// **************************
	//        swag loop
	// **************************

	vector<path> new_popul;

	for(int igen=0; igen<ngen; igen++)
	{
		cout << "GENERATION #" << igen+1 << '\r';
		
		//
		// check
		//
		for(int i=0; i<npaths; i++)
		{
			if(!is_permutation(cities.begin(), cities.end(), popul[i].begin()))
			{
				cout << "Invalid path! popul[" << i << "]" << endl;
				cout << popul[i] << endl;
				return 1;
			}
		}
		cout << "All paths are correct!" << '\r';
		
		//
		// evolve
		//
		sort(popul.begin(), popul.end()); //using overloaded  <
		
		if(igen%10==0 && igen!=0)
		{
			if(rank==0)
			{
				rankA = int(rnd.Rannyu(0,size));
				do rankB = int(rnd.Rannyu(0,size));
				while (rankA == rankB);
			}

			MPI_Bcast(&rankA, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&rankB, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
		
			ofstream Info(to_string(rank) + "info.dat", ios::app);
			Info << "rank " << rank << " gen " << igen << " rankA " << rankA << " rankB " << rankB << endl;
			
			MPI_Status statAx, statAy, statBx, statBy;
			vector<double> sendx, sendy, recvx, recvy;
			int tagAx = 0;
			int tagAy = 1;
			int tagBx = 2;
			int tagBy = 3;
			if(rank==rankA || rank==rankB)
			{
				for(int i=0; i<ncities; i++)
				{
					sendx.push_back(popul[0][i].x);
					sendy.push_back(popul[0][i].y);
				}

				recvx.resize(ncities);
				recvy.resize(ncities);
			}
			
			//exchanging x coordinates
			if(rank==rankA)
			{
				MPI_Send(sendx.data(), sendx.size(), MPI_DOUBLE, rankB, tagAx, MPI_COMM_WORLD);
				MPI_Recv(recvx.data(), recvx.size(), MPI_DOUBLE, rankB, tagBx, MPI_COMM_WORLD, &statBx);
			} 

			if(rank==rankB)
			{
			 	MPI_Send(sendx.data(), sendx.size(), MPI_DOUBLE, rankA, tagBx, MPI_COMM_WORLD);
                                MPI_Recv(recvx.data(), recvx.size(), MPI_DOUBLE, rankA, tagAx, MPI_COMM_WORLD, &statAx);
 
			}
			
			MPI_Barrier(MPI_COMM_WORLD);

			//exchanging y coordinates
			if(rank==rankA)
			{
				MPI_Send(sendy.data(), sendy.size(), MPI_DOUBLE, rankB, tagAy, MPI_COMM_WORLD);
				MPI_Recv(recvy.data(), recvy.size(), MPI_DOUBLE, rankB, tagBy, MPI_COMM_WORLD, &statBy);

			}

			if(rank==rankB)
			{
				MPI_Send(sendy.data(), sendy.size(), MPI_DOUBLE, rankA, tagBy, MPI_COMM_WORLD);
				MPI_Recv(recvy.data(), recvy.size(), MPI_DOUBLE, rankA, tagAy, MPI_COMM_WORLD, &statAy);
			}

			if(rank==rankA || rank==rankB)
			{
				// reconstructing path from recv_x and recv_y
				path r_path;
				for(int i=0; i<ncities; i++)
				{
					city r_city(recvx[i], recvy[i]);
					r_path.push_back(r_city);
				}

				popul[0] = r_path;
			}
			
			Info.close();
		}
			
		new_popul.clear();
		
		while(new_popul.size() != npaths)
		{
			int par1, par2, cut_pos;
			path child1(ncities), child2(ncities);

			//extract two different parents
			par1 = int(npaths* pow(rnd.Rannyu(), 2));
						   
			do par2 = int(npaths* pow(rnd.Rannyu(), 2));
			while (par1 == par2);
			
			// DEBUG
			//cout << "parents: " << par1 << " and " << par2 << endl;
			//cout << "bare childrens: " << child1 << endl << "and" << endl << child2 << endl;
			
			//reset childrens
			for(int i=0; i<ncities; i++)
			{
				child1[i].x = 0;
				child1[i].y = 0;
				child2[i].x = 0;
				child2[i].y = 0;
			}

			//crossing chromosomes
			if (rnd.Rannyu() < prob_sex)
			{
				cut_pos = ncities*rnd.Rannyu();

				// copy first part of chromosomes
				copy(
					popul[par1].begin(),
					popul[par1].begin() + cut_pos,
					child1.begin()
				);
				
				copy(
					popul[par2].begin(),
					popul[par2].begin() + cut_pos,
					child2.begin()
				);
				
				// copy second part in "the order of the other chromosome"
				copy_if(
					popul[par1].begin(),
					popul[par1].end(),
					child2.begin() + cut_pos,
					[&child2, cut_pos]( city n)
					{ return std::find(child2.begin(), child2.begin() + cut_pos, n) == child2.begin() + cut_pos; } //using overloaded ==
				);
				
				copy_if(
					popul[par2].begin(),
					popul[par2].end(),
					child1.begin() + cut_pos,
					[&child1, cut_pos]( city n)
					{ return std::find(child1.begin(), child1.begin() + cut_pos, n) == child1.begin() + cut_pos; }
				);
			}
			
			else
			{
				copy(popul[par1].begin(),
					 popul[par1].end(),
					 child1.begin()
					 );
				
				copy(popul[par2].begin(),
				popul[par2].end(),
				child2.begin()
				);
			}
			
			//DEBUG: printing childrens
		/*	cout << "want to add two children" << endl;
			if(! (is_permutation(cities.begin(), cities.end(), child1.begin()) &&
			   is_permutation(cities.begin(), cities.end(), child2.begin())) )
			{
				cout << "error!" << endl << "par1:" << endl << popul[par1] << endl << "par2:" << endl << popul[par2] << endl << "cut_pos = " << cut_pos << endl << "child1:" << endl << child1 << endl << "child2:" << endl << child2 << endl;
			}
			else cout << "they are ok" << endl;
		*/
			
			new_popul.push_back(child1);
			new_popul.push_back(child2);
			
			child1.clear();
			child2.clear();
			
		//	cout << "added! new_popul.size() = " << new_popul.size() << endl << endl;
		}
		
		//elitary
		new_popul.pop_back();
		new_popul.push_back(popul[0]);
		popul = new_popul;
		
		//DEBUG
		//cout << "new population: " << endl << new_popul;
		
		//
		// mutate
		//
		for(int ipath=0; ipath<npaths; ipath++) {
		   
		   if (rnd.Rannyu() < prob_swap) // swap
		   {
			 int i = int(rnd.Rannyu(1, ncities));
			 int j;
		   
			 do j = int(rnd.Rannyu(1, ncities));
			 while(j==i);
			 
			 swap(popul[ipath][i], popul[ipath][j]);
			 
		   }
		   
		   if (rnd.Rannyu() < prob_shift)
		   {
			 int start = int(rnd.Rannyu(1,ncities-1));
		 int width = int( rnd.Rannyu(1, (ncities-start)/2) );
			 int jump = int(rnd.Rannyu(width+1, (ncities-start)/2) );
			 for( int k=start; k<start+width; k++){
			   swap(popul[ipath][k], popul[ipath][k+jump]);
			 }
		   }
		 
		   if (rnd.Rannyu() < prob_reverse)
		   {
			 int start = int(rnd.Rannyu(1, ncities -1 ));
			 int end = 1 + int(rnd.Rannyu(start+1, ncities));
			 //cout << start << " " << end-1 << endl;
			 reverse( popul[ipath].begin() + start, popul[ipath].begin() + end);
		   }
		 }
	
		
		// statistics
		best_path = *min_element(popul.begin(), popul.end());
		sort(popul.begin(), popul.end()); //using overloaded  <
		for(int i=0; i<npaths/2; i++)
		{
			sup_ave_len += Length(popul[i]);
		}
		sup_ave_len = sup_ave_len/ (npaths/2);
		of << setw(13) << Length(best_path) << setw(13) << sup_ave_len << endl;
		
	}

	of.close();
	
	if (rank == 0) cout << mode + to_string(rank) + "_best_path.dat" << endl;
	ofstream best(mode + to_string(rank) + "_best_path.dat");
	best << best_path;
	best.close();
	

	rnd.SaveSeed();
	MPI::Finalize();
	return 0;
}



// /////////////
//  functions //
// /////////////

bool operator==(const city & city1, const city & city2)
{
	if(city1.x == city2.x && city1.y == city2.y)
		return true;
	
	else return false;
}


std::ostream &operator<< (ostream &os, const path & m_path)
{
	for(int i=0; i<m_path.size(); i++)
	{
		os  << setw(13) << m_path[i].x << setw(13) << m_path[i].y << '\n';
	}
	
	os << setw(13) << m_path[0].x << setw(13) << m_path[0].y << '\n';
	//os << "len = " << Length(m_path);

	return os;
}

std::ostream &operator<< (ostream & os, const vector<path> & vec)
{
	for(int j=0; j<vec.size(); j++)
	{
		os << vec[j] << endl << endl;
	}
	return os;
}

double Length(const path & mpath)
{
	double len=0;
	int size = mpath.size();
	for(int i=1; i<size; i++)
	{
		len += (mpath[i] - mpath[i-1])*(mpath[i] - mpath[i-1]);
	}
	
	len += (mpath[0] - mpath[size-1])*(mpath[0] - mpath[size-1]);

	return len;
}

city operator-(const city & city1, const city & city2)
{
	city new_city(city1.x - city2.x, city1.y - city2.y);
	return new_city;
}

double operator*(const city & city1, const city & city2)
{
	return city1.x*city2.x + city1.y*city2.y;
}

bool operator < (const path & path1, const path & path2)
{
	return Length(path1) < Length(path2);
}

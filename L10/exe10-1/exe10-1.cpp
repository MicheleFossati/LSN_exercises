#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
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
	int ncities, nstep, attempted=0, accepted=0, temp_size;
	string mode;
	double length, temperature, cool_rate, accept_thres, acceptance;
	
	path actual_path, new_path, best_path;
		
	// loading random gen
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open())
	{
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

	// loading params
	input.open("input.dat");
	input >> ncities;
	input >> mode;
	input >> length;
	input >> temp_size;
	input >> cool_rate;
	input >> accept_thres;
	
	
	input.close();
	
	// generating cities
	if (mode == "CIRCLE"){
		for(int i=0; i<ncities; i++) {
			double theta = 2*M_PI*rnd.Rannyu();
			city m_city( length*cos(theta), length*sin(theta) );
			actual_path.push_back(m_city);
		}
	}

	else if (mode == "SQUARE"){
		for (int i=0; i<ncities; i++)
		{
			city m_city(length*rnd.Rannyu(), length*rnd.Rannyu());
			actual_path.push_back(m_city);
		}
	}

	else {
		cout << "Insert valid mode CIRCLE or SQUARE" << endl;
		return 1;
	}
	
	best_path = actual_path;
	
	//cout << "STARTING PATH:" << endl;
	//cout << actual_path << endl << endl;
	
	temperature = 2*Length(actual_path);
	
	ofstream of(mode + "_length.dat");
	// **************************
	//        swag loop
	// **************************

	nstep = 1;
	do
	{
		new_path = actual_path;
		switch(int(rnd.Rannyu()*3))
		{
			case (0): // swap
			{
				int i = int(rnd.Rannyu(1, ncities));
				int j;
			  
				do j = int(rnd.Rannyu(1, ncities));
				while(j==i);
				
				swap(new_path[i], new_path[j]);
				break;
			}
							
			case (1):
			{
				int start = int(rnd.Rannyu(1,ncities-1));
				int width = int( rnd.Rannyu(1, (ncities-start)/2) );
				int jump = int(rnd.Rannyu(width+1, (ncities-start)/2) );
				for( int k=start; k<start+width; k++){
				  swap(new_path[k], new_path[k+jump]);
				}
				break;
			}
				
			case(2):
			{
				int start = int(rnd.Rannyu(1, ncities -1 ));
				int end = 1 + int(rnd.Rannyu(start+1, ncities));
				reverse( new_path.begin() + start, new_path.begin() + end);
			}
		}
		
		double len_diff = Length(new_path) - Length(actual_path);
		if (len_diff < 0)
		{
			accepted++;
			actual_path = new_path;
		}
		
		else if (rnd.Rannyu() < exp( - len_diff / temperature))
		{
			accepted++;
			actual_path = new_path;
		}

		attempted++;
		acceptance = double(accepted) / double(attempted);
		
		if(nstep % temp_size == 0)
		{
			of << Length(actual_path) << endl;
			cout << "Temp = " << temperature << "  acceptance = " << acceptance << endl;
			
			attempted = 0;
			accepted = 0;
			temperature = cool_rate*temperature;
					
		}

		nstep++;
	//	cout << "nstep = " << nstep << endl;
	//	cout << "Temp = " << temperature << "  acceptance = " << acceptance << endl;

	}
	while (attempted < temp_size-1 || acceptance > accept_thres );
	//while (nstep /1000 < 100);


	of.close();
	
	of.open(mode + "_best_path.dat");
	of << actual_path << endl;
	of.close();
				 

	rnd.SaveSeed();
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

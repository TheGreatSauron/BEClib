#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <complex>
#include <valarray>
#include <fstream>
#include <mpi.h>
#include <time.h>

#include <armadillo>
#define QICLIB_DONT_USE_NLOPT
#include <QIClib>
//#include "Bose_Hubbard_RK45.h"
#include <BEClib>

//using namespace std;
//using namespace arma;
//using namespace qic;
//using namespace bec;

// =============================================================================
// Function for mpi processor distribution
void para_range(int n1, int n2, int &isize, int &irank, int &istart, int &iend) 
	{
	int iwork1, iwork2;
	iwork1 = (n2 - n1) / isize;
	iwork2 = (( n2 - n1) % isize);
	istart = irank * iwork1 + n1 + std::min(irank, iwork2);
	iend = istart + iwork1;
	if (iwork2 > irank) iend = iend + 1;
	}

// =============================================================================
int main(int argc, char *argv[])
{
    // Parameters
	int L = 60;
	int len_args = 2;
	double dt = 1e-3;
	double acc = 1e-6;

    // Initial state
	std::valarray<comp> v(2*L*L);
	for (int i = 0; i < int(v.size()); i++)
	{
		v[i] = double(std::rand() % 1000) / 1000;
	}

	double args[len_args][6] = {
		{ 1.0 * L * L, 1.0, 0.0, 0.1, 0.1, 0.2 },
		{ 1.0 * L * L, 1.0, 0.0, 0.1, 0.1, 0.5 } };
		
	// MPI directives
	int size, rank, i_start, i_end;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	para_range(0, len_args, size, rank, i_start, i_end);
	
	for (int n = i_start; n < i_end; n++)
	{
		bec::RK45_2D_Spin rk(dt, acc, v, L, args[n]);

		rk.groundState();

		std::ofstream file;
		std::string filename = ("Data/G_L(" + std::to_string(L) +
			                    ")_n(" + std::to_string(args[n][0] / (L * L)) +
			                    ")_U(" + std::to_string(args[n][2]) +
			                    ")_tz(" + std::to_string(args[n][3]) +
			                    ")_tso(" + std::to_string(args[n][4]) +
			                    ")_Mz(" + std::to_string(args[n][5]) + ").txt");

		std::cout << filename << '\n';
		file.open(filename);
        	file << std::scientific;

        	for (int i = 0; i < rk.getVector().size(); i++)
		{
			comp z = rk.getVector()[i];
			if (z.imag() >= 0)
			{
				file << "(" << z.real() << "+" << z.imag() << "j) "<< '\n';
			}
			else
			{
				file << "(" << z.real() << "-" << std::abs(z.imag()) << "j) "<< '\n';
			}
		}


		file.close();
		
// =============================================================================
// Code to calculate Bogoliubov excitations and non-condensate densities

        bec::noncondensate(rk.getVector())  // This function does not exist as of now
        
// =============================================================================
// Code to re-calculate condensate density again using non-condensate densities

        bec::condensate(non-condensate density)  // This function does not exist as of now

// =============================================================================
// Code to check convergence etc.
        
// =============================================================================
	}
	MPI_Finalize();

	return 0;
}

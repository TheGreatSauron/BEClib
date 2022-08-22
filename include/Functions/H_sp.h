#ifndef H_SP_H
#define H_SP_H

#include <armadillo>
#include <string>

namespace bec
{

// Define complex type as comp (equivalent to arma::cx_double)
typedef std::complex<double> comp;
// Define imaginary unit macro
#ifndef i1
#define i1 comp(0.0, 1.0)
#endif

// Define Pauli matrices
static arma::Mat<comp> sig0(2, 2, arma::fill::eye);
static arma::Mat<comp> sig1 = { {comp(0.0), comp(1.0)}, {comp(1.0), comp(0.0)}};
static arma::Mat<comp> sig2 = { {comp(0.0), -i1}, {i1, comp(0.0)} };
static arma::Mat<comp> sig3 = { {comp(1.0), comp(0.0)}, {comp(0.0), comp(-1.0)}};

inline arma::Mat<comp> H_sp(std::string Lat, double tso, double tz, double Mz, double thop, double k1, double k2, double k01, double k02)
{
	double h0, h1, h2, h3;

	if (Lat == "Sq")
	{
		h1 = tso * (cos(k1) - cos(k2));
		h2 = tso * (cos(k1 + k2) - cos(k1 - k2));
		h3 = Mz + 2 * tz * (cos(k1) + cos(k2));
		h0 = -2 * thop * (cos(k1 - k01) + cos(k2 - k02));
	}
	else if (Lat == "Tri")
	{
		h1 = tso * (2 * cos(k1 + k2) - cos(k1) - cos(k2));
		h2 = -sqrt(3) * tso * (cos(k1) - cos(k2));
		h3 = Mz + 2 * tz * (cos(k1) + cos(k2) + cos(k1 + k2));
		h0 = -2 * thop * (cos(k1 - k01) + cos(k2 - k02) + cos(k1 + k2 - k01 - k02));
	}
	else return arma::Mat<comp>(2, 2, arma::fill::zeros);

	arma::Mat<comp> Hk = h1 * sig1 + h2 * sig2 + h3 * sig3 + h0 * sig0;
	return Hk;
}

}

#endif // H_SP_H
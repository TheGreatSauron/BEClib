#pragma once

#include <complex>
#include <vector>
#include <string>
#include <valarray>

typedef std::complex<double> comp;

class RK45 {
public:

	RK45(double step_size, double accuracy, std::valarray<comp> initial);

	virtual std::valarray<comp> func(std::valarray<comp> y1);

	int full_step(comp factor = comp(1.0, 0.0));

	std::valarray<comp> getVector() const;
	double getStepSize() const;
	double getNorm() const;
	double getTime() const;
	double getAcc() const;
	double getRMS() const;

	void setAcc(double accuracy);
	void setTime(double time);

	std::string printNorms() const;
	std::string printVector() const;
	std::string printPhases() const;

protected:

	std::valarray<comp> y;

	void setStepSize(double step_size);

private:

	double t;
	double h;
	double acc;

	static double B[5][5];
	static double CH[6];
	static double CT[6];
};
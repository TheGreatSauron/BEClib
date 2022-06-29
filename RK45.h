#ifndef RK45_H
#define RK45_H

#include <complex>
#include <vector>
#include <string>
#include <valarray>

// Define complex type as comp
typedef std::complex<double> comp;
// Define imaginary unit macro
#define i1 comp(0.0, 1.0)

// Class which handles base RK45 functions, inherited for particular5 problems and further applications
class RK45 {
public:
	
	// Standard constructor
	RK45(double step_size, double accuracy, std::valarray<comp> initial);

	// The defined function which is to be integrated, 
	// returns a valarray of the function for all sites in the vector
	virtual std::valarray<comp> func(std::valarray<comp> y1);

	// Takes one full step of integration using RK45 method,
	// factor argument is multiplied to the derivative function,
	// returns the number of steps taken including failed steps
	int full_step(comp factor = comp(1.0, 0.0));

	// Returns the current vector of values
	std::valarray<comp> getVector() const;
	// Returns the step size h
	double getStepSize() const;
	// Returns the total norm of the whole vector
	double getNorm() const;
	// Returns the current time
	double getTime() const;
	// Returns the accuracy
	double getAcc() const;
	// Returns the RMS [sqrt(<x^2> - <x>^2)] value for the vector,
	// note that this only works if the vector is one-dimensional in real space
	double getRMS() const;

	// Sets a new accuracy
	void setAcc(double accuracy);
	// Sets a new time
	void setTime(double time);

	// Returns a string of the norms of all terms in the vector
	std::string printNorms() const;
	// Returns a string of all terms in the vector
	std::string printVector() const;
	// Returns a string of the phases of all terms in the vector
	std::string printPhases() const;

protected:

	// The vector of all  current values being integrated
	std::valarray<comp> y;

	// Sets a new step_size,
	// protected since it can have negative effects if done in the middle of an integration.
	void setStepSize(double step_size);

private:

	// The time
	double t;
	// The step size
	double h;
	// The accuracy
	double acc;

	// Constants required for the integration
	static double B[5][5];
	static double CH[6];
	static double CT[6];
};

#endif // RK45_H
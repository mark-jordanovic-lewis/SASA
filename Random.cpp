#include "math.h"

#include "sasa_RandomGenerator.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

Random::Random(int aSeed) : seed(aSeed) {};

long double Random::Uniform()
{
	seed ^= MASK;
	int k = seed/IQ;
	seed = IA*(seed - k * IQ) - IR * k;
	if (seed<0)
		seed += IM;
	long double ans = AM * seed;
	seed ^= MASK;

	return ans;
};

long double Random::Normal()
{
	static bool iset = true;
	static long double gset = 0;
	long double v1, v2;
	long double r = 2.0;

	if(iset) {
		do {
			v1= 2.0*Uniform()-1.0;
			v2= 2.0*Uniform()-1.0;
			r = v1*v1 + v2*v2;
		} while (r>=1);

		long double fac = sqrt(-2.0*log(r)/r);

		gset = v1 * fac;
		iset = false;
		return v2 * fac;
	}
	else {
		iset = true;
		return gset;
	}

};

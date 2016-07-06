class Random
{
private:
	int seed;

public:
	Random(int aSeed);
	long double Uniform();
	long double Normal();
};

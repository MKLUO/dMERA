#include <iostream>
#include <fstream>
#include <vector>

#include <cmath>
#include <random>

#include "Dmera.h"

const int SITES = 10;

const double DELTA = 0.0;
const double DISORDER = 5.0;

const int EPOCH = 20;

int main()
{
	std::uniform_real_distribution<double> jdis(0., 1.);
	std::random_device rd;
	std::default_random_engine re = std::default_random_engine(rd());

	std::vector<double> j;

    for (int i = 0; i < SITES; ++i)
		j.push_back(std::pow(jdis(re), DISORDER));

	j = {1, 2, 3, 1, 12, 3, 2, 5, 7, 9};

	Dmera case1(j, DELTA);

	case1.check();

	for (int i = 0; i < EPOCH; ++i)
	{
		std::cout << "EPOCH\t" << i + 1 << " : ";
		case1.VarUpdate();
	}

	for (int i = 1; i < SITES/2 + 1; ++i)
		std::cout << "\nEntropy of length " << i << " : " << case1.AverageEntropy(i);

	std::cout << std::endl;

	return 0;
}

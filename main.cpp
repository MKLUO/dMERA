#include <iostream>
#include <fstream>
#include <vector>

#include <cmath>
#include <random>

#include "Dmera.h"

const int SITES = 30;

const int SAMPLES = 20;

//const double DELTA = 0.0;		//Delta-z
const double DISORDER = 1.0;	//Delta-J

std::vector<double> xx_case_entropy(double, int);

int main()
{
	std::vector<double> totalEntropys(SITES/2, 0.);

    for (int i = 0; i < SAMPLES; ++i)
	{
		std::cout << "Sample " << i + 1 << ",\t" << std::flush;
		std::vector<double> result(xx_case_entropy(DISORDER, SITES));

		for (int j = 0; j < SITES/2; ++j)
			totalEntropys[j] += result[j];
	}
	for (int j = 0; j < SITES/2; ++j)
		totalEntropys[j] /= SAMPLES;

	for (int j = 0; j < SITES/2; ++j)
		std::cout << "Entropy of disorder = " << DISORDER << ", " << SITES << " sites, " << SAMPLES << " samples:\t" << totalEntropys[j] << std::endl;

	return 0;
}

std::vector<double> xx_case_entropy(double disorder, int sites)
{
	std::uniform_real_distribution<double> jdis(0., 1.);
	std::random_device rd;
	std::default_random_engine re = std::default_random_engine(rd());

	std::vector<double> j;
    for (int i = 0; i < sites; ++i)
		j.push_back(std::pow(jdis(re), disorder));

	Dmera case1(j, 0.);

	//case1.check();
	case1.VarUpdateForEpochs();

	std::vector<double> entropys;
	for (int i = 1; i < sites/2 + 1; ++i)
	{
		std::cout << "|" << std::flush;
		entropys.push_back(case1.AverageEntropy(i));
	}
	std::cout << std::endl;

}

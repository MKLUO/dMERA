#include <iostream>
#include <fstream>
#include <vector>

#include <random>

#include "Dmera.h"

const int SITES = 2;
const int EPOCH = 500;

int main()
{
	
	std::uniform_real_distribution<double> jdis(0., 10.);
	std::random_device rd;
	std::default_random_engine re = std::default_random_engine(rd());

	std::vector<double> j;

    for (int i = 0; i < SITES; ++i)
		j.push_back(jdis(re));

	Dmera case1(j, 0.01);

	std::cout << case1.summary(false) << std::endl;

	std::ofstream of("energy.gp");
/*
    for (int i = 0; i < EPOCH; ++i)
	{
		std::cout << "EPOCH: " << i + 1 << std::endl;
		of << case1.VarUpdate() << std::endl;
	}
*/
	of.close();

	case1.check();

	return 0;
}

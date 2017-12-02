#include <iostream>
#include <fstream>
#include <vector>

#include <random>

#include "Dmera.h"

int main()
{
	
	std::uniform_real_distribution<double> jdis(0., 10.);
	std::random_device rd;
	std::default_random_engine re = std::default_random_engine(rd());

	std::vector<double> j;

    for (int i = 0; i < SITES; ++i)
		j.push_back(jdis(re));

	Dmera case1(j, 0.01);


	for (int i = 0; i < EPOCH; ++i)
	{
		std::cout << "EPOCH\t" << i + 1 << std::endl;
		case1.VarUpdate();
	}

	case1.check();
 



	return 0;
}

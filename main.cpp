#include <iostream>
#include <fstream>
#include <vector>

#include <random>

#include "Dmera.h"

int main()
{
	// TODO: disorder J distribution	
	std::uniform_real_distribution<double> jdis(0., 10.);
	std::random_device rd;
	std::default_random_engine re = std::default_random_engine(rd());

	std::vector<double> j;

    for (int i = 0; i < SITES; ++i)
		j.push_back(jdis(re));


	//std::vector<double> j(SITES, 1.0);
//	std::vector<double> j = {0.1, 0.1, 4.0, 0.1, 3.0, 0.1, 2.0, 0.1};


	Dmera case1(j, DELTA);


	for (int i = 0; i < EPOCH; ++i)
	{
		std::cout << "EPOCH\t" << i + 1 << std::endl;
		case1.VarUpdate();
	}

	case1.check();
 



	return 0;
}

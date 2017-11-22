#include <iostream>
#include <vector>

#include "Dmera.h"

int main()
{
	std::vector<double> j({1., 0.1, 22., 0.1, 5., 0.2, 33., 0.3});
	Dmera case1(j, 0.01);

	std::cout << case1.summary() << std::endl;

	return 0;
}

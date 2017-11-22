#include <vector>

#include "Dmera.h"

int main()
{
	std::vector<double> j({1., 0.1, 2., 0.1, 5., 0.2, 3., 0.3});
	Dmera case1(j, 0.01);
	
	return 0;
}

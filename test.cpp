#include <iostream>
#include <fstream>
#include <vector>

#include <random>

#include "Dmera.h"

int main()
{
	std::vector<uni10::Complex> elem = 
						{ 1., 1., 0., 0.,
					    2., 1., 0., 0.,
		   				0., 0., 10., 0.,
					    0., 0., 0., 1.	};
	uni10::Matrix val(4, 4);
	
	val.setElem(elem);	

	for (int i = 0; i < 16; ++i)
	{
		std::cout << val(i);
	}


	return 0;

}

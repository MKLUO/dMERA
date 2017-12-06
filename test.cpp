#include <iostream>
#include <fstream>
#include <vector>

#include <random>

#include "Dmera.h"

int main()
{
	uni10::Complex val[] = 
		{	-1., 0., 0., 0.,
		    0., 1., 0., 0.,
		    0., 0., -10., 0.,
		    0., 0., 0., -1.	};

	uni10::Matrix M(4, 4);
	M.setElem(val);

	uni10::Bond b_in(uni10::BD_IN, 2);
	uni10::Bond b_out(uni10::BD_OUT, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_in);
	B.push_back(b_out);
	B.push_back(b_out);

	uni10::UniTensor T(uni10::CTYPE, B);
	T.setRawElem(M);
/*
	uni10::Matrix M2 = T.getRawElem();

	for (int i : {0, 5, 10, 15})
		M2(i) = double((0. < M2(i).real()) - (M2(i).real() < 0.));

	T.setRawElem(M2);

	std::cout << T;

*/
	double ES;
	T = Dmera::eigenshift(T, ES);
	std::cout << T;

	double E;
    uni10::UniTensor Tm = Dmera::svdSolveMinimal(T, E);

	std::cout << Tm;

	std::cout << ((T.getRawElem()) * (Tm.getRawElem())).trace(uni10::CTYPE) << std::endl;

	return 0;

}

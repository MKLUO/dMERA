#include <vector>
#include <random>
#include <utility>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

#include <iostream>

#include <uni10.hpp>

#include "Dmera.h"

Dmera::Dmera(std::vector<double> js, double delta): width(js.size())
{                                  
	//iteratively contract nodes

	while (js.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < js.size(); ++i)
			if (js[i] > js[idx])
				idx = i;

		blocks.push_back(new Block(idx));

		int left	= (idx == 0)? (js.size() - 1): (idx - 1);
		int right	= (idx == js.size() - 1)? 0: (idx + 1);

		js[left] = Dmera::effective_j(js[idx], js[left], js[right], delta);

		js.erase(js.begin() + idx);
		js.erase(js.begin() + right);			
	}

	//UniTensor initialize

	for (Blcok* b : blocks)
	{
		b->t_l = Dmera::Random_Unitary();
		b->t_r = Dmera::Random_Unitary();
		b->t_u = Dmera::Random_Unitary();
	}

    BuildNetworkForms();

    /*
	den_mat_descend[0]	= new uni10::Network("NetworkSheets/dmd0");
	den_mat_descend[1]	= new uni10::Network("NetworkSheets/dmd1");
	den_mat_descend[2]	= new uni10::Network("NetworkSheets/dmd2");
	den_mat_descend[3]	= new uni10::Network("NetworkSheets/dmd3");
	den_mat_descend[4]	= new uni10::Network("NetworkSheets/dmd4");

	eff_ham_ascend[0]	= new uni10::Network("NetworkSheets/eha0");
	eff_ham_ascend[1]	= new uni10::Network("NetworkSheets/eha1");
	eff_ham_ascend[2]	= new uni10::Network("NetworkSheets/eha2");
	eff_ham_ascend[3]	= new uni10::Network("NetworkSheets/eha3");
	eff_ham_ascend[4]	= new uni10::Network("NetworkSheets/eha4");

	energy_u[0]			= new uni10::Network("NetworkSheets/eu0"); 
	energy_u[1]			= new uni10::Network("NetworkSheets/eu1"); 
	energy_u[2]			= new uni10::Network("NetworkSheets/eu2"); 

	energy_l[0]			= new uni10::Network("NetworkSheets/el0"); 
	energy_l[1]			= new uni10::Network("NetworkSheets/el1"); 
	energy_l[2]			= new uni10::Network("NetworkSheets/el2"); 
	energy_l[3]			= new uni10::Network("NetworkSheets/el3"); 

	energy_r[0]			= new uni10::Network("NetworkSheets/er0"); 
	energy_r[1]			= new uni10::Network("NetworkSheets/er1"); 
	energy_r[2]			= new uni10::Network("NetworkSheets/er2"); 
	energy_r[3]			= new uni10::Network("NetworkSheets/er3"); 

	svd_restore			= new uni10::Network("NetworkSheets/svd"); 
	*/
}

Dmera::VarUpdate()
{
	// Build descending density matrix

	// TODO: push 2 singlet dm at first

	std::vector<uni10::UniTensor> dm = dm_init;

	for (int i = blocks.size() - 1; i >= 0; --i)
	{
		Block* b = blocks[i];
		int idx = b->get_idx();
		bool boundry = (idx == (dm.size() + 1));
		
		if (boundry) // boundry contraction
			for (int i = 0; i < 3; ++i)
				b->dm[i] = dm[dm.size() + i - 2]; 
		else
			for (int i = 0; i < 3; ++i)
				b->dm[i] = dm[Dmera::index(idx + i - 2, dm.size())];  

		if (boundry) // boundry contraction
		{
			dm.insert(dm.begin(), 0);
			dm.push_back(0);
		}
		else
			dm.insert(dm.begin + idx, 2, 0);

		for (int i = 0; i < 5; ++i)
			dm[Dmera::index(idx + i, dm.size())] = network["dmd"][i].launch(b);

	}

	// Variational update & Build ascending effective hamiltonian

	std::vector<uni10::UniTensor> eh = eh_init;

	for (int i = 0; i < blocks.size(); ++i)
	{
		Block* b = blocks[i];
		int idx = b->get_idx();
		bool boundry = (idx == (eh.size() - 1));

		for (int i = 0; i < 5; ++i)
			b->eh[i] = eh[Dmera::index(idx + i - 2, eh.size())];


		if (boundry)
		{
			eh.erase(eh.begin());
			eh.pop_back();
		}
		else
			eh.erase(eh.begin() + idx, eh.begin() + idx + 2);



		if (boundry)
		{
			for (int i = 0; i < 3; ++i)
				eh[Dmera::index(i - 2, eh.size())] = network["eha"][i].launch(b);
		}
		else
			for (int i = 0; i < 3; ++i)
				eh[Dmera::index(idx + i - 2, eh.size())] = network["eha"][i].launch(b);

}

double Dmera::effective_j(double j, double j_l, double j_r, double delta) { return j_l * j_r / j / (1. + delta); }

int Dmera::index(int idx, int size)
{
	if (idx < 0) return idx + size;
	else if (idx >= size) return idx - size;
	else return idx;
}

uni10::UniTensor Dmera::Random_Unitary()
{
	uni10::Bond b_in(uni10::BD_IN, 2);
	uni10::Bond b_out(uni10::BD_OUT, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_in);
	B.push_back(b_out);
	B.push_back(b_out);

	uni10::UniTensor T(RCTYPE, B);

	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_int_distribution<int> unif(1, 50);

	for (int i = 0; i < unif(dre); ++i)
		T.orthoRand();

	return T;
}

uni10::UniTensor Dmera::Singlet()
{

	uni10::Complex m[] = {	0.,				1. / sqrt(2.),
							-1. / sqrt(2.),	0., };

	uni10::Bond b_in(uni10::BD_IN, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_in);

	uni10::UniTensor T(RCTYPE, B);

	T.setRawElem(m);

	return T;
}

uni10::UniTensor Dmera::Identity()
{

	uni10::Complex m[] = {	1.,	0., 0., 0.,
							0., 1., 0., 0.,
							0., 0., 1., 0.,
							0., 0., 0., 1. };

	uni10::Bond b_in(uni10::BD_IN, 2);
	uni10::Bond b_out(uni10::BD_OUT, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_in);
	B.push_back(b_out);
	B.push_back(b_out);

	uni10::UniTensor T(RCTYPE, B);

	T.setRawElem(m);

	return T;
}

uni10::UniTensor Dmera::Identity2()
{

	uni10::Complex m[] = {	1.,	0.,
							0., 1. };

	uni10::Bond b_in(uni10::BD_IN, 2);
	uni10::Bond b_out(uni10::BD_OUT, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_out);

	uni10::UniTensor T(RCTYPE, B);

	T.setRawElem(m);

	return T;
}

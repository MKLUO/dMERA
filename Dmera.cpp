#include <array>
#include <vector>
#include <random>
#include <utility>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <uni10.hpp>

#include "Dmera.h"

const std::string TEMP_FNAME("NetworkSheets/_temp");
const std::string SVD_FNAME("NetworkSheets/SVD");
const std::string DMS_0_FNAME("NetworkSheets/DMS_0");
const std::string DMS_1_FNAME("NetworkSheets/DMS_1");
const std::string PNAME("NetworkSheets/");

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

		int left	= Dmera::index(idx - 1, js.size());
		int right	= Dmera::index(idx + 1, js.size());

		js[left] = Dmera::effective_j(js[idx], js[left], js[right], delta);

		if (idx == js.size() - 1)
		{
			js.erase(js.begin());
			js.pop_back();
		}
		else
			js.erase(js.begin() + idx, js.begin() + idx + 2);
	}

	BuildNetworkForms();

}

Dmera::Block::Block(int idx_): idx(idx_)
{
	t["l"] = Dmera::Random_Unitary();
	t["r"] = Dmera::Random_Unitary();
	t["u"] = Dmera::Random_Unitary();
	t["s"] = Dmera::Singlet();

	// TODO: dm, eh initialize
}

uni10::UniTensor Dmera::Block::tensor(std::string name) const
{
	if		(name == "l")	return t["l"];
	else if	(name == "r")	return t["r"];
	else if	(name == "u")	return t["u"];
	else if	(name == "s")	return t["s"];

	else if	(name == "dm0")	return dm[0];
	else if	(name == "dm1") return dm[1];
	else if	(name == "dm2") return dm[2];

	else if	(name == "eh0") return eh[0];
	else if	(name == "eh1") return eh[1];
	else if	(name == "eh2") return eh[2];
	else if	(name == "eh3") return eh[3];
	else if	(name == "eh4") return eh[4];
}

void Dmera::BuildNetworkForms()
{
	network["dmd"].resize(5);
	network["eha"].resize(5);
	network["enl"].resize(4);
	network["enu"].resize(3);
	network["enr"].resize(4);

	network["dmd"][0] = NetworkForm("DM", 0, 0, 0);
	network["dmd"][1] = NetworkForm("DM", 1, 0, 1);
	network["dmd"][2] = NetworkForm("DM", 2, 0, 1);
	network["dmd"][3] = NetworkForm("DM", 3, 0, 1);
	network["dmd"][4] = NetworkForm("DM", 4, 0, 2);

	network["eha"][0] = NetworkForm("EH", 0, 0, 0);
	network["eha"][1] = NetworkForm("EH", 1, 0, 1);
	network["eha"][2] = NetworkForm("EH", 2, 0, 1);
	network["eha"][3] = NetworkForm("EH", 3, 0, 1);
	network["eha"][4] = NetworkForm("EH", 4, 0, 2);

	network["enl"][0] = NetworkForm("EN", 0, 0, 0);
	network["enl"][1] = NetworkForm("EN", 1, 0, 1);
	network["enl"][2] = NetworkForm("EN", 2, 0, 1);
	network["enl"][3] = NetworkForm("EN", 3, 0, 1);

	network["enu"][0] = NetworkForm("EN", 1, 1, 1);
	network["enu"][1] = NetworkForm("EN", 2, 1, 1);
	network["enu"][2] = NetworkForm("EN", 3, 1, 1);

	network["enr"][0] = NetworkForm("EN", 1, 2, 1);
	network["enr"][1] = NetworkForm("EN", 2, 2, 1);
	network["enr"][2] = NetworkForm("EN", 3, 2, 1);
	network["enr"][3] = NetworkForm("EN", 4, 2, 2);
}   

Dmera::NetworkForm::NetworkForm(std::string type, int idx1, int idx2, int idx3)
{
	class Nodes
	{
		public:
			void append(std::string name, int idx)
			{
				TensorList[name] = {up[idx], up[idx + 1], bond, bond + 1};
				up[idx]		= bond;
				up[idx + 1]	= bond + 1;
				bond += 2;
			}
			void appendT(std::string name, int idx)
			{
				TensorList[name] = {bond, bond + 1, down[idx], down[idx + 1]};
				down[idx]		= bond;
				down[idx + 1]	= bond + 1;
				bond += 2;
			}           
			void appendS(std::string name, int idx)
			{
				TensorList[name] = {up[idx], up[idx + 1], 0, 0};
				up[idx]		= 0;
				up[idx + 1]	= 0;
			}            
			void appendST(std::string name, int idx)
			{
				TensorList[name] = {0, 0, down[idx], down[idx + 1]};
				down[idx]		= 0;
				down[idx + 1]	= 0;
			}
			void appendF(std::string name, int idx)
			{
				int idx_next;
				switch (idx)
				{
					case (0): idx_next = 1; break;
					case (1): idx_next = 4; break;
					case (4): idx_next = 5; break;
				}
				TensorList[name] = {up[idx], up[idx_next], down[idx], down[idx_next]};
			}

		private:
			std::map<std::string, std::array<int, 4>> TensorList;

			std::array<int, 6> up	= {1, 2, 3, 4, 5, 6};
			std::array<int, 6> down	= {1, 2, 3, 4, 5, 6};
			int bond = 7;
	};

	Nodes nodes;

	if (type == "DM")	nodes.append("_", idx1);
	else				nodes.append("eh" + std::string(idx1), idx1);

	if ((type == "EN") && (idx2 == 1))	nodes.append("_", 2);			
	else								nodes.append("u", 2);

	if ((type == "EN") && (idx2 == 0))	nodes.append("_", 1);			
	else								nodes.append("l", 1);

	if ((type == "EN") && (idx2 == 2))	nodes.append("_", 3);			
	else								nodes.append("r", 3);

	nodes.appendT("ut", 2);
	nodes.appendT("lt", 1);
	nodes.appendT("rt", 3);

	nodes.appendS("s", 2);
	nodes.appendST("st", 2);

	if (type == "EH")
	{
		if (idx3 == 0)	nodes.appendF("_", 0);
		if (idx3 == 1)	nodes.appendF("_", 1);
		if (idx3 == 2)	nodes.appendF("_", 4);
	}
	else 
	{
		if (idx3 == 0)	nodes.appendF("dm" + std::string(idx3), 0);
		if (idx3 == 1)	nodes.appendF("dm" + std::string(idx3), 1);
		if (idx3 == 2)	nodes.appendF("dm" + std::string(idx3), 4);
	}

	std::map<std::string, std::array<int, 4>> tensorList = nodes.getList();

	std::string fname = PNAME + std::string(type) + std::string(idx1) + std::string(idx2) + std::string(idx3);
	std::ofstream of(fname);

	for (auto t : tensorList)
	{
		std::string name = t.first;
		std::array<int, 4> bonds = t.second;

		if (name == "_")
		{
			bonds_tout = bonds;
			continue;
		}

		if (bonds[0] == 0)
			of << name << " : ; " << bonds[2] << " " << bonds[3] << std::endl;
		else if (bonds[2] == 0)  
			of << name << " : " << bonds[0] << " " << bonds[1] << " ; " << std::endl;
		else
			of << name << " : " << bonds[0] << " " << bonds[1] << " ; " << bonds[2] << " " << bonds[3] << std::endl;

		symbols.push_back(name);
	}

	of << "TOUT : " << bonds_tout[2] << " " << bonds_tout[3] << " ; " << bonds_tout[0] << " " << bonds_tout[1] << std::endl;

	of.close();

	network = new uni10::Network(fname);
}

uni10::UniTensor Dmera::NetworkForm::launch(Block* b)
{
	for ( auto name : symbols )
	{
		if		(name == "lt")	network->putTensorD(name,	b->tensor("l"));
		else if	(name == "rt")	network->putTensorD(name,	b->tensor("r"));
		else if	(name == "ut")	network->putTensorD(name,	b->tensor("u"));
		else if	(name == "st")	network->putTensorD(name,	b->tensor("s"));
		else                    network->putTensor(name,	b->tensor(name));
	}

	return network->launch();
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
				b->dm[i] = dm[Dmera::index(i - 2, dm.size())]; 
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

		// Var. Update

		uni10::UniTensor t_l, t_u, t_r;

		t_l	= network["enl"][0].launch(b)
			+ network["enl"][1].launch(b) 
			+ network["enl"][2].launch(b) 
			+ network["enl"][3].launch(b);

		t_u	= network["enu"][0].launch(b)
			+ network["enu"][1].launch(b) 
			+ network["enu"][2].launch(b);

		t_r	= network["enr"][0].launch(b)
			+ network["enr"][1].launch(b) 
			+ network["enr"][2].launch(b) 
			+ network["enr"][3].launch(b);

		t_l = Dmera::eigenshift(t_l);
		t_u = Dmera::eigenshift(t_u);
		t_r = Dmera::eigenshift(t_r);

		b->update(	Dmera::svdSolveMinimal(t_l), 
				Dmera::svdSolveMinimal(t_u), 
				Dmera::svdSolveMinimal(t_r));


		// Ascend Hamiltonian

		int l, c, r;

		if (boundry)
		{
			l = eh.size() - 2;
			c = eh.size() - 1;
			r = 0;
		}
		else
		{
			l = Dmera::index(idx - 2, eh.size());
			c = Dmera::index(idx - 1, eh.size());
			r = Dmera::index(idx    , eh.size());
		}

		eh[l]	= network["eha"][0].launch(b);
		eh[c]	= network["eha"][1].launch(b) 
			+ network["eha"][2].launch(b)
			+ network["eha"][3].launch(b);
		eh[r]	= network["eha"][4].launch(b); 

	}

}

double Dmera::effective_j(double j, double j_l, double j_r, double delta) { return j_l * j_r / j / (1. + delta); }

int Dmera::index(int idx, int size)
{
	if (idx < 0)			return idx + size;
	else if (idx >= size)	return idx - size;
	else					return idx;
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

uni10::UniTensor Dmera::DmSinglet(int type)
{
	uni10::Network* N;
	if (type == 0)
		N = new uni10::Network(DMS_0_FNAME); 
	else if (type == 1)
		N = new uni10::Network(DMS_1_FNAME);

	N->putTensor("s", Dmera::Singlet());
	N->putTensorT("st", Dmera::Singlet());

	return N->launch();
}

uni10::UniTensor Dmera::eigenshift(uni10::UniTensor input)
{
	uni10::Matrix input_m = input.getRawElem();
	double data[4][4];

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			data[i][j] = input_m[i][j];

	gsl_matrix_view m	= gsl_matrix_view_array(data, 4, 4);
	gsl_vector *eval	= gsl_vector_alloc(4);
	gsl_matrix *evec	= gsl_matrix_alloc(4, 4);

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(4);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);

	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

	double eval_0 = gsl_vector_get(eval, 0);

	return input - eval_0 * Dmera::Identity();	
}

uni10::UniTensor Dmera::svdSolveMinimal(uni10::UniTensor input)
{
	std::vector<uni10::Matrix> M_svd;
	int labels[] = {-1, -2, -3, -4};
	int label_groups[] = {2, 2};

	std::vector<uni10::UniTensor> T_svd = input.hosvd(labels, label_groups, 2, M_svd);

	uni10::Network SVD(SVD_FNAME);
	UP.putTensor("1", T_svd[1]);
	UP.putTensorT("2", T_svd[0]);

	return (-1) * SVD.launch();
}

#include <array>
#include <vector>
#include <random>
#include <utility>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

#include <cmath>

#include <fstream>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <uni10.hpp>

#include "Dmera.h"

Dmera::Dmera(std::vector<double> js, double delta): width(js.size()) 
{                                  
	// DM & EH init

	dm_init.resize(2);
	dm_init[0] = Dmera::DmSinglet(0);
	dm_init[1] = Dmera::DmSinglet(1);

	eh_init.resize(width);
	for (int i = 0; i < width; ++i)
		eh_init[i] = TwoSiteHam(js[i], delta);

	es.resize(width);

	//iteratively contract nodes (in SDRG order)

	//std::cout << "SDRG order:";

	while (js.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < js.size(); ++i)
			if (js[i] > js[idx])
				idx = i;

		//std::cout << " " << idx;

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

	//std::cout << std::endl;

	// Build Network Forms
	BuildNetworkForms();

	// Build Full Network
	BuildFullNetwork();

	of = std::ofstream("_energy");
}

Dmera::~Dmera()
{
	for (auto b : blocks) delete b;
	delete FN;
}

Dmera::Block::Block(int idx_): idx(idx_)
{
	t["l"] = Dmera::Random_Unitary();
	t["r"] = Dmera::Random_Unitary();
	t["u"] = Dmera::Random_Unitary();
	t["s"] = Dmera::Singlet();
}

int Dmera::Block::get_idx() const { return idx; }

void Dmera::Block::update(std::string name, uni10::UniTensor new_t)
{
	t[name] = new_t;
}

uni10::UniTensor Dmera::Block::tensor(std::string name) const
{
	if		(name == "l")	return t.at("l");
	else if	(name == "r")	return t.at("r");
	else if	(name == "u")	return t.at("u");
	else if	(name == "s")	return t.at("s");

	else if	(name == "dm0")	return dm[0];
	else if	(name == "dm1") return dm[1];
	else if	(name == "dm2") return dm[2];

	else if	(name == "eh0") return eh[0];
	else if	(name == "eh1") return eh[1];
	else if	(name == "eh2") return eh[2];
	else if	(name == "eh3") return eh[3];
	else if	(name == "eh4") return eh[4];
	
	else if	(name == "ehs0") return ehs[0];
	else if	(name == "ehs1") return ehs[1];
	else if	(name == "ehs2") return ehs[2];
	else if	(name == "ehs3") return ehs[3];
	else if	(name == "ehs4") return ehs[4];
}

void Dmera::BuildNetworkForms()
{
	network["dmd"].resize(5);
	network["eha"].resize(5);
	network["enl"].resize(4);
	network["enu"].resize(3);
	network["enr"].resize(4);

	network["dmdi"].resize(5);
	network["ehai"].resize(5);
	network["enli"].resize(4);
	network["enui"].resize(3);
	network["enri"].resize(4);

	network["dmd"][0] = NetworkForm("DM", 0, 0, 0);
	network["dmd"][1] = NetworkForm("DM", 1, 0, 1);
	network["dmd"][2] = NetworkForm("DM", 2, 0, 1);
	network["dmd"][3] = NetworkForm("DM", 3, 0, 1);
	network["dmd"][4] = NetworkForm("DM", 4, 0, 2);

	network["dmdi"][0] = NetworkForm("DMI", 0, 0, 0);
	network["dmdi"][1] = network["dmd"][1];
	network["dmdi"][2] = network["dmd"][2];
	network["dmdi"][3] = network["dmd"][3];
	network["dmdi"][4] = network["dmdi"][0];

	network["eha"][0] = NetworkForm("EH", 0, 0, 0);
	network["eha"][1] = NetworkForm("EH", 1, 0, 1);
	network["eha"][2] = NetworkForm("EH", 2, 0, 1);
	network["eha"][3] = NetworkForm("EH", 3, 0, 1);
	network["eha"][4] = NetworkForm("EH", 4, 0, 2);

	network["ehai"][0] = network["eha"][0];
	network["ehai"][1] = network["eha"][1];
	network["ehai"][2] = network["eha"][2];
	network["ehai"][3] = network["eha"][3];
	network["ehai"][4] = network["eha"][4];

	network["enl"][0] = NetworkForm("EN", 0, 0, 0);
	network["enl"][1] = NetworkForm("EN", 1, 0, 1);
	network["enl"][2] = NetworkForm("EN", 2, 0, 1);
	network["enl"][3] = NetworkForm("EN", 3, 0, 1);

	network["enli"][0] = NetworkForm("ENI", 0, 0, 0);
	network["enli"][1] = network["enl"][1];
	network["enli"][2] = network["enl"][2];
	network["enli"][3] = network["enl"][3];

	network["enu"][0] = NetworkForm("EN", 1, 1, 1);
	network["enu"][1] = NetworkForm("EN", 2, 1, 1);
	network["enu"][2] = NetworkForm("EN", 3, 1, 1);

	network["enui"][0] = network["enu"][0];
	network["enui"][1] = network["enu"][1];
	network["enui"][2] = network["enu"][2];

	network["enr"][0] = NetworkForm("EN", 1, 2, 1);
	network["enr"][1] = NetworkForm("EN", 2, 2, 1);
	network["enr"][2] = NetworkForm("EN", 3, 2, 1);
	network["enr"][3] = NetworkForm("EN", 4, 2, 2);

	network["enri"][0] = network["enr"][0];
	network["enri"][1] = network["enr"][1];
	network["enri"][2] = network["enr"][2];
	network["enri"][3] = NetworkForm("ENI", 4, 2, 2);
}   

Dmera::NetworkForm::NetworkForm() { }

Dmera::NetworkForm::NetworkForm(std::string type, int idx1, int idx2, int idx3)
{
	std::string fname = PNAME + std::string(type) + std::to_string(idx1) + std::to_string(idx2) + std::to_string(idx3);

	network = new uni10::Network(fname);
    symbols = network->get_symbols();
}

uni10::UniTensor Dmera::NetworkForm::launch(Block* b) const
{
	// Insert tensors from 'b' into 'network' according to 'symbols'
	// TODO: Using putTensorT here, thus only allowing RTYPE

	for (auto name : symbols)
	{
		if		(name == "lt")	network->putTensorT(name,	b->tensor("l"));
		else if	(name == "rt")	network->putTensorT(name,	b->tensor("r"));
		else if	(name == "ut")	network->putTensorT(name,	b->tensor("u"));
		else if	(name == "st")	network->putTensorT(name,	b->tensor("s"));
		else                    network->putTensor (name,	b->tensor(name));
	}

	return network->launch();
}

void Dmera::BuildFullNetwork()
{
	FN = new Network(width);

	for (int i = 0; i < blocks.size(); ++i)
	{
		auto b = blocks[i];
	  	int idx = b->get_idx();
		FN->putTensor(idx		, blocks[i]->tensor("u"), "u");
		FN->putTensor(idx - 1	, blocks[i]->tensor("l"), "l");
		FN->putTensor(idx + 1 	, blocks[i]->tensor("r"), "r");
		FN->putTensor(idx	 	, blocks[i]->tensor("s"), "s");

		FN->coarse(idx);
	}		

	FN->putTensor(0, blocks[0]->tensor("s"), "s");
}

Dmera::Network::Network(int size_)
{
	size = size_;	
	maxDepth = 0;

	for (int i = 0; i < size_; ++i)
		nodes.push_back(node(i));

	port.resize(size_);
	porth.resize(size_);
}

Dmera::Network::~Network()
{
	for (auto t : tensors)
		delete(t);
}

Dmera::Network::node::node(int pos_)
{
	pos = pos_;
	T = 0;
	lr = "";
}

Dmera::Network::Tensor::Tensor(int lPos_, int rPos_, std::string type_, uni10::UniTensor data_)
{
	lPos = lPos_;
	rPos = rPos_;
	type = type_;
	data = data_;

	lParent = 0;
	rParent = 0;

	flag = false;
}

void Dmera::Network::Tensor::setLParent(Tensor* T) { lParent = T; }
void Dmera::Network::Tensor::setRParent(Tensor* T) { rParent = T; }

Dmera::Network::Tensor* Dmera::Network::Tensor::getLParent() { return lParent; }
Dmera::Network::Tensor* Dmera::Network::Tensor::getRParent() { return rParent; }

void Dmera::Network::putTensor(int idx, uni10::UniTensor t, std::string type)
{
	int idx1 = Dmera::index(idx    , nodes.size());
	int idx2 = Dmera::index(idx + 1, nodes.size());

	int pos1 = nodes[idx1].pos;
	int pos2 = nodes[idx2].pos;

	Tensor* NT = new Tensor(pos1, pos2, type, t);
	tensors.push_back(NT);

	if (nodes[idx1].T != 0)
	{
		if (nodes[idx1].lr == "l")
			nodes[idx1].T->setLParent(NT);
		else
		    nodes[idx1].T->setRParent(NT);
	} else {
     	port[pos1] = NT;
	}
		
	if (nodes[idx2].T != 0)
	{
		if (nodes[idx2].lr == "l")
			nodes[idx2].T->setLParent(NT);
		else
			nodes[idx2].T->setRParent(NT);
	} else {
     	port[pos2] = NT;
	}

	if (porth[pos1] == 0) porth[pos1] = NT;

	int d1, d2;

	if (nodes[idx1].T == 0) d1 = 0;
	else d1 = nodes[idx1].T->depth;

	if (nodes[idx2].T == 0) d2 = 0;
	else d2 = nodes[idx2].T->depth;

	NT->depth = (d1 > d2)? (d1 + 1): (d2 + 1);

	if (NT->depth > maxDepth) maxDepth = NT->depth;

	nodes[idx1].T = NT;
	nodes[idx1].lr = "l";
	nodes[idx2].T = NT;
	nodes[idx2].lr = "r";
}

void Dmera::Network::coarse(int idx)
{
	if (idx == nodes.size() - 1)
	{
		nodes.pop_back();
		nodes.erase(nodes.begin());

	} else {
     	nodes.erase(nodes.begin() + idx); 
     	nodes.erase(nodes.begin() + idx); 
	}
}

void Dmera::Network::resetFlag()
{
   for (auto t : tensors) t->unsetFlag();
}

void Dmera::Network::causalCone(int i)
{
	 Tensor* t = port[i];
	 t->flagSelfAndParents();
}

void Dmera::Network::causalConeH(int i)
{
	 Tensor* t = porth[i];
	 t->flagSelfAndParents();
}

bool Dmera::Network::Tensor::flaged() 
{	
	return flag;
}

void Dmera::Network::Tensor::unsetFlag()
{
	flag = false;
}

void Dmera::Network::Tensor::flagSelfAndParents()
{
	flag = true;
	if (type != "s")
	{
		lParent->flagSelfAndParents();
		rParent->flagSelfAndParents();
	}	
}

void Dmera::Network::printNetwork()
{
	std::string graph[size][maxDepth];
    for (int i = 0; i < size; ++i)
    	for (int j = 0; j < maxDepth; ++j)
			graph[i][j] = "   ";
		
    for (int i = 0; i < size; ++i)
		std::cout << " " << i << " ";
	std::cout << std::endl;

	std::string left, right, bar;
	for (auto t : tensors)
	{
		if (!t->flaged())
		{
			left	= " └─";
			right	= "─┘ ";
			bar		= "───";

		} else {
			left	= " ╚═";
			right	= "═╝ ";
			bar		= "═══";
		}

		int l = t->lPos;	
		int r = t->rPos;
		int d = t->depth - 1;

		if (l > r)
		{
			for (int i = 0; i < r; ++i)	    	graph[i][d] = bar;
			for (int i = l + 1; i < size; ++i)	graph[i][d] = bar;

		} else {
         	for (int i = l + 1; i < r; ++i)		graph[i][d] = bar; 
		}

		graph[l][d] = left;
		graph[r][d] = right;
	}

    for (int j = 0; j < maxDepth; ++j)
	{
    	for (int i = 0; i < size; ++i)
			std::cout << graph[i][j];

		std::cout << std::endl;
	}
}

Dmera::Network::UniNetworkAgent::UniNetworkAgent(int size_)
{
	for (int i = 0; i < size_; ++i)
	{
		upperLeg.push_back(i);
		lowerLeg.push_back(i);
		disjointed.push_back(false);
	}
	totalLegs = size_;
}

Dmera::Network::UniNetworkAgent::tensorInfo::tensorInfo(int index_, std::vector<int> inLegs_, std::vector<int> outLegs_, bool transpose_):
	index(index_),
	inLegs(inLegs_),
	outLegs(outLegs_),
	transpose(transpose_) {}


int Dmera::Network::UniNetworkAgent::newLeg()
{
	++totalLegs;
	return totalLegs - 1;
}

void Dmera::Network::UniNetworkAgent::setOperator(int pos, uni10::UniTensor t)
{
	if (disjointed[pos]) return;

	lowerLeg[pos] = newLeg();
	tensorInfos.push_back(tensorInfo(	tensorDatas.size(), 
										{lowerLeg[pos]}, 
										{upperLeg[pos]},
										false));
	tensorDatas.push_back(t);

    disjointed[pos] = true;
}

void Dmera::Network::UniNetworkAgent::putUnitaries(int pos1, int pos2, uni10::UniTensor t)
{
	int newUpperLeg1 = newLeg();
	int newUpperLeg2 = newLeg();
	tensorInfos.push_back(tensorInfo(   tensorDatas.size(),
										{upperLeg[pos1], upperLeg[pos2]},
										{newUpperLeg1, newUpperLeg2},
										false));
	upperLeg[pos1] = newUpperLeg1;
	upperLeg[pos2] = newUpperLeg2;
	tensorDatas.push_back(t);

	int newLowerLeg1 = newLeg();
	int newLowerLeg2 = newLeg();
	tensorInfos.push_back(tensorInfo(   tensorDatas.size(), 
										{newLowerLeg1, newLowerLeg2},
										{lowerLeg[pos1], lowerLeg[pos2]},
										true));
	lowerLeg[pos1] = newLowerLeg1;
	lowerLeg[pos2] = newLowerLeg2;
	tensorDatas.push_back(t);
}

void Dmera::Network::UniNetworkAgent::putSinglets(int pos1, int pos2)
{
	tensorInfos.push_back(tensorInfo(   tensorDatas.size(),
										{upperLeg[pos1], upperLeg[pos2]},
										{},
										false));
	tensorDatas.push_back(Dmera::Singlet());

	tensorInfos.push_back(tensorInfo(   tensorDatas.size(), 
										{},
										{lowerLeg[pos1], lowerLeg[pos2]},
										true));
	tensorDatas.push_back(Dmera::Singlet());
}

void Dmera::Network::UniNetworkAgent::disjoint(int pos)
{
	if (disjointed[pos]) return;

	lowerLeg[pos] = newLeg();
	legPair.push_back(std::pair<int, int>(upperLeg[pos], lowerLeg[pos]));

	disjointed[pos] = true;
}

uni10::UniTensor Dmera::Network::UniNetworkAgent::launch()
{
	std::ofstream of(TEMP_FNAME);

	for (auto ti : tensorInfos)
	{
		of << std::to_string(ti.index) << "\t:\t";
		for (auto ib : ti.inLegs) of << ib << "\t";
		of << ";\t";
		for (auto ob : ti.outLegs) of << ob << "\t";
		of << std::endl;
	}
	of << "TOUT\t:\t";
	for (auto lp : legPair) of << lp.first << "\t";
	of << ";\t";
	for (auto lp : legPair) of << lp.second << "\t";
	of << std::endl;

	of.close();
            
	uni10::Network network(TEMP_FNAME); 

	for (auto ti : tensorInfos)
		if (ti.transpose)
			network.putTensorT(std::to_string(ti.index), tensorDatas[ti.index]);
		else
			network.putTensor(std::to_string(ti.index), tensorDatas[ti.index]);

	return network.launch();
	
}

double Dmera::Network::Entropy(int pos1, int pos2)
{
	UniNetworkAgent una(size);
	resetFlag();
	causalConeH(pos1);
	causalConeH(pos2);

	for (auto t : tensors)
	{
		if (!t->flaged()) continue;

		if (pos2 > pos1)
		{
			if (((t->lPos) > pos1) && ((t->lPos) <= pos2)) una.disjoint(t->lPos);
			if (((t->rPos) > pos1) && ((t->rPos) <= pos2)) una.disjoint(t->rPos);
		} else {
			if (((t->lPos) > pos1) || ((t->lPos) <= pos2)) una.disjoint(t->lPos);
			if (((t->rPos) > pos1) || ((t->rPos) <= pos2)) una.disjoint(t->rPos);
		}
			

		if (t->type == "s")	una.putSinglets (t->lPos, t->rPos);
		else				una.putUnitaries(t->lPos, t->rPos, t->data); 
	}

	uni10::Matrix result = una.launch().getRawElem();
	uni10::Matrix result_diag = result.eigh().at(0);

	double E = 0.;

	int rows = result_diag.row();
//    std::cout << result;     

	for (int i = 0; i < rows; ++i)
	{   
		//uni10::Real e = result_diag.at(uni10::CTYPE, i * rows + i).real();
		uni10::Real e = result_diag.at(i, i);
		if (e > 1.e-9) E += - e * std::log(e) / log(2);
	}

//	std::cout << result;										 
	return E;
}

double Dmera::AverageEntropy(int length) const
{
	double E = 0.;
	for (int pos = 0; pos < width; ++pos)
	{
		int pos1 = pos;
		int pos2 = index(pos + length, width);
		E += FN->Entropy(pos1, pos2);
	}
	return E / double(width);
}

void Dmera::VarUpdateForEpochs()
{
	for (int i = 0; i < EPOCH; ++i)
	{
		//std::cout << "EPOCH\t" << i + 1 << " : ";
		VarUpdate();
	}
}

void Dmera::VarUpdate()
{
	// Build descending density matrix

	std::vector<uni10::UniTensor> dm = dm_init;

	for (int i = blocks.size() - 1; i >= 0; --i)
	{

		Block* b = blocks[i];
		const int idx = b->get_idx();
		const bool boundry = (idx == (dm.size() + 1));

		// Block obtains density matrices from global array
		if (boundry)
			for (int j = 0; j < 3; ++j)
				b->dm[j] = dm[Dmera::index(j - 2, dm.size())]; 
		else
			for (int j = 0; j < 3; ++j)
				b->dm[j] = dm[Dmera::index(idx + j - 2, dm.size())];  

		// Resize global den. mat. array
		if (boundry)
		{
			dm.insert(dm.begin(), 0);
			dm.push_back(0);
		}
		else
			dm.insert(dm.begin() + idx, 2, 0);

		// Descend den. mat.
		for (int j = 0; j < 5; ++j)
			if (i == blocks.size() - 1)
				dm[Dmera::index(idx + j - 2, dm.size())] = network["dmdi"][j].launch(b);
			else
				dm[Dmera::index(idx + j - 2, dm.size())] = network["dmd"][j].launch(b);

	}

	// Variational update & Build ascending effective hamiltonian

	std::vector<uni10::UniTensor> eh = eh_init;

	for (int i = 0; i < width; ++i)
		eh[i] = Dmera::eigenshift(eh[i], es[i]);
		

	for (int i = 0; i < blocks.size(); ++i)
	{
		Block* b = blocks[i];
		const int idx = b->get_idx();
		const bool boundry = (idx == (eh.size() - 1));

		// Block obtains effective hamiltonians from global array
		for (int j = 0; j < 5; ++j)
		{
			b->eh[j] = eh[Dmera::index(idx + j - 2, eh.size())];

			b->ehs[j] = b->eh[j];
		}


		// Resize global eff. ham. array
		if (boundry)
		{
			eh.erase(eh.begin());
			eh.pop_back();
		}
		else
			eh.erase(eh.begin() + idx, eh.begin() + idx + 2);

		// Var. update

		uni10::UniTensor t_l, t_u, t_r;
		double E; // energy recording

		for (int j = 0; j < VAR_TIME; ++j)
		{
			t_u	= network["enu"][0].launch(b)
				+ network["enu"][1].launch(b) 
				+ network["enu"][2].launch(b);

			b->update("u", Dmera::svdSolveMinimal(t_u, E));
		}
		for (int j = 0; j < VAR_TIME; ++j)
		{
			if (i == blocks.size() - 1)
				t_l	= network["enli"][0].launch(b)
					+ network["enli"][1].launch(b) 
					+ network["enli"][2].launch(b) 
					+ network["enli"][3].launch(b);
			else
				t_l	= network["enl"][0].launch(b)
					+ network["enl"][1].launch(b) 
					+ network["enl"][2].launch(b) 
					+ network["enl"][3].launch(b);

			b->update("l", Dmera::svdSolveMinimal(t_l, E));
		}
		for (int j = 0; j < VAR_TIME; ++j)
		{
			if (i == blocks.size() - 1)
				t_r	= network["enri"][0].launch(b)
					+ network["enri"][1].launch(b) 
					+ network["enri"][2].launch(b) 
					+ network["enri"][3].launch(b);
			else
				t_r	= network["enr"][0].launch(b)
					+ network["enr"][1].launch(b) 
					+ network["enr"][2].launch(b) 
					+ network["enr"][3].launch(b);

			b->update("r", Dmera::svdSolveMinimal(t_r, E));
		}
		
		// Ascend eff ham.

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

	double E = 0;

	uni10::Network MUL(MUL_FNAME);

	for (int i : {0, 1})
	{
   		MUL.putTensor("1", eh[i]);
		MUL.putTensor("2", dm_init[i]);

		uni10::UniTensor t = MUL.launch();
		E += t.trace().real();
	}
	for (double ess : es)
		E += ess;

	of << E << std::endl;
	//std::cout << E << std::endl;
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

	T.orthoRand(RCTYPE);

	return T;
}

uni10::UniTensor Dmera::Singlet()
{

	uni10::Complex m[] = 
	{	0.,				1. / sqrt(2.),
		-1. / sqrt(2.),	0., 			};

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

	uni10::Complex m[] = {	
		1.,	0., 0., 0.,
		0., 1., 0., 0.,
		0., 0., 1., 0.,
		0., 0., 0., 1. 
	};

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

uni10::UniTensor Dmera::DmSinglet(int type)
{
	uni10::Network* N;
	if (type == 0)
		N = new uni10::Network(DMS_0_FNAME); 
	else if (type == 1)
		N = new uni10::Network(DMS_1_FNAME);

	N->putTensor("s", Dmera::Singlet());
	N->putTensorT("st", Dmera::Singlet());

	uni10::UniTensor T = N->launch();
	delete N;

	return T;
}

uni10::UniTensor Dmera::TwoSiteHam(double J, double delta)
{
	uni10::Bond b_in(uni10::BD_IN, 2);
	uni10::Bond b_out(uni10::BD_OUT, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_in);
	B.push_back(b_out);
	B.push_back(b_out);

	uni10::UniTensor T(RCTYPE, B);

	// XXZ Ham.

	uni10::Real M[16] = 
	{	delta,	0.,		0.,		0.,
		0.,		-delta,	2.,		0.,
		0.,		2.,		-delta,	0.,
		0.,		0.,		0.,		delta	};

	T.setRawElem(M);

	return 0.25 * J * T;
}

uni10::UniTensor Dmera::eigenshift(uni10::UniTensor input, double& e)
{
	/*
	   uni10::Matrix input_m = input.getRawElem();
	   double data[16];

	   for (int i = 0; i < 4; ++i)
	   for (int j = 0; j < 4; ++j)
	   data[4 * i + j] = input_m[4 * i + j];

	   gsl_matrix_view m	= gsl_matrix_view_array(data, 4, 4);
	   gsl_vector *eval	= gsl_vector_alloc(4);
	   gsl_matrix *evec	= gsl_matrix_alloc(4, 4);

	   gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(4);
	   gsl_eigen_symmv(&m.matrix, eval, evec, w);
	   gsl_eigen_symmv_free(w);

	   gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

	   double eval_0 = gsl_vector_get(eval, 0);

	   return input + (-1) * eval_0 * Dmera::Identity();	
	 */

	std::vector<uni10::Matrix> m = input.getRawElem().eigh();
	double max_ev = m[0].max();

	return input + (-1) * max_ev * Dmera::Identity();
}

uni10::UniTensor Dmera::svdSolveMinimal(uni10::UniTensor input, double& E)
{
	std::vector<uni10::Matrix> M_svd;
	int labels[] = {-1, -2, -3, -4};
	int label_groups[] = {2, 2};

	input.setLabel(labels);
	std::vector<uni10::UniTensor> T_svd = input.hosvd(labels, label_groups, 2, M_svd);

	uni10::Matrix MS = T_svd[2].getRawElem();

	for (int i : {0, 5, 10, 15})
		if (MS(i).real() == 0.) MS(i) = 1.;
		else MS(i) = double((0. < MS(i).real()) - (MS(i).real() < 0.));
   
	uni10::UniTensor TS = T_svd[2];
	TS.setRawElem(MS);

	uni10::Network SVD(SVD_FNAME);
	SVD.putTensor("1", T_svd[1]);
	SVD.putTensor("I", TS);
	SVD.putTensorT("2", T_svd[0]);

	E = (-1) * (T_svd[2].getRawElem() * MS).trace(uni10::CTYPE).real();	
	//E = (-1) * T_svd[2].getRawElem().trace(uni10::RTYPE);	

	return (-1) * SVD.launch();
}

void Dmera::check() const
{
	FN->printNetwork();
}

void print(const uni10::UniTensor& t) { std::cout << t; }

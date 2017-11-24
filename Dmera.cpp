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

enum class Dmera::Tensor::Type
{
	Unitary,
	Singlet
};

void Dmera::Tensor::set_type(Dmera::Tensor::Type type_) { type = type_; }

Dmera::Tensor::Tensor(Bond* b1, Bond* b2): in1(b1), in2(b2), _idx1(b1->get_idx()), _idx2(b2->get_idx())
{
	type = Type::Unitary;
	out1 = 0;
	out2 = 0;

	int depth1, depth2;

	if (b1->get_in() == 0)
		depth1 = 0;
	else depth1 = b1->get_in()->get_depth() + 1;

	if (b2->get_in() == 0)
		depth2 = 0;
	else depth2 = b2->get_in()->get_depth() + 1;

	depth = (depth1 > depth2)? depth1: depth2;
}

void Dmera::Tensor::set_out(Bond* b1, Bond* b2)
{
	out1 = b1;
	out2 = b2;
}

bool Dmera::Tensor::open() const
{
	switch (type)
	{
		case (Type::Unitary): return ((in1 == 0) || (in2 == 0) || (out1 == 0) || (out2 == 0));
		case (Type::Singlet): return ((in1 == 0) || (in2 == 0));
	}
}

Dmera::Bond* Dmera::Tensor::get_in1() const { return in1; }

Dmera::Bond* Dmera::Tensor::get_in2() const { return in2; }

Dmera::Bond* Dmera::Tensor::get_out1() const { return out1; }

Dmera::Bond* Dmera::Tensor::get_out2() const { return out2; }

int Dmera::Tensor::get_depth() const { return depth; } 

int Dmera::Tensor::get_idx1() const { return _idx1; }

int Dmera::Tensor::get_idx2() const { return _idx2; }

void Dmera::Tensor::putTensor(uni10::UniTensor T_) { T = T_; }

uni10::UniTensor Dmera::Tensor::getTensor() const { return T; }

Dmera::Bond::Bond(Tensor* t, int idx): in(t), _idx(idx)
{
	out = 0;
}

void Dmera::Bond::set_in(Tensor* t)
{
	in = t;
}

void Dmera::Bond::set_out(Tensor* t)
{
	out = t;
}

bool Dmera::Bond::open() const
{
	return ((in == 0) || (out == 0));
}

Dmera::Tensor* Dmera::Bond::get_in() const { return in; }

Dmera::Tensor* Dmera::Bond::get_out() const { return out; }

int Dmera::Bond::get_idx() const { return _idx; }

Dmera::Sdrg_Node::Sdrg_Node(double j_, int idx_, Bond* b): _j(j_), _idx(idx_), _bond(b) {}

double Dmera::Sdrg_Node::j() const { return _j; }

Dmera::Bond* Dmera::Sdrg_Node::bond() const { return _bond; }

void Dmera::Sdrg_Node::set_j(double j_) { _j = j_; }

void Dmera::Sdrg_Node::set_bond(Bond* b) { _bond = b; }

Dmera::Dmera(std::vector<double> js, double delta): width(js.size())
{
	//build original sdrg-list

	int idx = 0;
	std::vector<Sdrg_Node*> nodes;
	for (double j : js)
	{
		Bond* b = new Bond(0, idx++);
		nodes.push_back(new Sdrg_Node(j, nodes.size(), b));
		//bonds.push_back(b);
		in_bonds.push_back(b);

		Bond* bo = new Bond(0, -1);
		op_bonds.push_back(bo);

		Tensor* to = new Tensor(bo, bo);
		operators.push_back(to);
	}

	for (int i = 0; i < width; ++i)
	{
		int r = (i == (width - 1))? 0: (i + 1);
		operators[i]->set_out(in_bonds[i], op_bonds[r]);
		in_bonds[i]->set_in(operators[i]);
		op_bonds[r]->set_in(operators[i]);
		op_bonds[r]->set_out(operators[r]);
	}

	//iteratively contract nodes

	while (nodes.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < nodes.size(); ++i)
			if (nodes[i]->j() > nodes[idx]->j())
				idx = i;

		Sdrg_Node* const node_1 = nodes[((idx - 1) < 0)? (idx + nodes.size() - 1): (idx - 1)];
		Sdrg_Node* const node_2 = nodes[idx];
		Sdrg_Node* const node_3 = nodes[((idx + 1) >= nodes.size())? (idx - nodes.size() + 1): (idx + 1)];
		Sdrg_Node* const node_4 = nodes[((idx + 2) >= nodes.size())? (idx - nodes.size() + 2): (idx + 2)];

		//construct new tensors and bonds

		Tensor* t_d		= Append_Tensor(node_2->bond(), node_3->bond());

		std::pair<Bond*, Bond*> bond_d = Append_Bonds(t_d);
		Bond* bd_d1		= bond_d.first;
		Bond* bd_d2		= bond_d.second;

		Tensor* t_u1	= Append_Tensor(node_1->bond(), bd_d1);
		Tensor* t_u2	= Append_Tensor(bd_d2, node_4->bond());

		std::pair<Bond*, Bond*> bond_u1 = Append_Bonds(t_u1);
		Bond* bd_u1		= bond_u1.first;
		Bond* bd_u2		= bond_u1.second;

		std::pair<Bond*, Bond*> bond_u2 = Append_Bonds(t_u2);
		Bond* bd_u3		= bond_u2.first;
		Bond* bd_u4		= bond_u2.second;

		Tensor* t_s		= Append_Tensor(bd_u2, bd_u3);
		t_s->set_type(Tensor::Type::Singlet);

		//node update

		node_1->set_bond(bd_u1);
		node_4->set_bond(bd_u4);

		node_1->set_j(Dmera::effective_j(node_2->j(), node_1->j(), node_3->j(), delta));

		delete node_2;
		delete node_3;

		if (idx == nodes.size() - 1)
		{
			nodes.erase(nodes.begin());
			nodes.erase(nodes.end() - 1);			
		} 
		else 
		{
			nodes.erase(nodes.begin() + idx, nodes.begin() + idx + 2);
		}
	}

	Tensor* t_s = Append_Tensor(nodes[0]->bond(), nodes[1]->bond());
	t_s->set_type(Tensor::Type::Singlet);


	//Tensors and bonds created.
    

	//Create Environment Network & Info
	for (Tensor* t : tensors)
		if (t->get_type() == Tensor::Type::Unitary)
			envs[t] = TensorEnv(t);
	
	//UniTensor initialize
	for (Tensor* t : tensors)
		switch (t->get_type())
		{
			case (Tensor::Type::Unitary): 
				t->putTensor(Dmera::Random_Unitary());
				break;
			case (Tensor::Type::Singlet):
				t->putTensor(Dmera::Singlet());
				break;
		}

	for (Tensor* t : operators)
		t->putTensor(Dmera::Random_Unitary());

	svd_restore = new uni10::Network("NetworkSheets/svd_restore");
}

Dmera::Tensor* Dmera::Append_Tensor(Bond* b1, Bond* b2)
{
	Tensor* t = new Tensor(b1, b2);
	b1->set_out(t);
	b2->set_out(t);

	tensors.push_back(t);
	return t;
}

std::pair<Dmera::Bond*, Dmera::Bond*> Dmera::Append_Bonds(Tensor* t)
{
	Bond* b1 = new Bond(t, t->get_idx1());
	Bond* b2 = new Bond(t, t->get_idx2());
	t->set_out(b1, b2);

	bonds.push_back(b1);
	bonds.push_back(b2);
	return std::pair<Bond*, Bond*>(b1, b2);
}

double Dmera::effective_j(double j, double j_l, double j_r, double delta)
{
	return j_l * j_r / j / (1. + delta);
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

	uni10::UniTensor T(uni10::CTYPE, B);

	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_int_distribution<int> unif(1, 50);

	for (int i = 0; i < unif(dre); ++i)
		T.orthoRand();

	return T;
}

uni10::UniTensor Dmera::Singlet()
{

	uni10::Complex m[] = {	0.,		1.,
							-1.,	0., };

	uni10::Bond b_in(uni10::BD_IN, 2);

	std::vector<uni10::Bond> B;
	B.push_back(b_in);
	B.push_back(b_in);

	uni10::UniTensor T(uni10::CTYPE, B);

    T.setRawElem(m);

	return T;
}

Dmera::Tensor::Type Dmera::Tensor::get_type() const { return type; }

std::string Dmera::Tensor::get_type_s() const
{
	switch (type)
	{
		case (Type::Unitary): return "Unitary";
		case (Type::Singlet): return "Singlet";
	}
}

std::string Dmera::Tensor::summary() const
{
	std::stringstream ss;

	ss << "Type: " << get_type_s() << " | Depth: " << get_depth() << " | " << get_idx1() << " ~ " << get_idx2() << std::endl;

	return ss.str();
}

std::string Dmera::summary(bool flagonly) const
{
	std::stringstream ss;

	ss << "System Width: " << width << std::endl;

	//Tensor summary

	for (Tensor* tensor : tensors)
	{
		ss << tensor->summary();
	}
                   
	ss << std::endl;

	//Diagram

	std::vector<std::vector<int>> diagram;
                 
	for (Tensor* tensor : tensors)
	{
		int d = tensor->get_depth();
		if (d + 1 > diagram.size())
			diagram.resize(d + 1);

		int x1 = tensor->get_idx1();
		int x2 = tensor->get_idx2();

		diagram[d].resize(width);

		int offset = tensor->flaged()? 5: (flagonly? -99: 0);

		diagram[d][x1] = 1 + offset;
		diagram[d][x2] = 3 + offset;
		if (x2 > x1)
		{
			for (int i = x1 + 1; i < x2; ++i)
				diagram[d][i] = 2 + offset;
		} 
		else
		{
			for (int i = x1 + 1; i < width; ++i)
				diagram[d][i] = 2 + offset;
			for (int i = 0; i < x2; ++i)
				diagram[d][i] = 2 + offset;
		}    	
	}

	for (int i = 0; i < width; ++i)
		ss << std::setw(3) << i;
	ss << std::endl;

	for (auto row : diagram)
	{
		ss << " ";
		for (int pixel : row)

			switch (pixel)
			{
				case (1):	ss << " └─"; break;
				case (2):	ss << "───"; break;
				case (3):	ss << "─┘ "; break;

				case (6):	ss << " ╘═"; break;
				case (7):	ss << "═══"; break;
				case (8):	ss << "═╛ "; break;

				default:	ss << "   "; break;
			}
		ss << std::endl;
	}              
	return ss.str();
}

// Flags

void Dmera::Tensor::flag() { _flag = true; }
void Dmera::Tensor::unflag() { _flag = false; }

bool Dmera::Tensor::flaged() const { return _flag; }

void Dmera::Tensor::flag_caus()
{
	_flag = true;

	if (type == Type::Unitary)
	{
		out1->flag_caus();
		out2->flag_caus();
	}
}

void Dmera::Bond::flag() { _flag = true; }
void Dmera::Bond::unflag() { _flag = false; }

bool Dmera::Bond::flaged() const { return _flag; }

void Dmera::Bond::flag_caus()
{
	_flag = true;

	out->flag_caus();
}

void Dmera::reset_flags()
{
	for (auto t : tensors) t->unflag();
	for (auto b : bonds ) b->unflag();
	for (auto b : in_bonds ) b->unflag();
}

Dmera::TensorInfo::TensorInfo(Tensor* t_, std::string name_, int in1_, int in2_, int out1_, int out2_, bool trans_):
	t(t_), name(name_), in1(in1_), in2(in2_), out1(out1_), out2(out2_), trans(trans_) {}

std::string Dmera::TensorInfo::bonds() const
{
	std::stringstream ss;
	ss << in1 << " " << in2 << " | " << out1 << " " << out2;

	return ss.str();
}

std::vector<Dmera::TensorInfo> Dmera::TensorEnv(Tensor* t0)
{
	reset_flags();

	t0->flag();

	std::vector<TensorInfo> info;
    std::map<Bond*, int> bond_symbol_up;
    std::map<Bond*, int> bond_symbol_down;

	int bond_idx = 0;

 	const int offset = bonds.size() + in_bonds.size();

	std::vector<Bond*> all_bonds = op_bonds;
	all_bonds.insert(all_bonds.end(), in_bonds.begin(), in_bonds.end());
	all_bonds.insert(all_bonds.end(), bonds.begin(), bonds.end());

	for (Bond* b : all_bonds)
	{        
		bond_symbol_up[b] = ++bond_idx;
		bond_symbol_down[b] = bond_symbol_up[b] + offset;
	}

	bond_symbol_up[t0->get_in1()]	= -3;
	bond_symbol_up[t0->get_in2()]	= -4;
	bond_symbol_up[t0->get_out1()]	= -1;
	bond_symbol_up[t0->get_out2()]	= -2;

	bond_symbol_up[0] = 0;
	bond_symbol_down[0] = 0;


	int tensor_idx = 0;
	for (Tensor* t : tensors)
	{
		if (!(t->flaged()))
		{
			info.push_back(TensorInfo(	t, 
										std::to_string(++tensor_idx), 
										bond_symbol_up[t->get_in1()], 
										bond_symbol_up[t->get_in2()], 
										bond_symbol_up[t->get_out1()], 
										bond_symbol_up[t->get_out2()],
										false ));
		}								

		info.push_back(TensorInfo(	t, 
									std::to_string(++tensor_idx), 
									bond_symbol_down[t->get_out1()], 
									bond_symbol_down[t->get_out2()], 
									bond_symbol_down[t->get_in1()], 
									bond_symbol_down[t->get_in2()],
									true ));

	}
		
	for (Tensor* t : operators)
		info.push_back(TensorInfo(	t, 
									std::to_string(++tensor_idx), 
									bond_symbol_up[t->get_in1()], 
									bond_symbol_down[t->get_out1()], 
									bond_symbol_up[t->get_out1()], 
									bond_symbol_up[t->get_out2()],
									false ));

//	for (auto i : info)
//		std::cout << i.t->summary() << " Bonds: " << i.bonds() << std::endl;
                                                       
  
	std::ofstream of("NetworkSheets/_temp");

	for (auto i : info)
		switch (i.t->get_type())
		{
			case(Tensor::Type::Unitary):
				of << i.name << ": " << i.in1 << " " << i.in2 << " ; " << i.out1 << " " << i.out2 << std::endl;
				break;

			case(Tensor::Type::Singlet):
				if (i.in1 == 0)
					of << i.name << ": " << " ; " << i.out1 << " " << i.out2 << std::endl;
				else
					of << i.name << ": " << i.in1 << " " << i.in2 << " ;" << std::endl;
				break;
		}

	of << "TOUT: -1 -2; -3 -4" << std::endl;
	of.close();		

	network[t0] = new uni10::Network("NetworkSheets/_temp");

	reset_flags();

	return info;
}
		
// Var-Update

void Dmera::VarUpdate()
{
	for (Tensor* t : tensors)
		if (t->get_type() == Tensor::Type::Unitary)
		{
			VarUpdateTensor(t);
			std::cout << "|" << std::flush;
		}
	
	std::cout << std::endl;
}

void Dmera::VarUpdateTensor(Tensor* t)
{
	for (auto info : envs[t])
	{
		if (!info.trans)
			network[t]->putTensor(info.name, info.t->getTensor());
		else
			network[t]->putTensorT(info.name, info.t->getTensor());
	}

	uni10::UniTensor eff_env = network[t]->launch();

	std::vector<uni10::Matrix> M_svd;
	int labels[] = {-1, -2, -3, -4};
	int label_groups[] = {2, 2};

	std::vector<uni10::UniTensor> eff_env_svd = eff_env.hosvd(labels, label_groups, 2, M_svd);

	svd_restore->putTensor("1", eff_env_svd[1]);
	svd_restore->putTensorT("2", eff_env_svd[0]);

	t->putTensor(svd_restore->launch() * -1);
}

void Dmera::check() const
{
}

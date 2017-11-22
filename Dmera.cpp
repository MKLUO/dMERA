#include <vector>
#include <utility>
//#include <iostream>
#include <sstream>

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

Dmera::Bond::Bond(Tensor* t, int idx): in(t), _idx(idx)
{
	out = 0;
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
	}

	//iteratively contract nodes

	while (nodes.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < nodes.size(); ++i)
			if (nodes[i]->j() > nodes[idx]->j())
				idx = i;

		//		for (int i = 0; i < nodes.size(); ++i) cout << nodes[i]->j() << " ";
		//		cout << endl << idx << " " << nodes.size() << endl;

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

std::string Dmera::Tensor::get_type() const
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

	ss << "Type: " << get_type() << " | Depth: " << get_depth() << " | " << get_idx1() << " ~ " << get_idx2() << std::endl;

	return ss.str();
}

std::string Dmera::summary() const
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

		diagram[d][x1] = 1;
		diagram[d][x2] = 3;
		if (x2 > x1)
		{
			for (int i = x1 + 1; i < x2; ++i)
				diagram[d][i] = 2;
		} 
		else
		{
			for (int i = x1 + 1; i < width; ++i)
				diagram[d][i] = 2;
			for (int i = 0; i < x2; ++i)
				diagram[d][i] = 2;
		}    	
	}

	for (int i = 0; i < width; ++i)
		ss << " " << i << " ";
	ss << std::endl;

	for (auto row : diagram)
	{
		for (int pixel : row)

			switch (pixel)
			{
				case (1):	ss << " ╘═"; break;
				case (2):	ss << "═══"; break;
				case (3):	ss << "═╛ "; break;
				default:	ss << "   "; break;
			}
		ss << std::endl;
	}              
	return ss.str();
}

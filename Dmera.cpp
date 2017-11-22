#include <vector>
#include <utility>
#include <iostream>

using std::cout;
using std::endl;

#include "Dmera.h"

enum class Dmera::Tensor::Type
{
	Unitary,
	Singlet
};

void Dmera::Tensor::set_type(Dmera::Tensor::Type type_) { type = type_; }

Dmera::Tensor::Tensor(Bond* b1, Bond* b2): in1(b1), in2(b2)
{
	type = Type::Unitary;
	out1 = 0;
	out2 = 0;
}

void Dmera::Tensor::set_out(Bond* b1, Bond* b2)
{
	out1 = b1;
	out2 = b2;
}

bool Dmera::Tensor::open() const
{
	return ((in1 == 0) || (in2 == 0) || (out1 == 0) || (out2 == 0));
}

Dmera::Bond::Bond(Tensor* t): in(t)
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

Dmera::Sdrg_Node::Sdrg_Node(double j_, int idx_, Bond* b): _j(j_), _idx(idx_), _bond(b) {}

double Dmera::Sdrg_Node::j() const { return _j; }

Dmera::Bond* Dmera::Sdrg_Node::bond() const { return _bond; }

void Dmera::Sdrg_Node::set_j(double j_) { _j = j_; }

void Dmera::Sdrg_Node::set_bond(Bond* b) { _bond = b; }

Dmera::Dmera(std::vector<double> js, double delta)
{
	//build original sdrg-list

	std::vector<Sdrg_Node*> nodes;
	for (double j : js)
		nodes.push_back(new Sdrg_Node(j, nodes.size(), new Bond(0)));

	//iteratively contract nodes

	while (nodes.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < nodes.size(); ++i)
			if (nodes[i]->j() > nodes[idx]->j())
				idx = i;

		for (int i = 0; i < nodes.size(); ++i) cout << nodes[i]->j() << " ";
		cout << endl << idx << " " << nodes.size() << endl;

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
	Bond* b1 = new Bond(t);
	Bond* b2 = new Bond(t);
	t->set_out(b1, b2);

	bonds.push_back(b1);
	bonds.push_back(b2);
	return std::pair<Bond*, Bond*>(b1, b2);
}

double Dmera::effective_j(double j, double j_l, double j_r, double delta)
{
	return j_l * j_r / j / (1. + delta);
}

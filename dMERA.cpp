
Dmera::Tensor::Tensor(Bond* b1, Bond* b2): in1(b1), in2(b2)
{
	type = Type::Unitary;
	out1 = 0;
	out2 = 0;
}

enum class Dmera::Tensor::Type
{
	Unitary,
	Singlet
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

Dmera::Dmera(Phys_Cond phys_cond)
{
	//build original sdrg-list

	//iteratively contract nodes

	while (nodes.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < nodes.size(); ++i)
			if (nodes[i]->j() > node[idx]->j())
				idx = i;
		
        Sdrg_Node* const node_1 = nodes[((idx - 1) < 0)? (idx + nodes.size()): idx];
        Sdrg_Node* const node_2 = nodes[idx];
        Sdrg_Node* const node_3 = nodes[((idx + 1) >= nodes.size())? (idx - nodes.size()): idx];
        Sdrg_Node* const node_4 = nodes[((idx + 2) >= nodes.size())? (idx - nodes.size()): idx];
                     
		//construct new tensors and bonds

		Tensor* t_d		= Append_Tensor(node_2->bond(), node_3->bond());

        std::pair<Bond*, Bond*> bonds = Append_Bonds(t_d);
		Bond* bd_d1		= bonds.first;
        Bond* bd_d2		= bonds.second;

        Tensor* t_u1	= Append_Tensor(node_1->bond(), bd_d1);
        Tensor* t_u2	= Append_Tensor(bd_d2, node_4->bond());

		std::pair<Bond*, Bond*> bonds = Append_Bonds(t_u1);
		Bond* bd_u1		= bonds.first;
		Bond* bd_u2		= bonds.second;

		std::pair<Bond*, Bond*> bonds = Append_Bonds(t_u2);
		Bond* bd_u3		= bonds.first;
		Bond* bd_u4		= bonds.second;

		Tensor* t_s		= Append_Tensor(bd_u2, bd_u3);
		t_s->set_type(Tensor::Type::Singlet);

		//node update

		node_1->set_bond(bd_u1);
		node_4->set_bond(bd_u4);

		node_1->set_next(node_4);
		node_4->set_prev(node_1);

		node_1->set_j(node_2->j(), node_1->j(), node_3->j());
	
		delete node_2;
		delete node_3;

		nodes.erase(nodes.begin() + idx, (idx + 1 == nodes.size())? nodes.begin(): (nodes.begin() + idx + 1));
	}

	Tensor* t_s = Append_Tensor(nodes[0]->bond(), nodes[1]->bond());
	t_s->set_type(Tensor::Type::Singlet);
}
                                 
Tensor* Dmera::Append_Tensor(Bond* b1, Bond* b2)
{
	Tensor* t = new Tensor(b1, b2);
	b1->set_out(t);
	b2->set_out(t);

	tensors.push_back(t);
	return t;
}

std::pair<Bond*, Bond*> Dmera::Append_Bonds(Tensor* t)
{
	Bond* b1 = new Bond(t);
	Bond* b2 = new Bond(t);
	t->set_out(b1, b2);

	bonds.push_back(b1);
	bonds.push_back(b2);
	return std::pair<Bond*, Bond*>(b1, b2);
}


Dmera::Dmera(Phys_Cond phys_cond)
{
	//build original sdrg-list

	//iteratively contract nodes

	while (nodes.size() >= 4)
	{
		int idx = 0;
		for (int i = 0; i < nodes.size(); ++i)
			if (nodes[i]->j > node[idx]->j)
				idx = i;
		
        Sdrg_Node* const node_1 = nodes[((idx - 1) < 0)? (idx + nodes.size()): idx];
        Sdrg_Node* const node_2 = nodes[idx];
        Sdrg_Node* const node_3 = nodes[((idx + 1) >= nodes.size())? (idx - nodes.size()): idx];
        Sdrg_Node* const node_4 = nodes[((idx + 2) >= nodes.size())? (idx - nodes.size()): idx];

		Tensor_Disentangle* t_d = new Tensor_Disentangle(node_2->bond(), node_3->bond());
		node_2->bond()->set_child(t_d);
		node_3->bond()->set_child(t_d);

		Bond* bd_d1 = new Bond(t_d);
        Bond* bd_d2 = new Bond(t_d);

        Tensor_Unitary* t_u1 = new Tensor_Unitary(node_1->bond(), bd_d1);
        Tensor_Unitary* t_u2 = new Tensor_Unitary(bd_d2, node_4->bond());

        Bond* bd_u1 = new Bond(t_u1);
        Bond* bd_u2 = new Bond(t_u1);
        Bond* bd_u3 = new Bond(t_u2);
        Bond* bd_u4 = new Bond(t_u2);

		node_1->set_bond(bd_u1);
		node_4->set_bond(bd_u4);

		Tensor_Singlet* t_s = new Tensor_Singlet(bd_u2, bd_u3);

		node_1->set_next(node_4);
		node_4->set_prev(node_1);

		delete node_2;
		delete node_3;

		nodes.erase(nodes.begin() + idx);
		nodes.erase(nodes.begin() + idx + 1);
	}
}

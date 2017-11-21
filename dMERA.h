
class Dmera
{
public:
	Dmera(Phys_Cond);

private:
	class Tensor;
	class Tensor_Disentangle	: Tensor;
	class Tensor_Unitary		: Tensor;
	class Tensor_Singlet		: Tensor;

	class Bond;

	class Sdrg_Node
	{
		double j;
		const int idx;
		const Bond* bond;
	}

	std::vector<Tensor*>	tensors;
	std::vector<Bond*>		bonds;
};


class Dmera
{
public:
	Dmera(Phys_Cond);

private:
	class Tensor
	{
	public:
		Tensor(Bond*, Bond*);
		enum class Type;
		void set_out(Bond*, Bond*);
		bool open() const;
	private:
		Type type;
		Bond* in1;
		Bond* in2;
		Bond* out1;
		Bond* out2;
	};

	class Bond;
	{
	public:
		Bond(Tensor*);
		void set_out(Tensor*);
		bool open() const;
	private:
		Tensor* in;
		Tensor* out;
	};

	class Sdrg_Node
	{
	public:
		Sdrg_Node(double, int, Bond*);
	private:
		double j;
		const int idx;
		const Bond* bond;
	};

	std::vector<Tensor*>	tensors;
	std::vector<Bond*>		bonds;
};

class 

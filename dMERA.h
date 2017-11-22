#ifndef DMERA_H
#define DMERA_H

#include <vector>
#include <utility>

class Dmera
{
	public:
		Dmera(std::vector<double>, double);

	private:
		class Tensor;
		class Bond;
		class Srdg_Node;

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
				double j() const;
				Bond* bond() const;

				void set_j(double);
				void set_bond(Bond*);

			private:
				double j;
				const int idx;
				Bond* bond;
		};

		Tensor* Append_Tensor(Bond*, Bond*);
		std::pair<Bond*, Bond*> Append_Bonds(Tensor*);

		static double effective_j(double, double, double, double);

		std::vector<Tensor*>	tensors;
		std::vector<Bond*>		bonds;
};

#endif //ifndef DMERA_H

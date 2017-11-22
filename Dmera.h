#ifndef DMERA_H
#define DMERA_H

#include <string>
#include <vector>
#include <utility>

class Dmera
{
	public:
		Dmera(std::vector<double>, double);

		std::string summary() const;

	private:
		class Tensor;
		class Bond;
		class Sdrg_Node;

		class Tensor
		{
			public:
				Tensor(Bond*, Bond*);
				enum class Type;
				void set_type(Type);
				void set_out(Bond*, Bond*);
				bool open() const;

				Bond* get_in1() const;
				Bond* get_in2() const;
				Bond* get_out1() const;
				Bond* get_out2() const;

				int get_depth() const;

				std::string get_type() const;
				std::string summary() const;

			private:
				Type type;
				int depth;
				Bond* in1;
				Bond* in2;
				Bond* out1;
				Bond* out2;
		};

		class Bond
		{
			public:
				Bond(Tensor*);
				void set_out(Tensor*);
				bool open() const;

				Tensor* get_in() const;
				Tensor* get_out() const;

			private:
				Tensor* in;
				Tensor* out;
		};

		class Sdrg_Node
		{
			public:
				Sdrg_Node(double, int, Bond*);
				double j() const;
				int idx() const;
				Bond* bond() const;

				void set_j(double);
				void set_bond(Bond*);

			private:
				double _j;
				const int _idx;
				Bond* _bond;
		};

		Tensor* Append_Tensor(Bond*, Bond*);
		std::pair<Bond*, Bond*> Append_Bonds(Tensor*);

		static double effective_j(double, double, double, double);

		std::vector<Tensor*>	tensors;
		std::vector<Bond*>		bonds;
		std::vector<Bond*>		in_bonds;
};

#endif //ifndef DMERA_H

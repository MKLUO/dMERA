#ifndef DMERA_H
#define DMERA_H

#include <string>
#include <vector>
#include <utility>

#include <uni10.hpp>

class Dmera
{
		class Tensor;
		class Bond;
		class Sdrg_Node;

	public:
		Dmera(std::vector<double>, double);

		std::string summary(bool) const;

		// Variational Updates

		void VarUpdate();

		void check() const;

	private:
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

				int get_idx1() const;
				int get_idx2() const;

				Type get_type() const;
				std::string get_type_s() const;
				std::string summary() const;

				void flag();
				void unflag();
				bool flaged() const;
				void flag_caus();

				void putTensor(uni10::UniTensor);
				uni10::UniTensor getTensor() const;

			private:
            	uni10::UniTensor T;

				Type type;
				int depth;
				Bond* in1;
				Bond* in2;
				Bond* out1;
				Bond* out2;

				const int _idx1, _idx2;

				bool _flag;
		};

		class Bond
		{
			public:
				Bond(Tensor*, int);
				void set_in(Tensor*);
				void set_out(Tensor*);
				bool open() const;

				Tensor* get_in() const;
				Tensor* get_out() const;

				int get_idx() const;

				void flag();
				void unflag();
				bool flaged() const;
				void flag_caus();

			private:
				Tensor* in;
				Tensor* out;

				const int _idx;

				bool _flag;
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
		
		struct TensorInfo
		{
			TensorInfo(Tensor*, std::string, int, int, int, int, bool);
			std::string bonds() const;			
			Tensor* t;
			std::string name;
			int in1, in2, out1, out2;
			bool trans;
		};

		Tensor* Append_Tensor(Bond*, Bond*);
		std::pair<Bond*, Bond*> Append_Bonds(Tensor*);

		static double effective_j(double, double, double, double);

		static uni10::UniTensor Random_Unitary();

		static uni10::UniTensor Singlet();

		uni10::Network* svd_restore;

		void reset_flags();

		std::vector<TensorInfo> TensorEnv(Tensor*);
		void VarUpdateTensor(Tensor*);

		const int width;

		std::vector<Tensor*>	tensors;
		std::vector<Tensor*>	operators;
		std::vector<Bond*>		bonds;
		std::vector<Bond*>		in_bonds;
		std::vector<Bond*>		op_bonds;

		std::map<Tensor*, std::vector<TensorInfo>> envs;
		std::map<Tensor*, uni10::Network*> network;
};

#endif //ifndef DMERA_H

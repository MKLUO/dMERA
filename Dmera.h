#ifndef DMERA_H
#define DMERA_H

#include <string>
#include <vector>
#include <utility>

#include <uni10.hpp>

const auto RCTYPE = uni10::RTYPE;

class Dmera
{
	public:
		Dmera(std::vector<double>, double);

		std::string summary(bool) const;

		double VarUpdate();

		void check() const;

	private:
    	class Block
		{
			public:
				int get_idx() const;

			private:
				const int idx;

				uni10::UniTensor t_l;
				uni10::UniTensor t_u;
				uni10::UniTensor t_r;

				std::array<uni10::UniTensor, 3> dm; // density matrix
				std::array<uni10::UniTensor, 5> eh; // effective hamiltonian
		};


		class NetworkForm
		{
			public:
				uni10::UniTensor launch(Block*) const;

			private:
				const int n;
				std::vector<std::string> symbols; //symbol used in block
				uni10::Network network;
		};


		// Utils
		static double effective_j(double, double, double, double);
		static int index(int, int);
		static uni10::UniTensor Random_Unitary();
		static uni10::UniTensor Singlet();
		static uni10::UniTensor Identity();
		static uni10::UniTensor Identity2();

		const int width;

		std::vector<Block*> blocks;
		std::map<NetworkForm> network;

		std::vector<uni10::UniTensor> dm_init;
		std::vector<uni10::UniTensor> eh_init;
};

#endif //ifndef DMERA_H

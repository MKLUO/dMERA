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

//		std::string summary(bool) const;

		double VarUpdate();

		void check() const;

	private:
    	class Block
		{
			public:
				Block(int);
				int get_idx() const;
				uni10::UniTensor tensor(std::string) const;

				void update(uni10::UniTensor, uni10::UniTensor, uni10::UniTensor);

				std::map<std::string, uni10::UniTensor> t;
				std::array<uni10::UniTensor, 3> dm; // density matrix
				std::array<uni10::UniTensor, 5> eh; // effective hamiltonian
			private:
				const int idx;
		};

		class NetworkForm
		{
			public:
				NetworkForm(std::string, int, int, int);
				uni10::UniTensor launch(Block*) const;

			private:
				std::vector<std::string> symbols; //symbol used in block
				uni10::Network* network;
		};

		void BuildMetworkForms();

		// Utils
		static double effective_j(double, double, double, double);
		static int index(int, int);
		static uni10::UniTensor Random_Unitary();
		static uni10::UniTensor Singlet();
		static uni10::UniTensor Identity();
		static uni10::UniTensor TwoSiteHam(double, double);

		static uni10::UniTensor eigenshift(uni10::UniTensor);
		static uni10::UniTensor svdSolveMinimal(uni10::UniTensor);
        
		const int width;

		std::vector<Block*> blocks;
		std::map<std::string, std::vector<NetworkForm>> network;

		std::array<uni10::UniTensor, 3> dm_init;
		std::vector<uni10::UniTensor> eh_init;
};

#endif //ifndef DMERA_H

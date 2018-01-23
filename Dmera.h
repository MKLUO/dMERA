#ifndef DMERA_H
#define DMERA_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <utility>

#include <uni10.hpp>

const auto RCTYPE = uni10::RTYPE;

const std::string SVD_FNAME("NetworkSheets/SVD");
const std::string MUL_FNAME("NetworkSheets/MUL");
const std::string DMS_0_FNAME("NetworkSheets/DMS_0");
const std::string DMS_1_FNAME("NetworkSheets/DMS_1");
const std::string PNAME("NetworkSheets/");

const bool FROM_FILE = true;

const int VAR_TIME = 5;

class Dmera
{
	public:
		Dmera(std::vector<double>, double);

//		std::string summary(bool) const;

		void VarUpdate();

		void check() const;

		static uni10::UniTensor eigenshift(uni10::UniTensor, double&);
		static uni10::UniTensor svdSolveMinimal(uni10::UniTensor, double&);

	private:
    	class Block
		{
			public:
				Block(int);
				int get_idx() const;
				uni10::UniTensor tensor(std::string) const;

				void update(std::string, uni10::UniTensor);

				std::map<std::string, uni10::UniTensor> t;
				std::array<uni10::UniTensor, 3> dm; // density matrix
				std::array<uni10::UniTensor, 5> eh; // effective hamiltonian
				std::array<uni10::UniTensor, 5> ehs; // effective hamiltonian - shift
			private:
				const int idx;
		};

		class NetworkForm
		{
			public:
				NetworkForm();
				NetworkForm(std::string, int, int, int);
				uni10::UniTensor launch(Block*) const;

			private:
				std::vector<std::string> symbols; //symbol used in block
				uni10::Network* network;
		};

		class Network
		{
			class Tensor
			{
				public:
					Tensor(int, int, int, std::string);
					Tensor* getLParent();
					Tensor* getRParent();
					Tensor* getLChild();
					Tensor* getRChild();

					void setLParent(Tensor*);
					void setRParent(Tensor*);
					void setLChild(Tensor*);
					void setRChild(Tensor*);

					int lPos, rPos, depth;

				private:
					Tensor* lParent;
					Tensor* rParent;
					Tensor* lChild;
					Tensor* rChild;

					int i; // i-th block in blocks
					std::string type;

					bool flag;
			};
			typedef struct node
			{
				node(int);
				int pos;
				Tensor* T;
				std::string lr;
			} node;

			public:
				Network(int);
				void putTensor(int, int, std::string);
				void coarse(int);

				void printNetwork(bool);

			private:
				int size, maxDepth;
				std::vector<node> nodes;
				std::vector<Tensor*> tensors;
		};

		void BuildNetworkForms();
		void BuildFullNetwork();

		// Utils
		static double effective_j(double, double, double, double);
		static int index(int, int);
		static uni10::UniTensor Random_Unitary();
		static uni10::UniTensor Singlet();
		static uni10::UniTensor Identity();
		static uni10::UniTensor DmSinglet(int);
		static uni10::UniTensor TwoSiteHam(double, double);
        
        // Datas

		const int width;

		std::vector<Block*> blocks;
		std::map<std::string, std::vector<NetworkForm>> network;

		Network* FN; // Full Network

		std::vector<uni10::UniTensor> dm_init;
		std::vector<uni10::UniTensor> eh_init;

		std::vector<double> es;

		std::ofstream of;
};

// for debug
void print(const uni10::UniTensor&);

#endif //ifndef DMERA_H

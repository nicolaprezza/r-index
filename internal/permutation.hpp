 /*
  * permutation: represents a permutation of the numbers in [1..r] using r * log r * (1 + 1/t), for some integer t > 0 specified at construction time.
  *
  * supports computing the permutation in O(1) time and the inverse permutation in O(t) time.
  *
  * from the paper:
  *
  * J. I. Munro, R. Raman, V. Raman, and S. S. Rao. Succinct representations of permutations and functions. Theor. Comput. Sci., 438:74–88, 2012.
  *
  */


#ifndef INTERNAL_PERMUTATION_HPP_
#define INTERNAL_PERMUTATION_HPP_

#include <definitions.hpp>

using namespace std;
using namespace sdsl;

#ifndef ulint
typedef uint64_t ulint;
#endif

#ifndef uint
typedef uint32_t uint;
#endif

namespace ri{

class permutation{

public:

	/*
	 * empty constructor. Initialize permutation with length 0.
	 */
	permutation(){}

	/*
	 * constructor. build permutation given in input the permutation P and the sampling
	 * factor t
	 */
	permutation(vector<ulint> Perm, ulint t = 6){

		ulint r = Perm.size();
		this->t = t;

		assert(r>0);

		uint w = r == 1 ? 1 :  (64-__builtin_clzll(r-1)); //log2 r

		sampled = bit_vector(r); //mark sampled positions
		vector<bool> visited(r,false); //has entry been visited during inversion?

		P = int_vector<>(r,0,w);

		for(ulint i=0;i<r;++i)
			P[i] = Perm[i];

		//scan permutation's entries
		for(ulint i=0;i<r;++i){

			//if not visited, navigate cycle
			if(not visited[i]){

				auto start = r;
				auto current = i;
				ulint k = 0;

				while(current != start){

					visited[current] = true;

					//sample
					if(k % t == 0){

						sampled[current] = true;

					}

					//in this case, we just started navigation
					//and current is the first element
					if(start==r) start = current;

					//next element
					current = P[current];
					k++;

				}

			}

		}

		rank1 = bit_vector::rank_1_type(&sampled);
		visited = vector<bool>(r,false);

		//init vector of samples
		S = int_vector<>(rank1(r),0,w);

		//scan permutation's entries and populate S
		for(ulint i=0;i<r;++i){

			//if not visited, navigate cycle
			if(not visited[i]){

				ulint k = 0;
				auto buf = vector<ulint>(t,r);

				auto start = r;
				auto current = i;

				while(current != start){

					visited[current] = true;

					//if sampled position and sample exists
					if(k % t == 0){

						assert(sampled[current]);

						S[rank1(current)] = buf[k%t];

					}

					buf[k%t] = current;

					//in this case, we just started navigation
					//and current is the first element
					if(start==r) start = current;

					//next element
					current = P[current];
					k++;

				}

				//t more cycles
				for(ulint j=0;j<t;++j){

					if(sampled[current]){

						S[rank1(current)] = buf[k%t];

					}

					buf[k%t] = current;

					k++;
					current = P[current];

				}

			}

		}

	}

	/*
	 * 	 returns image of i
	 */
	ulint map(ulint i){

		assert(i<size());
		return P[i];

	}
	/*
	 * 	 returns image of i
	 */
	ulint operator[](ulint i){

		assert(i<size());
		return P[i];

	}

	/*
	* returns: size of the bitvector
	*/
	ulint size(){return P.size();}

	/*
	 * return the inverse of i
	 */
	ulint inv(ulint i){

		auto orig = i;

		//have we gone backward yet?
		bool bw = false;

		for(ulint j=0;j<t;++j){

			bool smp = sampled[i];

			i = bw or (not smp) ? P[i] : S[rank1(i)];

			bw = bw or smp;

		}

		assert(map(i)==orig);

		return i;

	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint size=0;

		out.write((char*)&t, sizeof(t));
		size += sizeof(t);
		assert(size>0);

		size += sampled.serialize(out);
		assert(size>0);

		size += P.serialize(out);
		assert(size>0);

		size += S.serialize(out);
		assert(size>0);

		return size;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&t, sizeof(t));
		sampled.load(in);
		P.load(in);
		S.load(in);

		rank1 = bit_vector::rank_1_type(&sampled);

	}

private:

	ulint t = 0;

	bit_vector sampled; //sampled elements
	bit_vector::rank_1_type rank1;

	int_vector<> P; //the permutation
	int_vector<> S; //the samples


};

}


#endif /* INTERNAL_PERMUTATION_HPP_ */

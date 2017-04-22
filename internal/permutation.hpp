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

#include <vector>

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
	permutation(vector<ulint> Perm, ulint t = 8){

		ulint r = Perm.size();

		assert(r>0);

		uint w = r == 1 ? 1 :  (64-__builtin_clzll(r-1)); //log2 r

		bv = bit_vector(r);

		P = int_vector<>(r,0,w);

		for(ulint i=0;i<r;++i)
			P[i] = Perm[i];

		//TODO

		rank1 = bit_vector::rank_1_type(&bv);


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

		return 0;//TODO

	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint size=0;

		size += bv.serialize(out);
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

		bv.load(in);
		P.load(in);
		S.load(in);

		rank1 = bit_vector::rank_1_type(&bv);

	}

private:

	bit_vector bv; //sampled elements
	bit_vector::rank_1_type rank1;

	int_vector<> P; //the permutation
	int_vector<> S; //the samples


};

}


#endif /* INTERNAL_PERMUTATION_HPP_ */

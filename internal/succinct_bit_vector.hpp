 /*
  * succinct_bit_vector: an wrapper on biy_vector of the sdsl library, with support for rank and select
  */


#ifndef INTERNAL_SUCCINCT_BIT_VECTOR_HPP_
#define INTERNAL_SUCCINCT_BIT_VECTOR_HPP_

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

class succinct_bit_vector{

public:

	/*
	 * empty constructor. Initialize bitvector with length 0.
	 */
	succinct_bit_vector(){}

	/*
	 * constructor. build bitvector given a vector of bools
	 */
	succinct_bit_vector(vector<bool> b){

		bv = bit_vector(b.size());

		for(uint64_t i=0;i<b.size();++i)
			bv[i] = b[i];

		rank1 = bit_vector::rank_1_type(&bv);
		select1 = bit_vector::select_1_type(&bv);

	}

	succinct_bit_vector & operator= (const succinct_bit_vector & other) {

		bv = bit_vector(other.bv);
		rank1 = bit_vector::rank_1_type(&bv);
		select1 = bit_vector::select_1_type(&bv);

	    return *this;
	}

	/*
	 * not implemented
	 * argument: a boolean b
	 * behavior: append b at the end of the bitvector.
	 */
	//void push_back(bool b){}

	/*
	 * argument: position i in the bitvector
	 * returns: bit in position i
	 * only access! the bitvector is static.
	 */
	bool operator[](ulint i){

		assert(i<size());
		return bv[i];

	}

	/*
	 * argument: position i in the bitvector
	 * returns: number of bits equal to 1 before position i excluded
	 */
	ulint rank(ulint i){

		assert(i<=size());
		return rank1(i);

	}

	/*
	 * argument: integer i>=0
	 * returns: position of the i-th one in the bitvector. i starts from 0!
	 */
	ulint select(ulint i){

		assert(i<number_of_1());
		return select1(i+1);//in sd_vector, i starts from 1

	}

	/*
	* returns: size of the bitvector
	*/
	ulint size(){return bv.size();}

	/*
	 * returns: number of 1s in the bitvector
	 */
	ulint number_of_1(){return rank1(size()); }

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint size=0;

		size += bv.serialize(out);

		assert(size>0);

		return size;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		bv.load(in);
		rank1 = bit_vector::rank_1_type(&bv);
		select1 = bit_vector::select_1_type(&bv);

	}

private:

	bit_vector bv;
	bit_vector::rank_1_type rank1;
	bit_vector::select_1_type select1;

};

}


#endif /* INTERNAL_SUCCINCT_BIT_VECTOR_HPP_ */

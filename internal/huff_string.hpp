/*
 * huff_string.hpp
 *
 *  Created on: May 18, 2015
 *      Author: nicola
 *
 *  Huffman-compressed string with access/rank/select. The class is a wrapper on sdsl::wt_huff, with a simpler constructor
 */

#ifndef HUFF_STRING_HPP_
#define HUFF_STRING_HPP_

#include <sdsl/wavelet_trees.hpp>

using namespace sdsl;
using namespace std;

namespace ri{

class huff_string{

public:

	huff_string(){}

	huff_string(string &s){

		s.push_back(0);
		construct_im(wt, s.c_str(), 1);

		assert(wt.size()==s.size()-1);

	}

	uchar operator[](ulint i){

		assert(i<wt.size());
		return wt[i];

	}

	ulint size(){
		return wt.size();
	}

	ulint rank(ulint i, uchar c){

		assert(i<=wt.size());
		return wt.rank(i,c);

	}

	/*
	 * position of i-th character c. i starts from 0!
	 */
	ulint select(ulint i, uchar c){

		return wt.select(i+1,c);

	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		return wt.serialize(out);

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		wt.load(in);

	}

private:

	//wt_gmr<> wt;

	wt_huff<> wt;

};

}

#endif /* HUFF_STRING_HPP_ */

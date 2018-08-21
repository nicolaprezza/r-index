/*
 * rle_string.hpp
 *
 *  Created on: May 18, 2015
 *      Author: nicola
 *
 *  A run-length encoded string with rank/access functionalities.
 *
 *
 *  space of the structure: R * (H0 + log(n/R) + log(n/R)/B ) (1+o(1)) bits, n being text length,
 *  R number of runs, B block length, and H0 zero-order entropy of the run heads.
 *
 *  Time for all operations: O( B*(log(n/R)+H0) )
 *
 *  From the paper
 *
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza and Mathieu Raffinot.
 *  Flexible Indexing of Repetitive Collections. Computability in Europe (CiE) 2017)
 *
 */

#ifndef RLE_STRING_HPP_
#define RLE_STRING_HPP_

#include "definitions.hpp"
#include "huff_string.hpp"
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

namespace ri{

template<
	class sparse_bitvector_t = sparse_sd_vector, 	//predecessor structure storing run length
	class string_t	= huff_string 					//run heads
>
class rle_string{

public:

	rle_string(){}

	/*
	 * constructor: build structure on the input string
	 * \param input the input string without 0x0 bytes in it.
	 * \param B block size. The main sparse bitvector has R/B bits set (R being number of runs)
	 *
	 */
	rle_string(string &input, ulint B = 2){

		assert(not contains0(input));

		this->B = B;
		n = input.size();
		R = 0;

		auto runs_per_letter_bv = vector<vector<bool> >(256);

		//runs in main bitvector
		vector<bool> runs_bv;
		string run_heads_s;

		uchar last_c = input[0];

		for(ulint i=1;i<input.size();++i){

			if(uchar(input[i]) != last_c){

				run_heads_s.push_back(last_c);
				runs_per_letter_bv[last_c].push_back(true);

				last_c = input[i];

				//push back a bit set only at the end of a block
				runs_bv.push_back(R%B==B-1);

				R++;

			}else{

				runs_bv.push_back(false);
				runs_per_letter_bv[last_c].push_back(false);

			}

		}

		run_heads_s.push_back(last_c);
		runs_per_letter_bv[last_c].push_back(true);
		runs_bv.push_back(false);
		R++;

		assert(run_heads_s.size()==R);
		assert(R==count_runs(input));

		//cout << "runs in BWT(input) = " << count_runs(input) << endl;
		//cout << "runs in rle bwt = " << R << endl << endl;

		//now compact structures

		assert(runs_bv.size()==input.size());

		ulint t = 0;
		for(ulint i=0;i<256;++i)
			t += runs_per_letter_bv[i].size();

		assert(t==input.size());

		runs = sparse_bitvector_t(runs_bv);

		//a fast direct array: char -> bitvector.
		runs_per_letter = vector<sparse_bitvector_t>(256);
		for(ulint i=0;i<256;++i)
			runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);

		run_heads = string_t(run_heads_s);

		assert(run_heads.size()==R);


	}

	uchar operator[](ulint i){

		assert(i<n);
		return run_heads[run_of(i).first];

	}

	/*
	 * position of i-th character c. i starts from 0!
	 */
	ulint select(ulint i, uchar c){

		assert(i<runs_per_letter[c].size());

		//i-th c is inside j-th c-run (j starts from 0)
		assert(i<runs_per_letter[c].size());
		ulint j = runs_per_letter[c].rank(i);

		//starting position of i-th c inside its run
		assert(j==0 || i >= runs_per_letter[c].select(j-1) + 1);
		ulint before = (j==0 ? i : i - (runs_per_letter[c].select(j-1) + 1));

		//position in run_heads
		ulint r = run_heads.select(j,c);

		//k = number of bits before position of interest in the main string
		//here, k is initialized looking at the sampled runs
		assert(r/B==0 || r/B-1<runs.number_of_1());
		ulint k = (r/B==0?0 : runs.select(r/B-1)+1);

		//now add remaining run lengths to k
		for( ulint t = (r/B)*B; t<r; ++t ){

			k += run_at(t);

		}

		return k + before;

	}

	/*
	 * number of c before position i
	 */
	ulint rank(ulint i, uchar c){

		assert(i<=n);

		//letter does not exist in the text
		if(runs_per_letter[c].size()==0) return 0;

		if(i==n) return runs_per_letter[c].size();

		ulint last_block = runs.rank(i);
		ulint current_run = last_block*B;

		//current position in the string: the first of a block
		ulint pos = 0;
		if(last_block>0)
			pos = runs.select(last_block-1)+1;

		assert(pos <= i);

		ulint dist = i-pos;

		//otherwise, scan at most B runs
		while(pos < i){

			pos += run_at(current_run);
			current_run++;

			if(pos<=i) dist = i-pos;

		}

		if(pos>i) current_run--;

		//position i is inside run current_run
		assert(current_run<R);

		//number of c runs before the current run
		ulint rk = run_heads.rank(current_run,c);

		//number of c before i in the current run
		ulint tail = (run_heads[current_run]==c)*dist;

		//in this case, either there are no c before position i
		//or the current run is the first of the kind ccc...cc
		if(rk==0) return tail;

		return runs_per_letter[c].select(rk-1)+1+tail;

	}

	/*
	 * text position i is inside this run
	 */
	ulint run_of_position(ulint i){

		assert(i<n);

		ulint last_block = runs.rank(i);
		ulint current_run = last_block*B;

		//current position in the string: the first of a block
		ulint pos = 0;
		if(last_block>0)
			pos = runs.select(last_block-1)+1;

		assert(pos <= i);

		ulint dist = i-pos;

		//otherwise, scan at most B runs
		while(pos < i){

			pos += run_at(current_run);
			current_run++;

			if(pos<=i) dist = i-pos;

		}

		if(pos>i) current_run--;

		//position i is inside run current_run
		assert(current_run<R);

		return current_run;

	}
	//break range: given a range <l',r'> on the string and a character c, this function
	//breaks <l',r'> in maximal sub-ranges containing character c.
	//for simplicity and efficiency, we assume that characters at range extremities are both 'c'
	//thanks to the encoding (run-length), this function is quite efficient: O(|result|) ranks and selects
	vector<range_t> break_range(range_t rn,uchar c){

		auto l = rn.first;
		auto r = rn.second;

		assert(l<=r);
		assert(r<size());

		assert(operator[](l)==c);
		assert(operator[](r)==c);

		//retrieve runs that contain positions l and r
		auto run_l = run_of(l);
		auto run_r = run_of(r);

		//in this case rn contains only character c: do not break
		if(run_l.first==run_r.first) return {rn};

		vector<range_t> result;

		//first range: from l to the end of the run containing position l
		result.push_back({l,run_l.second});

		//rank of c's of interest in run_heads
		ulint rank_l = run_heads.rank(run_l.first,c);
		ulint rank_r = run_heads.rank(run_r.first,c);

		//now retrieve run bounds of all c-runs of interest
		for(ulint j = rank_l+1;j<rank_r;++j){

			result.push_back(run_range(run_heads.select(j,c)));

		}

		//now last (possibly incomplete) run

		auto range = run_range(run_heads.select(rank_r,c));
		result.push_back({range.first,r});

		return result;

	}

	ulint size(){return n;}

	/*
	 * return inclusive range of j-th run in the string
	 */
	pair<ulint,ulint> run_range(ulint j){

		assert(j<run_heads.size());

		ulint this_block = j/B;
		ulint current_run = this_block*B;
		ulint pos = (this_block==0?0:runs.select(this_block-1)+1);

		while(current_run < j){

 			pos += run_at(current_run);
			current_run++;

		}

		assert(current_run == j);

		return {pos,pos+run_at(j)-1};

	}

	//length of i-th run
	ulint run_at(ulint i){

		assert(i<R);
		uchar c = run_heads[i];

		return runs_per_letter[c].gapAt(run_heads.rank(i,c));

	}

	ulint number_of_runs(){return R;}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&R,sizeof(R));
		out.write((char*)&B,sizeof(B));

		w_bytes += sizeof(n) + sizeof(R) + sizeof(B);

		if(n==0) return w_bytes;

		w_bytes += runs.serialize(out);

		for(ulint i=0;i<256;++i)
			w_bytes += runs_per_letter[i].serialize(out);

		w_bytes += run_heads.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));
		in.read((char*)&R,sizeof(R));
		in.read((char*)&B,sizeof(B));

		if(n==0) return;

		runs.load(in);

		runs_per_letter = vector<sparse_bitvector_t>(256);

		for(ulint i=0;i<256;++i)
			runs_per_letter[i].load(in);

		run_heads.load(in);

	}

	string toString(){

		string s;

		for(ulint i=0;i<size();++i)
			s.push_back(operator[](i));

		return s;

	}


	ulint print_space(){

		ulint tot_bytes = 0;

		{
			std::ofstream out("/dev/null");
			auto bytesize = runs.serialize(out);

			tot_bytes += bytesize;

			cout << "main runs bitvector: " << bytesize << " Bytes" <<endl;

		}

		{
			std::ofstream out("/dev/null");

			ulint bytesize=0;

			for(auto r:runs_per_letter) bytesize += r.serialize(out);

			tot_bytes += bytesize;

			cout << "runs-per-letter bitvectors: " << bytesize << " Bytes" <<endl;

		}

		{
			std::ofstream out("/dev/null");

			ulint bytesize = run_heads.serialize(out);

			tot_bytes += bytesize;

			cout << "run heads: " << bytesize << " Bytes" <<endl;

		}

		return tot_bytes;

	}

	/*
	 * input: inclusive range rn, character c
	 *
	 * return the position j that is closest to rn.first,
	 * such that character in position j is c and that is
	 * adjacent to a position j' inside rn that contains a
	 * character != c
	 *
	 * rn must contain c and at least another character d!=c
	 *
	 */
	ulint closest_run_break(range_t rn, uchar c){

		/*
		 * case 1: range begins with a c-run: return last position of the run
		 */
		if(operator[](rn.first)==c){

			ulint i = run_of_position(rn.first);

			ulint j = run_range(i).second;

			//j must be inside rn, i.e. rn must not contain only c
			//j must not be last position of rn: this would imply
			//that rn contain only c
			assert(j<rn.second);

			return j;

		}else{

			//case 2: first c-run starts in the middle of the range

			//rank i of first c in the range
			ulint i = rank(rn.first,c);

			assert(i<rank(size(),c));

			//map from rank space to string position:
			//i now is the first position inside the range that contains c
			i = select(i,c);

			assert(operator[](i)==c);
			assert(i<=rn.second);

			return i;

		}

	}

private:

	static ulint count_runs(string &s){

		ulint runs=1;

		for(ulint i=1;i<s.size();++i){

			if(s[i]!=s[i-1]) runs++;

		}

		return runs;

	}

	//<j=run of position i, last position of j-th run>
	pair<ulint,ulint> run_of(ulint i){

		ulint last_block = runs.rank(i);
		ulint current_run = last_block*B;

		//current position in the string: the first of a block
		ulint pos = 0;
		if(last_block>0)
			pos = runs.select(last_block-1)+1;

		assert(pos <= i);

		while(pos < i){

 			pos += run_at(current_run);
			current_run++;

		}

		assert(pos >= i);

		if(pos>i){

			current_run--;

		}else{//pos==i

			pos += run_at(current_run);

		}

		assert(pos>0);
		assert(current_run<R);

		return {current_run,pos-1};

	}

	bool contains0(string &s){

		for(auto c : s)
			if(c==0) return true;

		return false;

	}

	//block size: bitvector 'runs' has R/B bits set (R being number of runs)
	ulint B=0;

	sparse_bitvector_t runs;

	//for each letter, its runs stored contiguously
	vector<sparse_bitvector_t> runs_per_letter;

	//store run heads in a compressed string supporting access/rank
	string_t run_heads;

	//text length and number of runs
	ulint n=0;
	ulint R=0;

};

typedef rle_string<sparse_sd_vector> rle_string_sd;
typedef rle_string<sparse_hyb_vector> rle_string_hyb;

}

#endif /* RLE_STRING_HPP_ */

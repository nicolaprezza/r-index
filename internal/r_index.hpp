/*
 * r_index.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Main class implementing the r-index
 *
 */

#ifndef R_INDEX_H_
#define R_INDEX_H_

#include <definitions.hpp>
#include <rle_string.hpp>
#include <sparse_sd_vector.hpp>
#include <permutation.hpp>

using namespace sdsl;

namespace ri{

class r_index{

public:

	using rle_string_t	= rle_string_sd;	//run-length encoded string
	using triple = std::tuple<range_t, ulint, ulint>;

	r_index(){}

	/*
	 * Build index
	 */
	r_index(string &input){

		assert(not contains_reserved_chars(input));

		cout << "Text length = " << input.size() << endl << endl;

		cout << "(1/3) Building BWT (libdivsufsort) ... " << flush;

		//build run-length encoded BWT
		{

			string bwt_s = build_bwt(input);

			cout << "done.\n(2/3) RLE encoding BWT ... " << flush;

			bwt = rle_string_t(bwt_s);

			//build F column
			F = vector<ulint>(256,0);
			for(uchar c : bwt_s)
				F[c]++;

			for(ulint i=255;i>0;--i)
				F[i] = F[i-1];

			F[0] = 0;

			for(ulint i=1;i<256;++i)
				F[i] += F[i-1];

			for(ulint i=0;i<bwt_s.size();++i)
				if(bwt_s[i]==TERMINATOR)
					terminator_position = i;

			assert(input.size()+1 == bwt.size());

		}

		cout << "done. " << endl<<endl;

		auto r = bwt.number_of_runs();

		cout << "Number of BWT equal-letter runs: r = " << r << endl;
		cout << "Rate n/r = " << double(bwt.size())/r << endl;
		cout << "log2(r) = " << log2(double(r)) << endl;
		cout << "log2(n/r) = " << log2(double(bwt.size())/r) << endl << endl;

		cout << "(3/3) Building U, D factorizations and UD, DR permutations ..." << flush;

		auto rank_up = bit_vector::rank_1_type(&sampled_up);
		auto rank_down = bit_vector::rank_1_type(&sampled_down);
		auto rank_down_bwt = bit_vector::rank_1_type(&sampled_down_bwt);

		{

			vector<ulint> UD_vec(r-1,r);
			vector<ulint> DR_vec(r-1,r);

			//positions in vectors sa_up_samples and sa_down_samples
			ulint u = 0;
			ulint d = 0;

			//scan BWT positions
			for(ulint k=0;k<bwt.size();++k){

				if(sampled_up_bwt[k]){

					assert(k>0);
					assert(bwt[k-1]!=bwt[k]);

					auto u_sample = rank_up(sa_up_samples[u++]);
					assert(d>0);
					auto d_sample = rank_down(sa_down_samples[d-1]);

					assert(u_sample<r-1);
					assert(d_sample<r-1);

					UD_vec[u_sample] = d_sample;

				}

				if(sampled_down_bwt[k]){

					assert(k<bwt.size()-1);
					assert(bwt[k]!=bwt[k+1]);

					//rank of run containing position k
					auto r_pos = rank_down_bwt(k);
					auto d_sample = rank_down(sa_down_samples[d++]);

					assert(d_sample<r-1);
					assert(r_pos<r-1);

					DR_vec[d_sample] = r_pos;

				}

			}

			//check that we filled all permutation's positions
			assert(not_contains(UD_vec,r));
			assert(not_contains(DR_vec,r));

			//build permutation structures
			UD = permutation(UD_vec);
			DR = permutation(DR_vec);

		}

		//build gap-encoded bitvectors
		U = sparse_sd_vector(sampled_up);
		D = sparse_sd_vector(sampled_down);

		assert(D.rank(D.size())==r-1);
		assert(U.rank(U.size())==r-1);


		cout << " done. " << endl<<endl;

	}


	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//inclusive range
		return {0,bwt_size()-1};

	}

	uchar operator[](ulint i ){
		return bwt[i];
	}

	/*
	 * \param r inclusive range of a string w
	 * \param c character
	 * \return inclusive range of cw
	 */
	range_t LF(range_t rn, uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		//number of c before the interval
		ulint c_before = bwt.rank(rn.first,c);

		//number of c inside the interval rn
		ulint c_inside = bwt.rank(rn.second+1,c) - c_before;

		//if there are no c in the interval, return empty range
		if(c_inside==0) return {1,0};

		ulint l = F[c] + c_before;

		return {l,l+c_inside-1};

	}

	//backward navigation of the BWT
	ulint LF(ulint  i){

		auto c = bwt[i];
		return F[c] + bwt.rank(i,c);

	}

	//forward navigation of the BWT
	ulint FL(ulint  i){

		//i-th character in first BWT column
		auto c = F_at(i);

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	//forward navigation of the BWT, where for efficiency we give c=F[i] as input
	ulint FL(ulint  i, uchar c){

		//i-th character in first BWT column
		assert(c == F_at(i));

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	/*
	 * access column F at position i
	 */
	uchar F_at(ulint i){

		ulint c = (upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
		assert(c<256);
		assert(i>=F[c]);

		return uchar(c);

	}

	/*
	 * Return BWT range of character c
	 */
	range_t get_char_range(uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		ulint l = F[c];
		ulint r = bwt_size()-1;

		if(c<255)
			r = F[c+1]-1;

		return {l,r};

	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){

		auto range = full_range();
		ulint m = P.size();

		for(ulint i=0;i<m and range.second>=range.first;++i)
			range = LF(range,P[m-i-1]);

		return range;

	}

	/*
	 * Retrieve range of input pattern, plus locate an occurrence <j, k> of the pattern, where j is text pos and k is corresponding bwt pos
	 */
	triple count_and_get_occ(string &P){

		auto range = full_range();
		ulint j = 0;
		ulint k = 0;

		return {range, j, k};

	}

	/*
	 * iterator locate(string &P){
	 *
	 * 	return iterator to iterate over all occurrences without storing them
	 * 	in memory
	 *
	 * }
	 */

	/*
	 * get number of runs in the BWT (terminator character included)
	 */
	ulint number_of_runs(){
		return bwt.number_of_runs();
	}

	/*
	 * get terminator (0x1) position in the BWT
	 */
	ulint get_terminator_position(){
		return terminator_position;
	}

	/*
	 * get BWT in string format
	 */
	string get_bwt(){
		return bwt.toString();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint w_bytes = 0;

		assert(F.size()>0);
		assert(bwt.size()>0);

		out.write((char*)&terminator_position,sizeof(terminator_position));
		out.write((char*)F.data(),256*sizeof(ulint));

		w_bytes += sizeof(terminator_position) + 256*sizeof(ulint);

		w_bytes += bwt.serialize(out);

		w_bytes += U.serialize(out);
		w_bytes += D.serialize(out);
		w_bytes += UD.serialize(out);
		w_bytes += DR.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&terminator_position,sizeof(terminator_position));

		F = vector<ulint>(256);
		in.read((char*)F.data(),256*sizeof(ulint));

		bwt.load(in);

		U.load(in);
		D.load(in);
		UD.load(in);
		DR.load(in);

	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".ri" will be automatically added
	 */
	void save_to_file(string path_prefix){

		string path = string(path_prefix).append(".ri");

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * load the structure from the path specified.
	 * \param path: full file name
	 */
	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}

	ulint text_size(){
		return bwt.size()-1;
	}

	ulint bwt_size(){
		return bwt.size();
	}

	uchar get_terminator(){
		return TERMINATOR;
	}

	string get_extension(){
		return string(".ri");
	}

	ulint print_space(){

		cout << "Number of runs = " << bwt.number_of_runs() << endl<<endl;

		ulint tot_bytes = bwt.print_space();

		cout << "\nTOT BWT space: " << tot_bytes << " Bytes" <<endl<<endl;

		return tot_bytes;

	}

private:

	/*
	 * temporary vectors used to speed-up construction
	 */

	//BWT positions associated to first position in each BWT run. In total r-1 entries
	//because we do not store up sample for the first BWT position
	vector<ulint> sa_up_samples;
	//BWT positions associated to first position in each BWT run
	//in total, r-1 entries: we do not store sample for last BWT position
	vector<ulint> sa_down_samples;

	//bitvectors of length n marking, respectively, first and last
	//BWT character in each run
	bit_vector sampled_up; //mark text positions associated with first position of a run (except first run)
	bit_vector sampled_down; //mark text positions associated with last position of a run (except last run)
	bit_vector sampled_up_bwt; //mark first position of each run on bwt (except first run)
	bit_vector sampled_down_bwt; //mark last position of each run on bwt (except last run)

	bool not_contains(vector<ulint> &V, ulint x){

		ulint r=0;

		for(auto y:V){

			if(y==x){

				cout << "failed at run " << r << " / " << V.size() << endl;
				return false;

			}

			r++;

		}

		return true;

	}

	uint8_t bitsize(uint64_t x){

		if(x==0) return 1;
		return 64 - __builtin_clzll(x);

	}

	/*
	 * check if range rn on column F contains a
	 * single character
	 */
	bool uniq_char(range_t rn){

		for(ulint i=0;i<256;++i){

			ulint l = F[i];
			ulint r = (i==255?bwt_size():F[i+1]);

			if(rn.first >= l and rn.second < r ) return true;

		}

		return false;

	}

	/*
	 * builds BWT of input string using SE_SAIS algorithm
	 * uses 0x1 character as terminator
	 *
	 */
	string build_bwt(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = LIBDIVSUFSORT;
	    construct_sa<8>(cc);

	    //now build BWT from SA
	    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cc));

	    {

	        for (ulint i=0; i<sa.size(); i++){
	            auto x = sa[i];

	            assert(x<=text.size());

	            if ( x > 0 )
	            	bwt_s.push_back((uchar)text[x-1]);
	            else
	            	bwt_s.push_back(TERMINATOR);

	        }

	    }

	    sampled_up = bit_vector(bwt_s.size());
	    sampled_down = bit_vector(bwt_s.size());
	    sampled_up_bwt = bit_vector(bwt_s.size());
	    sampled_down_bwt = bit_vector(bwt_s.size());

	    /*
	     * scan BWT and suffix array
	     */
	    for(ulint k=0;k<bwt_s.size();++k){

    		//text position associated with bwt position k
    		ulint j = sa[k]>0 ? sa[k]-1 : bwt_s.size()-1;

	    	//position k is the first in its run
	    	if(k>0 and bwt_s[k]!=bwt_s[k-1]){

	    		sampled_up_bwt[k] = true;
	    		sampled_up[j] = true;
	    		sa_up_samples.push_back( j );

	    	}

	    	//position k is the last in its run
			if(k<bwt_s.size()-1 and bwt_s[k]!=bwt_s[k+1]){

	    		sampled_down_bwt[k] = true;
	    		sampled_down[j] = true;

	    		sa_down_samples.push_back( j );

			}

	    }

	    sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
	    sdsl::remove(cache_file_name(conf::KEY_SA, cc));

	    return bwt_s;

	}

	static bool contains_reserved_chars(string &s){

		for(auto c : s)
			if(c == 0 or c == 1)
				return true;

		return false;

	}

	static const uchar TERMINATOR = 1;


	/*
	 * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
	 */

	//F column of the BWT (vector of 256 elements)
	vector<ulint> F;
	//L column of the BWT, run-length compressed
	rle_string_t bwt;
	ulint terminator_position = 0;

	/*
	 * Invertible permutations (i.e. efficient map and inverse):
	 *
	 * UD maps text positions associated to first (up) position of each BWT run to
	 * the last (down) position of the previous run. Size of the permutation is r-1 because
	 * the last run does not have runs after it.
	 *
	 * DR  maps text positions associated to Last (down) position of each BWT run to
	 * the corresponding BWT run (i.e. last position of that run)
	 *
	 * total: 2r log r * (1+epsilon) bits
	 *
	 */

	permutation UD; //Up samples to Down samples.
	permutation DR; //Down samples to BWT runs (last position of each run). Last run is excluded!

	/*
	 * gap-encoded bitvectors marking with a bit set sampled positions on the text
	 *
	 * U : marks text positions that are the first in their BWT run, except for the first run.
	 * D : marks text positions that are the last in their BWT run, except for the last run.
	 *
	 * r-1 bits set in each bitvector. Overall size = 2r log(n/r) bits
	 *
	 */

	sparse_sd_vector U;
	sparse_sd_vector D;

	/*
	 * overall: UD, DR, U, and D take r log n bits (plus low-order terms)
	 */

};

}

#endif /* R_INDEX_H_ */

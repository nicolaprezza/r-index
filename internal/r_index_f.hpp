/*
 * r_index_f.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Fast version of the r-index: O(r log (n/r)) words of space, O(1) locate time per occurrence
 *
 */

#ifndef R_INDEX_F_H_
#define R_INDEX_F_H_

#include <definitions.hpp>
#include <rle_string.hpp>
#include <sparse_sd_vector.hpp>

using namespace sdsl;

namespace ri{

class r_index_f{

public:

	using rle_string_t	= rle_string_sd;	//run-length encoded string
	using triple = std::tuple<range_t, ulint, ulint>;

	r_index_f(){}

	/*
	 * Build index of string input with sampling k and specified BWT algorithm
	 *
	 * - we store k samples before and after each BWT run. If k=0, k is set
	 *   as default to k = log(n/r)
	 *
	 * - if sais = true, SE-SAIS is used to build BWT. Otherwise, DIVSUFSORT is used
	 *
	 */
	r_index_f(string &input, ulint k = 0, bool sais = true){

		this->sais = sais;
		this->k = k;

		if(contains_reserved_chars(input)){

			cout << "Error: input string contains one of the reserved characters 0x0, 0x1" << endl;
			exit(1);

		}

		cout << "Text length = " << input.size() << endl << endl;

		cout << "(1/3) Building BWT ";
		if(sais) cout << "(SE-SAIS) ... " << flush;
		else cout << "(DIVSUFSORT) ... " << flush;

		//build run-length encoded BWT
		{

			string bwt_s = build_structures(input);

			k = this->k;

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



		cout << "Number of BWT equal-letter runs: r = " << r << endl;
		cout << "Rate n/r = " << double(bwt.size())/r << endl;
		cout << "log2(r) = " << log2(double(r)) << endl;
		cout << "log2(n/r) = " << log2(double(bwt.size())/r) << endl;
		cout << "Storing " << k << " SA samples before and after each BWT run break" << endl << endl;



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
	 * Return number of occurrences of P in the text
	 */
	ulint occ(string &P){

		auto rn = count(P);

		return (rn.second-rn.first)+1;

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
	 * locate all occurrences of P and return them in an array
	 * (space consuming if result is big).
	 */
	vector<ulint> locate_all(string& P){

		vector<ulint> OCC;

		//TODO
		return OCC;

	}


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

		//TODO

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

		//TODO

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

	ulint print_space(){

		cout << "Number of runs = " << bwt.number_of_runs() << endl<<endl;

		ulint tot_bytes = bwt.print_space();

		cout << "\nTOT BWT space: " << tot_bytes << " Bytes" <<endl<<endl;

		return tot_bytes;

	}

private:

	class pred_compare{

	public:

		pred_compare(int_vector<> * SSA){

			this->SSA=SSA;

		}

		bool operator()(ulint i, ulint j){

			assert(i<SSA->size());
			assert(j<SSA->size());

			return (*SSA)[i] < (*SSA)[j];

		}

	private:

		int_vector<> * SSA;

	};

	/*
	 * Retrieve range of input pattern, plus locate an occurrence <j, k> of the pattern, where j is text pos and k is corresponding bwt F-position, i.e. SA[k] = j
	 *
	 * returns <range, j,k>
	 *
	 */
	triple count_and_get_occ(string &P){

		//TODO

	}

	/*
	 * input: inclusive BWT range rn and character c
	 *
	 * find a BWT position k inside rn such that bwt[k]=c and k is sampled
	 * (i.e. k is the first or last position of its run and k is not first or last
	 * BWT position). Return k
	 *
	 */
	ulint get_sampled_position(range_t rn, uchar c){

		//check that rn contains at least 2 characters
		assert(not is_unary(rn));

		//check that BWT range contains character c
		assert(LF(rn,c) != range_t(1,0));

		/*
		 * Heuristic: retrieving SA[k-1] from SA[k] requires using permutation UD ( O(1/epsilon) time),
		 * while retrieving SA[k+1] from SA[k] requires using permutation DU ( O(1) time). It follows that is more
		 * efficient to extract SA entries from top to bottom of the range => better to locate top entries in the
		 * BWT range (so we minimize number of bottom-up steps). We therefore consider these 2 cases (rather than
		 * the symmetric ones):
		 *
		 * 1. bwt[rn] = ccc..cccx.... return position of last c in the first c-run
		 * 2. bwt[rn] = xx...xxcccc...cccc...yy... return position of first c in first c-run
		 */

		ulint k = bwt.closest_run_break(rn,c);

		return k;

	}

	/*
	 * input: inclusive BWT range [L..R]
	 *
	 * output: true iff BWT[L..R] is a unary string
	 *
	 */
	bool is_unary(range_t rn){

		return bwt.run_of_position(rn.first) == bwt.run_of_position(rn.second);

	}

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
	 * builds data structures and returns the BWT (in string format)
	 *
	 */
	string build_structures(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = sais ? SE_SAIS : LIBDIVSUFSORT;
	    construct_sa<8>(cc);

	    //now build BWT from SA
	    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cc));

	    r = 1;

	    //build BWT from SA and count runs
		for (ulint i=0; i<sa.size(); i++){
			auto x = sa[i];

			assert(x<=text.size());

			if ( x > 0 )
				bwt_s.push_back((uchar)text[x-1]);
			else
				bwt_s.push_back(TERMINATOR);

			if(i>0 && bwt_s[i]!=bwt_s[i-1]) r++;

		}

		n = bwt_s.size();

		assert(r>0);
		assert(n>0);

		//if sampling factor = 0, set it to log2(n/r)
		if(k==0){

			k = bitsize(n/r);

		}

	    //mark up ad down sampled BWT positions
	    auto sampled_position_up = vector<bool>(sa.size());
	    auto sampled_position_down = vector<bool>(sa.size());

	    //fill sampled_position_down
		for (ulint i=0; i<bwt_s.size()-1; i++){

			//if position i immediately precedes a run break
			if(bwt_s[i]!=bwt_s[i+1]){

				//sample k positions backward (including i). Stop as soon as an already marked position is found
				ulint j = 0;
				while(i>=j && not sampled_position_down[i-j] && j<k){

					sampled_position_down[i-j] = true;
					j++;

				}

			}

		}

	    //fill sampled_position_up
		for (ulint i=bwt_s.size()-1; i>0; i--){

			//if position i immediately follows a run break
			if(bwt_s[i]!=bwt_s[i-1]){

				//sample k positions forward (including i). Stop as soon as an already marked position is found
				ulint j = 0;
				while(i+j<n && not sampled_position_up[i+j] && j<k){

					sampled_position_up[i+j] = true;
					j++;

				}

			}

		}

		//count number of samples of all kinds
		ulint number_of_U_samples = 0;
		ulint number_of_D_samples = 0;
		ulint number_of_samples = 0; //not the same as number_of_U_samples + number_of_D_samples because there could be samples in common

		for(ulint i=0; i<bwt_s.size(); i++) number_of_U_samples += sampled_position_up[i];
		for(ulint i=0; i<bwt_s.size(); i++) number_of_D_samples += sampled_position_down[i];
		for(ulint i=0; i<bwt_s.size(); i++) number_of_samples += (sampled_position_up[i] or sampled_position_down[i]);

		logn = bitsize(n);

		SSA = int_vector<>(number_of_samples,0,logn);

		{

			//mark SSA samples that immediately precede a run break
			auto border_sample_vec = bit_vector(number_of_samples);
			ulint ssa_pos = 0;

			//populate SSA and border_sample_vec
			for (ulint i=0; i<bwt_s.size(); i++){

				//if position sampled
				if(sampled_position_up[i] or sampled_position_down[i]){

					assert(ssa_pos < SSA.size());

					//careful! SSA stores text positions associated to BWT positions
					SSA[ssa_pos] = sa[i] > 0 ? sa[i]-1 : n-1;

					//if position i immediately precedes a run break
					if(i<bwt_s.size()-1 && bwt_s[i]!=bwt_s[i+1]){

						assert(ssa_pos < border_sample_vec.size());
						border_sample_vec[ssa_pos] = true;

					}

					ssa_pos++;

				}

			}

			border_sample = bit_vector(border_sample_vec);
			border_sample_select = bit_vector::select_1_type(&border_sample);

		}

	    sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
	    sdsl::remove(cache_file_name(conf::KEY_SA, cc));

		ulint ssa_pos = 0;
		ulint U_pos = 0;
		ulint D_pos = 0;

		ulint log_nsamples = bitsize(ulint(SSA.size()));

		U = int_vector<>(number_of_U_samples,0,log_nsamples);
		D = int_vector<>(number_of_D_samples,0,log_nsamples);

		for (ulint i=0; i<bwt_s.size(); i++){

			if(sampled_position_up[i]){

				assert(ssa_pos < SSA.size());
				assert(U_pos<U.size());
				U[U_pos++] = ssa_pos;

			}

			if(sampled_position_down[i]){

				assert(ssa_pos < SSA.size());
				assert(D_pos<D.size());
				D[D_pos++] = ssa_pos;

			}

			if(sampled_position_up[i] or sampled_position_down[i]) ssa_pos++;

		}

		//sort U and D in text order
		std::sort(D.begin(), D.end(), pred_compare(&SSA));
		std::sort(U.begin(), U.end(), pred_compare(&SSA));

		//build suport for Elias-Fano-like predecessor queries

		logU = bitsize(ulint(U.size()));
		logD = bitsize(ulint(D.size()));

		assert(logU<=logn);
		assert(logD<=logn);

		if(logU==logn) logU--;
		if(logD==logn) logD--;

		assert(logU>0);
		assert(logD>0);

		buckets_U = bit_vector((ulint(1)<<logU));
		buckets_D = bit_vector((ulint(1)<<logD));

		partition_U = bit_vector(U.size());
		partition_D = bit_vector(D.size());

		partition_U[0] = true;
		partition_U[1] = true;

		ulint j = 0;//current position in U
		ulint last_pref = 0;

		for(auto u : U){

			assert(u<SSA.size());
			ulint x = SSA[u];//text position

			//x-prefix of logU bits
			ulint pref = x >> (logn-logU);

			assert(pref<buckets_U.size());

			buckets_U[pref] = true;

			if(j>0 and pref != last_pref){

				partition_U[j] = true;

			}

			last_pref = pref;
			j++;

		}

		j = 0;//current position in D

		for(auto d : D){

			assert(d<SSA.size());
			ulint x = SSA[d];//text position

			//x-prefix of logD bits
			ulint pref = x >> (logn-logD);

			assert(pref<buckets_D.size());

			buckets_D[pref] = true;

			if(j>0 and pref != last_pref){

				partition_D[j] = true;

			}

			last_pref = pref;
			j++;

		}

		return bwt_s;

	}

	static bool contains_reserved_chars(string &s){

		for(auto c : s)
			if(c == 0 or c == 1)
				return true;

		return false;

	}

	static const uchar TERMINATOR = 1;

	bool sais = true;

	/*
	 * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
	 */

	//F column of the BWT (vector of 256 elements)
	vector<ulint> F;
	//L column of the BWT, run-length compressed
	rle_string_t bwt;
	ulint terminator_position = 0;

	ulint k = 0;//sampling factor
	ulint r = 0;//number of BWT runs
	ulint n = 0;//BWT size
	ulint logn = 0;//number of bits required to store a text position

	//The compressed suffix array sampling. These positions correspond to BWT positions (not F-column positions)
	int_vector<> SSA;

	//pointers to SSA positions in text order
	int_vector<> U; //Up-samples
	int_vector<> D; //Down-samples

	//Elias-Fano support for predecessor on U (log n/|U| time)
	//buckets_U[x] = true iff there are SSA elements
	//pointed by U-elements starting with integer x (of bit size log|U|)
	//total space: O(|U|) bits
	ulint logU = 0;
	bit_vector buckets_U;
	bit_vector::rank_1_type buckets_U_rank;
	bit_vector partition_U;
	bit_vector::select_1_type partition_U_select;

	//Elias-Fano support for predecessor on D (log n/|D| time)
	ulint logD = 0;
	bit_vector buckets_D;
	bit_vector::rank_1_type buckets_D_rank;
	bit_vector partition_D;
	bit_vector::select_1_type partition_D_select;

	//mark with bit 1 samples in SSA that immediately precede a BWT run break
	bit_vector border_sample;
	bit_vector::select_1_type border_sample_select;

};

}

#endif /* R_INDEX_F_H_ */

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
	r_index(string &input, bool sais = true){

		this->sais = sais;

		assert(not contains_reserved_chars(input));

		cout << "Text length = " << input.size() << endl << endl;

		cout << "(1/3) Building BWT ";
		if(sais) cout << "(SE-SAIS) ... " << flush;
		else cout << "(DIVSUFSORT) ... " << flush;

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

		cout << "(3/3) Building U, D factorizations and DU, RD permutations ..." << flush;

		auto rank_up = bit_vector::rank_1_type(&sampled_up);
		auto rank_down = bit_vector::rank_1_type(&sampled_down);
		auto rank_down_bwt = bit_vector::rank_1_type(&sampled_down_bwt);

		{

			vector<ulint> DU_vec(r-1,r);
			vector<ulint> RD_vec(r-1,r);

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

					DU_vec[d_sample] = u_sample;

				}

				if(sampled_down_bwt[k]){

					assert(k<bwt.size()-1);
					assert(bwt[k]!=bwt[k+1]);

					//rank of run containing position k
					auto r_pos = rank_down_bwt(k);
					auto d_sample = rank_down(sa_down_samples[d++]);

					assert(d_sample<r-1);
					assert(r_pos<r-1);

					RD_vec[r_pos] = d_sample;

				}

			}

			//check that we filled all permutation's positions
			assert(not_contains(DU_vec,r));
			assert(not_contains(RD_vec,r));

			//build permutation structures
			DU = permutation(DU_vec);
			RD = permutation(RD_vec);

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

	/*
	 * input: position j on the text corresponding to some
	 * position k on the F-column of the BWT (i.e. SA[k]=j)
	 *
	 * output: SA[k+1]
	 */
	ulint next_SA(ulint j){

		auto rn = D.rank(j);

		auto pred_d_rank = rn == 0 ? D.rank(D.size())-1 : rn-1;//predecessor in rank space
		auto pred_d = D.select(pred_d_rank);//text position

		auto delta = pred_d == bwt.size()-1 ? j+1 : j-pred_d;

		auto pred_u = U.select(DU[pred_d_rank]);

		return (pred_u + delta) % bwt.size();

	}

	/*
	 * input: position j on the text corresponding to some
	 * position k on the F-column of the BWT (i.e. SA[k]=j)
	 *
	 * output: SA[k-1]
	 */
	ulint prev_SA(ulint j){

		//in this case, k = 0: no previous SA sample
		assert(j != bwt.size()-1);

		auto rn = U.rank(j);

		auto pred_u_rank = rn == 0 ? U.rank(U.size())-1 : rn-1;//predecessor in rank space
		auto pred_u = U.select(pred_u_rank);//text position

		auto delta = pred_u == bwt.size()-1 ? j+1 : j-pred_u;

		auto pred_d = D.select(DU.inv(pred_u_rank));

		return (pred_d + delta) % bwt.size();

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

		triple res = count_and_get_occ(P);

		range_t rn = std::get<0>(res);
		ulint k = std::get<2>(res);//SA position
		ulint j = std::get<1>(res);//text position: SA[k]

		assert(k >= rn.first and k <= rn.second);

		OCC.push_back(j);

		ulint before = k - rn.first;//occurrences to extract before k (k excluded)
		ulint after = rn.second - k;//occurrences to extract after k (k excluded)

		ulint j1 = j;

		for(ulint i = 0; i<before;++i){

			j1 = prev_SA(j1);
			OCC.push_back(j1);

		}

		j1 = j;

		for(ulint i = 0; i<after;++i){

			j1 = next_SA(j1);
			OCC.push_back(j1);

		}

		assert(OCC.size() == (rn.second-rn.first)+1);

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

		w_bytes += U.serialize(out);
		w_bytes += D.serialize(out);
		w_bytes += DU.serialize(out);
		w_bytes += RD.serialize(out);

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
		DU.load(in);
		RD.load(in);

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

	/*
	 * Retrieve range of input pattern, plus locate an occurrence <j, k> of the pattern, where j is text pos and k is corresponding bwt F-position, i.e. SA[k] = j
	 *
	 * returns <range, j,k>
	 *
	 */
	triple count_and_get_occ(string &P){

		//coordinates of an occurrence of current pattern suffix
		ulint j = 0;
		ulint k = 0;

		range_t range = full_range();
		range_t range1;

		ulint m = P.size();

		for(ulint i=0;i<m and range.second>=range.first;++i){

			uchar c = P[m-i-1];

			range1 = LF(range,c);

			//if suffix can be left-extended with char
			if(range1 != range_t(1,0)){

				if(is_unary(range)){

					assert(i>0);

					assert(j>0);

					//j-1 is an occurrence of the new pattern suffix
					j = j-1;

					//k is the new SA position corresponding to j
					k = LF(k);

				}else{

					//find a brand new occurrence
					k = get_sampled_position(range, c); //careful: k is an L-position
					j = k==0 ? 0 : bwt_sample(k);

					//map from L to F column
					k = LF(k);

				}

			}

			range = range1;

		}

		return triple(range, j, k);

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

	/*
	 * input: position k on the BWT
	 *
	 * k must be a sampled position, i.e. first or last in its run AND
	 * not the first or last BWT positions
	 *
	 * return: text position j associated with BWT position k
	 *
	 */
	ulint bwt_sample(ulint k){

		assert(k>0);
		assert(k<bwt.size()-1);
		assert(bwt[k] != bwt[k-1] or bwt[k] != bwt[k+1]);

		//is k the last position of its run? (recall that we sample down positions)
		bool down = bwt[k] != bwt[k+1];

		//run of BWT position k or k-1 (depending on who is a down position)
		ulint run = bwt.run_of_position(down ? k : k-1);

		//locate down position
		ulint rd_run = RD[run];

		//if k is a down position, map from D-rank space to text. Otherwise,
		//map from D-rank space to U-rank space (i.e. move down in the BWT)
		//and then pass from U-rank space to text position
		ulint sample = down ? D.select(rd_run) : U.select(DU[rd_run]);

		return sample;

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
	 * builds BWT of input string
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

	    construct_config::byte_algo_sa = sais ? SE_SAIS : LIBDIVSUFSORT;
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

	bool sais = true;

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
	 * DU maps (D-ranks of) text positions associated to last (down) position of each BWT run to
	 * the (U-ranks of) first (up) position of the next run. Size of the permutations is r-1 because
	 * the first/last run does not have runs before/after it.
	 *
	 * RD  maps BWT runs to text (D-ranks of) text positions associated to the Last (down) position
	 * of the corresponding BWT run (i.e. D-rank of text position associated to last position of that run)
	 *
	 * total: 2r log r * (1+epsilon) bits
	 *
	 */

	permutation DU; //Down samples to Up samples (all in rank space on D and U)
	permutation RD; //BWT runs to down samples on text in rank space (i.e. bitvector D). Last run is excluded!

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
	 * overall: DU, RD, U, and D take r log n bits (plus low-order terms)
	 */

};

}

#endif /* R_INDEX_H_ */

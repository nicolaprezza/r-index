#ifndef SRLBWT_H_
#define SRLBWT_H_

#include <definitions.hpp>
#include <sparse_sd_vector.hpp>
#include <rle_string.hpp>

using namespace sdsl;

namespace ri {

class s_rlbwt{

public:

	using rle_string_t	= rle_string_sd;	//run-length encoded string

	s_rlbwt(){}

	/*
	 * Build a run-length BWT from the input string
	 * input must not contain 0x0 and 0x1 characters
	 */
	s_rlbwt(string &input, ulint sa_rate = 128){

		assert(not contains_reserved_chars(input));
		assert(sa_rate>0);

		sr = sa_rate;

		//build run-length encoded BWT
		{

			string bwt_s = build_bwt(input);
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

		}

		ulint log2n = bitsize(uint64_t(bwt.size()));

		assert(log2n>0);

		//mark sampled BWT positions

		auto n = bwt.size();
		ulint sa_samples = n/sr + (n%sr != 0);

		assert(sa_samples>0);

		{

			vector<bool> sampled_BWT_pos_bools(bwt.size(),false);

			//position corresponding to last character
			ulint bwt_pos = terminator_position;
			ulint text_pos = bwt.size()-1;

			//now backward-navigate BWT

			while(text_pos>0){

				if(text_pos % sr == 0)
					sampled_BWT_pos_bools[bwt_pos] = true;

				bwt_pos = LF(bwt_pos);
				text_pos--;

			}

			assert(text_pos==0);

			//here text_pos=0
			sampled_BWT_pos_bools[bwt_pos] = true;

			//convert to gap-encoded vector
			sampled_BWT_pos = sparse_sd_vector(sampled_BWT_pos_bools);

			assert(sampled_BWT_pos.rank(sampled_BWT_pos.size())==sa_samples);

		}

		//build and populate SA sampling
		{

			SA = int_vector<>(sa_samples,0,log2n);

			//position corresponding to last character
			ulint bwt_pos = terminator_position;
			ulint text_pos = bwt.size()-1;

			//now backward-navigate BWT

			while(text_pos>0){

				if(text_pos % sr == 0){

					SA[ sampled_BWT_pos.rank(bwt_pos) ] = text_pos;

				}

				bwt_pos = LF(bwt_pos);
				text_pos--;

			}

			//here text_pos=0
			SA[ sampled_BWT_pos.rank(bwt_pos) ] = text_pos;

		}


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
	 * \param rn inclusive F-range containing a single character (say, a)
	 * \return 	disjoint inclusive ranges on column L whose union are all a's contained
	 * 			in the input range
	 *
	 * 	complexity: |result| * FL_cost
	 * 				where result is the output vector and FL_cost
	 * 				is the cost of FL mapping (log n/R with elias-fano RLE string)
	 *
	 */
	vector<range_t> FL(range_t rn){

		//F-range rn must contain only one character
		assert(uniq_char(rn));

		//unique character in F-range rn
		uchar c = F_at(rn.first);
		auto l = rn.first;
		auto r = rn.second;

		//break range: given a range <l',r'> on BWT and a character c, this function
		//breaks <l',r'> in maximal sub-ranges containing character c.
		return bwt.break_range({FL(l,c),FL(r,c)},c);

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
	 * Return BWT range of pattern P. If reverse = true, search P reversed
	 */
	range_t count(string &P, bool reverse = false){

		auto range = full_range();
		ulint m = P.size();

		for(ulint i=0;i<m and range.second>=range.first;++i)
			range = LF(range,(reverse?P[i]:P[m-i-1]) );

		return range;

	}

	/*
	 * locate F position i, i.e. return i-th suffix array entry
	 */
	ulint locate(ulint F_pos){

		return locate_L(FL(F_pos));

	}

	/*
	 * locate all occurrences in the inclusive BWT range
	 */
	vector<ulint> locate(range_t bwt_range){

		vector<ulint> occ;

		for(ulint f_pos = bwt_range.first; f_pos <= bwt_range.second; ++f_pos)
			occ.push_back(locate(f_pos));

		return occ;

	}

	/*
	 * locate all occurrences of pattern P
	 */
	vector<ulint> locate(string &P){

		return locate(count(P));

	}

	/*
	 * locate L position i (i.e. i is on the BWT)
	 */
	ulint locate_L(ulint L_pos){

		if(sampled_BWT_pos[L_pos]) return SA[sampled_BWT_pos.rank(L_pos)];

		return 1 + locate_L(LF(L_pos));

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


		out.write((char*)&sr,sizeof(ulint));
		w_bytes += sizeof(sr);

		w_bytes += SA.serialize(out);

		//marks sampled BWT positions
		w_bytes += sampled_BWT_pos.serialize(out);

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

		in.read((char*)&sr,sizeof(ulint));

		SA.load(in);

		sampled_BWT_pos.load(in);

	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".rlbwt" will be automatically added
	 */
	void save_to_file(string path_prefix){

		string path = string(path_prefix).append(".srlbwt");

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * load the structure from the path specified.
	 * \param path_prefix prefix of the index files. suffix ".rlbwt" will be automatically added
	 */
	void load_from_file(string path_prefix){

		string path = string(path_prefix).append(".srlbwt");
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
		return string(".rlbwt");
	}

	ulint print_space(){

		cout << "Number of runs = " << bwt.number_of_runs() << endl<<endl;

		ulint tot_bytes = bwt.print_space();

		cout << "\nTOT BWT space: " << tot_bytes << " Bytes" <<endl<<endl;

		ulint sa_space = 0;

		{
			std::ofstream out("/dev/null");
			auto bytesize = SA.serialize(out);

			sa_space += bytesize;
			tot_bytes += bytesize;

			cout << "SA sampling: " << bytesize << " Bytes" <<endl;

		}

		{
			std::ofstream out("/dev/null");
			auto bytesize = sampled_BWT_pos.serialize(out);

			sa_space += bytesize;
			tot_bytes += bytesize;

			cout << "marked sampled positions: " << bytesize << " Bytes" <<endl;

		}

		cout << "\nTOT SA sampling space: " << sa_space << " Bytes" <<endl<<endl;


		return tot_bytes;

	}

private:

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
	static string build_bwt(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = SE_SAIS;
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

	rle_string_t bwt;

	ulint terminator_position = 0;

	//F column of the BWT (vector of 256 elements)
	vector<ulint> F;

	//sample rate
	ulint sr = 0;

	//suffix array sampling
	int_vector<> SA;

	//marks sampled BWT positions
	sparse_sd_vector sampled_BWT_pos;

};

}

#endif /* SRLBWT_H_ */

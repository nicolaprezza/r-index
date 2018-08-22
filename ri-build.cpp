#include <iostream>

#include "internal/r_index.hpp"
#include "utils.hpp"
#include "internal/r_index.hpp"

using namespace ri;
using namespace std;

string out_basename=string();
string input_file=string();
int sa_rate = 512;
bool sais=true;
ulint T = 0;//Build fast index with SA rate = T
bool fast = false;//build fast index
bool hyb = false; //use hybrid bitvectors instead of sd_vectors?

void help(){
	cout << "ri-build: builds the r-index. Extension .ri is automatically added to output index file" << endl << endl;
	cout << "Usage: ri-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>        use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
	//cout << "   -h                   use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. Important: "<<endl;
	//cout << "                        if the index is built with -h, need to specify -h also when locating and counting (ri-count/ri-locte). "<<endl;
	//cout << "   -fast                build fast index (O(occ)-time locate, O(r log(n/r)) words of space). By default, "<<endl;
	//cout << "                        small index is built (O(occ*log(n/r))-time locate, O(r) words of space)"<<endl;
	//cout << "   -sa_rate <T>         T>0. if used, build the fast index (see option -fast) storing T SA samples before and after each"<<endl;
	//cout << "                        BWT equal-letter run. O(r*T) words of space, O(occ(log(n/r)/T) + log(n/r))-time locate. "<<endl;
	cout << "   -divsufsort          use divsufsort algorithm to build the BWT (fast, 7.5n Bytes of RAM). By default,"<<endl;
	cout << "                        SE-SAIS is used (about 4 time slower than divsufsort, 4n Bytes of RAM)."<<endl;
	cout << "   <input_file_name>    input text file." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_basename = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-divsufsort")==0){

		sais = false;

	}/*else if(s.compare("-h")==0){

		hyb=true;

	}/*else if(s.compare("-fast")==0){

		fast=true;

	}else if(s.compare("-T")==0){

		T = atoi(argv[ptr]);

		if(T<=0){
			cout << "Error: parameter T must be T>0" << endl;
			help();
		}

		ptr++;
		fast=true;

	}*/else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

	//parse options

    out_basename=string();
    input_file=string();
	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr]);

	if(out_basename.compare("")==0)
		out_basename = string(input_file);

	string idx_file = out_basename;
	idx_file.append(".ri");


	cout << "Building r-index of input file " << input_file << endl;
	cout << "Index will be saved to " << idx_file << endl;

	string input;

	{

		std::ifstream fs(input_file);
		std::stringstream buffer;
		buffer << fs.rdbuf();

		input = buffer.str();

	}

	string path = string(out_basename).append(".ri");
	std::ofstream out(path);

	//save flag storing whether index is fast or small
	out.write((char*)&fast,sizeof(fast));


	if(hyb){

		//auto idx = r_index<sparse_hyb_vector,rle_string_hyb>(input,sais);
		//idx.serialize(out);

	}else{

		auto idx = r_index<>(input,sais);
		idx.serialize(out);

	}


	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;


	out.close();

}

#include <iostream>

#include "internal/r_index.hpp"
#include "internal/utils.hpp"

using namespace ri;
using namespace std;

bool hyb=false;

void help(){
	cout << "ri-space: breakdown of index space usage" << endl;
	cout << "Usage:       ri-space <index>" << endl;
	//cout << "   -h        use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "             if the index was built with -h options enabled."<<endl;
	cout << "   <index>   index file (with extension .ri)" << endl;
	exit(0);
}


void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	/*if(s.compare("-h")==0){

		hyb=true;

	}else*/{

		cout << "Error: unknown option " << s << endl;
		help();

	}

}

int main(int argc, char** argv){

	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	if(hyb){

		r_index<sparse_hyb_vector,rle_string_hyb> idx;
		idx.load_from_file(argv[ptr]);
		auto space = idx.print_space();
		cout << "\nTOT space: " << space << " Bytes" <<endl;

	}else{

		r_index<> idx;
		idx.load_from_file(argv[ptr]);
		auto space = idx.print_space();
		cout << "\nTOT space: " << space << " Bytes" <<endl;

	}




}

#include <iostream>
#include "internal/r_index.hpp"
#include "internal/utils.hpp"

using namespace ri;
using namespace std;

void help(){
	cout << "ri-space: breakdown of index space usage" << endl;
	cout << "Usage:       ri-space <index>" << endl;
	cout << "   <index>   index file (with extension .ri)" << endl;
	exit(0);
}


int main(int argc, char** argv){

	if(argc != 2)
		help();

	r_index idx;
	idx.load_from_file(argv[1]);

	auto space = idx.print_space();

	cout << "\nTOT s-rlbwt space: " << space << " Bytes" <<endl;


}

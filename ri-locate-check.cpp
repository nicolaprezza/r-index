#include <iostream>
#include <r_index.hpp>
#include <utils.hpp>

using namespace ri;
using namespace std;

void help(){
	cout << "ri-locate-check: check correctness of r-index by verifying occurrences on text" << endl;
	cout << "of https://github.com/nicolaprezza/s-rlbwt" << endl << endl;

	cout << "Usage: ri-locate-check <idx> <text> <patterns>" << endl;
	cout << "   <idx>        r-index file name" << endl;
	cout << "   <text>       indexed text file name" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void search(string r_idx_file, string input_file, string patterns){

    r_index r_idx;

    ifstream ifs1(input_file);
    stringstream ss;
    ss << ifs1.rdbuf();//read the file
    string text = ss.str();
    ifs1.close();

    cout << "Loading r-index ... " << flush;
	r_idx.load_from_file(r_idx_file);
	cout << "done." << endl;

	cout << "searching patterns ... " << endl;
	ifstream ifs(patterns);

	//read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
	string header;
	std::getline(ifs, header);

	ulint n = get_number_of_patterns(header);
	ulint m = get_patterns_length(header);

	uint last_perc = 0;

	//extract patterns from file and search them in the indexes
	for(ulint i=0;i<n;++i){

		uint perc = (100*i)/n;
		if(perc>last_perc+4){
			cout << perc << "% done ..." << endl;
			last_perc=perc;
		}


		string p = string();

		for(ulint j=0;j<m;++j){
			char c;
			ifs.get(c);
			p+=c;
		}


		vector<ulint> OCC = r_idx.locate_all(p);

		if(OCC.size() != r_idx.occ(p)){

			cout << "Error: wrong number of located occurrences: " << OCC.size() << "/" << r_idx.occ(p) << endl;
			exit(0);

		}


		for(auto o:OCC){

			if(text.substr(o,p.size()).compare(p) != 0){

				cout << "Error: wrong occurrence: " << o << endl;
				exit(0);

			}

		}

	}

	cout << "Patterns located correctly!" << endl;

	ifs.close();

}

int main(int argc, char** argv){

	if(argc != 4)
		help();

	search(argv[1],argv[2],argv[3]);

}

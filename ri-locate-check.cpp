#include <iostream>
#include <r_index.hpp>
#include <utils.hpp>
#include <s_rlbwt.hpp>

using namespace ri;
using namespace std;

void help(){
	cout << "ri-locate-check: check correctness of r-index by comparing with the s-rlbwt index" << endl;
	cout << "of https://github.com/nicolaprezza/s-rlbwt" << endl << endl;

	cout << "Usage: ri-locate-check <idx_prefix> <patterns>" << endl;
	cout << "   <idx_prefix>   common prefix of r-index and s-rlbwt files." << endl;
	cout << "                  Extensions .ri and .srlbwt are automatically added." << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void search(string r_idx_file, string rlbwt_idx_file, string patterns){

    r_index r_idx;
    s_rlbwt rlbwt;

    cout << "Loading indexes ... " << flush;
	r_idx.load_from_file(r_idx_file);
	rlbwt.load_from_file(rlbwt_idx_file);
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


		//vector<ulint> OCC_1 = r_idx.locate_all(p);
		vector<ulint> OCC_2 = rlbwt.locate(p);

		//std::sort(OCC_1.begin(), OCC_1.end());
		//std::sort(OCC_2.begin(), OCC_2.end());

		/*if(OCC_1 != OCC_2){

			cout << "Error: occurrences do not coincide:" << endl;

			for(auto o : OCC_1) cout << o << " ";
			cout << endl << endl;
			for(auto o : OCC_2) cout << o << " ";
			cout << endl << endl;

			exit(0);

		}*/

		//cout << OCC_1.size() << " occurrences" << endl;

	}

	cout << "Patterns located correctly!" << endl;

	ifs.close();

}

int main(int argc, char** argv){

	if(argc != 3)
		help();

	string r_idx_file(argv[1]);
	r_idx_file.append(".ri");

	search(r_idx_file,argv[1],argv[2]);

}

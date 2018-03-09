#include <iostream>

#include "internal/r_index.hpp"

#include "internal/utils.hpp"

using namespace ri;
using namespace std;

string check = string();//check occurrences on this text

bool hyb=false;

void help(){
	cout << "ri-count: number of occurrences of the input patterns." << endl << endl;

	cout << "Usage: ri-count <index> <patterns>" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   <index>      index file (with extension .ri)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
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


template<class idx_t>
void count(std::ifstream& in, string patterns){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    string text;
    bool c = false;

    if(check.compare(string()) != 0){

    	c = true;

		ifstream ifs1(check);
		stringstream ss;
		ss << ifs1.rdbuf();//read the file
		text = ss.str();

    }

    auto t1 = high_resolution_clock::now();

    idx_t idx;

	idx.load(in);

	auto t2 = high_resolution_clock::now();

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

	ulint occ_tot=0;

	//extract patterns from file and search them in the index
	for(ulint i=0;i<n;++i){

		uint perc = (100*i)/n;
		if(perc>last_perc){
			cout << perc << "% done ..." << endl;
			last_perc=perc;
		}

		string p = string();

		for(ulint j=0;j<m;++j){
			char c;
			ifs.get(c);
			p+=c;
		}

		occ_tot += idx.occ(p);

	}

	double occ_avg = (double)occ_tot / n;

	cout << endl << occ_avg << " average occurrences per pattern" << endl;

	ifs.close();

	auto t3 = high_resolution_clock::now();

	//printRSSstat();

	uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Load time : " << load << " milliseconds" << endl;

	uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	cout << "number of patterns n = " << n << endl;
	cout << "pattern length m = " << m << endl;
	cout << "total number of occurrences  occ_t = " << occ_tot << endl;

	cout << "Total time : " << search << " milliseconds" << endl;
	cout << "Search time : " << (double)search/n << " milliseconds/pattern (total: " << n << " patterns)" << endl;
	cout << "Search time : " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;

}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	string idx_file(argv[ptr]);
	string patt_file(argv[ptr+1]);

	std::ifstream in(idx_file);

	bool fast;

	//fast or small index?
	in.read((char*)&fast,sizeof(fast));

	cout << "Loading r-index" << endl;

	if(hyb){

		count<r_index<sparse_hyb_vector,rle_string_hyb> >(in, patt_file);

	}else{

		count<r_index<> >(in, patt_file);

	}

	in.close();

}

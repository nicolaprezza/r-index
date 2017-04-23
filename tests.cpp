#include <iostream>
#include "r_index.hpp"
#include "utils.hpp"
#include <permutation.hpp>

using namespace ri;
using namespace std;

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    srand(time(NULL));

    auto n = 10000000;

    //vector<ulint> p {3,2,4,5,6,0,1,8,7};

	vector<ulint> p(n);
	for(ulint i=0;i<n;++i) p[i] = i;
	std::random_shuffle(p.begin(),p.end());

	permutation P(p);


	/*for(ulint i =0;i<n;++i){

		cout << i << " -> " << P.inv(i) << endl;

	}

	exit(0);*/


    auto t1 = high_resolution_clock::now();

    ulint xx=0;

	for(ulint i =0;i<n;++i){

		auto r = rand()%n;
		//xx |= r;

		xx |= P.map(r);

	}

    auto t2 = high_resolution_clock::now();

	for(ulint i =0;i<n;++i){

		auto r = rand()%n;
		//xx |= r;

		xx |= P.inv(r);

		//if(P[inv] != r) {cout << "error" << endl;exit(0);}

	}

	cout << xx << endl;

	auto t3 = high_resolution_clock::now();

	ulint map = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	ulint inv = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

	cout << "Map time : " << map << endl;
	cout << "Inv time : " << inv << endl;

}

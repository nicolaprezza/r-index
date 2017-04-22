r-index: the run-length BWT index
===============
Author: Nicola Prezza (nicola.prezza@gmail.com)
Joint work with Travis Gagie

### Brief description

WORK STILL IN PROGRESS!

The r-index is the first full-text index of size O(r), r being the number of BWT runs of the input text (of size n), supporting fast (almost optimal) locate of pattern occurrences. The r-index employs a novel suffix array sampling of size 2r; in classical FM-indexes, this sampling would result in a locate time of Omega(n/r) per occurrence. The r-index, on the other hand, reduces this time to O(log(n/r)).

Let s be the alphabet size and fix a constant eps>0. The r-index offers the following tradeoffs:

- Space: r * ( log s + log(n/r) + (2+eps)*log n ) bits
- Count time: O( (m/eps) * (log (n/r) + log s) )
- Locate time: O(log(n/r)) per occurrence (after count)

Another positive feature is that, during locate, we use only O(1) space on top of the index to locate all occ pattern occurrences. This is not true, e.g., in Lempel-Ziv indexes, which in the worst case use O(occ) space. 

### Download

To clone the repository, use the --recursive option:

> git clone --recursive http://github.com/nicolaprezza/r-index

### Compile

The library has been tested under linux using gcc 5.4.0. You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

We use cmake to generate the Makefile. Create a build folder in the main r-index folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  ri-build input

This command will create the r-index of the text file "input" and will store it as "input.ri". Use option -o to specify a different basename for the index file. 

Run

> ri-count index.ri patterns

to count number of occurrences of the patterns, where <patterns> is a file containing the patterns in pizza&chili format (http://pizzachili.dcc.uchile.cl/experiments.html)

Run

> ri-locate index patterns

to locate all occurrences of the patterns.

r-index: the run-length BWT index
===============
Author: Nicola Prezza (nicola.prezza@gmail.com)
Joint work with Travis Gagie

### Brief description

WORK STILL IN PROGRESS!

The r-index is the first full-text index of size O(r), r being the number of BWT runs of the input text (of size n). Unlikely other BWT-based indexes, the r-index stores only 2r suffix array samples. Despite the reduced suffix array sampling, the r-index is able to locate pattern occurrences in near-optimal time O(m + occ*log(n/r)): this is much faster than classic FM indexes using a regularly sampled suffix array to support locate. 

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

>  ri-build <input>

This command will create the r-index of the text file "input" and will store it as "input.ri". Use option -o to specify a different basename for the index file. 

Run

> ri-count <index.ri> <patterns>

to count number of occurrences of the patterns, where <patterns> is a file containing the patterns in pizza&chili format (http://pizzachili.dcc.uchile.cl/experiments.html)

Run

> ri-locate <index> <patterns>

to locate all occurrences of the patterns.

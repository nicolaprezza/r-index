# NB: WORK STILL IN PROGRESS! the library does not work yet

r-index: the run-length BWT index
===============
Author: Nicola Prezza (nicola.prezza@gmail.com)
Joint work with Travis Gagie

### Brief description

The r-index is the first full-text index of size O(r), r being the number of runs in the BWT of the input text (of size n), supporting fast (almost optimal) locate of pattern occurrences. The r-index employs a novel suffix array sampling of size 2r; in standard FM-indexes, this sampling would result in a locate time of Omega(n/r) per occurrence. The r-index, on the other hand, reduces this time to O(log(n/r)).

Let s be the alphabet size and fix a constant eps>0. The r-index offers the following tradeoffs:

- Space: r * ( log s + log(n/r) + (2+eps)*log n ) bits
- Count time: O( (m/eps) * (log (n/r) + log s) )
- Locate time: After count, O( log(n/r) + 1/eps ) time per occurrence 

Another positive feature is that, during locate, we use only O(1) space on top of the index to locate all occ pattern occurrences. This is not true, e.g., in Lempel-Ziv indexes, which in the worst case use O(occ) space. 







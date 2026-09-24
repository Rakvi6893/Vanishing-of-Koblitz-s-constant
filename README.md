# Congruence Obstructions in the Refined Koblitz Conjecture

This repository contains Magma code accompanying the paper *Congruence obstructions in the refined Koblitz conjecture* by Sung Min Lee, Jacob Mayle, and Rakvi. 

## Installation Instructions

1. Install an up-to-date version of [Magma](https://magma.maths.usyd.edu.au/magma/). We used Magma V2.28-21.
2. Download David Zywina's [Modular](https://github.com/davidzywina/Modular) and [OpenImage](https://github.com/davidzywina/OpenImage) repositories.
3. For certain files (for example, `Common.m`), you will need to modify the file path(s) for your machine to load `Modular` or `OpenImage`.

If you have any questions, comments, or suggestions, please contact us.

## A brief overview of files

1. Common.m contains data and helper routines for the congruence-obstruction computations.
2. GroupSearch.m contains finite group-theoretic search for primitive congruence obstructions.
3. Level 210.m contains computations for Proposition 24. This file restricts the possible mod 2, 3, 5, and 7 images of an elliptic curve with a primitive congruence obstruction modulo 210.
4. PrimObs.m contains a function to compute PrimObs(E) for a non-CM E/Q.
5. RulingOutj.m contains code to verify computational claims in Lemma 19.
6. congobstructions.m verifies that curves in Table 4 have primitive obstruction mod m (which is indicated on the leftmost column).
7. lemma rational points.m verifies computational claims of proposition 27.
8. multipleobstructions.m rules out simultaneous primitive congruence obstruction level possibilities.
9. primecongobstruction.m verifies computational claims of Proposition 10.
10. rational point computations.m computes number of rational points for curves 10.24.1.b.1, 10.24.1.b.2, 14.48.2.e.1, 14.48.2.e.2, 15.48.1.a.1 and 15.48.1.a.2.

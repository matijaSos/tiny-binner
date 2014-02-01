# Testing results

Tests were conducted on CARMA3 synthetic datasets consisting of 25 bacterias. Data is generated with different error models and read lengths.    
Correctness of LCA and MEGAN algorithms are compared.

* 454-400bp

Rank         |  True  | False | True  |  False
-------------|--------|-------|-------|--------
superkingdom |  19992 | 0     | 25000 | 0
phylum       |  19965 | 0     | 24990 | 1
class        |  19950 | 0     | 24989 | 2
order        |  19931 | 0     | 24984 | 8
family       |  19923 | 0     | 24980 | 18
genus        |  17981 | 0     | 24734 | 274
species      |  15851 | 0     | 23992 | 1281

* 454-256bp

Rank         |  True  | False | True  |  False
-------------|--------|-------|-------|--------
superkingdom |  21597 | 0     | 25000 | 0
phylum       |  21542 | 0     | 24991 | 1
class        |  21516 | 0     | 24989 | 3
order        |  21490 | 0     | 24982 | 12
family       |  21479 | 0     | 24978 | 24
genus        |  19513 | 0     | 24565 | 450
species      |  17218 | 0     | 23619 | 1831

* 454-80bp

Rank         |  True  | False | True  |  False
-------------|--------|-------|-------|--------
superkingdom |  21300 | 0     | 24908 | 0
phylum       |  21277 | 0     | 24890 | 7
class        |  21241 | 0     | 24887 | 17
order        |  21206 | 0     | 24876 | 37
family       |  21287 | 1     | 24870 | 63
genus        |  19330 | 3     | 24005 | 954
species      |  16627 | 5     | 22574 | 3288

* Illumina-80bp

Rank         |  True  | False | True  |  False
-------------|--------|-------|-------|--------
superkingdom |  22749 | 0     | 24908 | 0
phylum       |  22726 | 0     | 24890 | 5
class        |  22689 | 0     | 24887 | 15
order        |  22641 | 1     | 24876 | 40
family       |  22620 | 2     | 24870 | 70
genus        |  20678 | 3     | 24005 | 904
species      |  17768 | 3     | 22574 | 3159



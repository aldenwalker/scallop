# scallop

Copyright 2012, 2013, 2014 by Alden Walker

Released under the GPL. 

See INSTALL for installation instructions.


## Background

scallop computes with surfaces in free groups and free products of
cyclic groups.  The appropriate background material is contained 
in [1] (for the -local mode), [2] (for -train), and [3] (for -cyclic and -ball) 

The original version of scallop was written by Danny Calegari, implementing 
the algorithm described in [1], Chapter 4 to compute scl in free groups.  
The current version is a conglomeration of packages written by Alden Walker 
implementing pieces of [1], [2], and [3], all using the basic principle that 
surface maps into free groups factor through fatgraph maps, and we can 
usually produce these via linear programming.



## Installation

See the INSTALL file



## Executing scallop

scallop is run by choosing one of four main modes, and then 
providing more mode-specific arguments:

```
./scallop [-cyclic, -ball, -local, -train] <mode-specific arguments>
```

It is not necessary to explicitly give a mode, in which case 
scallop will default to the original scallop algorithm for free 
groups, which is equivalent to the options -local -pn, where n is 2*rank.

We discuss the four modes separately.  In all cases, giving no 
option-specific arguments will cause scallop to output a list of 
possible options and descriptions.  In an effort 
to produce a timeless document, we'll only discuss the options that 
are non-obvious or particularly important.  

In all modes, the option -v[n] is available, which gives verbose debug 
information.  The optional integer n says how much: default output (1), 
verbose output (2), really verbose (3), and so verbose you probably don't 
want it (4).

  ### A note on speed
  
  Scallop defaults to the nonrigorous original scallop algorithm because 
  for small rank it is much faster than the rigorous -cyclic.  At 
  rank approximately 3 or 4, -cyclic becomes faster.  Asymptotically, 
  -cyclic is polynomial time in the length, rank, and 
  orders of the finite factors.
  
  Gurobi is significantly faster than GLPK, so if possible, it's 
  worthwhile getting the free academic license.


  ### default (no mode specified) 
  
  ```
  ./scallop chain 
  ```

  Examples:
  
  ```
  $ ./scallop abAB
  scl( abAB ) = 1/2 = 0.5
  
  $ ./scallop AbaaBAAABabAAbAAABBBBaaababbabaBaabbABBA
  scl( AbaaBAAABabAAbAAABBBBaaababbabaBaabbABBA ) = 245/146 = 1.67808
  ```

  With no mode specified, scallop defaults to `-local`, and the default 
  `-local` behavior is to compute in a free group with the maximum fatgraph 
  valence set to 2*rank.  This is nonrigorous (but generically correct, and 
  experimentally accurate with probability more than 0.95 for words of 
  length up through 100ish ).  In some cases, this mode can fail:
  
  ```
  $ ./scallop abAAABBB aa bb
  No feasible solution found
  ```

  This is because a fatgraph bounding this chain must have valance at least 5.  
  If the default mode fails, simply run scallop in -cyclic mode to get a rigorous 
  answer:
  ```
  $./scallop -cyclic abAAABBB aa bb
  scl_{a*b}( 1abAAABBB 1aa 1bb ) = 3/4 = 0.75 
 ```

  ### `-cyclic`
  
  ```
  ./scallop -cyclic [options] [group] chain
  ```

  Examples: 
  
  ```
  $ ./scallop -cyclic a0b0 abAB
  scl_{a*b}( 1abAB ) = 1/2 = 0.5
  
  $ ./scallop -cyclic a0b0 abAAABBB aa bb
  scl_{a*b}( 1abAAABBB 1aa 1bb ) = 3/4 = 0.75
  
  $ ./scallop -cyclic a3b2 ab
  scl_{(a/3a)*(b/2b)}( 1ab ) = 1/12 = 0.0833333
  
  $ ./scallop -cyclic -mGUROBI a10b7 ab
  scl_{(a/10a)*(b/7b)}( 1ab ) = 53/140 = 0.378571
  
  $ ./scallop -cyclic -C a0b0 abAABB ab
  cl_{a*b}(1abAABB 1ab ) = 1 = 1
  
  $ ./scallop -cyclic -C a0b0c0 abcAABBCC abc
  cl_{a*b*c}(1abcAABBCC 1abc ) = 2 = 2

  $ ./scallop -cyclic -r 1,2,-1,-2 3,4,-3,-4 5,6,-5,-6 7,8,-7,-8 9,10,-9,-10 11,12,-11,-12 13,14,-13,-14 15,16,-15,-16 17,18,-17,-18 19,20,-19,-20 21,22,-21,-22 23,24,-23,-24 25,26,-25,-26 w2,27,28,-27,-28
  ```

  The 'raw' mode `-r` lets you compute in groups with more than 26 factors.

  The `-cyclic` mode computes scl and cl in free products of cyclic groups.  
  The group string lists the generators and their orders (0 means infinite).  
  The option -C causes scallop to compute commutator length.  This involves 
  an integer programming problem, which is far harder than the linear programming 
  problem used to compute scl.  
  The option -m[solver] lets the user choose GLPK (default), GUROBI (if compiled 
  with support), or EXLP (only available for free groups).  EXLP uses GMP for 
  exact solutions.   
  
  ### `-ball`
  
  ```
  ./scallop -ball [options] filename [group] chain1 , chain2
  ```

  Examples:
  
  ```
  $ ./scallop -ball test.eps abAB , abAABB ab
  Drew ball to file
  
  $ ./scallop -ball -mGUROBI test.eps a5b3 abAB , abAABB ab
  Drew ball to file
  ```

  The `-ball`` mode computes the scl norm ball in the subspace of B_1^H spanned 
  by the given chains.  It defaults to using `-mEXLP` as the solver to avoid 
  rounding errors.  However, to compute a ball in a group with finite cyclic 
  factors, the user must specify `-mGLPK` or `-mGUROBI`.  scallop will 
  write an eps file to whatever filename is specified.
  
  Theoretically, the algorithm can produce the ball in any dimension and 
  output the result as a polyhedron in CDD format (see http://www.inf.ethz.ch/personal/fukudak/cdd_home/).
  However, CDD output is not currently implemented, so the only 
  possibility is a 2d ball output to an eps file.
  
  Make sure to include the output file name!  Omitting it will 
  almost certainly cause scallop to crash (it will assume the first word in 
  the chain is the output file name, which will probably render the rest 
  not homologically trivial).  Also, the commas are necessary!
  
  
  
  ### `-local`
  
  ```
  ./scallop -local [options] chain
  ```

  Examples:
  
  ```
  $ ./scallop -local abAAABBB aa bb
  scl = 3/4 = 0.75
  
  $ ./scallop -local -f abAAABBB aa bb
  No feasible solution found
  
  $ ./scallop -local -p4 abAB
  scl = 1/2 = 0.5
  
  $ ./scallop -local -p4 -e abAB
  Feasible solution found
  
  $ ./scallop -local -p4 abcdABCD
  No feasible solution found
  
  $ ./scallop -local -p8 abcdABCD
  scl = 3/2 = 1.5
  
  $ ./scallop -local -ff1 abAB baBA a.ab.bAABB
  scl = 1/2 = 0.5
  
  $ ./scallop -local -ff1 abAB a.ab.bAABB
  No feasible solution found
  ```

  The -local option searchs for fatgraphs in a free group bounding the chain 
  which have the desired local properties.  The typical options are:
  -f requires that the result be Stallings folded
  -e only checks whether the surface exists, but does not find it (this is faster)
     (this is accomplished by optimizing the function 0 rather than -chi)
  -pn only searches for fatgraphs whose valence is at most n
  -ffn only searches for surfaces whose f-vertices are isolated (see [2]) and 
       partial^- doesn't touch itself, where n is the number of words in 
       partial^-, and the f-vertices in the rest are marked with periods, 
       as indicated above.
  -v[n]: verbose output (level n)
	-y: check if the chain is polygonal (overrides -f,-ff,-p)
       

  
  
  ### -train
  
  ```
  ./scallop -train [-sup, -scl] length chain
  ```

  Examples:
  
  ```
  $ ./scallop -train 3 abAABB ab
  scl( abA bAA AAB ABB BBa Bab aba bab ) = 1/2 = 0.5
  
  $ ./scallop -train -sup 3 abAABB ab
  sup_{Q_3} phi(C)/2D(phi) = inf_t(w + tE) = (t->1/3); 1/2 = 0.5

  $ ./scallop -train 4 abAABB ab
  scl( abAA bAAB AABB ABBa BBab BabA abab baba ) = 2/3 = 0.666667  
  
  $ ./scallop -train -sup -mGUROBI 4 abAABB ab
  sup_{Q_4} phi(C)/2D(phi) = inf_t(w + tE) = (t->523/5529); 2/3 = 0.666667

  $ ./scallop -train  -mGUROBI 4 abAAABBB aa bb
  scl( abAA bAAA AAAB AABB ABBB BBBa BBab BabA 2aaaa 2bbbb ) = 3/4 = 0.75
  
  $ ./scallop -train -sup -mGUROBI 4 abAAABBB aa bb
  sup_{Q_4} phi(C)/2D(phi) = inf_t(w + tE) = (t->35319/1304); 2/3 = 0.666667
  ```

  The `-train` mode computes with traintracks (see [2]).  If the length 
  specified is n, then it finds a minimal surface bounding the set 
  of length n words contained in the given chain.  It can also find 
  the supremum of phi(C)/2D(phi) over all counting quasimorphisms phi, 
  by giving the -sup option.  This linear programming problem is significantly 
  larger.  Note that this can certify an extremal quasimorphism.
  
  The -train mode has several other technical options for finding surfaces 
  that bound w - phi(w) for a collection of words w and counting quasi phi.  
  These ideas are discussed in [2].

  

## TODO / Troubleshooting

There are several pieces of scallop that aren't implemented.  A user 
encountering these issues *should* find that scallop explains itself and 
then exits, but it is possible that scallop will simply crash.

Make sure the options are correctly specified!  Often, a missing option 
will cause scallop to misinterpret what the chain is, which usually 
causes it to crash.

Missing pieces:
 - -cyclic should have the capability to output a fatgraph file 
 - -ball should be able to output an arbitrary dimensional ball
 - -ball should be able to write an eps or povray file for a 3d ball



## References

[1] D. Calegari. scl. MSJ Memoirs, 20. Mathematical Society of Japan, Tokyo, 2009.

[2] D. Calegari and A. Walker. Surface subgroups from Linear programming. preprint: arXiv:1212.2618

[3] A. Walker, Stable commutator length in free products of cyclic groups, Experimental Math 22 (2013), no. 3, 282-298.

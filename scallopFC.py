#!/usr/bin/python

"""this script takes an integer linear combination of words in Z/nZ * Z/mZ, 
generated by x and y, and computes scl of it using scallop by lifting the 
chain to the finite index free normal subgroup which is the kernel of the 
map to Z/nZ x Z/mZ.  It's generated by strings of the form 
x^iy^jx^{n-i}y^{m-j}, so each generator is identified with a pair (i,j), and 
generator (i,j) is given letter i*m + j

the generators can be denoted (i,j,k), where k is 0 or 1 depending on whether 
it's an inverse (1=inverse)

It takes whatever power of the chain is necessary so that it's in the kernel
"""

import fractions


def multiply_words(w1, w2=None):
  if type(w1) == str:
    i = 0
    w1L = len(w1)
    w2L = len(w2)
    while i < w1L and i < w2L and w1[w1L-i-1] == w2[i].swapcase():
      i += 1
    return w1[:w1L-i] + w2[i:]
  elif type(w1) == list:
    return reduce(multiply_words, w1)

def cyc_red(w):
  LW = len(w)
  if len(w) == 0 or len(w) == 1:
    return w
  else:
    i = 0
    while w[i] == w[LW-i-1].swapcase():
      i+=1
    return w[i:LW-i]

  
def inverse(w):
  return w[::-1].swapcase()


def gcd(a,b=None):
  if type(a) == int:
    t1 = max(a,b)
    t2 = min(a,b)
    while t2 != 0:
      t3 = t1 % t2
      t1 = t2
      t2 = t3
    return t1;
  elif type(a) == list:
    return reduce(gcd, a)


def lcm(a,b=None):
  if type(a) == int:
    return a*b/gcd(a,b)
  elif type(w) == list:
    return reduce(lcm, a)


def rewrite_word(w, g1, g2, n, m):
  """rewrite a word in the kernel in terms of gens.  no simplification."""
  #first, convert the word to a list of powers
  pows = [[w[0].lower(), 0]]
  wLen = len(w)
  i=0
  while i < wLen:
    if w[i].lower() == pows[-1][0]:
      if w[i].islower():
        pows[-1][1] += 1
      else:
        pows[-1][1] -= 1
    else:
      if w[i].islower():
        pows.append( [w[i], 1] )
      else:
        pows.append( [w[i].lower(), -1] )
    i += 1
  
  #now go through from the beginning and pull off the generators
  W = []  #denotes the gens, (i,j,I), where I=0 means no inverse, I=1 means inverse
  while len(pows) > 4:
    i,j,k,l = pows[0][1], pows[1][1], pows[2][1], pows[3][1]
    #the next gen is simply the first two characters
    if pows[0][0] == g1:
      W.append( (i, j, 0) )
    else:
      W.append( (j, i, 1) )
    #rewrite the beginning of the word so everything is still in the kernel
    if pows[0][0] == g1:  #this gen is a reguler gen
      if (i+k)%n == 0:    #the a powers are already ok
        if (j+l)%m == 0:  #and the b powers are too!
          pows = pows[4:]
        else:
          pows = [ [g2, ( l-(m-j) )%m] ] + pows[4:]
      else:
        pows = [ [g2, j], [g1, (k-(n-i))%n] ] + pows[3:]
    
    else:                 #this gen is an inverse
      if (i+k)%m == 0:    #the b powers are ok
        if (j+l)%n == 0:  #the a powers are too!
          pows = pows[4:]
        else:
          pows = [ [g1, (l-(n-j))%n] ] + pows[4:]
      else:
        pows = [ [g1, j], [g2, (k-(m-i))%m] ] + pows[3:]
  #end of the while loop
  #there's one remaining; throw it on
  i,j = pows[0][1], pows[1][1]
  if pows[0][0] == g1:
    W.append( (i, j, 0) )
  else:
    W.append( (j, i, 1) )
  if len(pows) != 4:
    print "Word not in the kernel?"
  return W
  
  
#the canonical form is cyclically reduced and lexicographically first
def cyclic_canonical_form(w):
  """return the canonical form of the cyclic word w"""
  

def powers_needed_to_kernelize(C, g1, g2, n, m):
  """return a list of the powers needed to bring each word into the kernel."""
  

def reduce_chain(C)
  """returns [p/q, D], where C = p/q * D, and D is in canonical form."""
  
  
  
  
#if C is a list of words, then that's the chain.  it also accepts 
#[weights, words]
def lift_chain(C, g1, g2, n, m)
  """lift a chain to one in the kernel.  Returns [p/q,c], where 
  c is the chain, and p/q  is the multiplier (i.e. C is p/q times the 
  projection of c)"""
  if type(C[0]) == str:
    new_C = [[1 for i in xrange(len(C))], C]
  else:
    new_C = C
  #first, get a list of the powers needed to get a word into the kernel
  powers_needed = powers_needed_to_kernelize(new_C[1], g1, g2, n, m)
  P = lcm(powers_needed)
  multiplier = fractions.Fraction(1, P)
  c = [ new_C[0], [word_power(w, P) for w in new_C[1]] ]
  c[1] = [rewrite_word(w, g1, g2, n, m) for w in c[1]]
  new_multiplier, c = reduce_chain(c)
  multiplier *= new_multiplier
  
  return [multiplier, c]
  














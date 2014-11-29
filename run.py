#!/usr/bin/env python

from fractions import Fraction
from collections import deque
from itertools import repeat

def reptend(numerator, denominator):
	'''
		returns a tuple (K,R) where

		K is the constant part of the quotient. 
		K is a list whose elements are indexed according to their base in the decimal system.
		So the constant part of the quotient can be reconstituted this way:

		n = len(K) - 1
		k = K[0]*10^0 + K[1]*10^-1 + K[2]*10^-2 + ... + K[n]*10^-n

		R is the repetitive part of the quotient, a.k.a. its cyclic number.
		R is a list whose elements are indexed according to their base in the decimal system,
		multiplied by 10^-n where n is the length of K.

		m = len(R) - 1
		r*10^(n+1) = R[0]*10^0 + R[1]*10^-1 + ... + R[m]*10^m
		
		Example:

		reptend(257, 5*7) => ([7,3], [4,2,8,5,7,1])

		k = 7*1 + 3*0.1 = 7.3
		r = 0.01*( 4*1 + 2*0.1 + 8*0.01 + 5*0.001 + 7*0.0001 + 1*0.00001) = 0.042871

	'''
	if denominator == 1:
		return ([numerator],[])

	remainders = set()
	quotient_first = dict()
	quotients = []
	r = numerator

	while True:
		q,r = divmod(r, denominator)
		quotients.append(q)
		if r not in remainders:
			remainders.add(r)
			quotient_first[r] = len(quotients)
		else:
			q_list = list(quotients)
			first = quotient_first[r]
			return (q_list[0:first], q_list[first:]) # infinite ~ division never terminates
		if r == 0:
			return (quotients, []) # finite ~ division terminates
		if r < denominator:
			r *= 10

def inverse_reptend(K, R):
	'''
		Given a list of quotient and a reptend, returns the corresponding (numerator, denominator)

		Knowns:
		n = len(K) - 1
		m = len(R) -1

		k = S[i=0..n]( K[i]*1e(-i) )
		r = S[i=0..m]( R[i]*1e(-n-1-i) )

		Unknowns:
		a: numerator
		b: denominator

		a - k*b = remainder
		remainder = b*r + 1e(-m-1)*remainder
		remainder(1 - 1e(-m-1)) = b*r
		remainder = b*r/(1 - 1e(-m-1))

		a - k*b = b*r/(1 - 1e(-m-1))
		a = b * [k + r/(1-1e(-m-1))]

		a/b = k + r/(1 - 1e(-m-1))
		a/b = [ k(1 - 1e(-m-1)) + r ] / [ 1 - 1e(-m-1) ]

		# we want integers on both numerator and denominator, so we need to multiply both by 1e(n+m+1).
		# That will allow us to use the GCD algorithm to simplify a/b down to the smallest (a,b) possible

		a/b = [ k(1e(n+m+1) - 1e(n)) + r*1e(n+m+1) ] / [ 1e(n+m+1) - 1e(n) ]
		a/b = [ U + V ] / W

		U = k(1e(n+m+1) - 1e(n)) = S[i=0..n]( K[i]*1e(-i)*(1e(n+m+1) -1e(n)) )
		U = S[i=0..n]( K[i]*1e(n-i)*(1e(m+1) - 1) )

		V = r*1e(n+m+1)
		V = S[i=0..m]( R[i]*1e(m-i) )

		W = 1e(n+m+1) - 1e(n)

	'''

	if not K:
		K = [0]
	if not R:
		R = [0]

	n = len(K) - 1
	m = len(R) - 1
	U = 0
	V = 0
	W = 0

	for (i,v) in enumerate(K):
		U += v*10**(n-i)*( 10**(m+1) - 1)

	for (i,v) in enumerate(R):
		V += v*10**(m-i) 

	W = 10**(n+m+1) - 10**(n)

	f = Fraction( U+V, W)
	# using a Fraction reduces both num and deno by their GCD

	return f.numerator, f.denominator
		
def to_int(reptend):
	''' convert a reptend sequence into a big int
	example: [1,2,3] => 123
	'''
	i = 0
	for n in reptend:
		i *= 10
		i += n
	return i

def to_list(reptend):
	'''
		convert an integer reptend to a list of digits
	'''
	rest = reptend
	digits = deque()
	while rest:
		rest, digit = divmod(rest,10)
		digits.appendleft(digit)
	return tuple(digits)

def hash( reptend ):
	''' the characteristic number of this reptend, i.e. a number based on the histogram of digits
	for instance (1,2,1,2) -> (0,2,2,0,0,0,0,0,0,0)
	'''
	m = dict( [ (i,0) for i in range(10) ] )
	for i in reptend:
		m[i] += 1
	return tuple(m.values())


def randomness(reptend):
	''' counts the number of times a given digit shows up in a reptend.
	This is to verify that long reptends have digits equally distributed
	'''

	count = dict(zip(range(10), repeat(0,10)))
	for digit in reptend:
		count[digit] += 1
	for n,c in count.items():
		count[n] = c*100//len(reptend)
	return count
	
def first_full_reptend(N):
	'''
		returns the first full reptend primes in the range [2..N]
		Example:
			first_full_reptend(10) => [(7, 142857)]
	'''
	all = []
	for n in range(2,N+1):
		_,r = reptend(1,n)
		if len(r) == n-1:
			all.append((n, to_int(r)))
	return all

# some findings:
#
# the period of the multiple of two full reptends is the product of the two periods divided by 2.
# example:
# 7 and 17  -> period(7*17) = 6*16/2 = 48
# 47 and 29 => period(47*29) = 46*28/2 = 644

#!/usr/bin/env python

def reptend(numerator, denominator):
		""" returns the sequence of repeating digits of numerator/denominator """
		if n <= 0:
			return (False, ([],[]) )
		if n == 1:
			return (False, ([1],[]) )

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
				return (True, (q_list[0:first], q_list[first:])) # infinite ~ program never terminates
			if r == 0:
				return (False, (quotients, [])) # finite ~ program terminates
			if r < denominator:
				r *= 10

def del_leading_zeros(reptend):
	''' turn [0, 0, 0, 4, 5, 9] -> [4,5,9] '''
	if not reptend:
		return []
	first_non_zero = 0
	while reptend[first_non_zero] == 0:
		first_non_zero += 1
	return reptend[first_non_zero:]

def to_int(reptend):
	''' convert a reptend sequence into a big int
	example: [1,2,3] => 123
	'''
	i = 0
	for n in reptend:
		i *= 10
		i += n
	return i


def hash( reptend ):
	''' the characteristic number of this reptend, i.e. a number based on the histogram of digits
	for instance (1,2,1,2) -> (0,2,2,0,0,0,0,0,0,0)
	'''
	m = dict( [ (i,0) for i in range(10) ] )
	for i in reptend:
		m[i] += 1
	return tuple(m.values())
	
all = []
for n in range(258):
	(infinite, T) = reptend(1,n)
	if len(T[1]) == n-1:
		all.append((n, to_int(del_leading_zeros(T[1]))))

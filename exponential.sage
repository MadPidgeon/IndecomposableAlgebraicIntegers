load("basic.sage")

def real_poly_check( f, r, a, s=False ):
	result = [True,True]
	if a <= 0.0:
		result[1] = False
	if f(r) > 1.0:
		result = [False,False]
	if f(0)*exp(a*r**2) > 1.0:
		result[0] = False

	X = f.parent().gen()
	B = (f.derivative()-2*a*X*f)
	for x in B.real_roots():
		if 0.0 <= x and f(x)*exp(-a*(x**2-r**2) ) > 1.0:
			if x <= r:
				result[0] = False
			if x >= r:
				result[1] = False
	if s:
		return tuple(result)
	else:
		return result[0] and result[1]

def real_poly_exp_interval( f, r, L, R ):
	if f(r) > 1.0: # Unsalvageble
		return [0.0,-1.0]
	if L == None:
		L, R = -1.0, 1.0
		while real_poly_check( f, r, R, True ) != (False,True): R *= 2
		while real_poly_check( f, r, L, True ) != (True,False): L *= 2
	A, B = L, R
	C, D = L, R
	for _ in range(20):
		M1 = (A+B)/2
		if real_poly_check( f, r, M1, True )[1] == True:
			B = M1
		else:
			A = M1
		M2 = (C+D)/2
		if real_poly_check( f, r, M2, True )[0] == True:
			C = M2
		else:
			D = M2
		if A > D:
			return [0.0,-1.0]
	return [(A+B)/2,(C+D)/2]	

def poly_exp_interval( f, r, places ):
	interval = [None,None]
	P.<X> = RR[]
	for sigma in places:
		g = apply_coefficient_wise( f, lambda x : abs(sigma(x)), P )
		interval = real_poly_exp_interval( g, r, *interval )
		if interval[0] > interval[1]:
			return interval
	return interval

def poly_exp_any( f, r, places ):
	P.<X> = RR[]
	for sigma in places:
		g = apply_coefficient_wise( f, lambda x : abs(sigma(x)), P )
		interval = real_poly_exp_interval( g, r )
		if interval[0] > interval[1]:
			return False
	return True

def is_good( f, r, places = None ):
	P = f.parent()
	K = P.base_ring()
	if places == None:
		places = K.places()
	return poly_exp_any( f, r, places )

def is_extra_good( f, r, places = None ):
	P = f.parent()
	K = P.base_ring()
	if places == None:
		places = K.places()
	interval = poly_exp_interval( f, r, places )
	return ( interval[0] < interval[1], (interval[0]+interval[1])/2 )

class DecompositionException:
	def __init__( self, g ):
		self.polynomial = g

def proves_indec( f, places = None, verbose = 0 ):
	if f.is_constant():
		return False
	X = f.parent().gen()
	a = f.parent().base_ring().gen()
	good, interval = is_extra_good( f( X+a/2 ), sqrt(elm_q( a/2 )), places )
	if verbose > 0:
		print( interval )
	if not good:
		return False
	for g, m in f.factor():
		if relpoly_is_decomp(g):
			raise DecompositionException( g )
	return is_integral( f )

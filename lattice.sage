from sage.rings.complex_field import is_ComplexField
from fpylll import IntegerMatrix, LLL, GSO, CVP, Enumeration, BKZ
import time
load("basic.sage")
load("exponential.sage")

class OrderLattice:
	def __init__( self, order, precision = 5 ):
		self.R = order
		self.prec = precision
		self.K = self.R.number_field()
		self.places = self.K.places( prec=200 )
		self.basis = self.R.basis()
		self.n = len( self.basis )
		self.B = IntegerMatrix.from_matrix( [ self.encode( b ) for b in self.basis ] )
		self.G = GSO.Mat( self.B )
		self.G.update_gso()
	def encode( self, x ):
		v = []
		for sigma in self.places:
			y = sigma(x)*10**self.prec
			if is_ComplexField( sigma.codomain() ):
				y = RR(sqrt(2))*y
				v.append( y.real().integer_part() )
				v.append( y.imag().integer_part() )
			else:
				v.append( y.integer_part() )
		return v

def coeff_list( f, n ):
	L = list(f)
	if len(L) > n:
		raise ValueError("Polynomial degree out of bounds.")
	return L + [0]*(n-len(L))


MODE_GEOMETRIC = 0 
MODE_BINOMIAL = 1
class PolynomialLattice:
	def __init__( self, polynomial_ring, order, n, center=0, radius=1.0, mode=MODE_GEOMETRIC ):
		self.P = polynomial_ring
		self.R = order
		self.c = center
		self.r = radius
		self.mode = mode
		self.order_lattice = OrderLattice( self.R )
		self.degree = n
		self.Y = self.P.gen()
		self.X = self.Y + self.c
		self.basis = [ b*(self.X)**i for i in range(self.degree) for b in self.order_lattice.basis ]
		self.B = IntegerMatrix.from_matrix( [ self.encode( b ) for b in self.basis ] )
		self.G = GSO.Mat( self.B )
		self.G.update_gso()
	def encode( self, f ):
		if self.mode == MODE_GEOMETRIC:
			g = f( self.r.exact_rational()*self.Y )
			return [ coord for coeff in coeff_list( g, self.degree ) for coord in self.order_lattice.encode( coeff ) ]
		else:
			g = f( self.r.exact_rational()*self.Y )
			return [ coord for k, coeff in enumerate( coeff_list( g, self.degree ) ) for coord in self.order_lattice.encode( factorial(k)*coeff ) ]
		
	def close_elements( self, target, sol_max, verbose=0 ):
		try:
			t = self.G.from_canonical( vector( self.encode( target ) ) )
			if verbose > 0:
				print( t )
				print( self.B )
			E = Enumeration( self.G, nr_solutions=sol_max )
			results = E.enumerate( 0, self.B.nrows, 10**30, 0, target=t )
			results_as_polynomials = [ (s, sum([ b*int(v) for (b,v) in zip( self.basis, vector ) ]) ) for s, vector in results ]
			return [ ( f, f(self.Y-self.c), s ) for s, f in results_as_polynomials ]
		except:
			print( "No elements found!" )
			return []
	def short_elements( self, sol_max, verbose=0 ):
		if verbose > 0:
			print( self.B )

		# Basis reduction
		Binv = Matrix( RR, self.B ).inverse()
		A = IntegerMatrix( self.B )
		BKZ.reduction( A, BKZ.Param(20) )
		F = GSO.Mat(A)
		F.update_gso()
		if verbose > 0:
			print(A)

		# Enumeration
		E = Enumeration( F, nr_solutions=sol_max )
		results_wrt_A = E.enumerate( 0, A.nrows, 10**30, 0 )
		results_wrt_B = [ (s, [ v.round() for v in vector( A.multiply_left( vec ) ) * Binv ]) for s, vec in results_wrt_A ]
		results_as_polynomials = [ (s, sum([ b*int(v) for (b,v) in zip( self.basis, vec ) ]) ) for s, vec in results_wrt_B ]

		return [ ( f, f(self.Y-self.c), s ) for s, f in results_as_polynomials ]

def find_szego_polynomial( f, degree, sol_max=40, verbose=0 ):
	start_time = time.time()
	K.<a> = QQ.extension(f)
	c = a / 2
	r = RR(elm_q(a/2)^(1/2))
	P.<Y> = K[]
	X = Y + c
	R = K.OK()
	L = PolynomialLattice( P, R, degree, center=c, radius=r )
	target = Y**degree - X**degree
	if verbose > 0:
		print( "target =", target )

	candidates = L.close_elements( target, sol_max, verbose=verbose-1 )
	end_time = time.time()
	if verbose > 0:
		print( "dtime =", end_time-start_time ) 

	places = L.order_lattice.places
	for g, g_shift, _ in candidates:
		# The resulting close integral polynomial is X^degree + g
		# However, since Y is the variable of the polynomial ring, we write Y^degree + g_shift
		# If so desired, renaming Y to X then gives the right polynomial
		h = Y**degree + g_shift
		val = relpoly_decomp_quality_places( h, places )
		if verbose > 0:
			print( h )
			print( val, h.is_irreducible() )
		if val > 0:
			return (False,h)
	return (None,None)

def find_fekete_polynomial( f, degree, sol_max=40, verbose=0 ):
	start = time.time()
	K.<a> = QQ.extension(f)
	c = a / 2
	r = RR(elm_q(a/2)^(1/2))
	P.<Y> = K[]
	X = Y + c
	R = K.OK()
	L = PolynomialLattice( P, R, degree+1, center=c, radius=r )
	candidates = L.short_elements( sol_max, verbose=verbose-1 )
	middle = time.time()
	if verbose > 0:
		print( "dtime1 =", middle-start )

	places = L.order_lattice.places 
	for g, g_shift, s in candidates:
		h = g_shift
		if verbose > 0:
			print( h )
			print( s )
		try:
			if proves_indec( h, places=places ):
				return ( True, h )
		except DecompositionException as e:
			return ( False, e.polynomial )

	end = time.time()
	if verbose > 0:
		print( "dtime2 =", end-middle )
	return (None,None)

def find_fekete_polynomial_square( f, degree, sol_max=40, verbose=0 ): # Possibly incorrect
	start = time.time()
	K.<a> = QQ.extension(f)
	c = a^2 / 4
	r = elm_q(a/2)
	P.<Y> = K[]
	X = Y + c
	R = K.OK()
	L = PolynomialLattice( P, R, degree+1, center=c, radius=r )
	candidates = L.short_elements( sol_max, verbose=verbose-1 )
	middle = time.time()
	if verbose > 0:
		print( "dtime1 =", middle-start )

	places = L.order_lattice.places 
	for g, g_shift, s in candidates:
		h = g_shift(Y*(Y-a))
		if verbose > 0:
			print( h )
			print( s )
		try:
			if proves_indec( h, places=places ):
				return ( True, h )
		except DecompositionException as e:
			return ( False, e.polynomial )

	end = time.time()
	if verbose > 0:
		print( "dtime2 =", end-middle )
	return (None,None)

def find_fekete_polynomial_binomial( f, degree, sol_max=40, verbose=0 ): # Certainly incorrect
	start = time.time()
	K.<a> = QQ.extension(f)
	c = a^2 / 4
	r = elm_q(a/2)
	P.<Y> = K[]
	X = Y + c
	R = K.OK()
	L = PolynomialLattice( P, R, degree+1, center=c, radius=r, mode=MODE_BINOMIAL )
	candidates = L.short_elements( sol_max, verbose=verbose-1 )
	middle = time.time()
	if verbose > 0:
		print( "dtime1 =", middle-start )

	places = L.order_lattice.places 
	for g, g_shift, s in candidates:
		h = g_shift(Y*(Y-a))
		if verbose > 0:
			print( h )
			print( s )
		try:
			if proves_indec( h, places=places ):
				return ( True, h )
		except DecompositionException as e:
			return ( False, e.polynomial )

	end = time.time()
	if verbose > 0:
		print( "dtime2 =", end-middle )
	return (None,None)

class strategy_find_szego_polynomial:
	def __init__( self, deg_max=10, sol_max=60, verbose=0 ):
		self.deg_max = deg_max
		self.sol_max = sol_max
		self.verbose = verbose
	def __call__( self, f ):
		for deg in range( 1, self.deg_max+1 ):
			if self.verbose > 0:
				print( "--- degree", deg, "---" )
			(is_indec,g) = find_szego_polynomial( f, deg, self.sol_max, verbose=self.verbose-1 )
			if is_indec != None:
				return (f,(is_indec,g))
		return (f,(None,None))

class strategy_find_fekete_polynomial:
	def __init__( self, deg_max=10, sol_max=60, verbose=0 ):
		self.deg_max = deg_max
		self.sol_max = sol_max
		self.verbose = verbose
	def __call__( self, f ):
		for deg in range( 1, self.deg_max+1 ):
			if self.verbose > 0:
				print( "--- degree", deg, "---" )
			(is_indec,g) = find_fekete_polynomial( f, deg, self.sol_max, verbose=self.verbose-1 )
			if is_indec != None:
				return (f,(is_indec,g))
		return (f,(None,None))

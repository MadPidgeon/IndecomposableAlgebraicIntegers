from sage.rings.complex_field import is_ComplexField
from fpylll import IntegerMatrix, LLL, GSO, CVP, Enumeration, BKZ
import time
load("basic.sage")
load("exponential.sage")

class OrderLattice:
	def __init__( self, order, precision = 5 ):
		self.R = order
		self.precision = precision
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
			y = sigma(x)*10**self.precision
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
MODE_SKEW_GEOMETRIC = 2
class PolynomialLattice:
	def __init__( self, polynomial_ring, order, n, center=0, radius=1.0, precision = 8, mode=MODE_GEOMETRIC ):
		self.P = polynomial_ring
		self.R = order
		self.c = center
		self.r = radius
		self.mode = mode
		self.order_lattice = OrderLattice( self.R, precision=precision )
		self.degree = n
		self.Y = self.P.gen()
		self.X = self.Y + self.c
		self.basis = [ b*(self.X)**i for i in range(self.degree) for b in self.order_lattice.basis ]
		self.B = IntegerMatrix.from_matrix( [ self.encode( b ) for b in self.basis ] )
		self.G = GSO.Mat( self.B )
		self.G.update_gso()
		self.A = None
		self.F = None
		self.Binv = None

	def lazy_basis_reduction( self, param=20 ):
		if self.A == None:
			self.Binv = Matrix( RR, self.B ).inverse()
			self.A = IntegerMatrix( self.B )
			BKZ.reduction( self.A, BKZ.Param(param) )
			self.F = GSO.Mat( self.A )
			self.F.update_gso()

	def encode( self, f ):
		rq = self.r.exact_rational()
		g = self.P(f)( rq*self.Y )
		if self.mode == MODE_GEOMETRIC:
			return [ coord for coeff in coeff_list( g, self.degree ) for coord in self.order_lattice.encode( coeff ) ]
		#elif self.mode == MODE_SKEW_GEOMETRIC:
		#	return [ coord for k, coeff in enumerate( coeff_list( g, self.degree ) ) for coord in self.order_lattice.encode( coeff ) ]
		else:
			return [ coord for k, coeff in enumerate( coeff_list( g, self.degree ) ) for coord in self.order_lattice.encode( factorial(k)*coeff * self.degree^(self.degree-k) ) ]

	def decode_from_B( self, vec ):
		return sum([ b*int(v) for (b,v) in zip( self.basis, vec ) ])

	def decode_from_A( self, vec ):
		return self.decode_from_B([ v.round() for v in vector( self.A.multiply_left( vec ) ) * self.Binv ])
		
	def close_elements( self, target, sol_max, verbose=0 ):
		t = self.G.from_canonical( vector( self.encode( target ) ) )
		if verbose > 0:
			print( t )
			print( self.B )

		# Basis reduction
		self.lazy_basis_reduction()
		if verbose > 0:
			print( self.A )

		# Enumeration
		E = Enumeration( self.F, nr_solutions=sol_max )
		results = E.enumerate( 0, self.A.nrows, 10**30, 0, target=t )

		# Decode
		results_as_polynomials = [ ( s, self.decode_from_A( vec ) ) for s, vec in results ]
		return [ ( f, f(self.Y-self.c), s ) for s, f in results_as_polynomials ]
		
	def short_elements( self, sol_max, verbose=0 ):
		t = self.G.from_canonical( vector( self.encode( 0 ) ) )
		if verbose > 0:
			print( self.B )

		# Basis reduction
		self.lazy_basis_reduction()
		if verbose > 0:
			print( self.A )

		# Enumeration
		E = Enumeration( self.F, nr_solutions=sol_max+1 )
		results = E.enumerate( 0, self.A.nrows, 10**30, 0, target=t )[1:]

		# Decode
		results_as_polynomials = [ ( s, self.decode_from_A( vec ) ) for s, vec in results ]
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
	candidates = L.close_elements( 0, sol_max+1, verbose=verbose-1 )[1:]
	# candidates = L.short_elements( sol_max, verbose=verbose-1 )
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
		h = g_shift(Y*(a-Y))
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
		h = g_shift(Y*(a-Y))
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
	def __init__( self, deg_max, deg_min=1, sol_max=50, verbose=0 ):
		self.deg_min = deg_min
		self.deg_max = deg_max
		self.sol_max = sol_max
		self.verbose = verbose
	def __call__( self, f ):
		for deg in range( self.deg_min, self.deg_max+1 ):
			if self.verbose > 0:
				print( "--- degree", deg, "---" )
			(is_indec,g) = find_szego_polynomial( f, deg, self.sol_max, verbose=self.verbose-1 )
			if is_indec != None:
				return (f,(is_indec,g))
		return (f,(None,None))

class strategy_find_fekete_polynomial:
	def __init__( self, deg_max, deg_min=1, sol_max=50, verbose=0 ):
		self.deg_min = deg_min
		self.deg_max = deg_max
		self.sol_max = sol_max
		self.verbose = verbose
	def __call__( self, f ):
		for deg in range( self.deg_min, self.deg_max+1 ):
			if self.verbose > 0:
				print( "--- degree", deg, "---" )
			(is_indec,g) = find_fekete_polynomial_square( f, deg, self.sol_max, verbose=self.verbose-1 )
			if is_indec != None:
				return (f,(is_indec,g))
		return (f,(None,None))

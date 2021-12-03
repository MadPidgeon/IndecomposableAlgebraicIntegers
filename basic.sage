import itertools

### Relative polynomials
def relpoly_to_abspoly( f ):
	K = f.base_ring()
	L.<beta> = K.extension(f)
	return beta.absolute_minpoly()

### Polynomial operations
def apply_coefficient_wise( f, sigma, P = None ):
	if P == None:
		P = f.parent()
	return P(list(map(sigma,list(f))))

def downcast_polynomial( f ):
	L = f.base_ring()
	factors = [ apply_coefficient_wise(f,sigma) for sigma in L.automorphisms() ]
	return reduce(lambda f,g: f*g,factors,1).change_ring( L.base_ring() )

def is_integral( f ):
	for c in list(f):
		if not c.is_integral():
			return False
	return True

def is_monic_integral( f ):
	return f.is_monic() and is_integral(f)

### Compute q
def abspoly_q( f ):
	roots = f.roots( ring=ComplexField(global_precision), multiplicities=False )
	return sum([ abs(x)**2 for x in roots ]) / f.degree()

def elm_q( alpha ):
	return abspoly_q( alpha.absolute_minpoly() )

def relpoly_q( f ):
	K = f.base_ring()
	L.<beta> = K.extension( f )
	return elm_q( beta )

def elm_q_places( alpha, places ):
	return sum([ (1+is_ComplexField(sigma.codomain()))*abs(sigma(alpha))**2 for sigma in places]) / alpha.parent().absolute_degree()

def relpoly_q_places( f, places ): 
	val = 0
	for sigma in places:
		P.<Y> = sigma.codomain()[]
		g = apply_coefficient_wise( f, sigma, P )
		if is_ComplexField( sigma.codomain() ): # WARNING: If you type K.places(200) instead of K.places(prec=200), you wrong answers here
			val += 2*sum([ abs(x)**2 for x in g.complex_roots() ])
		else:
			val += sum([ abs(x)**2 for x in g.complex_roots() ])
	return val / ( places[0].domain().absolute_degree() * f.degree() )

### Decomposition quality
def elm_decomp_quality( alpha, beta ):
	if beta==0 or beta==alpha:
		return -1000000
	return elm_q( alpha ) - elm_q( beta ) - elm_q( alpha-beta )

def relpoly_decomp_quality( f ):
	K = f.base_ring()
	alpha = K.gen()
	if not f.is_irreducible():
		return -1000000
	L.<beta> = K.extension( f )
	return elm_decomp_quality( alpha, beta )

def abspoly_decomp_quality( f, g ):
	K.<alpha> = QQ.extension( f )
	P.<Y> = K[]
	g2 = g(Y).factor()[0][0]
	return relpoly_decomp_quality( g2 )

def relpoly_decomp_quality_places( f, places ):
	K = f.base_ring()
	alpha = K.gen()
	val = elm_q_places( alpha, places ) - relpoly_q_places( f, places ) - relpoly_q_places( f( alpha - f.parent().gen() ), places )
	if val >= 0:
		if not f.is_irreducible():
			return -1000000
	return val

### Check decomposition
def elm_is_decomp( alpha, beta ):
	return elm_decomp_quality( alpha, beta ) >= -10**(-50) and beta.is_integral()

def relpoly_is_decomp( f ):
	return relpoly_decomp_quality( f ) >= -10**(-50) and is_monic_integral(f)

def abspoly_is_decomp( f, g ):
	return abspoly_decomp_quality( f, g ) >= -10**(-50) and is_monic_integral(g)

### Polynomial candidates of bounded q
def maclaurin_bound( deg, q = 4 ):
	RRR = RealField(global_precision)
	P.<X> = RRR[]
	f = (X+(q+RRR(0.00000001))**.5)^deg
	return [ x.floor() for x in list(f) ]

def maclaurin_candidates( deg, P, q = 4 ):
	b = maclaurin_bound( deg, q )
	return [ P(l) for l in itertools.product(*([range(-r,r+1) for r in b[:-1]]+[range(1,2)])) ]

def inverse_euler_phi_max(n):
	c = 0
	b = n
	while b > 1:
		b = euler_phi(b)
		c = c+1
	m = 2*3**c
	while euler_phi(m) > n:
		m = m-1
	return m

def get_all_candidates( deg, X = None, q = 4, statistics = False ):
	P = None
	if X == None:
		P.<X> = QQ[]
	else:
		P = X.parent()
	L1 = maclaurin_candidates( deg, P, q )
	L2 = list(filter(lambda f: f.is_irreducible(), L1 ))
	L3 = list(filter(lambda f: abspoly_q(f) <= q, L2 ))
	if statistics:
		return (L3, (len(list(L1)), len(list(L2)), len(list(L3))) )
	return L3

def get_all_candidates_up_to_isometry( deg, X = None, q = 4, statistics = False ):
	(L,stats) = get_all_candidates( deg, X, q, True )
	s = []
	for n in range(2,inverse_euler_phi_max(deg*deg)):
		K = CyclotomicField( n )
		Q.<Y> = K[]
		roots = [1]+[K.gen()**k for k in range(n) if gcd(k,n)==1 ]
		orbits = dict([(frozenset([f(zeta*Y)*(zeta)**(-deg) for zeta in roots]),f) for f in L] )
		L = orbits.values()
		s.append(len(L))
	if statistics:
		return (L,stats,n,s)
	return L

def candidate_count_table( max_deg ):
	P.<X> = QQ[]
	for d in range( 1, max_deg+1 ):
		print( "degree", d )
		mlc = maclaurin_candidates( d, P )
		print( "maclaurin", len(mlc) )
		qboundc = [ f for f in mlc if abspoly_q_slow(f) <= 4.0 ]
		print( "qbound", len(qboundc) )
		irredc = [ f for f in qboundc if f.is_irreducible() ]
		print( "irreducible", len(irredc) )



import itertools
import pickle
from multiprocessing import Pool 
load("exponential.sage")

### Read and write results

# Datastructure is a dict with polynomials in X as key and a pair (is_indecomposable,witness) as value, 
# where is_decomposable is a boolean and witness is a polynomial, with a decomposition of as root when is_indecomposable is False.
# Values may be (None,None) for undetermined
def write_file( data, filename, overwrite=True ):
	if not overwrite:
		try:
			file = open( filename, 'r' )
			file.close()
			return False
		except:
			pass
	file = open( filename, 'wb' )
	pickle.dump( data, file )
	file.close()
	return True

def read_file( filename ):
	try:
		file = open( filename, 'rb' )
		data = pickle.load( file )
		file.close()
		return data
	except:
		return None

def create_file( deg, filename_prefix='polynomial_data', overwrite=False ):
	data = dict([(candidate,(None,None)) for candidate in get_all_candidates_up_to_isometry(deg)])
	filename = filename_prefix + '_deg_' + str(deg) + '.txt'
	if write_file( data, filename, overwrite ):
		return filename

### Applying strategies
def apply_strategy( data, strat, verbose=2 ):
	n = len(data)
	for m, f in enumerate(data):
		if verbose == 1:
			sys.stdout.write("\rApplying to polynomial %i/%i..." % (m+1,n))
			sys.stdout.flush()
		elif verbose > 1:
			sys.stdout.write("Applying to polynomial %i/%i: " % (m+1,n))
			print( f )
			sys.stdout.flush()
		if data[f][0] == None:
			try:
				data[f] = strat( f, verbose=verbose-1 )
				if data[f][0] != None:
					print( "FOUND POLYNOMIAL:", data[f][1] )
			except:
				pass
	if verbose>0:
		sys.stdout.write("\rDone!                             \n")
		sys.stdout.flush()
	return data

def apply_strategy_parallel( data, strat, workers = 2 ):
	flatten_data = [ f for (f,(a,b)) in data.items() if a == None ]
	n = len(flatten_data)
	rdata = []
	p = Pool( workers )
	for m, retval in enumerate( p.imap_unordered( strat, flatten_data, 1) ):
		rdata.append(retval)
		if retval[1][0] != None:
			print( "FOUND POLYNOMIAL:", retval[1][1] )
		sys.stderr.write("\rApplying to polynomial %i/%i..." % (m+1,n))
		sys.stderr.flush()
	sys.stderr.write("\rDone!                             \n")
	sys.stderr.flush()
	return dict(rdata)

### Validate data
def validate_data( data, verbose = 0 ):
	P.<X> = QQ[]
	for f, (is_indecomposable,g) in data.items():
		if is_indecomposable == True:
			if g == None:
				if abspoly_q( f(X) ) >= 2:
					return f
			else:
				if not proves_indec( g, verbose=verbose-1 ):
					return f
		elif is_indecomposable == False:
			if not abspoly_is_decomp( f(X), g(X) ):
				return f
		else:
			if verbose>0:
				print( f, "skipped" )
	return None

def quality_data( data ):
	return dict([ (f,abspoly_decomp_quality(f,g)) for (f,(is_indecomposable,g)) in data.items() if is_indecomposable == False ])

def split_data( data ):
	return ( dict([ x for x in data.items() if x[1][0] == None ]),
		dict([ x for x in data.items() if x[1][0] == True ]),
		dict([ x for x in data.items() if x[1][0] == False ])
	  )

#def todo_data( data ):
#	return (len(filter(lambda x:x[0]==None,data.values())),len(data))

def todo_data( data ):
	return dict([(f,(u,g)) for (f,(u,g)) in data.items() if u == None])

def q_data( data ):
	return dict([(f,abspoly_q(f)) for f in data])

def diff_data( data1, data2, count = False ):
	if set(data1) != set(data2):
		return None
	diff = dict([(f,(data1[f],data2[f])) for f in data1 if data1[f][0] != data2[f][0] ])
	if count:
		return (len([f for f in data1 if data1[f][0] != None and data2[f][0] == None ]),
				len([f for f in data1 if data1[f][0] == None and data2[f][0] != None ]))
	else:
		return diff

def q_sort( data ):
	return sorted([ (abspoly_q(f).numerical_approx(digits=5),(f,g)) for f,g in data.items() ])

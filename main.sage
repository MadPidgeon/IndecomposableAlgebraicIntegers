load("lattice.sage")
load("bulk.sage")

## Load data
data = read_file("polynomial_data_deg_3.txt")
if data == None:
	print("Data file cannot be found. Creating new data file...")
	data = read_file( create_file(3) )
	print("Done")
todo = todo_data(data)

## Set environment
P.<X> = QQ[]

## Print how to
print(
"The variable data contains a dictionary with all (up to isometry) irreducible polynomials over the integers of degree 3 and square norm up to 4 as keys. "
"The value at a polynomial f is a pair (is_indec,g), where is_indec is either True (if f is indecomposable), False (if f is decomposable) or None (if it is not yet known whether f is indecomposable), and g is a polynomial over the ring of integers of Q[X]/f. "
"If is_indec is True, then g is exponentially bounded, proving f is indecomposable. "
"If is_indec is False, then g is the minimal polynomial of a decomposition of f. "
"The variable todo contains the subdictionary of data of all entries where is_indec is None. "
)

print()
print(
"The function validate_data validates whether a subdictionary of data is correct. "
"The function apply_strategy (resp. apply_strategy_parallel) apply to a subdictionary of data a given function that attempts to compute a correct value for each key. "
"Available strategies include strategy_find_szego_polynomial and strategy_find_fekete_polynomial from lattice.sage. "
)
print("Example:")
print("sage: results = apply_strategy(todo,strategy_find_szego_polynomial(10))")
print("sage: data = {**data,**results}")
print("sage: write_file(data,\"polynomial_data_deg_3.txt\")")


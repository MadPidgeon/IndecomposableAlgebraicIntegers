load("main.sage")
f = X^3-2*X^2-X+1
K.<a> = QQ.extension(f)
R = K.OK()
Q.<Y> = K[]
L = PolynomialLattice(Q,R,6,center=a/2,radius=RR(elm_q(a/2)^(1/2)))
print( L.close_elements(Y+a/2,10) )
print( L.short_elements(10) )
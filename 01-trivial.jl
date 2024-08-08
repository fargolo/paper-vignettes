using Pkg; 
Pkg.add(["AnalyticComb", "SymPy"])

using AnalyticComb, SymPy
# Defines z as a symbolic variable

z = SymPy.symbols("z")
# Creates the generating function in W
W = SEQ(2z)
#   1   
#-------
#1 - 2z
# Obtains the Taylor series expansion
series(W)
#1 + 2z + 4z^2  + 8z^3  + 16z^4  + 32z^5  + O(z)^6
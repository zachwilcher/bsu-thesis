from sympy import *

# utility function to write a sympy expression to a file in latex form
def write_to_file(filename, expr):
    with open(filename, 'w') as f:
        f.write(latex(expr))

#
# get the expression for the number of 3-matchings in the pulsar graph
#
E  = symbols('E')


n1, m, n2 = symbols('n_1 m n_2')

# number of edges in pulsar graph
E_expr = n1 + n1 + m + m + n2 + n2
E_expr = expand(E_expr)

three_matchings = binomial(n1 + n2, 3) \
    + (n1 + n2) * binomial(n1 + n2 - 1, 2) \
    + 2 * m * binomial(n1 + n2, 2) \
    + n1 * n2 * (n1 + n2 - 2) \
    + m * (n1 + n2) * (n1 + n2 - 1) \
    + m * (m - 1) * (n1 + n2)

two_matchings = binomial(n1 + n2, 2) \
    + (n1 + n2) * (n1 + n2 - 1) \
    + 2 * m * (n1 + n2) \
    + (n1 + m) * (n2 + m) - m

one_matchings = E
one_matchings = expand(simplify(one_matchings.subs([(E, E_expr)])))


zero_matchings = 1

#
# compute the number k-cubes in the 3-point configuration space of the pulsar graph
#

# number of vertices in the graph
V = n1 + n1 + 1 + m + 1 + n2 + n2

zero_cubes = zero_matchings * binomial(V - 2 * 0, 3 - 0)
one_cubes = one_matchings * binomial(V - 2 * 1, 3 - 1)
two_cubes = two_matchings * binomial(V - 2 * 2, 3 - 2)
three_cubes = three_matchings * binomial(V - 2 * 3, 3 - 3)

euler_characteristic = expand(simplify(zero_cubes - one_cubes + two_cubes - three_cubes))

write_to_file("pulsar_graph_euler.txt", euler_characteristic)

#
# Subtract \Chi(UConf_3(pulsar-graph)) - \Chi(UConf_3(theta-graph))
#

# formula from Justin's paper
theta_euler = m * (m - 2) * (m - 7) / 6

write_to_file("theta_euler.txt", expand(simplify(theta_euler)))

difference = expand(simplify(euler_characteristic - theta_euler))

write_to_file("difference.txt", difference)

n2_zero = expand(simplify(difference.subs([(n2, 0)])))

write_to_file("n2_zero.txt", n2_zero)

print(f"Euler characteristic of UConf_3(Pulsar graph(5,3,4)):{ euler_characteristic.subs([(n1, 5), (m, 3), (n2, 4)])}")

counted_generators = n1 * (m - 1) + n1 * (n1 - 1) / 2 + n1 * (n1 - 1) * (m - 1) / 2 + n1 * (n1 - 1) * (n1 - 2) / 6 * 2
write_to_file("counted_generators.txt", expand(simplify(counted_generators)))
counted_generators += n1 * (m - 1) + n1 * (n1 - 1)
write_to_file("counted_generators2.txt", expand(simplify(counted_generators)))

counted_generators += counted_generators.subs([(n1, n2)])
write_to_file("counted_generators3.txt", expand(simplify(counted_generators)))
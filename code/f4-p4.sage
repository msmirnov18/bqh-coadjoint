# Computations for the homogenous variety F4/P4 (Section 7 of the paper)

reset()

import subprocess
subprocess.call('git clone https://github.com/msmirnov18/littlewood-richardson-rule.git', shell = True)

load('./littlewood-richardson-rule/littlewood-richardson.sage')

# Computing 4-point Gromov--Witten invariants for Proposition 7.3.
# This computation takes place in the classical cohomology ring of the homogenous variety G(3,6) = A5/P3.
print('')
print('Computing 4-point Gromov--Witten invariants for Proposition 7.3.')
print('This computation takes place in the classical cohomology ring of the homogenous variety D5/P4')
print('')

# Defining the cohomology ring of D5/P4
F = CohomologyPartialFlagVariety('D5', (4,), QQ)

# define relabelling from E6 to D5
R = dict({6: 1, 5: 2, 4: 3, 3: 4, 2: 5})
Rinv = dict((b, a) for a,b in R.items())


def to_reduced_words(input):

    return dict((tuple(term.leading_support().reduced_word()), term.leading_coefficient()) for term in input.terms())


# Defining shatbar
shatbar = F.module.monomial(F.weyl_group.from_reduced_word([R[2],R[4],R[3]]))
print('Defining the cohomology class shatbar:', shatbar)
print('')

# defining the hyperplane class h
hyperplane_class = F.module.monomial(F.weyl_group.from_reduced_word([4]))

# Computing shatbar^3 * h
print('Computing shatbar^3 * h:')
print('shatbar^3 * h = ', F.cup_product(hyperplane_class, F.cup_product(shatbar,F.cup_product(shatbar,shatbar))))
print('')

# Computing the Poincare dual of shatbar * h
print('Computing the Poincare dual of shatbar * h:')
pd_shatbar_h = F.poincare_dual(F.cup_product(shatbar, hyperplane_class))
pd_shatbar_h_reduced_words = to_reduced_words(pd_shatbar_h)
print(pd_shatbar_h)
print('')





# Root notation for the lift of Poicaré dual of shatbar * h from D5/P4 to E6/P1 and subsequent restriction to F4/P4.
# This computation takes place in the classical cohomology ring of the homogenous variety F4/P4.
print('')
print('Root notation for the lift of Poicaré dual of shatbar * h from D5/P4 to E6/P1 and subsequent restriction to F4/P4.')
print('This computation takes place in the classical cohomology ring of the homogenous variety F4/P4.')
print('')

# define relabelling from E6 to F4
S = dict({1: 4, 6: 4, 3: 3, 5: 3, 4: 2, 2: 1})

# Defining the cohomology ring of F4/P4
H = CohomologyPartialFlagVariety('F4', (4,), QQ)
print(H)
print('')


def lift_restrict_from_reduced_words(input):

    lift_to_E6_P1 = dict((tuple([Rinv[i] for i in word] + [1]), coeff) for (word,coeff) in input.items())

    restriction_to_F4_P4 = sum(coeff*H.module.monomial(H.weyl_group.from_reduced_word([S[i] for i in word])) for (word,coeff) in lift_to_E6_P1.items())

    return restriction_to_F4_P4


print('The restriction of the lift of the Poicaré dual of shatbar * h is:')
print(lift_restrict_from_reduced_words(pd_shatbar_h_reduced_words))
print('In the root notation we have:')
print(H.to_roots(lift_restrict_from_reduced_words(pd_shatbar_h_reduced_words)))
print('')




# Computing auxiliary cohomology classes appearing in the proof of Proposition 7.4.
# This computation takes place in the classical cohomology ring of the homogenous variety F4/P4.


# defining the hyperplane class h
h = H.module.monomial(H.weyl_group.from_reduced_word([4]))

# defining the cohomology class s
s = H.module.monomial(H.weyl_group.from_reduced_word([1,2,3,4]))

def raise_to_power(element, power):

    assert element in H.module, ''
    assert power >= 0, ''

    if power == 0:
        return H.module.monomial(H.weyl_group.from_reduced_word([]))

    else:
        return H.cup_product(element, raise_to_power(element, power - 1))

# Computing the class \sigma appearing in the proof of Proposition 7.4.
print('Computing the class \sigma appearing in the proof of Proposition 7.4:')
sigma = -2*raise_to_power(h, 7) + 6*H.cup_product(s, raise_to_power(h, 3))
print(sigma)
print('In the root notation we have:')
print(H.to_roots(sigma))

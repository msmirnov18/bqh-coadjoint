# Computations for the homogenous variety E8/P8 (Section 6.3 of the paper)

reset()

import subprocess
subprocess.call('git clone https://github.com/msmirnov18/littlewood-richardson-rule.git', shell = True)

load('./littlewood-richardson-rule/littlewood-richardson.sage')

# Computations for Lemma 6.19.
print('')
print('Computing cup-products appearing in Lemma 6.19.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E8/P8.')
print('')

# Defining the cohomology ring of E8/P8
H = CohomologyPartialFlagVariety('E8', (8,), QQ)
print(H)
print('')

# Defining the cohomology class s
s = H.module.monomial(H.weyl_group.from_reduced_word([2, 4, 5, 6, 7, 8]))
print('Defining the cohomology class s:', s)
print('')

# Defining the cohomology class t
t = H.module.monomial(H.weyl_group.from_reduced_word([6, 5, 4, 3, 2, 4, 5, 6, 7, 8]))
print('Defining the cohomology class t:', t)
print('')

# Computing t^2
print('Computing t^2:')
t_squared = H.cup_product(t,t)
print('t^2 =', t_squared)
print('')
print('In the root notation we have:')
print('t^2 =', H.to_roots(t_squared))
print('')

# Computing s^2
print('Computing s^2:')
s_squared = H.cup_product(s,s)
print('s^2 =', s_squared)
print('')
print('In the root notation we have:')
print('s^2 =', H.to_roots(s_squared))
print('')

# Computing s^3
print('Computing s^3:')
s_cubed = H.cup_product(s_squared,s)
print('s^3 =', s_cubed)
print('')
print('In the root notation we have:')
print('s^3 =', H.to_roots(s_cubed))
print('')

# Computing s^4
print('Computing s^4:')
s_4 = H.cup_product(s_cubed,s)
print('s^4 =', s_4)
print('')
print('In the root notation we have:')
print('s^4 =', H.to_roots(s_4))
print('')




# Computing 4-point Gromov--Witten invariants for Proposition 6.20.
# This computation takes place in the classical cohomology ring of the homogenous variety E7/P7.
print('')
print('Computing 4-point GW invariants for Proposition 6.20.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E7/P7.')
print('')

# Defining the cohomology ring of E7/P7
F = CohomologyPartialFlagVariety('E7', (7,), QQ)
print(F)
print('')

def to_reduced_words(input):

    return dict((tuple(term.leading_support().reduced_word()), term.leading_coefficient()) for term in input.terms())

# defining tbar
tbar = F.module.monomial(F.weyl_group.from_reduced_word([6, 5, 4, 3, 2, 4, 5, 6, 7]))
print('Defining the cohomology class tbar:', tbar)
print('')

# Computing tbar^3
print('Computing tbar^3:')
print('tbar^3 = ', F.cup_product(tbar,F.cup_product(tbar,tbar)))
print('')

# Computing the Poincare dual of tbar
print('Computing the Poincare dual of tbar:')
pd_tbar = F.poincare_dual(tbar)
pd_tbar_reduced_words = to_reduced_words(pd_tbar)
print(pd_tbar)
print('')

# defining sbar
sbar = F.module.monomial(F.weyl_group.from_reduced_word([2, 4, 5, 6, 7]))
print('Defining the cohomology class sbar:', sbar)
print('')

# Computing the Poincare dual of sbar
print('Computing the Poincare dual of sbar:')
pd_sbar = F.poincare_dual(sbar)
pd_sbar_reduced_words = to_reduced_words(pd_sbar)
print(pd_sbar)
print('')

# Computing sbar^2
print('Computing sbar^2:')
sbar_squared = F.cup_product(sbar,sbar)
print('sbar^2 = ', sbar_squared)
print('')

# Computing the Poincare dual of sbar^2
print('Computing the Poincare dual of sbar^2:')
pd_sbar_squared = F.poincare_dual(sbar_squared)
# pd_sbar_squared_reduced_words = [term.leading_support().reduced_word() for term in pd_sbar_squared.monomials()]
pd_sbar_squared_reduced_words = to_reduced_words(pd_sbar_squared)
print(pd_sbar_squared)
print('')



# Root notation for lifts from E7/P7 to E8/P8 of Poicaré duals of tbar and sbar.
# This computation takes place in the classical cohomology ring of the homogenous variety E8/P8.
print('')
print('Root notation for lifts from E7/P7 to E8/P8 of Poicaré duals of tbar, sbar, and sbar^2.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E8/P8.')
print('')

def lift_from_reduced_words(input):

    return sum(coeff*H.module.monomial(H.weyl_group.from_reduced_word([i for i in word] + [8])) for (word,coeff) in input.items())


print('The lift of the Poicaré dual of tbar is:')
print(lift_from_reduced_words(pd_tbar_reduced_words))
print('In the root notation we have:')
print(H.to_roots(lift_from_reduced_words(pd_tbar_reduced_words)))
print('')

print('The lift of the Poicaré dual of sbar is:')
print(lift_from_reduced_words(pd_sbar_reduced_words))
print('In the root notation we have:')
print(H.to_roots(lift_from_reduced_words(pd_sbar_reduced_words)))
print('')

print('The lift of the Poicaré dual of sbar^2 is:')
print(lift_from_reduced_words(pd_sbar_squared_reduced_words))
print('In the root notation we have:')
print(H.to_roots(lift_from_reduced_words(pd_sbar_squared_reduced_words)))
print('')




# Computing auxiliary cohomology classes appearing in the proof of Proposition 6.21.
# The answer is given directly in the root notation.
# This computation takes place in the classical cohomology ring of the homogenous variety E8/P8.
print('')
print('Computing auxiliary cohomology classes appearing in the proof of Proposition 6.21.')
print('The answer is given directly in the root notation.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E8/P8.')
print('')

# Computing the inverse image of 3t^2 under the multiplication by the hyperplane class.
print('Inverse image of 3t^2 under the multiplication by the hyperplane class:')
print(H.to_roots(H.inverse_lefschetz_map(3*t_squared)))
print('')

# Computing the inverse image of 5s^4 under the multiplication by the hyperplane class.
print('Inverse image of 5s^4 under the multiplication by the hyperplane class:')
print(H.to_roots(H.inverse_lefschetz_map(5*s_4)))
print('')

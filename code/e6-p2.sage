# Computations for the homogenous variety E6/P2 (Section 6.1 of the paper)

reset()

import subprocess
subprocess.call('git clone https://github.com/msmirnov18/littlewood-richardson-rule.git', shell = True)

load('./littlewood-richardson-rule/littlewood-richardson.sage')

# Computations for Lemma 6.5
print('')
print('Computing cup-products appearing in Lemma 6.5.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E6/P2.')
print('')

# Defining the cohomology ring of E6/P2
H = CohomologyPartialFlagVariety('E6', (2,), QQ)
print(H)
print('')

# Defining the cohomology class s
s = H.module.monomial(H.weyl_group.from_reduced_word([3,4,2]))
print('Defining the cohomology class s:', s)
print('')

# Defining the cohomology class t
t = H.module.monomial(H.weyl_group.from_reduced_word([1,3,4,2]))
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




# Computing 4-point Gromov--Witten invariants for Proposition 6.6.
# This computation takes place in the classical cohomology ring of the homogenous variety G(3,6) = A5/P3.
print('')
print('Computing 4-point GW invariants for Proposition 6.6.')
print('This computation takes place in the classical cohomology ring of the homogenous variety G(3,6) = A5/P3.')
print('')

# Defining the cohomology ring of G(3,6) = A5/P3
F = CohomologyPartialFlagVariety('A5', (3,), QQ)
print(F)
print('')

# Define relabelling from E6 to A5
R = dict({6: 5, 5: 4, 4: 3, 3: 2, 1: 1})
Rinv = dict((b, a) for a,b in R.items())

def to_reduced_words(input):

    return dict((tuple(term.leading_support().reduced_word()), term.leading_coefficient()) for term in input.terms())

# Defining tbar
tbar = F.module.monomial(F.weyl_group.from_reduced_word([R[1],R[3],R[4]]))
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

# Defining sbar
sbar = F.module.monomial(F.weyl_group.from_reduced_word([R[3],R[4]]))
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
pd_sbar_squared_reduced_words = to_reduced_words(pd_sbar_squared)
print(pd_sbar_squared)
print('')


# Root notation for lifts from A5/P3 to E6/P2 of Poicaré duals of tbar and sbar.
# This computation takes place in the classical cohomology ring of the homogenous variety E6/P2.
print('')
print('Root notation for lifts from A5/P3 to E6/P2 of Poicaré duals of tbar and sbar.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E6/P2.')
print('')

def lift_from_reduced_words(input):

    return sum(coeff*H.module.monomial(H.weyl_group.from_reduced_word([Rinv[i] for i in word] + [2])) for (word,coeff) in input.items())


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



# Computing auxiliary cohomology classes appearing in the proof of Proposition 6.7.
# The answer is given directly in the root notation.
# This computation takes place in the classical cohomology ring of the homogenous variety E6/P2.
print('')
print('Computing auxiliary cohomology classes appearing in the proof of Proposition 6.7.')
print('The answer is given directly in the root notation.')
print('This computation takes place in the classical cohomology ring of the homogenous variety E6/P2.')
print('')

# Computing the inverse image of 6t^2 under the multiplication by the hyperplane class.
print('Inverse image of 6t^2 under the multiplication by the hyperplane class:')
print(H.to_roots(H.inverse_lefschetz_map(6*t_squared)))
print('')

# Computing the inverse image of 2s^3 under the multiplication by the hyperplane class.
print('Inverse image of 2s^3 under the multiplication by the hyperplane class:')
print(H.to_roots(H.inverse_lefschetz_map(2*s_cubed)))
print('')

# SageNP: Newman-Penrose calculations for SageMath. The code is based on SageManifolds.

####################################################
# SageNP: Newman-Penrose Calculations for SageMath # 
####################################################
# Code by:                                         #
# Tolga Birkandan (Corr.: birkandant@itu.edu.tr)   #
# Emir Baysazan                                    #
# Pelin Ozturk                                     #
# Special thanks to Eric Gourgoulhon               #
# Based on SageMath (SageManifolds)                #
####################################################
r"""
####################################################

The class "SageNP" includes functions for some 
calculations defined in the Newman-Penrose formalism.

###########
REFERENCE:
The reference for all definitions and calculations:

H. Stephani, D. Kramer, M. MacCallum, 
C. Hoenselaers, and E. Herlt, 
"Exact Solutions of Einstein’s Field Equations", 
2nd ed. Cambridge: Cambridge University Press, 2003.
###########

BASIC DEFINITIONS AND NOTATION:
* We will use the Metric signature: (- + + +)                      

* For the null-tetrad vector names,
- The ref. book uses 
    (k,l,m,mbar)

- However, in the code we will use
    (l,n,m,mbar)
  like the rest of literature.

- Therefore on should set 
    k->l, l->n in the ref. book

- Products of the vectors are given by:
    l*n = -1, m*mbar = 1, all others zero.  
  
- The metric is found using the covariant null-tetrad
  vectors as:
    g = -2*l*n + 2*m*mbar                            
  and,
    g = [0  1  0  0]                                 
        [1  0  0  0]                                 
        [0  0  0 -1]                                 
        [0  0 -1  0]

- Please check the reference book for the details 
  and further definitions.
###########

INSTRUCTIONS WITH EXAMPLES:

- Import the class:
    from SageNP import *
    
- Define your manifold:
    MyManifold = Manifold(4 , 'MyManifold', r'\mathcal{Man}')

- Define your coordinates:
    MyCoordinates.<t,r,th,ph> = MyManifold.chart(r't r th:\theta ph:\phi')

- Define the metric functions (if needed):
    var('M')
    Delta=function('Delta',imag_part_func=0)(r)
    Delta=r^2-2*M*r
    
- Enter null tetrad elements:
    lvec=[1,-(r^2)/Delta,0,0]
    nvec=[Delta/(2*r^2),1/2,0,0]
    mvec=[0,0,(-r/sqrt(2)),(-I*r/sqrt(2))*sin(th)]
    mbarvec=[0,0,(-r/sqrt(2)),(I*r/sqrt(2))*sin(th)]

    *Here, the element ordering is the same as the coordinate ordering.
     (The first element is the t element, the second is the r element, etc.)

- Define an object of the class:
    schw=SageNP(MyManifold,MyCoordinates,lvec,nvec,mvec,mbarvec,'covariant')
  
  *Here, our null-tetrad vectors lvec, nvec, mvec and mbarvec 
   are covariant. Thus we used the keyword 'covariant'.
   If they were contravariant, 
   then we should use the keyword 'contravariant'.

- Once the object is defined, the code calculates the metric 
  and displays it on the screen. It is recommended that you check your metric.

###########

FUNCTIONS:

- All page and equation numbers belong to the reference book.


- test_nulltetrad(): Checks the products of the vectors l*n = -1, m*mbar = 1, all others zero.


- Spin coefficients (Page 75-76, Eq.(7.2))
    
    - calculate_spincoefficients(): Calculates the spin coefficients.
    
    - show_spincoefficients(): Displays the spin coefficients
    
    - All spin coefficients are available under their names:
      kappaNP, kappabarNP, tauNP, taubarNP, sigmaNP, sigmabarNP,
      rhoNP, rhobarNP, piNP, pibarNP, nuNP, nubarNP, muNP, mubarNP,
      lambdaNP, lambdabarNP, epsilonNP, epsilonbarNP, gammaNP, gammabarNP,
      betaNP, betabarNP, alphaNP, alphabarNP

      
- Directional Derivatives (Page 43, Eq.(3.82))

    - DlNP(X): Given X, calculates the D derivative (l direction).

    - DeltanNP(X): Given X, calculates the Delta derivative (n direction)
    
    - deltamNP(X): Given X, calculates the delta derivative (m direction)
    
    - deltambarNP(X): Given X, calculates the deltabar derivative (mbar direction)


- Commutators (Page 77, Eq.(7.6)):
    
    - Deltan_Dl_commNP(X): Given X, calculates the [Delta,D] commutator.
    
    - deltam_Dl_commNP(X): Given X, calculates the [delta,D] commutator.
    
    - deltam_Deltan_commNP(X): Given X, calculates the [delta,Delta] commutator.
    
    - deltambar_deltam_commNP(X): Given X, calculates the [deltabar,delta] commutator.


- Weyl tensor components (Page 38, Eq.(3.59)):
    
    - calculate_Weyl(): Calculates the Weyl tensor components.
    
    - show_Weyl(): Displays the Weyl tensor components.
    
    - All Weyl tensor components are available under their names:
      Psi0NP, Psi1NP, Psi2NP, Psi3NP, Psi4NP


- Ricci components (Page 78, Eq.(7.10-7.15)):
    
    - calculate_Ricci(): Calculates the Ricci tensor components.
    
    - show_Ricci(): Displays the Ricci tensor components.
    
    - All Ricci tensor components are available under their names:
      Phi00NP, Phi01NP, Phi10NP, Phi02NP, Phi20NP, 
      Phi11NP, Phi12NP, Phi21NP, Phi22NP, LambdaNP


- Ricci (Newman-Penrose) Equations (Page 79, Eq.(7.21)):
    
    - All Newman-Penrose equations are defined as 
      0 = -(left hand side)+(right hand side) of the equations
    
    - calculate_NPeq(): Calculates the Newman-Penrose equations
    
    - show_NPeq(): Displays the Newman-Penrose equations
    
    - All Newman-Penrose equations are available under their names
      in the order they are given in the reference:
      NPeq1, NPeq2, NPeq3, NPeq4, NPeq5, NPeq6, NPeq7, NPeq8, NPeq9, NPeq10, 
      NPeq11, NPeq12, NPeq13, NPeq14, NPeq15, NPeq16, NPeq17, NPeq18,


- Bianchi identities (Page 81, Eq.(7.32)):
    
    - All Bianchi identities are defined as 
      0 = -(left hand side)+(right hand side) of the equations
    
    - calculate_Bianchi(): Calculates the Bianchi identities
    
    - show_Bianchi(): Displays the Bianchi identities
    
    - All Bianchi identities are available under their names
      in the order they are given in the reference:
      BI1, BI2, BI3, BI4, BI5, BI6, BI7, BI8, BI9, BI10, BI11


- Petrov invariants I, J, K, L, N (Kramer p.121, 9.6; p.54, 4.19):
  (also check diagram Fig. 9.1 on p. 122)
  
    - calculate_PetrovinvINP(): Calculates the Petrov invariant I
  
    - calculate_PetrovinvJNP(): Calculates the Petrov invariant J
    
    - calculate_PetrovinvKNP(): Calculates the Petrov invariant K
    
    - calculate_PetrovinvLNP(): Calculates the Petrov invariant L
    
    - calculate_PetrovinvNNP(): Calculates the Petrov invariant N

    - All Petrov invariants are available under their names:
      PetrovinvINP, PetrovinvJNP, PetrovinvKNP, PetrovinvLNP, PetrovinvNNP


- Petrov type of the spacetime:
    - Petrov_frominvariants(): Calculates the Petrov type using I, J, K, L, N.
    - Petrov_fromWeyl(): Calculates the Petrov type using the Weyl components


- calculate_allNP(): Runs the following functions:
    - calculate_spincoefficients()
    - calculate_Weyl()
    - calculate_Ricci()
    - calculate_NPeq()
    - calculate_Bianchi()
    - Petrov_frominvariants()
    - Petrov_fromWeyl()

- show_allNP(): Runs the following functions:
    - show_spincoefficients()
    - show_Weyl()
    - show_Ricci()
    - show_NPeq()
    - show_Bianchi()
    - Petrov_frominvariants()
    - Petrov_fromWeyl()

####################################################
SAMPLE WORKSHEET 1:

(SCHWARZSCHILD METRIC WITH COVARIANT NULL-TETRAD VECTORS)

from SageNP import *

########################################
# Define 4-dim. the manifold:
MyManifold = Manifold(4 , 'MyManifold', r'\mathcal{Man}')
########################################
MyCoordinates.<t,r,th,ph> = MyManifold.chart(r't r th:\theta ph:\phi')
########################################
# Define the metric functions
var('M')
Delta=function('Delta',imag_part_func=0)(r)
Delta=r^2-2*M*r
########################################
# Enter null tetrad elements
# These vectors define the Schwarzschild spacetime
lvec=[1,-(r^2)/Delta,0,0]
nvec=[Delta/(2*r^2),1/2,0,0]
mvec=[0,0,(-r/sqrt(2)),(-I*r/sqrt(2))*sin(th)]
mbarvec=[0,0,(-r/sqrt(2)),(I*r/sqrt(2))*sin(th)]
########################################

# Define the object "schw" of the class "SageNP":
schw=SageNP(MyManifold,MyCoordinates,lvec,nvec,mvec,mbarvec,'covariant')

# Test the null tetrad with vector products:
schw.test_nulltetrad()

# Calculate and display the spin coefficients:
schw.calculate_spincoefficients()
schw.show_spincoefficients()

# Display a single spin coefficient:
show(schw.gammaNP.expr())

# Calculate and display the Weyl tensor components:
schw.calculate_Weyl()
schw.show_Weyl()

# Display a single Weyl tensor component:
show(schw.Psi2NP.expr())

# Calculate and display the Ricci tensor components:
schw.calculate_Ricci()
schw.show_Ricci()

# Display a single Ricci tensor component:
show(schw.Phi00NP.expr())

# Calculate and display the Newman-Penrose equations:
schw.calculate_NPeq()
schw.show_NPeq()

# Display a single Newman-Penrose equation:
show(schw.NPeq8.expr())

# Calculate and display the Bianchi identities:
schw.calculate_Bianchi()
schw.show_Bianchi()

# Display a single Bianchi identity:
show(schw.BI7.expr())

# Find the Petrov class using the Petrov invariants:
schw.Petrov_frominvariants()

# Find the Petrov class using the Weyl tensor components:
schw.Petrov_fromWeyl()

# Calculate and display the directional derivatives 
# of the spin coefficient gammaNP:
show(schw.DlNP(schw.gammaNP).expr())
show(schw.DeltanNP(schw.gammaNP).expr())
show(schw.deltamNP(schw.gammaNP).expr())
show(schw.deltambarNP(schw.gammaNP).expr())

# Calculate and display the commutators 
# for the spin coefficient gammaNP:
show(schw.Deltan_Dl_commNP(schw.gammaNP).expr())
show(schw.deltam_Dl_commNP(schw.gammaNP).expr())
show(schw.deltam_Deltan_commNP(schw.gammaNP).expr())
show(schw.deltambar_deltam_commNP(schw.gammaNP).expr())

# Use these functions to calculate and display
# spin coefficients, Weyl components , Ricci components,
# Newman-Penrose equations, Bianchi identities,
# Petrov class using the Petrov invariants,
# Petrov class using the Weyl tensor components
schw.calculate_allNP()
schw.show_allNP()

####################################################
SAMPLE WORKSHEET 2:

(REISSNER-NORDSTROM METRIC WITH CONTRAVARIANT NULL-TETRAD VECTORS)

from SageNP import *

########################################
# Define 4-dim. the manifold:
MyManifold = Manifold(4 , 'MyManifold', r'\mathcal{Man}')
########################################
MyCoordinates.<t,r,th,ph> = MyManifold.chart(r't r th:\theta ph:\phi')
########################################
# Define the metric functions
var('M,Q')
Delta=function('Delta',imag_part_func=0)(r)
Delta=r^2-2*M*r+Q^2
########################################
# Enter null tetrad elements
# These vectors define the Reissner-Nordstrom spacetime
lveccont=[(r^2)/Delta,1,0,0]
nveccont=[1/2,-Delta/(2*r^2),0,0]
mveccont=[0,0,1/(r*sqrt(2)),I*csc(th)/(r*sqrt(2))]
mbarveccont=[0,0,1/(r*sqrt(2)),-I*csc(th)/(r*sqrt(2))]
########################################

# Define the object "reisnor" of the class "SageNP":
reisnor=SageNP(MyManifold,MyCoordinates,lveccont,nveccont,mveccont,mbarveccont,'contravariant')

# Test the null tetrad with vector products:
reisnor.test_nulltetrad()

# Calculate and display the spin coefficients:
reisnor.calculate_spincoefficients()
reisnor.show_spincoefficients()

# Display a single spin coefficient:
show(reisnor.gammaNP.expr())

# Calculate and display the Weyl tensor components:
reisnor.calculate_Weyl()
reisnor.show_Weyl()

# Display a single Weyl tensor component:
show(reisnor.Psi2NP.expr())

# Calculate and display the Ricci tensor components:
reisnor.calculate_Ricci()
reisnor.show_Ricci()

# Display a single Ricci tensor component:
show(reisnor.Phi00NP.expr())

# Calculate and display the Newman-Penrose equations:
reisnor.calculate_NPeq()
reisnor.show_NPeq()

# Display a single Newman-Penrose equation:
show(reisnor.NPeq8.expr())

# Calculate and display the Bianchi identities:
reisnor.calculate_Bianchi()
reisnor.show_Bianchi()

# Display a single Bianchi identity:
show(reisnor.BI7.expr())

# Find the Petrov class using the Petrov invariants:
reisnor.Petrov_frominvariants()

# Find the Petrov class using the Weyl tensor components:
reisnor.Petrov_fromWeyl()
####################################################

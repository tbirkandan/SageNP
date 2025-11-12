from sage.calculus.var import var
from sage.calculus.var import function
from sage.functions.other import conjugate
from sage.symbolic.ring import SR
from sage.symbolic.relation import solve
from sage.rings.integer import Integer
from sage.rings.imaginary_unit import I
from sage.repl.rich_output.pretty_print import show
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from multiprocessing import Pool, cpu_count
from sage.all import *


_sage_const_4 = Integer(4); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_10 = Integer(10); _sage_const_3 = Integer(3); _sage_const_5 = Integer(5); _sage_const_6 = Integer(6); _sage_const_7 = Integer(7); _sage_const_8 = Integer(8); _sage_const_9 = Integer(9); _sage_const_24 = Integer(24); _sage_const_12 = Integer(12); _sage_const_27 = Integer(27); _sage_const_54 = Integer(54)

###########################################################
# NewmanPenrose: Newman-Penrose Calculations for SageMath # 
###########################################################
# Code by:                                                #
# Tolga Birkandan (Corr.: birkandant@itu.edu.tr)          #
# Onur Arman                                              #
# Emir Baysazan                                           #
# Selinay Sude Binici                                     #
# Pelin Ozturk                                            #
# Efe Ozyurek                                             #
# Special thanks to Eric Gourgoulhon                      #
# Based on SageMath (SageManifolds)                       #
###########################################################
r"""
####################################################

The class "NewmanPenrose" introduces functions into SageMath
for some calculations defined in the Newman-Penrose formalism.

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
    from SageNP import NewmanPenrose
    
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
    schw=NewmanPenrose(MyManifold,MyCoordinates,lvec,nvec,mvec,mbarvec,'covariant')
  
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
    
    - The right-hand sides of the commutation relations are calculated.
    
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
    - Petrov_fromWeyl()

- show_allNP(): Runs the following functions:
    - show_spincoefficients()
    - show_Weyl()
    - show_Ricci()
    - show_NPeq()
    - show_Bianchi()

####################################################
SAMPLE WORKSHEET 1:

(SCHWARZSCHILD METRIC WITH COVARIANT NULL-TETRAD VECTORS)

from SageNP import NewmanPenrose

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

# Define the object "schw" of the class "NewmanPenrose":
schw=NewmanPenrose(MyManifold,MyCoordinates,lvec,nvec,mvec,mbarvec,'covariant')

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

from SageNP import NewmanPenrose

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

# Define the object "reisnor" of the class "NewmanPenrose":
reisnor=NewmanPenrose(MyManifold,MyCoordinates,lveccont,nveccont,mveccont,mbarveccont,'contravariant')

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
"""
####################################################
def compute_spin_component(args):
    self, component = args
    if component ==0:
        return "kappaNP", -self.nab(self.lNP)['_ij']*(self.mupNP*self.lupNP)['^ij']
    elif component ==1:    
        return "tauNP", -self.nab(self.lNP)['_ij']*(self.mupNP*self.nupNP)['^ij']
    elif component ==2:    
        return "sigmaNP", -self.nab(self.lNP)['_ij']*(self.mupNP*self.mupNP)['^ij']
    elif component ==3:    
        return "rhoNP", -self.nab(self.lNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij']
    elif component ==4:    
        return "piNP", self.nab(self.nNP)['_ij']*(self.mbarupNP*self.lupNP)['^ij']
    elif component ==5:    
        return "nuNP", self.nab(self.nNP)['_ij']*(self.mbarupNP*self.nupNP)['^ij']
    elif component ==6:    
        return "muNP", self.nab(self.nNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij']
    elif component ==7:    
        return "lambdaNP", self.nab(self.nNP)['_ij']*(self.mbarupNP*self.mbarupNP)['^ij']
    elif component ==8:    
        return "epsilonNP", -(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.lupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.lupNP)['^ij'])
    elif component ==9:    
        return "gammaNP", (_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.nupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.nupNP)['^ij']) 
    elif component ==10:    
        return "betaNP", -(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.mupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij'])
    elif component ==11:    
        return "alphaNP", (_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.mbarupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij']) 
    elif component ==12:    
        return "kappabarNP", -self.nab(self.lNP)['_ij']*(self.mbarupNP*self.lupNP)['^ij']
    elif component ==13:    
        return "taubarNP", -self.nab(self.lNP)['_ij']*(self.mbarupNP*self.nupNP)['^ij']
    elif component ==14:    
        return "sigmabarNP", -self.nab(self.lNP)['_ij']*(self.mbarupNP*self.mbarupNP)['^ij']
    elif component ==15:    
        return "rhobarNP", -self.nab(self.lNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij']
    elif component ==16:    
        return "pibarNP", self.nab(self.nNP)['_ij']*(self.mupNP*self.lupNP)['^ij']
    elif component ==17:    
        return "nubarNP", self.nab(self.nNP)['_ij']*(self.mupNP*self.nupNP)['^ij']
    elif component ==18:    
        return "mubarNP", self.nab(self.nNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij']
    elif component ==19:    
        return "lambdabarNP", self.nab(self.nNP)['_ij']*(self.mupNP*self.mupNP)['^ij']
    elif component ==20:    
        return "epsilonbarNP", -(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.lupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.lupNP)['^ij'])
    elif component ==21:    
        return "gammabarNP", (_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.nupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.nupNP)['^ij'])
    elif component ==22:    
        return "betabarNP", -(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.mbarupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij'])
    elif component ==23:    
        return "alphabarNP", (_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.mupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij'])

def compute_weyl_component(args):
    self, component = args
    if component ==0:
        return "Psi0NP", self.weyl_C_allindicesdown['_{pqrs}']*(self.lupNP*self.mupNP*self.lupNP*self.mupNP)['^{pqrs}']
    elif component ==1: 
        return "Psi1NP", self.weyl_C_allindicesdown['_{pqrs}']*(self.lupNP*self.nupNP*self.lupNP*self.mupNP)['^{pqrs}']
    elif component ==2: 
        return "Psi2NP", self.weyl_C_allindicesdown['_{pqrs}']*(self.lupNP*self.mupNP*self.mbarupNP*self.nupNP)['^{pqrs}']
    elif component ==3: 
        return "Psi3NP", self.weyl_C_allindicesdown['_{pqrs}']*(self.lupNP*self.nupNP*self.mbarupNP*self.nupNP)['^{pqrs}']
    elif component ==4: 
        return "Psi4NP", self.weyl_C_allindicesdown['_{pqrs}']*(self.mbarupNP*self.nupNP*self.mbarupNP*self.nupNP)['^{pqrs}']

def compute_ricci_component(args):
    self, component = args
    if component ==0:
    	return "Phi00NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.lupNP)['^{ab}']
    elif component ==1: 
        return "Phi01NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.mupNP)['^{ab}']
    elif component ==2: 
        return "Phi10NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.mbarupNP)['^{ab}']
    elif component ==3: 
        return "Phi02NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.mupNP*self.mupNP)['^{ab}']
    elif component ==4: 
        return "Phi20NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.mbarupNP*self.mbarupNP)['^{ab}']
    elif component ==5: 
        return "Phi11NP", (_sage_const_1 /_sage_const_4 )*(self.ricci_R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.nupNP)['^{ab}']+self.ricci_R_traceless_allindicesdown['_{ab}']*(self.mupNP*self.mbarupNP)['^{ab}'])
    elif component ==6: 
        return "Phi12NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.nupNP*self.mupNP)['^{ab}']
    elif component ==7: 
        return "Phi21NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.nupNP*self.mbarupNP)['^{ab}']
    elif component ==8: 
        return "Phi22NP", (_sage_const_1 /_sage_const_2 )*self.ricci_R_traceless_allindicesdown['_{ab}']*(self.nupNP*self.nupNP)['^{ab}']
    elif component ==9: 
        return "LambdaNP", (_sage_const_1 /_sage_const_24 )*self.ricci_Ricciscalar

def compute_npeq_component(args):
    self, component = args
    if component ==0:
        return "NPeq1", -self.DlNP(self.rhoNP)+self.deltambarNP(self.kappaNP)+(self.rhoNP**_sage_const_2 +self.sigmaNP*self.sigmabarNP)+(self.epsilonNP+self.epsilonbarNP)*self.rhoNP-self.kappabarNP*self.tauNP-self.kappaNP*(_sage_const_3 *self.alphaNP+self.betabarNP-self.piNP)+self.Phi00NP
    elif component == 1: 
        return "NPeq2", -self.DlNP(self.sigmaNP)+self.deltamNP(self.kappaNP)+self.sigmaNP*(_sage_const_3 *self.epsilonNP-self.epsilonbarNP+self.rhoNP+self.rhobarNP)+self.kappaNP*(self.pibarNP-self.tauNP-_sage_const_3 *self.betaNP-self.alphabarNP)+self.Psi0NP
    elif component == 2: 
        return "NPeq3", -self.DlNP(self.tauNP)+self.DeltanNP(self.kappaNP)+(self.tauNP+self.pibarNP)*self.rhoNP+(self.taubarNP+self.piNP)*self.sigmaNP+(self.epsilonNP-self.epsilonbarNP)*self.tauNP-(_sage_const_3 *self.gammaNP+self.gammabarNP)*self.kappaNP+self.Psi1NP+self.Phi01NP
    elif component == 3: 
        return "NPeq4", -self.DlNP(self.alphaNP)+self.deltambarNP(self.epsilonNP)+(self.rhoNP+self.epsilonbarNP-_sage_const_2 *self.epsilonNP)*self.alphaNP+self.betaNP*self.sigmabarNP-self.betabarNP*self.epsilonNP-self.kappaNP*self.lambdaNP-self.kappabarNP*self.gammaNP+(self.epsilonNP+self.rhoNP)*self.piNP+self.Phi10NP
    elif component == 4: 
        return "NPeq5", -self.DlNP(self.betaNP) + self.deltamNP(self.epsilonNP) +(self.alphaNP + self.piNP)*self.sigmaNP + (self.rhobarNP - self.epsilonbarNP)*self.betaNP - (self.muNP + self.gammaNP)*self.kappaNP - (self.alphabarNP - self.pibarNP)*self.epsilonNP + self.Psi1NP
    elif component == 5: 
        return "NPeq6", -self.DlNP(self.gammaNP) + self.DeltanNP(self.epsilonNP) + (self.tauNP + self.pibarNP)*self.alphaNP + (self.taubarNP+self.piNP)*self.betaNP - (self.epsilonNP + self.epsilonbarNP)*self.gammaNP - (self.gammaNP + self.gammabarNP)*self.epsilonNP + self.tauNP*self.piNP - self.nuNP*self.kappaNP + self.Psi2NP + self.Phi11NP - self.LambdaNP
    elif component == 6: 
        return "NPeq7", -self.DlNP(self.lambdaNP) + self.deltambarNP(self.piNP) + (self.rhoNP*self.lambdaNP + self.sigmabarNP*self.muNP) + self.piNP**_sage_const_2  + (self.alphaNP - self.betabarNP)*self.piNP - self.nuNP*self.kappabarNP - (_sage_const_3 *self.epsilonNP - self.epsilonbarNP)*self.lambdaNP + self.Phi20NP
    elif component == 7: 
        return "NPeq8", -self.DlNP(self.muNP) + self.deltamNP(self.piNP) + (self.rhobarNP*self.muNP + self.sigmaNP*self.lambdaNP) +self.piNP*self.pibarNP - (self.epsilonNP + self.epsilonbarNP)*self.muNP - (self.alphabarNP - self.betaNP)*self.piNP - self.nuNP*self.kappaNP + self.Psi2NP + _sage_const_2 *self.LambdaNP
    elif component == 8: 
        return "NPeq9", -self.DlNP(self.nuNP) + self.DeltanNP(self.piNP) +(self.piNP + self.taubarNP)*self.muNP + (self.pibarNP + self.tauNP)*self.lambdaNP + (self.gammaNP - self.gammabarNP)*self.piNP - (_sage_const_3 *self.epsilonNP + self.epsilonbarNP)*self.nuNP + self.Psi3NP + self.Phi21NP
    elif component == 9: 
        return "NPeq10", -self.DeltanNP(self.lambdaNP) + self.deltambarNP(self.nuNP) - (self.muNP + self.mubarNP)*self.lambdaNP - (_sage_const_3 *self.gammaNP - self.gammabarNP)*self.lambdaNP + (_sage_const_3 *self.alphaNP + self.betabarNP + self.piNP - self.taubarNP)*self.nuNP - self.Psi4NP
    elif component == 10: 
        return "NPeq11", -self.deltamNP(self.rhoNP) + self.deltambarNP(self.sigmaNP) + (self.alphabarNP + self.betaNP)*self.rhoNP - (_sage_const_3 *self.alphaNP - self.betabarNP)*self.sigmaNP + (self.rhoNP - self.rhobarNP)*self.tauNP + (self.muNP - self.mubarNP)*self.kappaNP - self.Psi1NP + self.Phi01NP
    elif component == 11: 
        return "NPeq12", -self.deltamNP(self.alphaNP) + self.deltambarNP(self.betaNP) + (self.muNP*self.rhoNP - self.lambdaNP*self.sigmaNP) + self.alphaNP*self.alphabarNP + self.betaNP*self.betabarNP - _sage_const_2 *self.alphaNP*self.betaNP + (self.rhoNP-self.rhobarNP)*self.gammaNP + (self.muNP-self.mubarNP)*self.epsilonNP - self.Psi2NP + self.Phi11NP + self.LambdaNP
    elif component == 12: 
        return "NPeq13", -self.deltamNP(self.lambdaNP) + self.deltambarNP(self.muNP) + (self.rhoNP - self.rhobarNP)*self.nuNP + (self.muNP - self.mubarNP)*self.piNP + (self.alphaNP + self.betabarNP)*self.muNP + (self.alphabarNP - _sage_const_3 *self.betaNP)*self.lambdaNP - self.Psi3NP + self.Phi21NP
    elif component == 13: 
        return "NPeq14", -self.deltamNP(self.nuNP) + self.DeltanNP(self.muNP) +(self.muNP**_sage_const_2  + self.lambdaNP*self.lambdabarNP) + (self.gammaNP + self.gammabarNP)*self.muNP - self.nubarNP*self.piNP + (self.tauNP - _sage_const_3 *self.betaNP - self.alphabarNP)*self.nuNP + self.Phi22NP
    elif component == 14: 
        return "NPeq15", -self.deltamNP(self.gammaNP) + self.DeltanNP(self.betaNP) + (self.tauNP - self.alphabarNP - self.betaNP)*self.gammaNP + self.muNP*self.tauNP - self.sigmaNP* self.nuNP - self.epsilonNP*self.nubarNP - (self.gammaNP - self.gammabarNP - self.muNP)*self.betaNP + self.alphaNP*self.lambdabarNP + self.Phi12NP
    elif component == 15: 
        return "NPeq16", -self.deltamNP(self.tauNP) + self.DeltanNP(self.sigmaNP) + (self.muNP*self.sigmaNP + self.lambdabarNP*self.rhoNP) + (self.tauNP + self.betaNP - self.alphabarNP)*self.tauNP - (_sage_const_3 *self.gammaNP - self.gammabarNP)*self.sigmaNP - self.kappaNP*self.nubarNP + self.Phi02NP
    elif component == 16: 
        return "NPeq17", -self.DeltanNP(self.rhoNP) + self.deltambarNP(self.tauNP) - (self.rhoNP*self.mubarNP + self.sigmaNP*self.lambdaNP) + (self.betabarNP - self.alphaNP - self.taubarNP)*self.tauNP + (self.gammaNP + self.gammabarNP)*self.rhoNP + self.nuNP*self.kappaNP - self.Psi2NP - _sage_const_2 *self.LambdaNP
    elif component == 17: 
        return "NPeq18", -self.DeltanNP(self.alphaNP) + self.deltambarNP(self.gammaNP) + (self.rhoNP + self.epsilonNP)*self.nuNP - (self.tauNP + self.betaNP)*self.lambdaNP + (self.gammabarNP- self.mubarNP)*self.alphaNP + (self.betabarNP - self.taubarNP)*self.gammaNP - self.Psi3NP

def compute_bianchi_component(args):
    self, component = args
    if component ==0:	
        return "BI1", -(self.deltambarNP(self.Psi0NP)-self.DlNP(self.Psi1NP)+self.DlNP(self.Phi01NP)-self.deltamNP(self.Phi00NP))+(_sage_const_4 *self.alphaNP-self.piNP)*self.Psi0NP-_sage_const_2 *(_sage_const_2 *self.rhoNP+self.epsilonNP)*self.Psi1NP+_sage_const_3 *self.kappaNP*self.Psi2NP+(self.pibarNP-_sage_const_2 *self.alphabarNP-_sage_const_2 *self.betaNP)*self.Phi00NP+_sage_const_2 *(self.epsilonNP+self.rhobarNP)*self.Phi01NP+_sage_const_2 *self.sigmaNP*self.Phi10NP-_sage_const_2 *self.kappaNP*self.Phi11NP-self.kappabarNP*self.Phi02NP
    elif component == 1:
        return "BI2", -self.DeltanNP(self.Psi0NP) + self.deltamNP(self.Psi1NP) - self.DlNP(self.Phi02NP) + self.deltamNP(self.Phi01NP) + (_sage_const_4 *self.gammaNP-self.muNP)*self.Psi0NP - _sage_const_2 *(_sage_const_2 *self.tauNP + self.betaNP)*self.Psi1NP + _sage_const_3 *self.sigmaNP*self.Psi2NP + (_sage_const_2 *self.epsilonNP - _sage_const_2 *self.epsilonbarNP + self.rhobarNP)*self.Phi02NP + _sage_const_2 *(self.pibarNP - self.betaNP)*self.Phi01NP + _sage_const_2 *self.sigmaNP*self.Phi11NP - _sage_const_2 *self.kappaNP*self.Phi12NP - self.lambdabarNP*self.Phi00NP
    elif component == 2:
        return "BI3", -self.deltambarNP(self.Psi3NP) + self.DlNP(self.Psi4NP) - self.deltambarNP(self.Phi21NP) + self.DeltanNP(self.Phi20NP) + (_sage_const_4 *self.epsilonNP - self.rhoNP)*self.Psi4NP - _sage_const_2 *(_sage_const_2 *self.piNP + self.alphaNP)*self.Psi3NP + _sage_const_3 *self.lambdaNP*self.Psi2NP + (_sage_const_2 *self.gammaNP - _sage_const_2 *self.gammabarNP + self.mubarNP)*self.Phi20NP + _sage_const_2 *(self.taubarNP - self.alphaNP)*self.Phi21NP + _sage_const_2 *self.lambdaNP*self.Phi11NP - _sage_const_2 *self.nuNP*self.Phi10NP - self.sigmabarNP*self.Phi22NP
    elif component == 3:
        return "BI4", -self.DeltanNP(self.Psi3NP) + self.deltamNP(self.Psi4NP) - self.deltambarNP(self.Phi22NP) + self.DeltanNP(self.Phi21NP) + (_sage_const_4 *self.betaNP - self.tauNP)*self.Psi4NP - _sage_const_2 *(_sage_const_2 *self.muNP + self.gammaNP)*self.Psi3NP + _sage_const_3 *self.nuNP*self.Psi2NP + (self.taubarNP - _sage_const_2 *self.betabarNP - _sage_const_2 *self.alphaNP)*self.Phi22NP + _sage_const_2 *(self.gammaNP+self.mubarNP)*self.Phi21NP + _sage_const_2 *self.lambdaNP*self.Phi12NP - _sage_const_2 *self.nuNP*self.Phi11NP - self.nubarNP*self.Phi20NP
    elif component == 4:
        return "BI5", -self.DlNP(self.Psi2NP) + self.deltambarNP(self.Psi1NP) - self.DeltanNP(self.Phi00NP) + self.deltambarNP(self.Phi01NP) - _sage_const_2 *self.DlNP(self.LambdaNP) - self.lambdaNP*self.Psi0NP + _sage_const_2 *(self.piNP-self.alphaNP)*self.Psi1NP + _sage_const_3 *self.rhoNP*self.Psi2NP - _sage_const_2 *self.kappaNP*self.Psi3NP + (_sage_const_2 *self.gammaNP + _sage_const_2 *self.gammabarNP - self.mubarNP)*self.Phi00NP - _sage_const_2 *(self.taubarNP + self.alphaNP)*self.Phi01NP - _sage_const_2 *self.tauNP*self.Phi10NP + _sage_const_2 *self.rhoNP*self.Phi11NP + self.sigmabarNP*self.Phi02NP
    elif component == 5:
        return "BI6", -self.DeltanNP(self.Psi2NP) + self.deltamNP(self.Psi3NP) - self.DlNP(self.Phi22NP) + self.deltamNP(self.Phi21NP) - _sage_const_2 * self.DeltanNP(self.LambdaNP) + self.sigmaNP*self.Psi4NP + _sage_const_2 *(self.betaNP-self.tauNP)*self.Psi3NP - _sage_const_3 *self.muNP*self.Psi2NP + _sage_const_2 *self.nuNP*self.Psi1NP + (self.rhobarNP - _sage_const_2 *self.epsilonNP - _sage_const_2 *self.epsilonbarNP)*self.Phi22NP + _sage_const_2 *(self.pibarNP + self.betaNP)*self.Phi21NP + _sage_const_2 *self.piNP*self.Phi12NP - _sage_const_2 *self.muNP*self.Phi11NP - self.lambdabarNP*self.Phi20NP
    elif component == 6:
        return "BI7", -self.DlNP(self.Psi3NP) + self.deltambarNP(self.Psi2NP) + self.DlNP(self.Phi21NP) -self.deltamNP(self.Phi20NP) + _sage_const_2 *self.deltambarNP(self.LambdaNP) - self.kappaNP*self.Psi4NP + _sage_const_2 *(self.rhoNP-self.epsilonNP)*self.Psi3NP + _sage_const_3 *self.piNP*self.Psi2NP - _sage_const_2 *self.lambdaNP*self.Psi1NP + (_sage_const_2 *self.alphabarNP-_sage_const_2 *self.betaNP - self.pibarNP)*self.Phi20NP - _sage_const_2 *(self.rhobarNP-self.epsilonNP)*self.Phi21NP - _sage_const_2 *self.piNP*self.Phi11NP + _sage_const_2 *self.muNP*self.Phi10NP + self.kappabarNP*self.Phi22NP
    elif component == 7:
        return "BI8", -self.DeltanNP(self.Psi1NP) + self.deltamNP(self.Psi2NP) + self.DeltanNP(self.Phi01NP) - self.deltambarNP(self.Phi02NP) + _sage_const_2 *self.deltamNP(self.LambdaNP) + self.nuNP*self.Psi0NP + _sage_const_2 *(self.gammaNP-self.muNP)*self.Psi1NP - _sage_const_3 *self.tauNP*self.Psi2NP + _sage_const_2 *self.sigmaNP*self.Psi3NP + (self.taubarNP - _sage_const_2 *self.betabarNP + _sage_const_2 *self.alphaNP)*self.Phi02NP + _sage_const_2 *(self.mubarNP - self.gammaNP)*self.Phi01NP + _sage_const_2 *self.tauNP*self.Phi11NP-_sage_const_2 *self.rhoNP*self.Phi12NP-self.nubarNP*self.Phi00NP
    elif component == 8:
        return "BI9", -self.DlNP(self.Phi11NP) + self.deltamNP(self.Phi10NP) + self.deltambarNP(self.Phi01NP) - self.DeltanNP(self.Phi00NP) - _sage_const_3 *self.DlNP(self.LambdaNP) + (_sage_const_2 *self.gammaNP - self.muNP + _sage_const_2 *self.gammabarNP - self.mubarNP)*self.Phi00NP + (self.piNP - _sage_const_2 *self.alphaNP - _sage_const_2 *self.taubarNP)*self.Phi01NP + (self.pibarNP - _sage_const_2 *self.alphabarNP - _sage_const_2 *self.tauNP)*self.Phi10NP + _sage_const_2 *(self.rhoNP+self.rhobarNP)*self.Phi11NP + self.sigmabarNP*self.Phi02NP + self.sigmaNP*self.Phi20NP - self.kappabarNP*self.Phi12NP - self.kappaNP*self.Phi21NP
    elif component == 9:
        return "BI10", -self.DlNP(self.Phi12NP) + self.deltamNP(self.Phi11NP) + self.deltambarNP(self.Phi02NP) - self.DeltanNP(self.Phi01NP) - _sage_const_3 *self.deltamNP(self.LambdaNP) + (-_sage_const_2 *self.alphaNP + _sage_const_2 *self.betabarNP + self.piNP - self.taubarNP)*self.Phi02NP + (self.rhobarNP + _sage_const_2 *self.rhoNP - _sage_const_2 *self.epsilonbarNP)*self.Phi12NP + _sage_const_2 *(self.pibarNP - self.tauNP)*self.Phi11NP + (_sage_const_2 *self.gammaNP - _sage_const_2 *self.mubarNP - self.muNP)*self.Phi01NP + self.nubarNP*self.Phi00NP - self.lambdabarNP*self.Phi10NP + self.sigmaNP*self.Phi21NP - self.kappaNP*self.Phi22NP
    elif component == 10:
        return "BI11", -self.DlNP(self.Phi22NP) + self.deltamNP(self.Phi21NP) + self.deltambarNP(self.Phi12NP) - self.DeltanNP(self.Phi11NP) - _sage_const_3 *self.DeltanNP(self.LambdaNP) + (self.rhoNP + self.rhobarNP - _sage_const_2 *self.epsilonNP - _sage_const_2 *self.epsilonbarNP)*self.Phi22NP + (_sage_const_2 *self.betabarNP + _sage_const_2 *self.piNP - self.taubarNP)*self.Phi12NP + (_sage_const_2 *self.betaNP + _sage_const_2 *self.pibarNP - self.tauNP)*self.Phi21NP - _sage_const_2 *(self.muNP+self.mubarNP)*self.Phi11NP + self.nuNP*self.Phi01NP + self.nubarNP*self.Phi10NP - self.lambdabarNP*self.Phi20NP - self.lambdaNP*self.Phi02NP



class NewmanPenrose():
    ####################################################
    def __init__(self,Man,CO,lvecaux,nvecaux,mvecaux,mbarvecaux,vectype):
        self.Man=Man
        self.CO=CO
        self.spincoeffscalculated=_sage_const_0
        self.weylcalculated=_sage_const_0
        self.riccicalculated=_sage_const_0
        self.PetrovinvINPcalculated=_sage_const_0
        self.PetrovinvJNPcalculated=_sage_const_0
        self.PetrovinvKNPcalculated=_sage_const_0
        self.PetrovinvLNPcalculated=_sage_const_0
        self.PetrovinvNNPcalculated=_sage_const_0
        self.petrovtype = None
        self.weyl_C = None
        self.weyl_C_allindicesdown= None
        self.ricci_Ricciscalar= None
        self.ricci_R_traceless_allindicesdown= None
        #########################
        if vectype=='covariant':
            self.lNPvec=lvecaux
            self.nNPvec=nvecaux
            self.mNPvec=mvecaux
            self.mbarNPvec=mbarvecaux
            self.createmetric()
        #########################
        elif vectype=='contravariant':
            self.lNPvecup=lvecaux
            self.nNPvecup=nvecaux
            self.mNPvecup=mvecaux
            self.mbarNPvecup=mbarvecaux

            print("Inverting tetrad...")
            self.lNPvec=vector(SR, len(self.CO[:]))
            self.nNPvec=vector(SR, len(self.CO[:]))
            self.mNPvec=vector(SR, len(self.CO[:]))
            self.mbarNPvec=vector(SR, len(self.CO[:]))
            #########################   
            try:
                var('l1,l2,l3,l4,n1,n2,n3,n4,m1r,m2r,m3r,m4r,m1i,m2i,m3i,m4i')
                self.lNPvec=[l1,l2,l3,l4]
                self.nNPvec=[n1,n2,n3,n4]
                self.mNPvec=[m1r+I*m1i,m2r+I*m2i,m3r+I*m3i,m4r+I*m4i]
                self.mbarNPvec=[m1r-I*m1i,m2r-I*m2i,m3r-I*m3i,m4r-I*m4i]
                
                ll=_sage_const_0 
                nn=_sage_const_0 
                mm=_sage_const_0 
                mbarmbar=_sage_const_0 
            
                lm=_sage_const_0 
                lmbar=_sage_const_0 
                ml=_sage_const_0 
                mbarl=_sage_const_0 
            
                nm=_sage_const_0 
                nmbar=_sage_const_0 
                mn=_sage_const_0 
                mbarn=_sage_const_0 
            
                ln=_sage_const_1 
                nl=_sage_const_1 
                mmbar=-_sage_const_1 
                mbarm=-_sage_const_1 
            
                for i in range(_sage_const_4 ):
                    ll=ll+self.lNPvecup[i]*self.lNPvec[i]
                    nn=nn+self.nNPvecup[i]*self.nNPvec[i]
                    mm=mm+self.mNPvecup[i]*self.mNPvec[i]
                    mbarmbar=mbarmbar+self.mbarNPvecup[i]*self.mbarNPvec[i]
                    lm=lm+self.lNPvecup[i]*self.mNPvec[i]
                    lmbar=lmbar+self.lNPvecup[i]*self.mbarNPvec[i]
                    ml=ml+self.mNPvecup[i]*self.lNPvec[i]
                    mbarl=mbarl+self.mbarNPvecup[i]*self.lNPvec[i]
                    nm=nm+self.nNPvecup[i]*self.mNPvec[i]
                    nmbar=nmbar+self.nNPvecup[i]*self.mbarNPvec[i]
                    mn=mn+self.mNPvecup[i]*self.nNPvec[i]
                    mbarn=mbarn+self.mbarNPvecup[i]*self.nNPvec[i]
                    ln=ln+self.lNPvecup[i]*self.nNPvec[i]
                    nl=nl+self.nNPvecup[i]*self.lNPvec[i]
                    mmbar=mmbar+self.mNPvecup[i]*self.mbarNPvec[i]
                    mbarm=mbarm+self.mbarNPvecup[i]*self.mNPvec[i]
                
                covtetrad=solve([ll==_sage_const_0 ,nn==_sage_const_0 ,mm==_sage_const_0 ,mbarmbar==_sage_const_0 ,lm==_sage_const_0 ,lmbar==_sage_const_0 ,ml==_sage_const_0 ,mbarl==_sage_const_0 ,nm==_sage_const_0 ,nmbar==_sage_const_0 ,mn==_sage_const_0 ,mbarn==_sage_const_0 ,ln==_sage_const_0 ,nl==_sage_const_0 ,mmbar==_sage_const_0 ,mbarm==_sage_const_0 ],[l1,l2,l3,l4,n1,n2,n3,n4,m1r,m1i,m2r,m2i,m3r,m3i,m4r,m4i])
                for i in range(_sage_const_4 ):
                    self.lNPvec[i]=NewmanPenrose.simplify_fullfull(self.lNPvec[i].subs(covtetrad[_sage_const_0 ]))
                    self.nNPvec[i]=NewmanPenrose.simplify_fullfull(self.nNPvec[i].subs(covtetrad[_sage_const_0 ]))
                    self.mNPvec[i]=NewmanPenrose.simplify_fullfull(self.mNPvec[i].subs(covtetrad[_sage_const_0 ]))
                    self.mbarNPvec[i]=NewmanPenrose.simplify_fullfull(self.mbarNPvec[i].subs(covtetrad[_sage_const_0 ]))
                self.createmetric()
            #########################
            except:
                print("Inverting tetrad with another method (this can take a while)...")
                self.inversemetric=matrix(SR, len(self.CO[:]), len(self.CO[:]))
                metricforinversion=matrix(SR, len(self.CO[:]), len(self.CO[:]))
                #print("Calculating inverse metric...")
                for i in range(len(self.CO[:])):
                    for j in range(len(self.CO[:])):
                        self.inversemetric[i,j]=NewmanPenrose.simplify_fullfull((-self.lNPvecup[i]*self.nNPvecup[j]-self.nNPvecup[i]*self.lNPvecup[j]
                                            +self.mNPvecup[i]*self.mbarNPvecup[j]+self.mbarNPvecup[i]*self.mNPvecup[j]))
                metricforinversion=self.inversemetric.inverse()
                for i in range(len(self.CO[:])):
                    for j in range(len(self.CO[:])):
                        metricforinversion[i,j]=NewmanPenrose.simplify_fullfull(metricforinversion[i,j])
                        
                for i in range(len(CO[:])):
                    self.lNPvec[i]=0
                    self.nNPvec[i]=0
                    self.mNPvec[i]=0
                    self.mbarNPvec[i]=0
                    for j in range(len(self.CO[:])):
                        self.lNPvec[i]=(self.lNPvec[i]+metricforinversion[i,j]*self.lNPvecup[j])
                        self.nNPvec[i]=(self.nNPvec[i]+metricforinversion[i,j]*self.nNPvecup[j])
                        self.mNPvec[i]=(self.mNPvec[i]+metricforinversion[i,j]*self.mNPvecup[j])
                        self.mbarNPvec[i]=(self.mbarNPvec[i]+metricforinversion[i,j]*self.mbarNPvecup[j])
                self.createmetric()
        #########################                   
        else:
            print(r"Please check the type of the null tetrad ('covariant' or 'contravariant') and try again.")

    ####################################################
    # Let us define our own simplification routine:
    def simplify_fullfull(x):
        #This routine is used at some points in the functions 
        return x.expand().canonicalize_radical().simplify_rectform().simplify_full()
    
    ####################################################    
    def createmetric(self):
        self.lNP=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.lNPvec)
        self.nNP=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.nNPvec)
        self.mNP=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.mNPvec)
        self.mbarNP=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.mbarNPvec)
        self.g = self.Man.lorentzian_metric('g')

        #########################
        #print("Calculating metric...")
        for i in range(len(self.CO[:])):
            for j in range(len(self.CO[:])):
                self.g[i,j]=NewmanPenrose.simplify_fullfull(-self.lNP[i].expr()*self.nNP[j].expr()
                                                     -self.nNP[i].expr()*self.lNP[j].expr()
                                                     +self.mNP[i].expr()*self.mbarNP[j].expr()
                                                     +self.mbarNP[i].expr()*self.mNP[j].expr())
        #########################
        # Display the metric
        show('The metric:')
        show(self.g.display())
        #########################
        # Calculate the inverse metric:
        try:
            self.lNPup=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.lNPvecup)
            self.nNPup=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.nNPvecup)
            self.mNPup=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.mNPvecup)
            self.mbarNPup=self.Man.tensor_field(_sage_const_0 ,_sage_const_1 ,self.mbarNPvecup)
            self.inversegmat=[[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ]]
            self.inverseg=self.Man.tensor_field(_sage_const_0 ,_sage_const_2 ,self.inversegmat)
            for i in range(len(self.CO[:])):
                for j in range(len(self.CO[:])):
                    self.inverseg[i,j]=simplify_fullfull((-self.lNPup[i].expr()*self.nNPup[j].expr()
                                                          -self.nNPup[i].expr()*self.lNPup[j].expr()
                                                          +self.mNPup[i].expr()*self.mbarNPup[j].expr()
                                                          +self.mbarNPup[i].expr()*self.mNPup[j].expr()))    
        #########################
        except:
            self.inverseg = self.g.inverse()

        #########################
        try:
            self.lupNP=self.Man.tensor_field(_sage_const_1 ,_sage_const_0 ,self.lNPvecup)
            self.nupNP=self.Man.tensor_field(_sage_const_1 ,_sage_const_0 ,self.nNPvecup)
            self.mupNP=self.Man.tensor_field(_sage_const_1 ,_sage_const_0 ,self.mNPvecup)
            self.mbarupNP=self.Man.tensor_field(_sage_const_1 ,_sage_const_0 ,self.mbarNPvecup)
        #########################
        except:
            self.lupNP=self.lNP.up(self.g)
            self.nupNP=self.nNP.up(self.g)
            self.mupNP=self.mNP.up(self.g)
            self.mbarupNP=self.mbarNP.up(self.g)
        #########################
        #########################  
        self.nab = self.g.connection()
        #########################
    
    ####################################################
    def test_nulltetrad(self):
        # All must be zero
        print("Testing null tetrad...")
        testvectorNP=vector(SR,_sage_const_10 )
        testvectorNP[_sage_const_0 ]=NewmanPenrose.simplify_fullfull((self.lupNP['^a']*self.mNP['_a']).expr())
        testvectorNP[_sage_const_1 ]=NewmanPenrose.simplify_fullfull((self.lupNP['^a']*self.mbarNP['_a']).expr())
        testvectorNP[_sage_const_2 ]=NewmanPenrose.simplify_fullfull((self.nupNP['^a']*self.mNP['_a']).expr())
        testvectorNP[_sage_const_3 ]=NewmanPenrose.simplify_fullfull((self.nupNP['^a']*self.mbarNP['_a']).expr())
        testvectorNP[_sage_const_4 ]=NewmanPenrose.simplify_fullfull((self.lupNP['^a']*self.lNP['_a']).expr())
        testvectorNP[_sage_const_5 ]=NewmanPenrose.simplify_fullfull((self.nupNP['^a']*self.nNP['_a']).expr())
        testvectorNP[_sage_const_6 ]=NewmanPenrose.simplify_fullfull((self.mupNP['^a']*self.mNP['_a']).expr())
        testvectorNP[_sage_const_7 ]=NewmanPenrose.simplify_fullfull((self.mbarupNP['^a']*self.mbarNP['_a']).expr())
        testvectorNP[_sage_const_8 ]=NewmanPenrose.simplify_fullfull((self.lupNP['^a']*self.nNP['_a']).expr()+_sage_const_1 )
        testvectorNP[_sage_const_9 ]=NewmanPenrose.simplify_fullfull((self.mupNP['^a']*self.mbarNP['_a']).expr()-_sage_const_1 )
        if sum(testvectorNP)==_sage_const_0 :
            print("PASSED")
        else:
            print("FAILED! Please check your null tetrad elements!")
            for i in range(len(testvectorNP)):
                if not(testvectorNP[i]==_sage_const_0 ):
                    print("Failed test number: ",i)

    ####################################################
    # Spin coefficients:
    # Page 75-76, Eq.(7.2)    
    
    def calculate_spincoefficients(self):
        #########################
        self.kappaNP=self.Man.scalar_field(function('kappaNP')(*self.CO));
        self.kappabarNP=self.Man.scalar_field(function('kappabarNP')(*self.CO))
        self.tauNP=self.Man.scalar_field(function('tauNP')(*self.CO));
        self.taubarNP=self.Man.scalar_field(function('taubarNP')(*self.CO))
        self.sigmaNP=self.Man.scalar_field(function('sigmaNP')(*self.CO));
        self.sigmabarNP=self.Man.scalar_field(function('sigmabarNP')(*self.CO))
        self.rhoNP=self.Man.scalar_field(function('rhoNP')(*self.CO));
        self.rhobarNP=self.Man.scalar_field(function('rhobarNP')(*self.CO))
        self.piNP=self.Man.scalar_field(function('piNP')(*self.CO));
        self.pibarNP=self.Man.scalar_field(function('pibarNP')(*self.CO))
        self.nuNP=self.Man.scalar_field(function('nuNP')(*self.CO));
        self.nubarNP=self.Man.scalar_field(function('nubarNP')(*self.CO))
        self.muNP=self.Man.scalar_field(function('muNP')(*self.CO));
        self.mubarNP=self.Man.scalar_field(function('mubarNP')(*self.CO))
        self.lambdaNP=self.Man.scalar_field(function('lambdaNP')(*self.CO));
        self.lambdabarNP=self.Man.scalar_field(function('lambdabarNP')(*self.CO))
        self.epsilonNP=self.Man.scalar_field(function('epsilonNP')(*self.CO));
        self.epsilonbarNP=self.Man.scalar_field(function('epsilonbarNP')(*self.CO))
        self.gammaNP=self.Man.scalar_field(function('gammaNP')(*self.CO));
        self.gammabarNP=self.Man.scalar_field(function('gammabarNP')(*self.CO))
        self.betaNP=self.Man.scalar_field(function('betaNP')(*self.CO));
        self.betabarNP=self.Man.scalar_field(function('betabarNP')(*self.CO))
        self.alphaNP=self.Man.scalar_field(function('alphaNP')(*self.CO));
        self.alphabarNP=self.Man.scalar_field(function('alphabarNP')(*self.CO))
        #########################
        print("Calculating spin coefficients...")
        self.spincoeffscalculated=_sage_const_1 
        #########################
        self.kappaNP=-self.nab(self.lNP)['_ij']*(self.mupNP*self.lupNP)['^ij']
        self.tauNP=-self.nab(self.lNP)['_ij']*(self.mupNP*self.nupNP)['^ij']
        self.sigmaNP=-self.nab(self.lNP)['_ij']*(self.mupNP*self.mupNP)['^ij']
        self.rhoNP=-self.nab(self.lNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij']
        self.piNP=self.nab(self.nNP)['_ij']*(self.mbarupNP*self.lupNP)['^ij']
        self.nuNP=self.nab(self.nNP)['_ij']*(self.mbarupNP*self.nupNP)['^ij']
        self.muNP=self.nab(self.nNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij']
        self.lambdaNP=self.nab(self.nNP)['_ij']*(self.mbarupNP*self.mbarupNP)['^ij']
        self.epsilonNP=-(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.lupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.lupNP)['^ij'])
        self.gammaNP=(_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.nupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.nupNP)['^ij']) 
        self.betaNP=-(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.mupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij'])
        self.alphaNP=(_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.mbarupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij']) 
        self.kappabarNP=-self.nab(self.lNP)['_ij']*(self.mbarupNP*self.lupNP)['^ij']
        self.taubarNP=-self.nab(self.lNP)['_ij']*(self.mbarupNP*self.nupNP)['^ij']
        self.sigmabarNP=-self.nab(self.lNP)['_ij']*(self.mbarupNP*self.mbarupNP)['^ij']
        self.rhobarNP=-self.nab(self.lNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij']
        self.pibarNP=self.nab(self.nNP)['_ij']*(self.mupNP*self.lupNP)['^ij']
        self.nubarNP=self.nab(self.nNP)['_ij']*(self.mupNP*self.nupNP)['^ij']
        self.mubarNP=self.nab(self.nNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij']
        self.lambdabarNP=self.nab(self.nNP)['_ij']*(self.mupNP*self.mupNP)['^ij']
        self.epsilonbarNP=-(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.lupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.lupNP)['^ij'])
        self.gammabarNP=(_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.nupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.nupNP)['^ij'])
        self.betabarNP=-(_sage_const_1 /_sage_const_2 )*(self.nab(self.lNP)['_ij']*(self.nupNP*self.mbarupNP)['^ij']-self.nab(self.mbarNP)['_ij']*(self.mupNP*self.mbarupNP)['^ij'])
        self.alphabarNP=(_sage_const_1 /_sage_const_2 )*(self.nab(self.nNP)['_ij']*(self.lupNP*self.mupNP)['^ij']-self.nab(self.mNP)['_ij']*(self.mbarupNP*self.mupNP)['^ij'])



    #########################
    def show_spincoefficients(self):
        if self.spincoeffscalculated==_sage_const_1 :
            show("kappaNP=",self.kappaNP.expr())
            show("tauNP=",self.tauNP.expr())
            show("sigmaNP=",self.sigmaNP.expr())
            show("rhoNP=",self.rhoNP.expr())
            show("piNP=",self.piNP.expr())
            show("nuNP=",self.nuNP.expr())
            show("muNP=",self.muNP.expr())
            show("lambdaNP=",self.lambdaNP.expr())
            show("epsilonNP=",self.epsilonNP.expr())
            show("gammaNP=",self.gammaNP.expr())
            show("betaNP=",self.betaNP.expr())
            show("alphaNP=",self.alphaNP.expr())
        else:
            self.calculate_spincoefficients()
            self.show_spincoefficients()
   
    ####################################################
    # Directional Derivatives
    # Page 43, Eq.(3.82)
    #########################
    def DlNP(self,X):
        if X==_sage_const_0 :
            return _sage_const_0 
        else:
            return self.lupNP['^i']*self.nab(X)['_i']
    #########################
    def DeltanNP(self,X):
        if X==_sage_const_0 :
            return _sage_const_0 
        else:
            return self.nupNP['^i']*self.nab(X)['_i']
    #########################
    def deltamNP(self,X):
        if X==_sage_const_0 :
            return _sage_const_0 
        else:
            return self.mupNP['^i']*self.nab(X)['_i']
    #########################
    def deltambarNP(self,X):
        if X==_sage_const_0 :
            return _sage_const_0 
        else:
            return self.mbarupNP['^i']*self.nab(X)['_i']
  
    ####################################################
    # Commutators
    # Page 77, Eq.(7.6)
    # These define the right hand sides of the equations (7.6)
    #########################
    def Deltan_Dl_commNP(self,X):
        if self.spincoeffscalculated==_sage_const_0:
            self.calculate_spincoefficients()    
        comresult=(self.gammaNP+self.gammabarNP)*self.DlNP(X)+(self.epsilonNP+self.epsilonbarNP)*self.DeltanNP(X)-(self.tauNP+self.pibarNP)*self.deltambarNP(X)-(self.taubarNP+self.piNP)*self.deltamNP(X)
        return comresult
    #########################
    def deltam_Dl_commNP(self,X):
        if self.spincoeffscalculated==_sage_const_0:
            self.calculate_spincoefficients()
        comresult=(self.alphabarNP+self.betaNP-self.pibarNP)*self.DlNP(X)+self.kappaNP*self.DeltanNP(X)-self.sigmaNP*self.deltambarNP(X)-(self.rhobarNP+self.epsilonNP-self.epsilonbarNP)*self.deltamNP(X)
        return comresult
    #########################
    def deltam_Deltan_commNP(self,X):
        if self.spincoeffscalculated==_sage_const_0:
            self.calculate_spincoefficients()
        comresult=-self.nubarNP*self.DlNP(X)+(self.tauNP-self.alphabarNP-self.betaNP)*self.DeltanNP(X)+self.lambdabarNP*self.deltambarNP(X)+(self.muNP-self.gammaNP+self.gammabarNP)*self.deltamNP(X)
        return comresult
    #########################
    def deltambar_deltam_commNP(self,X):
        if self.spincoeffscalculated==_sage_const_0:
            self.calculate_spincoefficients()
        comresult=(self.mubarNP-self.muNP)*self.DlNP(X)+(self.rhobarNP-self.rhoNP)*self.DeltanNP(X)-(self.alphabarNP-self.betaNP)*self.deltambarNP(X)-(self.betabarNP-self.alphaNP)*self.deltamNP(X)
        return comresult
    
    ####################################################
    # Weyl tensor components:
    # Page 38, Eq.(3.59)
    def calculate_Weyl(self):
        self.Psi0NP=self.Man.scalar_field(function('Psi0NP')(*self.CO))
        self.Psi1NP=self.Man.scalar_field(function('Psi1NP')(*self.CO))
        self.Psi2NP=self.Man.scalar_field(function('Psi2NP')(*self.CO)) 
        self.Psi3NP=self.Man.scalar_field(function('Psi3NP')(*self.CO))
        self.Psi4NP=self.Man.scalar_field(function('Psi4NP')(*self.CO))
        #########################
        print("Calculating Weyl components...")
        self.weylcalculated=_sage_const_1 
        C = self.g.weyl()
        C_allindicesdown=C.down(self.g)
        self.Psi0NP=C_allindicesdown['_{pqrs}']*(self.lupNP*self.mupNP*self.lupNP*self.mupNP)['^{pqrs}']
        self.Psi1NP=C_allindicesdown['_{pqrs}']*(self.lupNP*self.nupNP*self.lupNP*self.mupNP)['^{pqrs}']
        self.Psi2NP=C_allindicesdown['_{pqrs}']*(self.lupNP*self.mupNP*self.mbarupNP*self.nupNP)['^{pqrs}']
        self.Psi3NP=C_allindicesdown['_{pqrs}']*(self.lupNP*self.nupNP*self.mbarupNP*self.nupNP)['^{pqrs}']
        self.Psi4NP=C_allindicesdown['_{pqrs}']*(self.mbarupNP*self.nupNP*self.mbarupNP*self.nupNP)['^{pqrs}']

            
    #########################    
    def show_Weyl(self):
        if self.weylcalculated==_sage_const_1 :
            show("Psi0NP=",self.Psi0NP.expr())
            show("Psi1NP=",self.Psi1NP.expr())
            show("Psi2NP=",self.Psi2NP.expr())
            show("Psi3NP=",self.Psi3NP.expr())
            show("Psi4NP=",self.Psi4NP.expr())
        else:
            self.calculate_Weyl()
            self.show_Weyl()
    
    ####################################################
    # Ricci components:
    # Page 78, Eq.(7.10-7.15)     
    def calculate_Ricci(self):
        self.Phi00NP=self.Man.scalar_field(function('Phi00NP')(*self.CO))
        self.Phi01NP=self.Man.scalar_field(function('Phi01NP')(*self.CO))
        self.Phi10NP=self.Man.scalar_field(function('Phi10NP')(*self.CO))
        self.Phi02NP=self.Man.scalar_field(function('Phi02NP')(*self.CO))
        self.Phi20NP=self.Man.scalar_field(function('Phi20NP')(*self.CO))
        self.Phi11NP=self.Man.scalar_field(function('Phi11NP')(*self.CO))
        self.Phi12NP=self.Man.scalar_field(function('Phi12NP')(*self.CO))
        self.Phi21NP=self.Man.scalar_field(function('Phi21NP')(*self.CO))
        self.Phi22NP=self.Man.scalar_field(function('Phi22NP')(*self.CO))
        self.LambdaNP=self.Man.scalar_field(function('LambdaNP')(*self.CO))
        #########################
        print("Calculating Ricci components...")
        self.riccicalculated=_sage_const_1 
        Ricciscalar=self.g.ricci_scalar()
        R_traceless_allindicesdown=self.g.ricci()['_{ab}']-(self.g.ricci_scalar()/_sage_const_4 )*self.g['_{ab}']
        self.Phi00NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.lupNP)['^{ab}']
        self.Phi01NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.mupNP)['^{ab}']
        self.Phi10NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.mbarupNP)['^{ab}']
        self.Phi02NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.mupNP*self.mupNP)['^{ab}']
        self.Phi20NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.mbarupNP*self.mbarupNP)['^{ab}']
        self.Phi11NP=(_sage_const_1 /_sage_const_4 )*(R_traceless_allindicesdown['_{ab}']*(self.lupNP*self.nupNP)['^{ab}']+R_traceless_allindicesdown['_{ab}']*(self.mupNP*self.mbarupNP)['^{ab}'])
        self.Phi12NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.nupNP*self.mupNP)['^{ab}']
        self.Phi21NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.nupNP*self.mbarupNP)['^{ab}']
        self.Phi22NP=(_sage_const_1 /_sage_const_2 )*R_traceless_allindicesdown['_{ab}']*(self.nupNP*self.nupNP)['^{ab}']
        self.LambdaNP=(_sage_const_1 /_sage_const_24 )*Ricciscalar # Stephani uses this as R and it is compatible with Lambda

            
    #########################    
    def show_Ricci(self):
        if self.riccicalculated==_sage_const_1 :
            show("Phi00NP=",self.Phi00NP.expr())
            show("Phi01NP=",self.Phi01NP.expr())
            show("Phi10NP=",self.Phi10NP.expr())
            show("Phi02NP=",self.Phi02NP.expr())
            show("Phi20NP=",self.Phi20NP.expr())
            show("Phi11NP=",self.Phi11NP.expr())
            show("Phi12NP=",self.Phi12NP.expr())
            show("Phi21NP=",self.Phi21NP.expr())
            show("Phi22NP=",self.Phi22NP.expr())
            show("LambdaNP=",self.LambdaNP.expr())
        else:
            self.calculate_Ricci()
            self.show_Ricci()

    ####################################################
    # Ricci (Newman-Penrose) Equations
    # Page 79, Eq.(7.21)
    def calculate_NPeq(self):
        if self.spincoeffscalculated==_sage_const_0:
            self.calculate_spincoefficients()
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        if self.riccicalculated==_sage_const_0:
            self.calculate_Ricci()
        #########################
        print("Calculating NP equations...")
        
        self.NPeq1=-self.DlNP(self.rhoNP)+self.deltambarNP(self.kappaNP)+(self.rhoNP**_sage_const_2 +self.sigmaNP*self.sigmabarNP)+(self.epsilonNP+self.epsilonbarNP)*self.rhoNP-self.kappabarNP*self.tauNP-self.kappaNP*(_sage_const_3 *self.alphaNP+self.betabarNP-self.piNP)+self.Phi00NP
        
        self.NPeq2=-self.DlNP(self.sigmaNP)+self.deltamNP(self.kappaNP)+self.sigmaNP*(_sage_const_3 *self.epsilonNP-self.epsilonbarNP+self.rhoNP+self.rhobarNP)+self.kappaNP*(self.pibarNP-self.tauNP-_sage_const_3 *self.betaNP-self.alphabarNP)+self.Psi0NP
        
        self.NPeq3=-self.DlNP(self.tauNP)+self.DeltanNP(self.kappaNP)+(self.tauNP+self.pibarNP)*self.rhoNP+(self.taubarNP+self.piNP)*self.sigmaNP+(self.epsilonNP-self.epsilonbarNP)*self.tauNP-(_sage_const_3 *self.gammaNP+self.gammabarNP)*self.kappaNP+self.Psi1NP+self.Phi01NP
        
        self.NPeq4=-self.DlNP(self.alphaNP)+self.deltambarNP(self.epsilonNP)+(self.rhoNP+self.epsilonbarNP-_sage_const_2 *self.epsilonNP)*self.alphaNP+self.betaNP*self.sigmabarNP-self.betabarNP*self.epsilonNP-self.kappaNP*self.lambdaNP-self.kappabarNP*self.gammaNP+(self.epsilonNP+self.rhoNP)*self.piNP+self.Phi10NP
        
        self.NPeq5=-self.DlNP(self.betaNP) + self.deltamNP(self.epsilonNP) +(self.alphaNP + self.piNP)*self.sigmaNP + (self.rhobarNP - self.epsilonbarNP)*self.betaNP - (self.muNP + self.gammaNP)*self.kappaNP - (self.alphabarNP - self.pibarNP)*self.epsilonNP + self.Psi1NP
        
        self.NPeq6=-self.DlNP(self.gammaNP) + self.DeltanNP(self.epsilonNP) + (self.tauNP + self.pibarNP)*self.alphaNP + (self.taubarNP+self.piNP)*self.betaNP - (self.epsilonNP + self.epsilonbarNP)*self.gammaNP - (self.gammaNP + self.gammabarNP)*self.epsilonNP + self.tauNP*self.piNP - self.nuNP*self.kappaNP + self.Psi2NP + self.Phi11NP - self.LambdaNP
        
        self.NPeq7=-self.DlNP(self.lambdaNP) + self.deltambarNP(self.piNP) + (self.rhoNP*self.lambdaNP + self.sigmabarNP*self.muNP) + self.piNP**_sage_const_2  + (self.alphaNP - self.betabarNP)*self.piNP - self.nuNP*self.kappabarNP - (_sage_const_3 *self.epsilonNP - self.epsilonbarNP)*self.lambdaNP + self.Phi20NP
        
        self.NPeq8=-self.DlNP(self.muNP) + self.deltamNP(self.piNP) + (self.rhobarNP*self.muNP + self.sigmaNP*self.lambdaNP) +self.piNP*self.pibarNP - (self.epsilonNP + self.epsilonbarNP)*self.muNP - (self.alphabarNP - self.betaNP)*self.piNP - self.nuNP*self.kappaNP + self.Psi2NP + _sage_const_2 *self.LambdaNP
        
        self.NPeq9=-self.DlNP(self.nuNP) + self.DeltanNP(self.piNP) +(self.piNP + self.taubarNP)*self.muNP + (self.pibarNP + self.tauNP)*self.lambdaNP + (self.gammaNP - self.gammabarNP)*self.piNP - (_sage_const_3 *self.epsilonNP + self.epsilonbarNP)*self.nuNP + self.Psi3NP + self.Phi21NP
        
        self.NPeq10=-self.DeltanNP(self.lambdaNP) + self.deltambarNP(self.nuNP) - (self.muNP + self.mubarNP)*self.lambdaNP - (_sage_const_3 *self.gammaNP - self.gammabarNP)*self.lambdaNP + (_sage_const_3 *self.alphaNP + self.betabarNP + self.piNP - self.taubarNP)*self.nuNP - self.Psi4NP
        
        self.NPeq11=-self.deltamNP(self.rhoNP) + self.deltambarNP(self.sigmaNP) + (self.alphabarNP + self.betaNP)*self.rhoNP - (_sage_const_3 *self.alphaNP - self.betabarNP)*self.sigmaNP + (self.rhoNP - self.rhobarNP)*self.tauNP + (self.muNP - self.mubarNP)*self.kappaNP - self.Psi1NP + self.Phi01NP
        
        self.NPeq12=-self.deltamNP(self.alphaNP) + self.deltambarNP(self.betaNP) + (self.muNP*self.rhoNP - self.lambdaNP*self.sigmaNP) + self.alphaNP*self.alphabarNP + self.betaNP*self.betabarNP - _sage_const_2 *self.alphaNP*self.betaNP + (self.rhoNP-self.rhobarNP)*self.gammaNP + (self.muNP-self.mubarNP)*self.epsilonNP - self.Psi2NP + self.Phi11NP + self.LambdaNP
        
        self.NPeq13=-self.deltamNP(self.lambdaNP) + self.deltambarNP(self.muNP) + (self.rhoNP - self.rhobarNP)*self.nuNP + (self.muNP - self.mubarNP)*self.piNP + (self.alphaNP + self.betabarNP)*self.muNP + (self.alphabarNP - _sage_const_3 *self.betaNP)*self.lambdaNP - self.Psi3NP + self.Phi21NP
        
        self.NPeq14=-self.deltamNP(self.nuNP) + self.DeltanNP(self.muNP) +(self.muNP**_sage_const_2  + self.lambdaNP*self.lambdabarNP) + (self.gammaNP + self.gammabarNP)*self.muNP - self.nubarNP*self.piNP + (self.tauNP - _sage_const_3 *self.betaNP - self.alphabarNP)*self.nuNP + self.Phi22NP
        
        self.NPeq15=-self.deltamNP(self.gammaNP) + self.DeltanNP(self.betaNP) + (self.tauNP - self.alphabarNP - self.betaNP)*self.gammaNP + self.muNP*self.tauNP - self.sigmaNP* self.nuNP - self.epsilonNP*self.nubarNP - (self.gammaNP - self.gammabarNP - self.muNP)*self.betaNP + self.alphaNP*self.lambdabarNP + self.Phi12NP
        
        self.NPeq16=-self.deltamNP(self.tauNP) + self.DeltanNP(self.sigmaNP) + (self.muNP*self.sigmaNP + self.lambdabarNP*self.rhoNP) + (self.tauNP + self.betaNP - self.alphabarNP)*self.tauNP - (_sage_const_3 *self.gammaNP - self.gammabarNP)*self.sigmaNP - self.kappaNP*self.nubarNP + self.Phi02NP
        
        self.NPeq17=-self.DeltanNP(self.rhoNP) + self.deltambarNP(self.tauNP) - (self.rhoNP*self.mubarNP + self.sigmaNP*self.lambdaNP) + (self.betabarNP - self.alphaNP - self.taubarNP)*self.tauNP + (self.gammaNP + self.gammabarNP)*self.rhoNP + self.nuNP*self.kappaNP - self.Psi2NP - _sage_const_2 *self.LambdaNP
        
        self.NPeq18=-self.DeltanNP(self.alphaNP) + self.deltambarNP(self.gammaNP) + (self.rhoNP + self.epsilonNP)*self.nuNP - (self.tauNP + self.betaNP)*self.lambdaNP + (self.gammabarNP- self.mubarNP)*self.alphaNP + (self.betabarNP - self.taubarNP)*self.gammaNP - self.Psi3NP

    
    #########################
    def show_NPeq(self):
        show("NPeq1=",self.NPeq1.expr())
        show("NPeq2=",self.NPeq2.expr())
        show("NPeq3=",self.NPeq3.expr())
        show("NPeq4=",self.NPeq4.expr())
        show("NPeq5=",self.NPeq5.expr())
        show("NPeq6=",self.NPeq6.expr())
        show("NPeq7=",self.NPeq7.expr())
        show("NPeq8=",self.NPeq8.expr())
        show("NPeq9=",self.NPeq9.expr())
        show("NPeq10=",self.NPeq10.expr())
        show("NPeq11=",self.NPeq11.expr())
        show("NPeq12=",self.NPeq12.expr())
        show("NPeq13=",self.NPeq13.expr())
        show("NPeq14=",self.NPeq14.expr())
        show("NPeq15=",self.NPeq15.expr())
        show("NPeq16=",self.NPeq16.expr())
        show("NPeq17=",self.NPeq17.expr())
        show("NPeq18=",self.NPeq18.expr())
    
    ####################################################
    # Bianchi identities
    # Page 81, Eq.(7.32)
    def calculate_Bianchi(self):
        if self.spincoeffscalculated==_sage_const_0:
            self.calculate_spincoefficients()
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        if self.riccicalculated==_sage_const_0:
            self.calculate_Ricci()
        #########################
        print("Calculating Bianchi identities...")
        self.BI1=-(self.deltambarNP(self.Psi0NP)-self.DlNP(self.Psi1NP)+self.DlNP(self.Phi01NP)-self.deltamNP(self.Phi00NP))+(_sage_const_4 *self.alphaNP-self.piNP)*self.Psi0NP-_sage_const_2 *(_sage_const_2 *self.rhoNP+self.epsilonNP)*self.Psi1NP+_sage_const_3 *self.kappaNP*self.Psi2NP+(self.pibarNP-_sage_const_2 *self.alphabarNP-_sage_const_2 *self.betaNP)*self.Phi00NP+_sage_const_2 *(self.epsilonNP+self.rhobarNP)*self.Phi01NP+_sage_const_2 *self.sigmaNP*self.Phi10NP-_sage_const_2 *self.kappaNP*self.Phi11NP-self.kappabarNP*self.Phi02NP
        
        self.BI2=-self.DeltanNP(self.Psi0NP) + self.deltamNP(self.Psi1NP) - self.DlNP(self.Phi02NP) + self.deltamNP(self.Phi01NP) + (_sage_const_4 *self.gammaNP-self.muNP)*self.Psi0NP - _sage_const_2 *(_sage_const_2 *self.tauNP + self.betaNP)*self.Psi1NP + _sage_const_3 *self.sigmaNP*self.Psi2NP + (_sage_const_2 *self.epsilonNP - _sage_const_2 *self.epsilonbarNP + self.rhobarNP)*self.Phi02NP + _sage_const_2 *(self.pibarNP - self.betaNP)*self.Phi01NP + _sage_const_2 *self.sigmaNP*self.Phi11NP - _sage_const_2 *self.kappaNP*self.Phi12NP - self.lambdabarNP*self.Phi00NP
        
        self.BI3=-self.deltambarNP(self.Psi3NP) + self.DlNP(self.Psi4NP) - self.deltambarNP(self.Phi21NP) + self.DeltanNP(self.Phi20NP) + (_sage_const_4 *self.epsilonNP - self.rhoNP)*self.Psi4NP - _sage_const_2 *(_sage_const_2 *self.piNP + self.alphaNP)*self.Psi3NP + _sage_const_3 *self.lambdaNP*self.Psi2NP + (_sage_const_2 *self.gammaNP - _sage_const_2 *self.gammabarNP + self.mubarNP)*self.Phi20NP + _sage_const_2 *(self.taubarNP - self.alphaNP)*self.Phi21NP + _sage_const_2 *self.lambdaNP*self.Phi11NP - _sage_const_2 *self.nuNP*self.Phi10NP - self.sigmabarNP*self.Phi22NP
        
        self.BI4=-self.DeltanNP(self.Psi3NP) + self.deltamNP(self.Psi4NP) - self.deltambarNP(self.Phi22NP) + self.DeltanNP(self.Phi21NP) + (_sage_const_4 *self.betaNP - self.tauNP)*self.Psi4NP - _sage_const_2 *(_sage_const_2 *self.muNP + self.gammaNP)*self.Psi3NP + _sage_const_3 *self.nuNP*self.Psi2NP + (self.taubarNP - _sage_const_2 *self.betabarNP - _sage_const_2 *self.alphaNP)*self.Phi22NP + _sage_const_2 *(self.gammaNP+self.mubarNP)*self.Phi21NP + _sage_const_2 *self.lambdaNP*self.Phi12NP - _sage_const_2 *self.nuNP*self.Phi11NP - self.nubarNP*self.Phi20NP
        
        self.BI5=-self.DlNP(self.Psi2NP) + self.deltambarNP(self.Psi1NP) - self.DeltanNP(self.Phi00NP) + self.deltambarNP(self.Phi01NP) - _sage_const_2 *self.DlNP(self.LambdaNP) - self.lambdaNP*self.Psi0NP + _sage_const_2 *(self.piNP-self.alphaNP)*self.Psi1NP + _sage_const_3 *self.rhoNP*self.Psi2NP - _sage_const_2 *self.kappaNP*self.Psi3NP + (_sage_const_2 *self.gammaNP + _sage_const_2 *self.gammabarNP - self.mubarNP)*self.Phi00NP - _sage_const_2 *(self.taubarNP + self.alphaNP)*self.Phi01NP - _sage_const_2 *self.tauNP*self.Phi10NP + _sage_const_2 *self.rhoNP*self.Phi11NP + self.sigmabarNP*self.Phi02NP
        
        self.BI6=-self.DeltanNP(self.Psi2NP) + self.deltamNP(self.Psi3NP) - self.DlNP(self.Phi22NP) + self.deltamNP(self.Phi21NP) - _sage_const_2 * self.DeltanNP(self.LambdaNP) + self.sigmaNP*self.Psi4NP + _sage_const_2 *(self.betaNP-self.tauNP)*self.Psi3NP - _sage_const_3 *self.muNP*self.Psi2NP + _sage_const_2 *self.nuNP*self.Psi1NP + (self.rhobarNP - _sage_const_2 *self.epsilonNP - _sage_const_2 *self.epsilonbarNP)*self.Phi22NP + _sage_const_2 *(self.pibarNP + self.betaNP)*self.Phi21NP + _sage_const_2 *self.piNP*self.Phi12NP - _sage_const_2 *self.muNP*self.Phi11NP - self.lambdabarNP*self.Phi20NP
        
        self.BI7=-self.DlNP(self.Psi3NP) + self.deltambarNP(self.Psi2NP) + self.DlNP(self.Phi21NP) -self.deltamNP(self.Phi20NP) + _sage_const_2 *self.deltambarNP(self.LambdaNP) - self.kappaNP*self.Psi4NP + _sage_const_2 *(self.rhoNP-self.epsilonNP)*self.Psi3NP + _sage_const_3 *self.piNP*self.Psi2NP - _sage_const_2 *self.lambdaNP*self.Psi1NP + (_sage_const_2 *self.alphabarNP-_sage_const_2 *self.betaNP - self.pibarNP)*self.Phi20NP - _sage_const_2 *(self.rhobarNP-self.epsilonNP)*self.Phi21NP - _sage_const_2 *self.piNP*self.Phi11NP + _sage_const_2 *self.muNP*self.Phi10NP + self.kappabarNP*self.Phi22NP
        
        self.BI8=-self.DeltanNP(self.Psi1NP) + self.deltamNP(self.Psi2NP) + self.DeltanNP(self.Phi01NP) - self.deltambarNP(self.Phi02NP) + _sage_const_2 *self.deltamNP(self.LambdaNP) + self.nuNP*self.Psi0NP + _sage_const_2 *(self.gammaNP-self.muNP)*self.Psi1NP - _sage_const_3 *self.tauNP*self.Psi2NP + _sage_const_2 *self.sigmaNP*self.Psi3NP + (self.taubarNP - _sage_const_2 *self.betabarNP + _sage_const_2 *self.alphaNP)*self.Phi02NP + _sage_const_2 *(self.mubarNP - self.gammaNP)*self.Phi01NP + _sage_const_2 *self.tauNP*self.Phi11NP-_sage_const_2 *self.rhoNP*self.Phi12NP-self.nubarNP*self.Phi00NP
        
        self.BI9=-self.DlNP(self.Phi11NP) + self.deltamNP(self.Phi10NP) + self.deltambarNP(self.Phi01NP) - self.DeltanNP(self.Phi00NP) - _sage_const_3 *self.DlNP(self.LambdaNP) + (_sage_const_2 *self.gammaNP - self.muNP + _sage_const_2 *self.gammabarNP - self.mubarNP)*self.Phi00NP + (self.piNP - _sage_const_2 *self.alphaNP - _sage_const_2 *self.taubarNP)*self.Phi01NP + (self.pibarNP - _sage_const_2 *self.alphabarNP - _sage_const_2 *self.tauNP)*self.Phi10NP + _sage_const_2 *(self.rhoNP+self.rhobarNP)*self.Phi11NP + self.sigmabarNP*self.Phi02NP + self.sigmaNP*self.Phi20NP - self.kappabarNP*self.Phi12NP - self.kappaNP*self.Phi21NP
        
        self.BI10=-self.DlNP(self.Phi12NP) + self.deltamNP(self.Phi11NP) + self.deltambarNP(self.Phi02NP) - self.DeltanNP(self.Phi01NP) - _sage_const_3 *self.deltamNP(self.LambdaNP) + (-_sage_const_2 *self.alphaNP + _sage_const_2 *self.betabarNP + self.piNP - self.taubarNP)*self.Phi02NP + (self.rhobarNP + _sage_const_2 *self.rhoNP - _sage_const_2 *self.epsilonbarNP)*self.Phi12NP + _sage_const_2 *(self.pibarNP - self.tauNP)*self.Phi11NP + (_sage_const_2 *self.gammaNP - _sage_const_2 *self.mubarNP - self.muNP)*self.Phi01NP + self.nubarNP*self.Phi00NP - self.lambdabarNP*self.Phi10NP + self.sigmaNP*self.Phi21NP - self.kappaNP*self.Phi22NP
        
        self.BI11=-self.DlNP(self.Phi22NP) + self.deltamNP(self.Phi21NP) + self.deltambarNP(self.Phi12NP) - self.DeltanNP(self.Phi11NP) - _sage_const_3 *self.DeltanNP(self.LambdaNP) + (self.rhoNP + self.rhobarNP - _sage_const_2 *self.epsilonNP - _sage_const_2 *self.epsilonbarNP)*self.Phi22NP + (_sage_const_2 *self.betabarNP + _sage_const_2 *self.piNP - self.taubarNP)*self.Phi12NP + (_sage_const_2 *self.betaNP + _sage_const_2 *self.pibarNP - self.tauNP)*self.Phi21NP - _sage_const_2 *(self.muNP+self.mubarNP)*self.Phi11NP + self.nuNP*self.Phi01NP + self.nubarNP*self.Phi10NP - self.lambdabarNP*self.Phi20NP - self.lambdaNP*self.Phi02NP


    
    #########################
    def show_Bianchi(self):
        show("BI1=",self.BI1.expr())
        show("BI2=",self.BI2.expr())
        show("BI3=",self.BI3.expr())
        show("BI4=",self.BI4.expr())
        show("BI5=",self.BI5.expr())
        show("BI6=",self.BI6.expr())
        show("BI7=",self.BI7.expr())
        show("BI8=",self.BI8.expr())
        show("BI9=",self.BI9.expr())
        show("BI10=",self.BI10.expr())
        show("BI11=",self.BI11.expr())

    ####################################################
    # Petrov Type
    # Petrov invariant I and J (Kramer p.54, 4.19)
    ######################### 
    def calculate_PetrovinvINP(self):
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        self.PetrovinvINPcalculated=_sage_const_1 
        self.PetrovinvINP=self.Psi0NP*self.Psi4NP-_sage_const_4 *self.Psi1NP*self.Psi3NP+_sage_const_3 *self.Psi2NP**_sage_const_2 
    #########################
    def calculate_PetrovinvJNP(self):
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        self.PetrovinvJNPcalculated=_sage_const_1 
        PetrovinvJmatrix=matrix([[self.Psi4NP,self.Psi3NP,self.Psi2NP],
                                 [self.Psi3NP,self.Psi2NP,self.Psi1NP],
                                 [self.Psi2NP,self.Psi1NP,self.Psi0NP]])
        self.PetrovinvJNP=PetrovinvJmatrix.determinant()
    #########################
    # Petrov invariant K, L, and N (Kramer p.121, 9.6)
    # also check diagram Fig. 9.1 on p. 122
    ######################### 
    def calculate_PetrovinvKNP(self):
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        self.PetrovinvKNPcalculated=_sage_const_1 
        self.PetrovinvKNP=self.Psi1NP*self.Psi4NP**_sage_const_2 -_sage_const_3 *self.Psi4NP*self.Psi3NP*self.Psi2NP+_sage_const_2 *self.Psi3NP**_sage_const_3 
    #########################    
    def calculate_PetrovinvLNP(self):
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        self.PetrovinvLNPcalculated=_sage_const_1 
        self.PetrovinvLNP=self.Psi2NP*self.Psi4NP-self.Psi3NP**_sage_const_2 
    ######################### 
    def calculate_PetrovinvNNP(self):
        if self.weylcalculated==_sage_const_0:
            self.calculate_Weyl()
        if self.PetrovinvLNPcalculated==_sage_const_0 :
            self.calculate_PetrovinvLNP()
        if self.PetrovinvINPcalculated==_sage_const_0 :
            self.calculate_PetrovinvINP()
        self.PetrovinvNNPcalculated=_sage_const_1 
        self.PetrovinvNNP=_sage_const_12 *self.PetrovinvLNP**_sage_const_2 -(self.Psi4NP**_sage_const_2 )*self.PetrovinvINP
    #########################    
    # Petrov type from the invariants
    def Petrov_frominvariants(self):
        if self.weylcalculated==_sage_const_0 :
            self.calculate_Weyl()
            self.Petrov_frominvariants()
        else:
            psi0=NewmanPenrose.simplify_fullfull(self.Psi0NP.expr())
            psi1=NewmanPenrose.simplify_fullfull(self.Psi1NP.expr())
            psi2=NewmanPenrose.simplify_fullfull(self.Psi2NP.expr())
            psi3=NewmanPenrose.simplify_fullfull(self.Psi3NP.expr())
            psi4=NewmanPenrose.simplify_fullfull(self.Psi4NP.expr())
            if psi0==_sage_const_0  and psi1==_sage_const_0  and psi2==_sage_const_0  and psi3==_sage_const_0  and psi4==_sage_const_0 :
                self.petrovtype = "O"
                print("Petrov Type O")
                return
            if self.PetrovinvINPcalculated==_sage_const_0 :
                self.calculate_PetrovinvINP()
            if self.PetrovinvJNPcalculated==_sage_const_0 :
                self.calculate_PetrovinvJNP()
            if self.PetrovinvKNPcalculated==_sage_const_0 :
                self.calculate_PetrovinvKNP()
            if self.PetrovinvLNPcalculated==_sage_const_0 :
                self.calculate_PetrovinvLNP()
            if self.PetrovinvNNPcalculated==_sage_const_0 :
                self.calculate_PetrovinvNNP()
            #########################
            print("Calculating Petrov Type...")
            LL=NewmanPenrose.simplify_fullfull(self.PetrovinvLNP.expr())
            JJ=NewmanPenrose.simplify_fullfull(self.PetrovinvJNP.expr())
            II=NewmanPenrose.simplify_fullfull(self.PetrovinvINP.expr())
            KK=NewmanPenrose.simplify_fullfull(self.PetrovinvKNP.expr())
            NN=NewmanPenrose.simplify_fullfull(self.PetrovinvNNP.expr())
            IIJJterm=NewmanPenrose.simplify_fullfull((II**_sage_const_3 -_sage_const_27 *JJ**_sage_const_2 ))
            if IIJJterm==_sage_const_0 :
                if II==_sage_const_0  and JJ==_sage_const_0 :
                    if KK==_sage_const_0  and LL==_sage_const_0 :
                        self.petrovtype = "N"
                        print("Petrov Type N")
                    else:
                        self.petrovtype = "III"
                        print("Petrov Type III")
                else:
                    if KK==_sage_const_0  and NN==_sage_const_0 :
                        self.petrovtype = "D"
                        print("Petrov Type D")
                    else:
                        self.petrovtype = "II"
                        print("Petrov Type II")
            else:
                self.petrovtype = "I"
                print("Petrov Type I")
            print("Attention: This procedure depends on the simplification of the structures. Therefore the Petrov type can be simpler.")
    
    #########################
    # Petrov type from the Weyl components
    def Petrov_fromWeyl(self):
        if self.weylcalculated==_sage_const_0 :
            self.calculate_Weyl()
            self.Petrov_fromWeyl()
        else:
            print("Calculating Petrov Type...")
            psi0=NewmanPenrose.simplify_fullfull(self.Psi0NP.expr())
            psi1=NewmanPenrose.simplify_fullfull(self.Psi1NP.expr())
            psi2=NewmanPenrose.simplify_fullfull(self.Psi2NP.expr())
            psi3=NewmanPenrose.simplify_fullfull(self.Psi3NP.expr())
            psi4=NewmanPenrose.simplify_fullfull(self.Psi4NP.expr())
            if psi0==_sage_const_0  and psi1==_sage_const_0  and psi2==_sage_const_0  and psi3==_sage_const_0  and psi4==_sage_const_0 :
                self.petrovtype = "O"
                print("Petrov Type O")
                return 
            E = var('E')
            
            P = psi0 - _sage_const_4*psi1*E + _sage_const_6*psi2*E**_sage_const_2 - _sage_const_4*psi3*E**_sage_const_3 + psi4*E**_sage_const_4

            try:
                degP = P.degree(E)
                dP = diff(P, E)
                G = gcd(P, dP)
    
                deg = G.degree(E) if G.has(E) else _sage_const_0
                
                if degP == _sage_const_4 or degP == _sage_const_3:
                    if deg == _sage_const_1:
                        self.petrovtype = "II"
                        print("Petrov Type II")
                        return
                    elif deg == _sage_const_2:
                        if gcd(G**_sage_const_2, P) == P:
                            self.petrovtype = "D"
                            print("Petrov Type D")
                            return
                        else:
                            self.petrovtype = "III"
                            print("Petrov Type III")
                            return
                    elif deg == _sage_const_3:
                        self.petrovtype = "N"
                        print("Petrov Type N")
                        return
                    else:
                        self.petrovtype = "I"
                        print("Petrov Type I")
                        return
                elif degP == _sage_const_2:
                    if deg == _sage_const_1:
                        self.petrovtype = "D"
                        print("Petrov Type D")
                        return
                    else:
                        self.petrovtype = "II"
                        print("Petrov Type II")
                        return
                elif degP == sage_const_1:
                    self.petrovtype = "III"
                    print("Petrov Type III")
                    return
                else:
                    self.petrovtype = " "
                    print("Petrov Type could not found")
            except:
                try:
                    self.Petrov_frominvariants()
                except:
                    self.petrovtype = " "
                    print("Error: Unable to calculate the Petrov type.")
                     
    ####################################################
    # Calculate everything about NP
    #########################
    def calculate_allNP(self):
        self.calculate_spincoefficients()
        self.calculate_Weyl()
        self.calculate_Ricci()
        self.calculate_NPeq()
        self.calculate_Bianchi()
        self.Petrov_fromWeyl()


    #########################
    def show_allNP(self):
        self.show_spincoefficients()
        self.show_Weyl()
        self.show_Ricci()
        try:
            self.show_NPeq()
        except:
            self.calculate_NPeq()
            self.show_NPeq()
        try:
            self.show_Bianchi()
        except:
            self.calculate_Bianchi()
            self.show_Bianchi()
        try:
            show("Petrov type " + self.petrovtype)
        except:
            self.Petrov_fromWeyl()
            show("Petrov type " + self.petrovtype)
    ####################################################
    # Type A Transformations defined in Carmeli and Kaye, Annals of Physics 99, 188 (1976)
    # Null rotation about l_mu
    def type_A_transformation(self,z):
        # Definitions:
        self.kappaNP_trA=self.Man.scalar_field(function('kappaNP_trA')(*self.CO));
        self.tauNP_trA=self.Man.scalar_field(function('tauNP_trA')(*self.CO));    
        self.sigmaNP_trA=self.Man.scalar_field(function('sigmaNP_trA')(*self.CO));
        self.rhoNP_trA=self.Man.scalar_field(function('rhoNP_trA')(*self.CO));
        self.piNP_trA=self.Man.scalar_field(function('piNP_trA')(*self.CO));
        self.nuNP_trA=self.Man.scalar_field(function('nuNP_trA')(*self.CO));
        self.muNP_trA=self.Man.scalar_field(function('muNP_trA')(*self.CO));
        self.lambdaNP_trA=self.Man.scalar_field(function('lambdaNP_trA')(*self.CO));
        self.epsilonNP_trA=self.Man.scalar_field(function('epsilonNP_trA')(*self.CO));
        self.gammaNP_trA=self.Man.scalar_field(function('gammaNP_trA')(*self.CO));
        self.betaNP_trA=self.Man.scalar_field(function('betaNP_trA')(*self.CO));
        self.alphaNP_trA=self.Man.scalar_field(function('alphaNP_trA')(*self.CO));
        
        self.Psi0NP_trA=self.Man.scalar_field(function('Psi0NP_trA')(*self.CO))
        self.Psi1NP_trA=self.Man.scalar_field(function('Psi1NP_trA')(*self.CO))
        self.Psi2NP_trA=self.Man.scalar_field(function('Psi2NP_trA')(*self.CO)) 
        self.Psi3NP_trA=self.Man.scalar_field(function('Psi3NP_trA')(*self.CO))
        self.Psi4NP_trA=self.Man.scalar_field(function('Psi4NP_trA')(*self.CO))

        self.Phi00NP_trA=self.Man.scalar_field(function('Phi00NP_trA')(*self.CO))
        self.Phi01NP_trA=self.Man.scalar_field(function('Phi01NP_trA')(*self.CO))
        self.Phi10NP_trA=self.Man.scalar_field(function('Phi10NP_trA')(*self.CO))
        self.Phi02NP_trA=self.Man.scalar_field(function('Phi02NP_trA')(*self.CO))
        self.Phi20NP_trA=self.Man.scalar_field(function('Phi20NP_trA')(*self.CO))
        self.Phi11NP_trA=self.Man.scalar_field(function('Phi11NP_trA')(*self.CO))
        self.Phi12NP_trA=self.Man.scalar_field(function('Phi12NP_trA')(*self.CO))
        self.Phi21NP_trA=self.Man.scalar_field(function('Phi21NP_trA')(*self.CO))
        self.Phi22NP_trA=self.Man.scalar_field(function('Phi22NP_trA')(*self.CO))
        self.LambdaNP_trA=self.Man.scalar_field(function('LambdaNP_trA')(*self.CO))
        
        # Transformations for the null tetrad
        self.lNP_trA = self.lNP[:] 
        self.mNP_trA = [self.mNP[i] + z*self.lNP[i] for i in range(4)]
        self.mbarNP_trA = [self.mbarNP[i] + conjugate(z)*self.lNP[i] for i in range(4)]
        self.nNP_trA = [self.nNP[i] + z*self.mbarNP[i] + conjugate(z)*self.mNP[i] + z*conjugate(z)*self.lNP[i] for i in range(4)]

        #Transformations for the spin coefficients
        self.rhoNP_trA = self.rhoNP.expr() + z*self.kappaNP.expr()
        self.sigmaNP_trA = self.sigmaNP.expr() + z*self.kappaNP.expr()
        self.tauNP_trA = self.tauNP.expr() + z*self.epsilonNP.expr() + z**2*self.kappaNP.expr()
        self.kappaNP_trA = self.kappaNP.expr()
        self.nuNP_trA = self.nuNP.expr() + z*(self.piNP.expr() + 2*self.tauNP.expr()) + z**2*(2*self.epsilonNP.expr() + self.rhoNP.expr()) \
                                 + z**3*self.kappaNP.expr()
        self.muNP_trA = self.muNP.expr() + 2*z*self.rhoNP.expr() + z*self.piNP.expr() + z**2*self.sigmaNP.expr() \
                                     + 2*z**2*self.epsilonNP.expr() + z**2*conjugate(z)*self.kappaNP.expr()
        self.lambdaNP_trA = self.lambdaNP.expr() + z*self.kappaNP.expr()
        self.piNP_trA = self.piNP.expr() + z*self.sigmaNP.expr() + z*conjugate(z)*self.kappaNP.expr()
        self.epsilonNP_trA = self.epsilonNP.expr() + z*self.kappaNP.expr()
        self.betaNP_trA = self.betaNP.expr() + z*self.sigmaNP.expr() + conjugate(z)*self.epsilonNP.expr() + z**2*self.kappaNP.expr()
        self.gammaNP_trA = self.gammaNP.expr() + z*(self.tauNP.expr() + self.betaNP.expr()) + conjugate(z)*(self.epsilonNP.expr() + self.rhoNP.expr()) + z**2*self.sigmaNP.expr() + z*conjugate(z)*self.epsilonNP.expr() + z**2*conjugate(z)*self.kappaNP.expr()
        self.alphaNP_trA = self.alphaNP.expr() + z*(self.rhoNP.expr() + self.epsilonNP.expr()) + z**2*self.kappaNP.expr()
        
        # Transformations for the Weyl tensor components 
        self.Psi0NP_trA = self.Psi0NP.expr()
        self.Psi1NP_trA = self.Psi1NP.expr() + z*self.Psi0NP.expr()
        self.Psi2NP_trA = self.Psi2NP.expr() + 2*z*self.Psi1NP.expr() + z**2*self.Psi0NP.expr()
        self.Psi3NP_trA = self.Psi3NP.expr() + 3*z*self.Psi2NP.expr() + 3*z**2*self.Psi1NP.expr() + z**3*self.Psi0NP.expr()
        self.Psi4NP_trA = self.Psi4NP.expr() + 4*z*self.Psi3NP.expr() + 6*z**2*self.Psi2NP.expr() + 4*z**3*self.Psi1NP.expr() + z**4*self.Psi0NP.expr()
      
        # Transformations for the Ricci tensor components
        self.Phi00NP_trA = self.Phi00NP.expr()
        self.Phi01NP_trA = self.Phi01NP.expr() + z*self.Phi00NP.expr()
        self.Phi02NP_trA = self.Phi02NP.expr() + 2*z*self.Phi01NP.expr() + z**2*self.Phi00NP.expr()
        self.Phi10NP_trA = self.Phi10NP.expr() + conjugate(z)*self.Phi00NP.expr()
        self.Phi11NP_trA = self.Phi11NP.expr() + z*self.Phi10NP.expr() + conjugate(z)*self.Phi01NP.expr() + z*conjugate(z)*self.Phi00NP.expr()
        self.Phi12NP_trA = self.Phi12NP.expr() + 2*z*self.Phi11NP.expr() + z**2*self.Phi10NP.expr() + conjugate(z)*self.Phi02NP.expr() \
                     + 2*z*conjugate(z)*self.Phi01NP.expr() + z**2*conjugate(z)*self.Phi00NP.expr()
        self.Phi20NP_trA = self.Phi20NP.expr() + 2*conjugate(z)*self.Phi10NP.expr() + conjugate(z)**2*self.Phi00NP.expr()
        self.Phi21NP_trA = self.Phi21NP.expr() + z*self.Phi20NP.expr() + 2*conjugate(z)*self.Phi11NP.expr() + 2*z*conjugate(z)*self.Phi10NP.expr() \
                     + conjugate(z)**2*self.Phi01NP.expr() + z*conjugate(z)**2*self.Phi00NP.expr()
        self.Phi22NP_trA = self.Phi22NP.expr() + 2*z*self.Phi21NP.expr() + z**2*self.Phi20NP.expr() + 2*conjugate(z)*self.Phi12NP.expr() \
                     + 4*z*conjugate(z)*self.Phi11NP.expr() + 2*z**2*conjugate(z)*self.Phi10NP.expr() + conjugate(z)**2*self.Phi02NP.expr() \
                     + 2*z*conjugate(z)**2*self.Phi01NP.expr() + z**2*conjugate(z)**2*self.Phi00NP.expr()
        
        self.LambdaNP_trA = self.LambdaNP.expr()
    #########################
    # Type B Transformations defined in Carmeli and Kaye, Annals of Physics 99, 188 (1976)
    # Boost in l^mu-n^mu plane and spacial rotation in m^mu-mbar^mu plane
    # z can be given as z = sqrt(A)*exp(I*theta/2)
    def type_B_transformation(self,z):
        # Definitions:
        self.kappaNP_trB=self.Man.scalar_field(function('kappaNP_trB')(*self.CO));
        self.tauNP_trB=self.Man.scalar_field(function('tauNP_trB')(*self.CO));    
        self.sigmaNP_trB=self.Man.scalar_field(function('sigmaNP_trB')(*self.CO));
        self.rhoNP_trB=self.Man.scalar_field(function('rhoNP_trB')(*self.CO));
        self.piNP_trB=self.Man.scalar_field(function('piNP_trB')(*self.CO));
        self.nuNP_trB=self.Man.scalar_field(function('nuNP_trB')(*self.CO));
        self.muNP_trB=self.Man.scalar_field(function('muNP_trB')(*self.CO));
        self.lambdaNP_trB=self.Man.scalar_field(function('lambdaNP_trB')(*self.CO));
        self.epsilonNP_trB=self.Man.scalar_field(function('epsilonNP_trB')(*self.CO));
        self.gammaNP_trB=self.Man.scalar_field(function('gammaNP_trB')(*self.CO));
        self.betaNP_trB=self.Man.scalar_field(function('betaNP_trB')(*self.CO));
        self.alphaNP_trB=self.Man.scalar_field(function('alphaNP_trB')(*self.CO));
        
        self.Psi0NP_trB=self.Man.scalar_field(function('Psi0NP_trB')(*self.CO))
        self.Psi1NP_trB=self.Man.scalar_field(function('Psi1NP_trB')(*self.CO))
        self.Psi2NP_trB=self.Man.scalar_field(function('Psi2NP_trB')(*self.CO)) 
        self.Psi3NP_trB=self.Man.scalar_field(function('Psi3NP_trB')(*self.CO))
        self.Psi4NP_trB=self.Man.scalar_field(function('Psi4NP_trB')(*self.CO))

        self.Phi00NP_trB=self.Man.scalar_field(function('Phi00NP_trB')(*self.CO))
        self.Phi01NP_trB=self.Man.scalar_field(function('Phi01NP_trB')(*self.CO))
        self.Phi10NP_trB=self.Man.scalar_field(function('Phi10NP_trB')(*self.CO))
        self.Phi02NP_trB=self.Man.scalar_field(function('Phi02NP_trB')(*self.CO))
        self.Phi20NP_trB=self.Man.scalar_field(function('Phi20NP_trB')(*self.CO))
        self.Phi11NP_trB=self.Man.scalar_field(function('Phi11NP_trB')(*self.CO))
        self.Phi12NP_trB=self.Man.scalar_field(function('Phi12NP_trB')(*self.CO))
        self.Phi21NP_trB=self.Man.scalar_field(function('Phi21NP_trB')(*self.CO))
        self.Phi22NP_trB=self.Man.scalar_field(function('Phi22NP_trB')(*self.CO))
        self.LambdaNP_trB=self.Man.scalar_field(function('LambdaNP_trB')(*self.CO))
        
        # Transformations for the null tetrad
        self.lNP_trB = [abs(z)**2 * x for x in self.lNP[:]]
        self.nNP_trB = [1/abs(z)**2 * x for x in self.nNP[:]]
        self.mNP_trB = [conjugate(z)/abs(z) * x for x in self.mNP[:]]
        self.mbarNP_trB = [z/abs(z) * x for x in self.mbarNP[:]]
        
        # Transformations for the spin coefficients
        self.rhoNP_trB = z*conjugate(z)*self.rhoNP.expr()
        self.sigmaNP_trB = z**2*self.sigmaNP.expr()
        self.tauNP_trB = conjugate(z)**2*self.tauNP.expr()
        self.kappaNP_trB = z**3*self.kappaNP.expr()
        self.nuNP_trB = z**(-3)*self.nuNP.expr()
        self.muNP_trB = z**(-1)*conjugate(z)**(-1)*self.muNP.expr()
        self.lambdaNP_trB = z**(-2)*self.lambdaNP.expr()
        self.piNP_trB = z**(-1)*conjugate(z)*self.piNP.expr()
        self.epsilonNP_trB = z*conjugate(z)*self.epsilonNP.expr()
        self.betaNP_trB = conjugate(z)/z*self.betaNP.expr()
        self.gammaNP_trB = z**(-1)*conjugate(z)**(-1)*self.gammaNP.expr()
        self.alphaNP_trB = z**(-1)*conjugate(z)*self.alphaNP.expr()
        
        # Transformations for the Weyl tensor components 
        self.Psi0NP_trB = z**4*self.Psi0NP.expr()
        self.Psi1NP_trB = z**3*self.Psi1NP.expr()
        self.Psi2NP_trB = z**2*self.Psi2NP.expr()
        self.Psi3NP_trB = z*self.Psi3NP.expr()
        self.Psi4NP_trB = self.Psi4NP.expr()
        
        # Transformations for the Ricci tensor components
        self.Phi00NP_trB = z**2*conjugate(z)**2*self.Phi00NP.expr()
        self.Phi01NP_trB = z**2*conjugate(z)*self.Phi01NP.expr()
        self.Phi02NP_trB = z**2*self.Phi02NP.expr()
        self.Phi10NP_trB = z*conjugate(z)**2*self.Phi10NP.expr()
        self.Phi11NP_trB = z*conjugate(z)*self.Phi11NP.expr()
        self.Phi12NP_trB = z*self.Phi12NP.expr()
        self.Phi20NP_trB = conjugate(z)**2*self.Phi20NP.expr()
        self.Phi21NP_trB = conjugate(z)*self.Phi21NP.expr()
        self.Phi22NP_trB = self.Phi22NP.expr()
        self.LambdaNP_trB = self.LambdaNP.expr()
    #########################
    # Type C Transformations defined in Carmeli and Kaye, Annals of Physics 99, 188 (1976)
    # Null rotation about n_mu
    # z can be given as z = sqrt(A)*exp(I*theta/2)
    def type_C_transformation(self,z):
        # Definitions:
        self.kappaNP_trC=self.Man.scalar_field(function('kappaNP_trC')(*self.CO));
        self.tauNP_trC=self.Man.scalar_field(function('tauNP_trC')(*self.CO));    
        self.sigmaNP_trC=self.Man.scalar_field(function('sigmaNP_trC')(*self.CO));
        self.rhoNP_trC=self.Man.scalar_field(function('rhoNP_trC')(*self.CO));
        self.piNP_trC=self.Man.scalar_field(function('piNP_trC')(*self.CO));
        self.nuNP_trC=self.Man.scalar_field(function('nuNP_trC')(*self.CO));
        self.muNP_trC=self.Man.scalar_field(function('muNP_trC')(*self.CO));
        self.lambdaNP_trC=self.Man.scalar_field(function('lambdaNP_trC')(*self.CO));
        self.epsilonNP_trC=self.Man.scalar_field(function('epsilonNP_trC')(*self.CO));
        self.gammaNP_trC=self.Man.scalar_field(function('gammaNP_trC')(*self.CO));
        self.betaNP_trC=self.Man.scalar_field(function('betaNP_trC')(*self.CO));
        self.alphaNP_trC=self.Man.scalar_field(function('alphaNP_trC')(*self.CO));
        
        self.Psi0NP_trC=self.Man.scalar_field(function('Psi0NP_trC')(*self.CO))
        self.Psi1NP_trC=self.Man.scalar_field(function('Psi1NP_trC')(*self.CO))
        self.Psi2NP_trC=self.Man.scalar_field(function('Psi2NP_trC')(*self.CO)) 
        self.Psi3NP_trC=self.Man.scalar_field(function('Psi3NP_trC')(*self.CO))
        self.Psi4NP_trC=self.Man.scalar_field(function('Psi4NP_trC')(*self.CO))

        self.Phi00NP_trC=self.Man.scalar_field(function('Phi00NP_trC')(*self.CO))
        self.Phi01NP_trC=self.Man.scalar_field(function('Phi01NP_trC')(*self.CO))
        self.Phi10NP_trC=self.Man.scalar_field(function('Phi10NP_trC')(*self.CO))
        self.Phi02NP_trC=self.Man.scalar_field(function('Phi02NP_trC')(*self.CO))
        self.Phi20NP_trC=self.Man.scalar_field(function('Phi20NP_trC')(*self.CO))
        self.Phi11NP_trC=self.Man.scalar_field(function('Phi11NP_trC')(*self.CO))
        self.Phi12NP_trC=self.Man.scalar_field(function('Phi12NP_trC')(*self.CO))
        self.Phi21NP_trC=self.Man.scalar_field(function('Phi21NP_trC')(*self.CO))
        self.Phi22NP_trC=self.Man.scalar_field(function('Phi22NP_trC')(*self.CO))
        self.LambdaNP_trC=self.Man.scalar_field(function('LambdaNP_trC')(*self.CO))
        
        # Transformations for the null tetrad
        self.lNP_trC = [self.lNP[i] + conjugate(z)*self.mNP[i] + z*self.mbarNP[i] + z*conjugate(z)*self.nNP[i] for i in range(4)]
        self.mNP_trC = [self.mNP[i] + z*self.nNP[i] for i in range(4)]
        self.mbarNP_trC = [self.mbarNP[i] + conjugate(z)*self.nNP[i] for i in range(4)]
        self.nNP_trC = self.nNP[:] # n' = n
        
        # Transformations for the spin coefficients
        self.muNP_trC = self.muNP.expr() + 2*z*self.alphaNP.expr() + 2*self.tauNP.expr() + z**2*self.lambdaNP.expr() \
                   + 2*z*conjugate(z)*self.gammaNP.expr() + z**2*conjugate(z)*self.nuNP.expr()
        self.alphaNP_trC = self.alphaNP.expr() + z*self.piNP.expr() + z*self.gammaNP.expr() + z**2*self.nuNP.expr()
        self.lambdaNP_trC = self.lambdaNP.expr() + z*self.nuNP.expr()
        self.kappaNP_trC = self.kappaNP.expr() + z*(self.sigmaNP.expr() + 2*self.epsilonNP.expr()) + z*conjugate(z)*self.betaNP.expr() \
                     + z**2*(2*self.alphaNP.expr() + self.piNP.expr()) + z**2*conjugate(z)*(2*self.gammaNP.expr() + self.muNP.expr()) \
                     + z**3*self.lambdaNP.expr() + z**3*conjugate(z)*self.nuNP.expr()
        self.epsilonNP_trC = self.epsilonNP.expr() + z*(self.piNP.expr() + self.alphaNP.expr()) + z*conjugate(z)*self.betaNP.expr() \
                       + z**2*self.lambdaNP.expr() + z*conjugate(z)*(self.muNP.expr() + self.gammaNP.expr()) \
                       + z**2*conjugate(z)*self.nuNP.expr()
        self.piNP_trC = self.piNP.expr() + 2*z*self.piNP.expr() + z**2*self.nuNP.expr()
        self.sigmaNP_trC = self.sigmaNP.expr() + z*(self.tauNP.expr() + 2*self.betaNP.expr()) \
                     + z**2*(self.muNP.expr() + 2*self.gammaNP.expr()) + z**3*self.nuNP.expr()
        self.rhoNP_trC = self.rhoNP.expr() + z*self.muNP.expr() + conjugate(z)*self.muNP.expr() + z*conjugate(z)*self.nuNP.expr()
        self.tauNP_trC = self.tauNP.expr() + 2*z*self.gammaNP.expr() + z**2*self.nuNP.expr()
        self.gammaNP_trC = self.gammaNP.expr() + z*self.nuNP.expr()
        self.nuNP_trC = self.nuNP.expr()
        self.betaNP_trC = self.betaNP.expr() + z*self.gammaNP.expr() + z*self.muNP.expr() + z**2*self.nuNP.expr()
        
        # Transformations for the Weyl tensor components 
        self.Psi0NP_trC = self.Psi0NP.expr() + 4*z*self.Psi1NP.expr() + 6*z**2*self.Psi2NP.expr() \
                     + 4*z**3*self.Psi3NP.expr() + z**4*self.Psi4NP.expr()
        self.Psi1NP_trC = self.Psi1NP.expr() + 3*z*self.Psi2NP.expr() + 3*z**2*self.Psi3NP.expr() + z**3*self.Psi4NP.expr()
        self.Psi2NP_trC = self.Psi2NP.expr() + 2*z*self.Psi3NP.expr() + z**2*self.Psi4NP.expr()
        self.Psi3NP_trC = self.Psi3NP.expr() + z*self.Psi4NP.expr()
        self.Psi4NP_trC = self.Psi4NP.expr()
        
        # Transformations for the Ricci tensor components
        self.Phi00NP_trC = self.Phi00NP.expr() + 2*z*self.Phi01NP.expr() + 2*conjugate(z)*self.Phi10NP.expr() \
                     + 4*z*conjugate(z)*self.Phi11NP.expr() + z**2*self.Phi02NP.expr() \
                     + conjugate(z)**2*self.Phi20NP.expr() + 2*z**2*conjugate(z)*self.Phi12NP.expr() \
                     + 2*z*conjugate(z)**2*self.Phi21NP.expr() + z**2*conjugate(z)**2*self.Phi22NP.expr()
        self.Phi01NP_trC = self.Phi01NP.expr() + 2*z*self.Phi11NP.expr() + z**2*self.Phi21NP.expr() \
                     + conjugate(z)*self.Phi02NP.expr() + 2*z*conjugate(z)*self.Phi12NP.expr() \
                     + z**2*conjugate(z)*self.Phi22NP.expr()
        self.Phi02NP_trC = self.Phi02NP.expr() + 2*z*self.Phi12NP.expr() + z**2*self.Phi22NP.expr()
        self.Phi10NP_trC = self.Phi10NP.expr() + 2*conjugate(z)*self.Phi11NP.expr() + conjugate(z)**2*self.Phi12NP.expr() \
                     + z*self.Phi20NP.expr() + 2*z*conjugate(z)*self.Phi21NP.expr() + z*conjugate(z)**2*self.Phi22NP.expr()
        self.Phi11NP_trC = self.Phi11NP.expr() + z*self.Phi21NP.expr() + conjugate(z)*self.Phi12NP.expr() \
                     + z*conjugate(z)*self.Phi22NP.expr()
        self.Phi12NP_trC = self.Phi12NP.expr() + z*self.Phi22NP.expr()
        self.Phi20NP_trC = self.Phi20NP.expr() + 2*conjugate(z)*self.Phi21NP.expr() + conjugate(z)**2*self.Phi22NP.expr() 
        self.Phi21NP_trC = self.Phi21NP.expr() + conjugate(z)*self.Phi22NP.expr()
        self.Phi22NP_trC = self.Phi22NP.expr()
        self.LambdaNP_trC = self.LambdaNP.expr()
    #########################
    def show_type_A_transformation(self):
        show("Results of the Type A Transformations:")
        show("lNP_trA=",self.lNP_trA)
        show("nNP_trA=",self.nNP_trA)
        show("mNP_trA=",self.mNP_trA)
        show("mbarNP_trA=",self.mbarNP_trA)
        
        show("kappaNP_trA=",self.kappaNP_trA)
        show("tauNP_trA=",self.tauNP_trA)
        show("sigmaNP_trA=",self.sigmaNP_trA)
        show("rhoNP_trA=",self.rhoNP_trA)
        show("piNP_trA=",self.piNP_trA)
        show("nuNP_trA=",self.nuNP_trA)
        show("muNP_trA=",self.muNP_trA)
        show("lambdaNP_trA=",self.lambdaNP_trA)
        show("epsilonNP_trA=",self.epsilonNP_trA)
        show("gammaNP_trA=",self.gammaNP_trA)
        show("betaNP_trA=",self.betaNP_trA)
        show("alphaNP_trA=",self.alphaNP_trA)

        show("Psi0NP_trA=",self.Psi0NP_trA)
        show("Psi1NP_trA=",self.Psi1NP_trA)
        show("Psi2NP_trA=",self.Psi2NP_trA)
        show("Psi3NP_trA=",self.Psi3NP_trA)
        show("Psi4NP_trA=",self.Psi4NP_trA)

        show("Phi00NP_trA=",self.Phi00NP_trA)
        show("Phi01NP_trA=",self.Phi01NP_trA)
        show("Phi02NP_trA=",self.Phi02NP_trA)
        show("Phi10NP_trA=",self.Phi10NP_trA)
        show("Phi11NP_trA=",self.Phi11NP_trA)
        show("Phi12NP_trA=",self.Phi12NP_trA)
        show("Phi20NP_trA=",self.Phi20NP_trA)
        show("Phi21NP_trA=",self.Phi21NP_trA)
        show("Phi22NP_trA=",self.Phi22NP_trA)
        show("LambdaNP_trA=",self.LambdaNP_trA)
    #########################
    def show_type_B_transformation(self):
        show("Results of the Type B Transformations:")
        show("lNP_trB=",self.lNP_trB)
        show("nNP_trB=",self.nNP_trB)
        show("mNP_trB=",self.mNP_trB)
        show("mbarNP_trB=",self.mbarNP_trB)
        
        show("kappaNP_trB=",self.kappaNP_trB)
        show("tauNP_trB=",self.tauNP_trB)
        show("sigmaNP_trB=",self.sigmaNP_trB)
        show("rhoNP_trB=",self.rhoNP_trB)
        show("piNP_trB=",self.piNP_trB)
        show("nuNP_trB=",self.nuNP_trB)
        show("muNP_trB=",self.muNP_trB)
        show("lambdaNP_trB=",self.lambdaNP_trB)
        show("epsilonNP_trB=",self.epsilonNP_trB)
        show("gammaNP_trB=",self.gammaNP_trB)
        show("betaNP_trB=",self.betaNP_trB)
        show("alphaNP_trB=",self.alphaNP_trB)

        show("Psi0NP_trB=",self.Psi0NP_trB)
        show("Psi1NP_trB=",self.Psi1NP_trB)
        show("Psi2NP_trB=",self.Psi2NP_trB)
        show("Psi3NP_trB=",self.Psi3NP_trB)
        show("Psi4NP_trB=",self.Psi4NP_trB)

        show("Phi00NP_trB=",self.Phi00NP_trB)
        show("Phi01NP_trB=",self.Phi01NP_trB)
        show("Phi02NP_trB=",self.Phi02NP_trB)
        show("Phi10NP_trB=",self.Phi10NP_trB)
        show("Phi11NP_trB=",self.Phi11NP_trB)
        show("Phi12NP_trB=",self.Phi12NP_trB)
        show("Phi20NP_trB=",self.Phi20NP_trB)
        show("Phi21NP_trB=",self.Phi21NP_trB)
        show("Phi22NP_trB=",self.Phi22NP_trB)
        show("LambdaNP_trB=",self.LambdaNP_trB)
    #########################
    def show_type_C_transformation(self):
        show("Results of the Type C Transformations:")
        show("lNP_trC=",self.lNP_trC)
        show("nNP_trC=",self.nNP_trC)
        show("mNP_trC=",self.mNP_trC)
        show("mbarNP_trC=",self.mbarNP_trC)
        
        show("kappaNP_trC=",self.kappaNP_trC)
        show("tauNP_trC=",self.tauNP_trC)
        show("sigmaNP_trC=",self.sigmaNP_trC)
        show("rhoNP_trC=",self.rhoNP_trC)
        show("piNP_trC=",self.piNP_trC)
        show("nuNP_trC=",self.nuNP_trC)
        show("muNP_trC=",self.muNP_trC)
        show("lambdaNP_trC=",self.lambdaNP_trC)
        show("epsilonNP_trC=",self.epsilonNP_trC)
        show("gammaNP_trC=",self.gammaNP_trC)
        show("betaNP_trC=",self.betaNP_trC)
        show("alphaNP_trC=",self.alphaNP_trC)

        show("Psi0NP_trC=",self.Psi0NP_trC)
        show("Psi1NP_trC=",self.Psi1NP_trC)
        show("Psi2NP_trC=",self.Psi2NP_trC)
        show("Psi3NP_trC=",self.Psi3NP_trC)
        show("Psi4NP_trC=",self.Psi4NP_trC)

        show("Phi00NP_trC=",self.Phi00NP_trC)
        show("Phi01NP_trC=",self.Phi01NP_trC)
        show("Phi02NP_trC=",self.Phi02NP_trC)
        show("Phi10NP_trC=",self.Phi10NP_trC)
        show("Phi11NP_trC=",self.Phi11NP_trC)
        show("Phi12NP_trC=",self.Phi12NP_trC)
        show("Phi20NP_trC=",self.Phi20NP_trC)
        show("Phi21NP_trC=",self.Phi21NP_trC)
        show("Phi22NP_trC=",self.Phi22NP_trC)
        show("LambdaNP_trC=",self.LambdaNP_trC)
    ####################################################
    # Klein-Gordon equation for a massive scalar field
    # as defined in G. Silva-Ortigoza, Rev. Mex. Fis. 4, 543 (1996).
    # We used an overall minus sign due to our definitions.
    # Set as metricalias.kgNP=0.
    # Example:
    # var('M2') # mass
    # Phi = manifoldalias.scalar_field(function('Phi')(t,r,th,ph))
    # metricalias.kleingordon(Phi,M2)
    # show(metricalias.kgNP.expr())
    def kleingordon(self,f,mass):
        Deltanf=self.DeltanNP(f)
        deltambarf=self.deltambarNP(f)
        deltamf=self.deltamNP(f)
        Dlf=self.DlNP(f)
        term1=self.DlNP(Deltanf)+(self.epsilonNP+self.epsilonbarNP-self.rhobarNP-self.rhoNP)*Deltanf
        term2=self.deltamNP(deltambarf)+(self.betaNP-self.alphabarNP+self.pibarNP-self.tauNP)*deltambarf
        term3=self.DeltanNP(Dlf)+(-self.gammaNP-self.gammabarNP+self.mubarNP+self.muNP)*Dlf
        term4=self.deltambarNP(deltamf)+(self.betabarNP-self.alphaNP+self.piNP-self.taubarNP)*deltamf
        self.kgNP=-term1+term2-term3+term4-mass*mass*f
    ####################################################
    # Dirac equation for a massive spinor
    # as defined in S. Chandrasekhar, "Mathematical Theory of Black Holes", 
    # Oxford Univ. Press, New York (1983), p.544.
    # We used an overall minus sign due to our definitions.
    # Set as metricalias.diracNP=0 where metricalias.diracNP is an array containing four equations.
    # Example:
    # var('M2') # mass
    # f1 = manifoldalias.scalar_field(function('f1')(t,r,th,ph))
    # f2 = manifoldalias.scalar_field(function('f2')(t,r,th,ph))
    # g1 = manifoldalias.scalar_field(function('g1')(t,r,th,ph))
    # g2 = manifoldalias.scalar_field(function('g2')(t,r,th,ph))
    # metricalias.dirac(f1,f2,g1,g2,M2)
    # show(metricalias.diracNP)
    # show(metricalias.diracNP[0].expr())
    # show(metricalias.diracNP[1].expr())
    # show(metricalias.diracNP[2].expr())
    # show(metricalias.diracNP[3].expr())
    def dirac(self,f1,f2,g1,g2,mass):
        dirac1=-(self.DlNP(f1)+(self.epsilonNP-self.rhoNP)*f1+self.deltambarNP(f2)+(self.piNP-self.alphaNP)*f2)-I*mass*g1
        dirac2=-(self.DeltanNP(f2)+(self.muNP-self.gammaNP)*f2+self.deltamNP(f1)+(self.betaNP-self.tauNP)*f1)-I*mass*g2
        dirac3=-(self.DlNP(g2)+(self.epsilonbarNP-self.rhobarNP)*g2-self.deltamNP(g1)-(self.pibarNP-self.alphabarNP)*g1)-I*mass*f2
        dirac4=-(self.DeltanNP(g1)+(self.mubarNP-self.gammabarNP)*g1-self.deltambarNP(g2)-(self.betabarNP-self.taubarNP)*g2)-I*mass*f1
        self.diracNP=[dirac1,dirac2,dirac3,dirac4]
    ####################################################


        
# SageNP: Newman-Penrose calculations for SageMath. 

This work is funded by TUBITAK 1001 Program, Grant Number 123R114.

The class **SageNP** includes functions for some calculations defined in the Newman-Penrose formalism. The code is based on SageManifolds.

# Coded by:                                        

- [Tolga Birkandan](https://web.itu.edu.tr/birkandant/) (Corr.: birkandant@itu.edu.tr)

- [Onur Arman](https://www.linkedin.com/in/onur-arman-709478337/)

- [Emir Baysazan](https://scholar.google.com/citations?user=kq9ia_oAAAAJ&hl=en)

- [Selinay Sude Binici](https://scholar.google.com/citations?user=OjiOpogAAAAJ&hl=tr)

- [Pelin Ozturk](https://www.linkedin.com/in/pelin-%C3%B6zt%C3%BCrk-3904572b2/?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=ios_app)

- Efe Özyürek

- **Special thanks to [Eric Gourgoulhon](https://luth.obspm.fr/~luthier/gourgoulhon/en/)**


# FILES:

- **SageNP.py**: Main file to import in SageMath.

- **[SageNP_Tutorial.ipynb](https://github.com/tbirkandan/SageNP/blob/main/Notebooks/SageNP_Tutorial.ipynb)**: Tutorial (ipynb file) - Definitions and calculations for the Schwarzschild (with covariant null-tetrad vectors) and Reissner-Nordstrom (with contravariant null-tetrad vectors) spacetimes.

- **[SageNP_Tutorial.pdf](https://github.com/tbirkandan/SageNP/blob/main/Notebooks/SageNP_Tutorial.pdf)**: Tutorial (PDF file)  


# REFERENCE:

Main reference for all definitions and calculations:

H. Stephani, D. Kramer, M. MacCallum, C. Hoenselaers, and E. Herlt, "Exact Solutions of Einstein’s Field Equations", 2nd ed. Cambridge: Cambridge University Press, 2003.


# BASIC DEFINITIONS AND NOTATION:

* We will use the Metric signature: (- + + +)                      

* For the null-tetrad vector names, the ref. book uses (k,l,m,mbar). However, in the code we will use (l,n,m,mbar) like the rest of the literature. Therefore, one should set k->l, l->n in the ref. book

- Products of the vectors are given by: l*n = -1, m*mbar = 1, all others zero.  
  
- The metric is found using the covariant null-tetrad vectors as:
  
    g = -2*l*n + 2*m*mbar
                            
  and,

    g = [[0  1  0  0], [1  0  0  0], [0  0  0 -1], [0  0 -1  0]]

- Please check the reference book for the details and further definitions.

# INSTRUCTIONS WITH EXAMPLES:

- **Import the class:**
  
    from SageNP import NewmanPenrose
    
- **Define your manifold:**
  
    MyManifold = Manifold(4 , 'MyManifold', r'\mathcal{Man}')

- **Define your coordinates:**
  
    MyCoordinates.<t,r,th,ph> = MyManifold.chart(r't r th:\theta ph:\phi')

- **Define the metric functions (if needed):**

    var('M')
    Delta=r^2-2*M*r
    
- **Enter null tetrad elements:**

    lvec=[1,-(r^2)/Delta,0,0]
  
    nvec=[Delta/(2*r^2),1/2,0,0]

    mvec=[0,0,(-r/sqrt(2)),(-I*r/sqrt(2))*sin(th)]

    mbarvec=[0,0,(-r/sqrt(2)),(I*r/sqrt(2))*sin(th)]

    - *Here, the element ordering is the same as the coordinate ordering. (The first element is the t element, the second is the r element, etc.)*

- **Define an object of the class:**
  
    schw=SageNP(MyManifold,MyCoordinates,lvec,nvec,mvec,mbarvec,'covariant')
  
  - *Here, our null-tetrad vectors lvec, nvec, mvec and mbarvec are covariant. Thus we used the keyword 'covariant'.*
  
  - *If they were contravariant, then we should use the keyword 'contravariant'.*

- **Once the object is defined, the code calculates the metric and displays it on the screen. It is recommended that you check your metric.**

# FUNCTIONS:

- *All page and equation numbers belong to the reference book.*


- **test_nulltetrad()**: Checks the products of the vectors l*n = -1, m*mbar = 1, all others zero.


- **Spin coefficients (Page 75-76, Eq.(7.2))**:
    
    - **calculate_spincoefficients()**: Calculates the spin coefficients.
    
    - **show_spincoefficients()**: Displays the spin coefficients
    
    - All spin coefficients are available under their names:
      
      **kappaNP, kappabarNP, tauNP, taubarNP, sigmaNP, sigmabarNP,
      rhoNP, rhobarNP, piNP, pibarNP, nuNP, nubarNP, muNP, mubarNP,
      lambdaNP, lambdabarNP, epsilonNP, epsilonbarNP, gammaNP, gammabarNP,
      betaNP, betabarNP, alphaNP, alphabarNP**

      
- **Directional derivatives (Page 43, Eq.(3.82))**:

    - **DlNP(X)**: Given X, calculates the D derivative (l direction).

    - **DeltanNP(X)**: Given X, calculates the Delta derivative (n direction)
    
    - **deltamNP(X)**: Given X, calculates the delta derivative (m direction)
    
    - **deltambarNP(X)**: Given X, calculates the deltabar derivative (mbar direction)


- **Commutators (Page 77, Eq.(7.6))**:

    - *The right-hand sides of the commutation relations are calculated.*
    
    - **Deltan_Dl_commNP(X)**: Given X, calculates the [Delta,D] commutator.
    
    - **deltam_Dl_commNP(X)**: Given X, calculates the [delta,D] commutator.
    
    - **deltam_Deltan_commNP(X)**: Given X, calculates the [delta,Delta] commutator.
    
    - **deltambar_deltam_commNP(X)**: Given X, calculates the [deltabar,delta] commutator.


- **Weyl tensor components (Page 38, Eq.(3.59))**:
    
    - **calculate_Weyl()**: Calculates the Weyl tensor components.
    
    - **show_Weyl()**: Displays the Weyl tensor components.
    
    - All Weyl tensor components are available under their names:
      
      **Psi0NP, Psi1NP, Psi2NP, Psi3NP, Psi4NP**


- **Ricci components (Page 78, Eq.(7.10-7.15))**:
    
    - **calculate_Ricci()**: Calculates the Ricci tensor components.
    
    - **show_Ricci()**: Displays the Ricci tensor components.
    
    - All Ricci tensor components are available under their names:
      
      **Phi00NP, Phi01NP, Phi10NP, Phi02NP, Phi20NP, 
      Phi11NP, Phi12NP, Phi21NP, Phi22NP, LambdaNP**


- **Ricci (Newman-Penrose) equations (Page 79, Eq.(7.21))**:
    
    - *All Newman-Penrose equations are defined as 0 = -(left hand side)+(right hand side) of the equations.*
    
    - **calculate_NPeq()**: Calculates the Newman-Penrose equations
    
    - **show_NPeq()**: Displays the Newman-Penrose equations
    
    - All Newman-Penrose equations are available under their names
      in the order they are given in the reference:
  
      **NPeq1, NPeq2, NPeq3, NPeq4, NPeq5, NPeq6, NPeq7, NPeq8, NPeq9, NPeq10, 
      NPeq11, NPeq12, NPeq13, NPeq14, NPeq15, NPeq16, NPeq17, NPeq18**


- **Bianchi identities (Page 81, Eq.(7.32))**:
    
    - *All Bianchi identities are defined as 0 = -(left hand side)+(right hand side) of the equations.*
    
    - **calculate_Bianchi()**: Calculates the Bianchi identities
    
    - **show_Bianchi()**: Displays the Bianchi identities
    
    - All Bianchi identities are available under their names
      in the order they are given in the reference:
      
      **BI1, BI2, BI3, BI4, BI5, BI6, BI7, BI8, BI9, BI10, BI11**


- **Petrov invariants I, J, K, L, N (Kramer p.121, 9.6; p.54, 4.19)**:

  *(also check diagram Fig. 9.1 on p. 122)*
  
    - **calculate_PetrovinvINP()**: Calculates the Petrov invariant I
  
    - **calculate_PetrovinvJNP()**: Calculates the Petrov invariant J
    
    - **calculate_PetrovinvKNP()**: Calculates the Petrov invariant K
    
    - **calculate_PetrovinvLNP()**: Calculates the Petrov invariant L
    
    - **calculate_PetrovinvNNP()**: Calculates the Petrov invariant N

    - All Petrov invariants are available under their names:
      
      **PetrovinvINP, PetrovinvJNP, PetrovinvKNP, PetrovinvLNP, PetrovinvNNP**


- **Petrov type of the spacetime**:
    - **Petrov_frominvariants()**: Calculates the Petrov type using I, J, K, L, N.
      
    - **Petrov_fromWeyl()**: Calculates the Petrov type using the Weyl components

	
- **Massive Klein-Gordon equation**:

	**kleingordon(Phi,M2)**: Calculates the Klein-Gordon equation for a massive scalar field where Phi is a scalar field on the manifold and M2 is the mass of the scalar field. (Ref.: G. Silva-Ortigoza, Rev. Mex. Fis. 4, 543 (1996))
	
	- The result is available under the name: **kgNP**


- **Massive Dirac equation**:

	**dirac(f1,f2,g1,g2,M2)**: Calculates the Dirac equation for a massive spinor field where f1,f2,g1,g2 are the components of the spinor field (defined as scalar fields on the manifold) and M2 is the mass of the spinor field. (Ref.: S. Chandrasekhar, "Mathematical Theory of Black Holes", Oxford Univ. Press, New York (1983), p.544.)
	
	- The result is available under the name: **diracNP**

	
- **SL(2,C) Transformations**: SL(2,C) transformations as defined in Carmeli and Kaye, Annals of Physics 99, 188 (1976).

	**type_A_transformation(z)**: Calculates the Type A transformations where z is a complex variable.

	**show_type_A_transformation()**: Shows the results of the Type A transformations.

	**type_B_transformation(z)**: Calculates the Type B transformations where z is a complex variable.

	**show_type_B_transformation()**: Shows the results of the Type B transformations.

	**type_C_transformation(z)**: Calculates the Type C transformations where z is a complex variable.

	**show_type_C_transformation()**: Shows the results of the Type C transformations.
	
	The results are available under their names:
	**lNP_trA, nNP_trA, mNP_trA, mbarNP_trA, kappaNP_trA, rhoNP_trA, etc., Psi0NP_trA, Psi1NP_trA, etc., Phi00NP_trA, Phi01NP_trA, etc.** and the same notation for **lNP_trB, nNP_trB, etc.** and **lNP_trC, nNP_trC, etc.** for other types.


- **calculate_allNP()**: Runs the following functions:
  
    - calculate_spincoefficients()

    - calculate_Weyl()

    - calculate_Ricci()

    - calculate_NPeq()

    - calculate_Bianchi()

    - Petrov_fromWeyl()

- **show_allNP()**: Runs the following functions:

    - show_spincoefficients()

    - show_Weyl()

    - show_Ricci()

    - show_NPeq()

    - show_Bianchi()

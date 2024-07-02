####################################################
# Newman-Penrose Calculations                      # 
####################################################
# Code by:                                         #
# Tolga Birkandan (Corr.: birkandant@itu.edu.tr)   #
# Emir Baysazan                                    #
# Pelin Ozturk                                     #
# Based on SageMath (SageManifolds)                #
####################################################
# Reference:                                       #
# Stephani H., Kramer D., MacCallum, M. A. H.,     #
# Hoenselaers C., Herlt E.,                        #
# "Exact Solutions of Einstein's Field Equations", #
# 2nd ed. (2003), Cambridge,                       #
# Cambridge University Press.                      #
####################################################
# NOTATION:                                        #
# Metric signature: (- + + +)                      #
# Ref. book uses (k,l,m,mbar)                      #
# The code uses  (l,n,m,mbar)                      #
# Therefore set k->l, l->n in the ref. book        #
# l*n = -1, m*mbar = 1                             #
# g = -2*l*n + 2*m*mbar                            #
# g = [0  1  0  0]                                 #
#     [1  0  0  0]                                 #
#     [0  0  0 -1]                                 #
#     [0  0 -1  0]                                 #
####################################################
##########################################################################
# For a wider screen:
#from IPython.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))
# Turn off Deprecation warnings:
#import warnings
#warnings.filterwarnings("ignore", category=DeprecationWarning)

#reset()
##########################################################################
# Define 4-dim. the manifold "Man":
Man = Manifold(4, 'Man', r'\mathcal{Man}')
##########################################################################
# Let us define our own simplification routine:
def simplify_fullfull(theinput):
    #This routine is used at some points in the functions 
    return theinput.expand().canonicalize_radical().simplify_rectform().simplify_full()
##########################################################################
# Inverting tetrad elements to covariant form
# Input (global): lNPvecup, nNPvecup,mNPvecup,mbarNPvecup
# Output (global): lNPvec,nNPvec,mNPvec,mbarNPvec
def invert_tetradNP():
    print("Inverting tetrad...")
    global lNPvec,nNPvec,mNPvec,mbarNPvec
    lNPvec=vector(SR, len(CO[:]))
    nNPvec=vector(SR, len(CO[:]))
    mNPvec=vector(SR, len(CO[:]))
    mbarNPvec=vector(SR, len(CO[:]))
    var('l1,l2,l3,l4,n1,n2,n3,n4,m1r,m2r,m3r,m4r,m1i,m2i,m3i,m4i')
    lNPvec=[l1,l2,l3,l4]
    nNPvec=[n1,n2,n3,n4]
    mNPvec=[m1r+I*m1i,m2r+I*m2i,m3r+I*m3i,m4r+I*m4i]
    mbarNPvec=[m1r-I*m1i,m2r-I*m2i,m3r-I*m3i,m4r-I*m4i]
    
    ll=0
    nn=0
    mm=0
    mbarmbar=0

    lm=0
    lmbar=0
    ml=0
    mbarl=0

    nm=0
    nmbar=0
    mn=0
    mbarn=0

    ln=1
    nl=1
    mmbar=-1
    mbarm=-1

    for i in range(4):
        ll=ll+lNPvecup[i]*lNPvec[i]
        nn=nn+nNPvecup[i]*nNPvec[i]
        mm=mm+mNPvecup[i]*mNPvec[i]
        mbarmbar=mbarmbar+mbarNPvecup[i]*mbarNPvec[i]
        lm=lm+lNPvecup[i]*mNPvec[i]
        lmbar=lmbar+lNPvecup[i]*mbarNPvec[i]
        ml=ml+mNPvecup[i]*lNPvec[i]
        mbarl=mbarl+mbarNPvecup[i]*lNPvec[i]
        nm=nm+nNPvecup[i]*mNPvec[i]
        nmbar=nmbar+nNPvecup[i]*mbarNPvec[i]
        mn=mn+mNPvecup[i]*nNPvec[i]
        mbarn=mbarn+mbarNPvecup[i]*nNPvec[i]
        ln=ln+lNPvecup[i]*nNPvec[i]
        nl=nl+nNPvecup[i]*lNPvec[i]
        mmbar=mmbar+mNPvecup[i]*mbarNPvec[i]
        mbarm=mbarm+mbarNPvecup[i]*mNPvec[i]
    
    covtetrad=solve([ll==0,nn==0,mm==0,mbarmbar==0,lm==0,lmbar==0,ml==0,mbarl==0,nm==0,nmbar==0,mn==0,mbarn==0,ln==0,nl==0,mmbar==0,mbarm==0],[l1,l2,l3,l4,n1,n2,n3,n4,m1r,m1i,m2r,m2i,m3r,m3i,m4r,m4i])
    for i in range(4):
        lNPvec[i]=simplify_fullfull(lNPvec[i].subs(covtetrad[0]))
        nNPvec[i]=simplify_fullfull(nNPvec[i].subs(covtetrad[0]))
        mNPvec[i]=simplify_fullfull(mNPvec[i].subs(covtetrad[0]))
        mbarNPvec[i]=simplify_fullfull(mbarNPvec[i].subs(covtetrad[0]))
 
    """
    #Using the inverse metric:
    inversemetric=matrix(SR, len(CO[:]), len(CO[:]))
    metricforinversion=matrix(SR, len(CO[:]), len(CO[:]))
    #print("Calculating inverse metric...")
    for i in range(len(CO[:])):
        for j in range(len(CO[:])):
            inversemetric[i,j]=simplify_fullfull((-lNPvecup[i]*nNPvecup[j]-nNPvecup[i]*lNPvecup[j]
                                +mNPvecup[i]*mbarNPvecup[j]+mbarNPvecup[i]*mNPvecup[j]))
    #print("Ters metrik hesabi bitti, metrige ters alma islemiyle geciliyor")
    metricforinversion=inversemetric.inverse()
    for i in range(len(CO[:])):
        for j in range(len(CO[:])):
            metricforinversion[i,j]=simplify_fullfull(metricforinversion[i,j])
            
    #print("Metrige gecildi, metrikten kovaryant tetrad olusturuluyor")
    for i in range(len(CO[:])):
        lNPvec[i]=0
        nNPvec[i]=0
        mNPvec[i]=0
        mbarNPvec[i]=0
        for j in range(len(CO[:])):
            lNPvec[i]=(lNPvec[i]+metricforinversion[i,j]*lNPvecup[j])
            nNPvec[i]=(nNPvec[i]+metricforinversion[i,j]*nNPvecup[j])
            mNPvec[i]=(mNPvec[i]+metricforinversion[i,j]*mNPvecup[j])
            mbarNPvec[i]=(mbarNPvec[i]+metricforinversion[i,j]*mbarNPvecup[j])
    """        
##########################################################################
file_name = str(inputfile)+'.sage'
load(file_name)
##########################################################################
##########################################################################
######                 DO  NOT CHANGE THE REST                      ######
##########################################################################
##########################################################################
# We need tensor fields for index calculations
lNP=Man.tensor_field(0,1,lNPvec)
nNP=Man.tensor_field(0,1,nNPvec)
mNP=Man.tensor_field(0,1,mNPvec)
mbarNP=Man.tensor_field(0,1,mbarNPvec)
##########################################################################
# Metric components:
# Metric "g" on manifold "Man":
g = Man.lorentzian_metric('g')

#print("Calculating metric...")
for i in range(len(CO[:])):
    for j in range(len(CO[:])):
        g[i,j]=simplify_fullfull((-lNP[i].expr()*nNP[j].expr()
                                  -nNP[i].expr()*lNP[j].expr()
                                  +mNP[i].expr()*mbarNP[j].expr()
                                  +mbarNP[i].expr()*mNP[j].expr()))
# Calculate the inverse metric:
try:
    lNPup=Man.tensor_field(0,1,lNPvecup)
    nNPup=Man.tensor_field(0,1,nNPvecup)
    mNPup=Man.tensor_field(0,1,mNPvecup)
    mbarNPup=Man.tensor_field(0,1,mbarNPvecup)
    inversegmat=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    inverseg=Man.tensor_field(0,2,inversegmat)
    for i in range(len(CO[:])):
        for j in range(len(CO[:])):
            inverseg[i,j]=simplify_fullfull((-lNPup[i].expr()*nNPup[j].expr()
                                             -nNPup[i].expr()*lNPup[j].expr()
                                             +mNPup[i].expr()*mbarNPup[j].expr()
                                             +mbarNPup[i].expr()*mNPup[j].expr()))    
except:
    inverseg = g.inverse()
#####################
# Display the metric
show('The metric:')
show(g.display())
#####################
##########################################################################
try:
    lupNP=Man.tensor_field(1,0,lNPvecup)
    nupNP=Man.tensor_field(1,0,nNPvecup)
    mupNP=Man.tensor_field(1,0,mNPvecup)
    mbarupNP=Man.tensor_field(1,0,mbarNPvecup)
except:
    lupNP=lNP.up(g)
    nupNP=nNP.up(g)
    mupNP=mNP.up(g)
    mbarupNP=mbarNP.up(g)
##########################################################################    
nab = g.connection()
##########################################################################
def test_nulltetrad():
    # All must be zero
    print("Testing null tetrad...")
    testvectorNP=vector(SR,10)
    testvectorNP[0]=simplify_fullfull((lupNP['^a']*mNP['_a']).expr())
    testvectorNP[1]=simplify_fullfull((lupNP['^a']*mbarNP['_a']).expr())
    testvectorNP[2]=simplify_fullfull((nupNP['^a']*mNP['_a']).expr())
    testvectorNP[3]=simplify_fullfull((nupNP['^a']*mbarNP['_a']).expr())
    testvectorNP[4]=simplify_fullfull((lupNP['^a']*lNP['_a']).expr())
    testvectorNP[5]=simplify_fullfull((nupNP['^a']*nNP['_a']).expr())
    testvectorNP[6]=simplify_fullfull((mupNP['^a']*mNP['_a']).expr())
    testvectorNP[7]=simplify_fullfull((mbarupNP['^a']*mbarNP['_a']).expr())
    testvectorNP[8]=simplify_fullfull((lupNP['^a']*nNP['_a']).expr()+1)
    testvectorNP[9]=simplify_fullfull((mupNP['^a']*mbarNP['_a']).expr()-1)
    if sum(testvectorNP)==0:
        print("PASSED")
    else:
        print("FAILED! Please check your null tetrad elements!")
        for i in range(len(testvectorNP)):
            if not(testvectorNP[i]==0):
                print("Failed test number: ",i)
##########################################################################
# Spin coefficients:
# Page 75-76, Eq.(7.2)
spincoeffscalculated=0
kappaNP=Man.scalar_field(function('kappaNP')(*CO));kappabarNP=Man.scalar_field(function('kappabarNP')(*CO))
tauNP=Man.scalar_field(function('tauNP')(*CO));taubarNP=Man.scalar_field(function('taubarNP')(*CO))
sigmaNP=Man.scalar_field(function('sigmaNP')(*CO));sigmabarNP=Man.scalar_field(function('sigmabarNP')(*CO))
rhoNP=Man.scalar_field(function('rhoNP')(*CO));rhobarNP=Man.scalar_field(function('rhobarNP')(*CO))
piNP=Man.scalar_field(function('piNP')(*CO));pibarNP=Man.scalar_field(function('pibarNP')(*CO))
nuNP=Man.scalar_field(function('nuNP')(*CO));nubarNP=Man.scalar_field(function('nubarNP')(*CO))
muNP=Man.scalar_field(function('muNP')(*CO));mubarNP=Man.scalar_field(function('mubarNP')(*CO))
lambdaNP=Man.scalar_field(function('lambdaNP')(*CO));lambdabarNP=Man.scalar_field(function('lambdabarNP')(*CO))
epsilonNP=Man.scalar_field(function('epsilonNP')(*CO));epsilonbarNP=Man.scalar_field(function('epsilonbarNP')(*CO))
gammaNP=Man.scalar_field(function('gammaNP')(*CO));gammabarNP=Man.scalar_field(function('gammabarNP')(*CO))
betaNP=Man.scalar_field(function('betaNP')(*CO));betabarNP=Man.scalar_field(function('betabarNP')(*CO))
alphaNP=Man.scalar_field(function('alphaNP')(*CO));alphabarNP=Man.scalar_field(function('alphabarNP')(*CO))

def calculate_spincoefficients():
    global spincoeffscalculated
    global kappaNP,tauNP,sigmaNP,rhoNP,piNP,nuNP,muNP,lambdaNP,epsilonNP,gammaNP,betaNP,alphaNP
    global kappabarNP,taubarNP,sigmabarNP,rhobarNP,pibarNP,nubarNP,mubarNP,lambdabarNP,epsilonbarNP,gammabarNP,betabarNP,alphabarNP
    print("Calculating spin coefficients...")
    spincoeffscalculated=1
    kappaNP=-nab(lNP)['_ij']*(mupNP*lupNP)['^ij']
    tauNP=-nab(lNP)['_ij']*(mupNP*nupNP)['^ij']
    sigmaNP=-nab(lNP)['_ij']*(mupNP*mupNP)['^ij']
    rhoNP=-nab(lNP)['_ij']*(mupNP*mbarupNP)['^ij']
    piNP=nab(nNP)['_ij']*(mbarupNP*lupNP)['^ij']
    nuNP=nab(nNP)['_ij']*(mbarupNP*nupNP)['^ij']
    muNP=nab(nNP)['_ij']*(mbarupNP*mupNP)['^ij']
    lambdaNP=nab(nNP)['_ij']*(mbarupNP*mbarupNP)['^ij']
    epsilonNP=-(1/2)*(nab(lNP)['_ij']*(nupNP*lupNP)['^ij']-nab(mNP)['_ij']*(mbarupNP*lupNP)['^ij'])
    gammaNP=(1/2)*(nab(nNP)['_ij']*(lupNP*nupNP)['^ij']-nab(mbarNP)['_ij']*(mupNP*nupNP)['^ij']) 
    betaNP=-(1/2)*(nab(lNP)['_ij']*(nupNP*mupNP)['^ij']-nab(mNP)['_ij']*(mbarupNP*mupNP)['^ij'])
    alphaNP=(1/2)*(nab(nNP)['_ij']*(lupNP*mbarupNP)['^ij']-nab(mbarNP)['_ij']*(mupNP*mbarupNP)['^ij']) 
    kappabarNP=-nab(lNP)['_ij']*(mbarupNP*lupNP)['^ij']
    taubarNP=-nab(lNP)['_ij']*(mbarupNP*nupNP)['^ij']
    sigmabarNP=-nab(lNP)['_ij']*(mbarupNP*mbarupNP)['^ij']
    rhobarNP=-nab(lNP)['_ij']*(mbarupNP*mupNP)['^ij']
    pibarNP=nab(nNP)['_ij']*(mupNP*lupNP)['^ij']
    nubarNP=nab(nNP)['_ij']*(mupNP*nupNP)['^ij']
    mubarNP=nab(nNP)['_ij']*(mupNP*mbarupNP)['^ij']
    lambdabarNP=nab(nNP)['_ij']*(mupNP*mupNP)['^ij']
    epsilonbarNP=-(1/2)*(nab(lNP)['_ij']*(nupNP*lupNP)['^ij']-nab(mbarNP)['_ij']*(mupNP*lupNP)['^ij'])
    gammabarNP=(1/2)*(nab(nNP)['_ij']*(lupNP*nupNP)['^ij']-nab(mNP)['_ij']*(mbarupNP*nupNP)['^ij'])
    betabarNP=-(1/2)*(nab(lNP)['_ij']*(nupNP*mbarupNP)['^ij']-nab(mbarNP)['_ij']*(mupNP*mbarupNP)['^ij'])
    alphabarNP=(1/2)*(nab(nNP)['_ij']*(lupNP*mupNP)['^ij']-nab(mNP)['_ij']*(mbarupNP*mupNP)['^ij'])

def show_spincoefficients():
    if spincoeffscalculated==1:
        show("kappaNP=",kappaNP.expr())
        show("tauNP=",tauNP.expr())
        show("sigmaNP=",sigmaNP.expr())
        show("rhoNP=",rhoNP.expr())
        show("piNP=",piNP.expr())
        show("nuNP=",nuNP.expr())
        show("muNP=",muNP.expr())
        show("lambdaNP=",lambdaNP.expr())
        show("epsilonNP=",epsilonNP.expr())
        show("gammaNP=",gammaNP.expr())
        show("betaNP=",betaNP.expr())
        show("alphaNP=",alphaNP.expr())
    else:
        calculate_spincoefficients()
        show_spincoefficients()
##########################################################################
# Directional Derivatives
# Page 43, Eq.(3.82)
def DlNP(X):
    if X==0:
        return 0
    else:
        return lupNP['^i']*nab(X)['_i']
def DeltanNP(X):
    if X==0:
        return 0
    else:
        return nupNP['^i']*nab(X)['_i']
def deltamNP(X):
    if X==0:
        return 0
    else:
        return mupNP['^i']*nab(X)['_i']
def deltambarNP(X):
    if X==0:
        return 0
    else:
        return mbarupNP['^i']*nab(X)['_i']
##########################################################################
# Commutators
# Page 77, Eq.(7.6)
def Deltan_Dl_commNP(X):
    comresult=(gammaNP+gammabarNP)*DlNP(X)+(epsilonNP+epsilonbarNP)*DeltanNP(X)-(tauNP+pibarNP)*deltambarNP(X)-(taubarNP+piNP)*deltamNP(X)
    return comresult
def deltam_Dl_commNP(X):
    comresult=(alphabarNP+betaNP-pibarNP)*DlNP(X)+kappaNP*DeltanNP(X)-sigmaNP*deltambarNP(X)-(rhobarNP+epsilonNP-epsilonbarNP)*deltamNP(X)
    return comresult
def deltam_Deltan_commNP(X):
    comresult=-nubarNP*DlNP(X)+(tauNP-alphabarNP-betaNP)*DeltanNP(X)+lambdabarNP*deltambarNP(X)+(muNP-gammaNP+gammabarNP)*deltamNP(X)
    return comresult
def deltambar_deltam_commNP(X):
    comresult=(mubarNP-muNP)*DlNP(X)+(rhobarNP-rhoNP)*DeltanNP(X)-(alphabarNP-betaNP)*deltambarNP(X)-(betabarNP-alphaNP)*deltamNP(X)
    return comresult
##########################################################################
# Weyl tensor components:
# Page 38, Eq.(3.59)
weylcalculated=0
Psi0NP=Man.scalar_field(function('Psi0NP')(*CO))
Psi1NP=Man.scalar_field(function('Psi1NP')(*CO))
Psi2NP=Man.scalar_field(function('Psi2NP')(*CO)) 
Psi3NP=Man.scalar_field(function('Psi3NP')(*CO))
Psi4NP=Man.scalar_field(function('Psi4NP')(*CO))

def calculate_Weyl():
    global weylcalculated,Psi0NP,Psi1NP,Psi2NP,Psi3NP,Psi4NP
    print("Calculating Weyl components...")
    weylcalculated=1
    C = g.weyl()
    C_allindicesdown=C.down(g)
    Psi0NP=C_allindicesdown['_{pqrs}']*(lupNP*mupNP*lupNP*mupNP)['^{pqrs}']
    Psi1NP=C_allindicesdown['_{pqrs}']*(lupNP*nupNP*lupNP*mupNP)['^{pqrs}']
    Psi2NP=C_allindicesdown['_{pqrs}']*(lupNP*mupNP*mbarupNP*nupNP)['^{pqrs}']
    Psi3NP=C_allindicesdown['_{pqrs}']*(lupNP*nupNP*mbarupNP*nupNP)['^{pqrs}']
    Psi4NP=C_allindicesdown['_{pqrs}']*(mbarupNP*nupNP*mbarupNP*nupNP)['^{pqrs}']
    
def show_Weyl():
    if weylcalculated==1:
        show("Psi0NP=",Psi0NP.expr())
        show("Psi1NP=",Psi1NP.expr())
        show("Psi2NP=",Psi2NP.expr())
        show("Psi3NP=",Psi3NP.expr())
        show("Psi4NP=",Psi4NP.expr())
    else:
        calculate_Weyl()
        show_Weyl()
##########################################################################
# Ricci components:
# Page 78, Eq.(7.10-7.15)
riccicalculated=0
Phi00NP=Man.scalar_field(function('Phi00NP')(*CO))
Phi01NP=Man.scalar_field(function('Phi01NP')(*CO))
Phi10NP=Man.scalar_field(function('Phi10NP')(*CO))
Phi02NP=Man.scalar_field(function('Phi02NP')(*CO))
Phi20NP=Man.scalar_field(function('Phi20NP')(*CO))
Phi11NP=Man.scalar_field(function('Phi11NP')(*CO))
Phi12NP=Man.scalar_field(function('Phi12NP')(*CO))
Phi21NP=Man.scalar_field(function('Phi21NP')(*CO))
Phi22NP=Man.scalar_field(function('Phi22NP')(*CO))
LambdaNP=Man.scalar_field(function('LambdaNP')(*CO))

def calculate_Ricci():
    global riccicalculated,Phi00NP,Phi01NP,Phi10NP,Phi02NP,Phi20NP,Phi11NP,Phi12NP,Phi21NP,Phi22NP,LambdaNP
    print("Calculating Ricci components...")
    riccicalculated=1
    Ricciscalar=g.ricci_scalar()
    R_traceless_allindicesdown=g.ricci()['_{ab}']-(g.ricci_scalar()/4)*g['_{ab}']
    Phi00NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(lupNP*lupNP)['^{ab}']
    Phi01NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(lupNP*mupNP)['^{ab}']
    Phi10NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(lupNP*mbarupNP)['^{ab}']
    Phi02NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(mupNP*mupNP)['^{ab}']
    Phi20NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(mbarupNP*mbarupNP)['^{ab}']
    Phi11NP=(1/4)*(R_traceless_allindicesdown['_{ab}']*(lupNP*nupNP)['^{ab}']+R_traceless_allindicesdown['_{ab}']*(mupNP*mbarupNP)['^{ab}'])
    Phi12NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(nupNP*mupNP)['^{ab}']
    Phi21NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(nupNP*mbarupNP)['^{ab}']
    Phi22NP=(1/2)*R_traceless_allindicesdown['_{ab}']*(nupNP*nupNP)['^{ab}']
    LambdaNP=(1/24)*Ricciscalar # Stephani tanimlamamis, R diye kullaniyor ama bu tanim denklemlerle uyumlu
    
def show_Ricci():
    if riccicalculated==1:
        show("Phi00NP=",Phi00NP.expr())
        show("Phi01NP=",Phi01NP.expr())
        show("Phi10NP=",Phi10NP.expr())
        show("Phi02NP=",Phi02NP.expr())
        show("Phi20NP=",Phi20NP.expr())
        show("Phi11NP=",Phi11NP.expr())
        show("Phi12NP=",Phi12NP.expr())
        show("Phi21NP=",Phi21NP.expr())
        show("Phi22NP=",Phi22NP.expr())
        show("LambdaNP=",LambdaNP.expr())
    else:
        calculate_Ricci()
        show_Ricci()
##########################################################################
# Ricci (Newman-Penrose) Equations
# Page 79, Eq.(7.21)
def calculate_NPeq():
    global NPeq1,NPeq2,NPeq3,NPeq4,NPeq5,NPeq6,NPeq7,NPeq8,NPeq9,NPeq10,NPeq11,NPeq12,NPeq13,NPeq14,NPeq15,NPeq16,NPeq17,NPeq18
    print("Calculating NP equations...")
    NPeq1=-DlNP(rhoNP)+deltambarNP(kappaNP)+(rhoNP^2+sigmaNP*sigmabarNP)+(epsilonNP+epsilonbarNP)*rhoNP-kappabarNP*tauNP-kappaNP*(3*alphaNP+betabarNP-piNP)+Phi00NP
    NPeq2=-DlNP(sigmaNP)+deltamNP(kappaNP)+sigmaNP*(3*epsilonNP-epsilonbarNP+rhoNP+rhobarNP)+kappaNP*(pibarNP-tauNP-3*betaNP-alphabarNP)+Psi0NP
    NPeq3=-DlNP(tauNP)+DeltanNP(kappaNP)+(tauNP+pibarNP)*rhoNP+(taubarNP+piNP)*sigmaNP+(epsilonNP-epsilonbarNP)*tauNP-(3*gammaNP+gammabarNP)*kappaNP+Psi1NP+Phi01NP
    NPeq4=-DlNP(alphaNP)+deltambarNP(epsilonNP)+(rhoNP+epsilonbarNP-2*epsilonNP)*alphaNP+betaNP*sigmabarNP-betabarNP*epsilonNP-kappaNP*lambdaNP-kappabarNP*gammaNP+(epsilonNP+rhoNP)*piNP+Phi10NP
    NPeq5=-DlNP(betaNP) + deltamNP(epsilonNP) +(alphaNP + piNP)*sigmaNP + (rhobarNP - epsilonbarNP)*betaNP - (muNP + gammaNP)*kappaNP - (alphabarNP - pibarNP)*epsilonNP + Psi1NP
    NPeq6=-DlNP(gammaNP) + DeltanNP(epsilonNP) + (tauNP + pibarNP)*alphaNP + (taubarNP+piNP)*betaNP - (epsilonNP + epsilonbarNP)*gammaNP - (gammaNP + gammabarNP)*epsilonNP + tauNP*piNP - nuNP*kappaNP + Psi2NP + Phi11NP - LambdaNP
    NPeq7=-DlNP(lambdaNP) + deltambarNP(piNP) + (rhoNP*lambdaNP + sigmabarNP*muNP) + piNP^2 + (alphaNP - betabarNP)*piNP - nuNP*kappabarNP - (3*epsilonNP - epsilonbarNP)*lambdaNP + Phi20NP
    NPeq8=-DlNP(muNP) + deltamNP(piNP) + (rhobarNP*muNP + sigmaNP*lambdaNP) +piNP*pibarNP - (epsilonNP + epsilonbarNP)*muNP - (alphabarNP - betaNP)*piNP - nuNP*kappaNP + Psi2NP + 2*LambdaNP
    NPeq9=-DlNP(nuNP) + DeltanNP(piNP) +(piNP + taubarNP)*muNP + (pibarNP + tauNP)*lambdaNP + (gammaNP - gammabarNP)*piNP - (3*epsilonNP + epsilonbarNP)*nuNP + Psi3NP + Phi21NP
    NPeq10=-DeltanNP(lambdaNP) + deltambarNP(nuNP) - (muNP + mubarNP)*lambdaNP - (3*gammaNP - gammabarNP)*lambdaNP + (3*alphaNP + betabarNP + piNP - taubarNP)*nuNP - Psi4NP
    NPeq11=-deltamNP(rhoNP) + deltambarNP(sigmaNP) + (alphabarNP + betaNP)*rhoNP - (3*alphaNP - betabarNP)*sigmaNP + (rhoNP - rhobarNP)*tauNP + (muNP - mubarNP)*kappaNP - Psi1NP + Phi01NP
    NPeq12=-deltamNP(alphaNP) + deltambarNP(betaNP) + (muNP*rhoNP - lambdaNP*sigmaNP) + alphaNP*alphabarNP + betaNP*betabarNP - 2*alphaNP*betaNP + (rhoNP-rhobarNP)*gammaNP + (muNP-mubarNP)*epsilonNP - Psi2NP + Phi11NP + LambdaNP
    NPeq13=-deltamNP(lambdaNP) + deltambarNP(muNP) + (rhoNP - rhobarNP)*nuNP + (muNP - mubarNP)*piNP + (alphaNP + betabarNP)*muNP + (alphabarNP - 3*betaNP)*lambdaNP - Psi3NP + Phi21NP
    NPeq14=-deltamNP(nuNP) + DeltanNP(muNP) +(muNP^2 + lambdaNP*lambdabarNP) + (gammaNP + gammabarNP)*muNP - nubarNP*piNP + (tauNP - 3*betaNP - alphabarNP)*nuNP + Phi22NP
    NPeq15=-deltamNP(gammaNP) + DeltanNP(betaNP) + (tauNP - alphabarNP - betaNP)*gammaNP + muNP*tauNP - sigmaNP* nuNP - epsilonNP*nubarNP - (gammaNP - gammabarNP - muNP)*betaNP + alphaNP*lambdabarNP + Phi12NP
    NPeq16=-deltamNP(tauNP) + DeltanNP(sigmaNP) + (muNP*sigmaNP + lambdabarNP*rhoNP) + (tauNP + betaNP - alphabarNP)*tauNP - (3*gammaNP - gammabarNP)*sigmaNP - kappaNP*nubarNP + Phi02NP
    NPeq17=-DeltanNP(rhoNP) + deltambarNP(tauNP) - (rhoNP*mubarNP + sigmaNP*lambdaNP) + (betabarNP - alphaNP - taubarNP)*tauNP + (gammaNP + gammabarNP)*rhoNP + nuNP*kappaNP - Psi2NP - 2*LambdaNP
    NPeq18=-DeltanNP(alphaNP) + deltambarNP(gammaNP) + (rhoNP + epsilonNP)*nuNP - (tauNP + betaNP)*lambdaNP + (gammabarNP- mubarNP)*alphaNP + (betabarNP - taubarNP)*gammaNP - Psi3NP

# Only the right-hand-sides of the NP equations
def calculate_NPeqrhs():
    global NPeq1rhs,NPeq2rhs,NPeq3rhs,NPeq4rhs,NPeq5rhs,NPeq6rhs,NPeq7rhs,NPeq8rhs,NPeq9rhs,NPeq10rhs,NPeq11rhs,NPeq12rhs,NPeq13rhs,NPeq14rhs,NPeq15rhs,NPeq16rhs,NPeq17rhs,NPeq18rhs
    print("Calculating right hand sides of the NP equations...")
    NPeq1rhs=(rhoNP^2+sigmaNP*sigmabarNP)+(epsilonNP+epsilonbarNP)*rhoNP-kappabarNP*tauNP-kappaNP*(3*alphaNP+betabarNP-piNP)+Phi00NP
    NPeq2rhs=sigmaNP*(3*epsilonNP-epsilonbarNP+rhoNP+rhobarNP)+kappaNP*(pibarNP-tauNP-3*betaNP-alphabarNP)+Psi0NP
    NPeq3rhs=(tauNP+pibarNP)*rhoNP+(taubarNP+piNP)*sigmaNP+(epsilonNP-epsilonbarNP)*tauNP-(3*gammaNP+gammabarNP)*kappaNP+Psi1NP+Phi01NP
    NPeq4rhs=(rhoNP+epsilonbarNP-2*epsilonNP)*alphaNP+betaNP*sigmabarNP-betabarNP*epsilonNP-kappaNP*lambdaNP-kappabarNP*gammaNP+(epsilonNP+rhoNP)*piNP+Phi10NP
    NPeq5rhs=(alphaNP + piNP)*sigmaNP + (rhobarNP - epsilonbarNP)*betaNP - (muNP + gammaNP)*kappaNP - (alphabarNP - pibarNP)*epsilonNP + Psi1NP
    NPeq6rhs=(tauNP + pibarNP)*alphaNP + (taubarNP+piNP)*betaNP - (epsilonNP + epsilonbarNP)*gammaNP - (gammaNP + gammabarNP)*epsilonNP + tauNP*piNP - nuNP*kappaNP + Psi2NP + Phi11NP - LambdaNP
    NPeq7rhs=(rhoNP*lambdaNP + sigmabarNP*muNP) + piNP^2 + (alphaNP - betabarNP)*piNP - nuNP*kappabarNP - (3*epsilonNP - epsilonbarNP)*lambdaNP + Phi20NP
    NPeq8rhs=(rhobarNP*muNP + sigmaNP*lambdaNP) +piNP*pibarNP - (epsilonNP + epsilonbarNP)*muNP - (alphabarNP - betaNP)*piNP - nuNP*kappaNP + Psi2NP + 2*LambdaNP
    NPeq9rhs=(piNP + taubarNP)*muNP + (pibarNP + tauNP)*lambdaNP + (gammaNP - gammabarNP)*piNP - (3*epsilonNP + epsilonbarNP)*nuNP + Psi3NP + Phi21NP
    NPeq10rhs=(muNP + mubarNP)*lambdaNP - (3*gammaNP - gammabarNP)*lambdaNP + (3*alphaNP + betabarNP + piNP - taubarNP)*nuNP - Psi4NP
    NPeq11rhs=(alphabarNP + betaNP)*rhoNP - (3*alphaNP - betabarNP)*sigmaNP + (rhoNP - rhobarNP)*tauNP + (muNP - mubarNP)*kappaNP - Psi1NP + Phi01NP
    NPeq12rhs=(muNP*rhoNP - lambdaNP*sigmaNP) + alphaNP*alphabarNP + betaNP*betabarNP - 2*alphaNP*betaNP + (rhoNP-rhobarNP)*gammaNP + (muNP-mubarNP)*epsilonNP - Psi2NP + Phi11NP + LambdaNP
    NPeq13rhs=(rhoNP - rhobarNP)*nuNP + (muNP - mubarNP)*piNP + (alphaNP + betabarNP)*muNP + (alphabarNP - 3*betaNP)*lambdaNP - Psi3NP + Phi21NP
    NPeq14rhs=(muNP^2 + lambdaNP*lambdabarNP) + (gammaNP + gammabarNP)*muNP - nubarNP*piNP + (tauNP - 3*betaNP - alphabarNP)*nuNP + Phi22NP
    NPeq15rhs=(tauNP - alphabarNP - betaNP)*gammaNP + muNP*tauNP - sigmaNP* nuNP - epsilonNP*nubarNP - (gammaNP - gammabarNP - muNP)*betaNP + alphaNP*lambdabarNP + Phi12NP
    NPeq16rhs=(muNP*sigmaNP + lambdabarNP*rhoNP) + (tauNP + betaNP - alphabarNP)*tauNP - (3*gammaNP - gammabarNP)*sigmaNP - kappaNP*nubarNP + Phi02NP
    NPeq17rhs=-(rhoNP*mubarNP + sigmaNP*lambdaNP) + (betabarNP - alphaNP - taubarNP)*tauNP + (gammaNP + gammabarNP)*rhoNP + nuNP*kappaNP - Psi2NP - 2*LambdaNP
    NPeq18rhs=(rhoNP + epsilonNP)*nuNP - (tauNP + betaNP)*lambdaNP + (gammabarNP- mubarNP)*alphaNP + (betabarNP - taubarNP)*gammaNP - Psi3NP
            
def show_NPeq():
    show("NPeq1=",NPeq1.expr())
    show("NPeq2=",NPeq2.expr())
    show("NPeq3=",NPeq3.expr())
    show("NPeq4=",NPeq4.expr())
    show("NPeq5=",NPeq5.expr())
    show("NPeq6=",NPeq6.expr())
    show("NPeq7=",NPeq7.expr())
    show("NPeq8=",NPeq8.expr())
    show("NPeq9=",NPeq9.expr())
    show("NPeq10=",NPeq10.expr())
    show("NPeq11=",NPeq11.expr())
    show("NPeq12=",NPeq12.expr())
    show("NPeq13=",NPeq13.expr())
    show("NPeq14=",NPeq14.expr())
    show("NPeq15=",NPeq15.expr())
    show("NPeq16=",NPeq16.expr())
    show("NPeq17=",NPeq17.expr())
    show("NPeq18=",NPeq18.expr())

def show_NPeqrhs():
    show("NPeq1rhs=",NPeq1rhs.expr())
    show("NPeq2rhs=",NPeq2rhs.expr())
    show("NPeq3rhs=",NPeq3rhs.expr())
    show("NPeq4rhs=",NPeq4rhs.expr())
    show("NPeq5rhs=",NPeq5rhs.expr())
    show("NPeq6rhs=",NPeq6rhs.expr())
    show("NPeq7rhs=",NPeq7rhs.expr())
    show("NPeq8rhs=",NPeq8rhs.expr())
    show("NPeq9rhs=",NPeq9rhs.expr())
    show("NPeq10rhs=",NPeq10rhs.expr())
    show("NPeq11rhs=",NPeq11rhs.expr())
    show("NPeq12rhs=",NPeq12rhs.expr())
    show("NPeq13rhs=",NPeq13rhs.expr())
    show("NPeq14rhs=",NPeq14rhs.expr())
    show("NPeq15rhs=",NPeq15rhs.expr())
    show("NPeq16rhs=",NPeq16rhs.expr())
    show("NPeq17rhs=",NPeq17rhs.expr())
    show("NPeq18rhs=",NPeq18rhs.expr())
##########################################################################
# Bianchi identities
# Page 81, Eq.(7.32)
def calculate_Bianchi():
    global BI1,BI2,BI3,BI4,BI5,BI6,BI7,BI8,BI9,BI10,BI11
    print("Calculating Bianchi identities...")
    BI1=-(deltambarNP(Psi0NP)-DlNP(Psi1NP)+DlNP(Phi01NP)-deltamNP(Phi00NP))+(4*alphaNP-piNP)*Psi0NP-2*(2*rhoNP+epsilonNP)*Psi1NP+3*kappaNP*Psi2NP+(pibarNP-2*alphabarNP-2*betaNP)*Phi00NP+2*(epsilonNP+rhobarNP)*Phi01NP+2*sigmaNP*Phi10NP-2*kappaNP*Phi11NP-kappabarNP*Phi02NP
    BI2=-DeltanNP(Psi0NP) + deltamNP(Psi1NP) - DlNP(Phi02NP) + deltamNP(Phi01NP) + (4*gammaNP-muNP)*Psi0NP - 2*(2*tauNP + betaNP)*Psi1NP + 3*sigmaNP*Psi2NP + (2*epsilonNP - 2*epsilonbarNP + rhobarNP)*Phi02NP + 2*(pibarNP - betaNP)*Phi01NP + 2*sigmaNP*Phi11NP - 2*kappaNP*Phi12NP - lambdabarNP*Phi00NP
    BI3=-deltambarNP(Psi3NP) + DlNP(Psi4NP) - deltambarNP(Phi21NP) + DeltanNP(Phi20NP) + (4*epsilonNP - rhoNP)*Psi4NP - 2*(2*piNP + alphaNP)*Psi3NP + 3*lambdaNP*Psi2NP + (2*gammaNP - 2*gammabarNP + mubarNP)*Phi20NP + 2*(taubarNP - alphaNP)*Phi21NP + 2*lambdaNP*Phi11NP - 2*nuNP*Phi10NP - sigmabarNP*Phi22NP
    BI4=-DeltanNP(Psi3NP) + deltamNP(Psi4NP) - deltambarNP(Phi22NP) + DeltanNP(Phi21NP) + (4*betaNP - tauNP)*Psi4NP - 2*(2*muNP + gammaNP)*Psi3NP + 3*nuNP*Psi2NP + (taubarNP - 2*betabarNP - 2*alphaNP)*Phi22NP + 2*(gammaNP+mubarNP)*Phi21NP + 2*lambdaNP*Phi12NP - 2*nuNP*Phi11NP - nubarNP*Phi20NP
    BI5=-DlNP(Psi2NP) + deltambarNP(Psi1NP) - DeltanNP(Phi00NP) + deltambarNP(Phi01NP) - 2*DlNP(LambdaNP) - lambdaNP*Psi0NP + 2*(piNP-alphaNP)*Psi1NP + 3*rhoNP*Psi2NP - 2*kappaNP*Psi3NP + (2*gammaNP + 2*gammabarNP - mubarNP)*Phi00NP - 2*(taubarNP + alphaNP)*Phi01NP - 2*tauNP*Phi10NP + 2*rhoNP*Phi11NP + sigmabarNP*Phi02NP
    BI6=-DeltanNP(Psi2NP) + deltamNP(Psi3NP) - DlNP(Phi22NP) + deltamNP(Phi21NP) - 2* DeltanNP(LambdaNP) + sigmaNP*Psi4NP + 2*(betaNP-tauNP)*Psi3NP - 3*muNP*Psi2NP + 2*nuNP*Psi1NP + (rhobarNP - 2*epsilonNP - 2*epsilonbarNP)*Phi22NP + 2*(pibarNP + betaNP)*Phi21NP + 2*piNP*Phi12NP - 2*muNP*Phi11NP - lambdabarNP*Phi20NP
    BI7=-DlNP(Psi3NP) + deltambarNP(Psi2NP) + DlNP(Phi21NP) -deltamNP(Phi20NP) + 2*deltambarNP(LambdaNP) - kappaNP*Psi4NP + 2*(rhoNP-epsilonNP)*Psi3NP + 3*piNP*Psi2NP - 2*lambdaNP*Psi1NP + (2*alphabarNP-2*betaNP - pibarNP)*Phi20NP - 2*(rhobarNP-epsilonNP)*Phi21NP - 2*piNP*Phi11NP + 2*muNP*Phi10NP + kappabarNP*Phi22NP
    BI8=-DeltanNP(Psi1NP) + deltamNP(Psi2NP) + DeltanNP(Phi01NP) - deltambarNP(Phi02NP) + 2*deltamNP(LambdaNP) + nuNP*Psi0NP + 2*(gammaNP-muNP)*Psi1NP - 3*tauNP*Psi2NP + 2*sigmaNP*Psi3NP + (taubarNP - 2*betabarNP + 2*alphaNP)*Phi02NP + 2*(mubarNP - gammaNP)*Phi01NP + 2*tauNP*Phi11NP-2*rhoNP*Phi12NP-nubarNP*Phi00NP
    BI9=-DlNP(Phi11NP) + deltamNP(Phi10NP) + deltambarNP(Phi01NP) - DeltanNP(Phi00NP) - 3*DlNP(LambdaNP) + (2*gammaNP - muNP + 2*gammabarNP - mubarNP)*Phi00NP + (piNP - 2*alphaNP - 2*taubarNP)*Phi01NP + (pibarNP - 2*alphabarNP - 2*tauNP)*Phi10NP + 2*(rhoNP+rhobarNP)*Phi11NP + sigmabarNP*Phi02NP + sigmaNP*Phi20NP - kappabarNP*Phi12NP - kappaNP*Phi21NP
    BI10=-DlNP(Phi12NP) + deltamNP(Phi11NP) + deltambarNP(Phi02NP) - DeltanNP(Phi01NP) - 3*deltamNP(LambdaNP) + (-2*alphaNP + 2*betabarNP + piNP - taubarNP)*Phi02NP + (rhobarNP + 2*rhoNP - 2*epsilonbarNP)*Phi12NP + 2*(pibarNP - tauNP)*Phi11NP + (2*gammaNP - 2*mubarNP - muNP)*Phi01NP + nubarNP*Phi00NP - lambdabarNP*Phi10NP + sigmaNP*Phi21NP - kappaNP*Phi22NP
    BI11=-DlNP(Phi22NP) + deltamNP(Phi21NP) + deltambarNP(Phi12NP) - DeltanNP(Phi11NP) - 3*DeltanNP(LambdaNP) + (rhoNP + rhobarNP - 2*epsilonNP - 2*epsilonbarNP)*Phi22NP + (2*betabarNP + 2*piNP - taubarNP)*Phi12NP + (2*betaNP + 2*pibarNP - tauNP)*Phi21NP - 2*(muNP+mubarNP)*Phi11NP + nuNP*Phi01NP + nubarNP*Phi10NP - lambdabarNP*Phi20NP - lambdaNP*Phi02NP

# Only the right-hand-sides of the Bianchi identities
def calculate_Bianchirhs():
    global BI1rhs,BI2rhs,BI3rhs,BI4rhs,BI5rhs,BI6rhs,BI7rhs,BI8rhs,BI9rhs,BI10rhs,BI11rhs
    print("Calculating the right had sides of the Bianchi identities...")
    BI1rhs=(4*alphaNP-piNP)*Psi0NP-2*(2*rhoNP+epsilonNP)*Psi1NP+3*kappaNP*Psi2NP+(pibarNP-2*alphabarNP-2*betaNP)*Phi00NP+2*(epsilonNP+rhobarNP)*Phi01NP+2*sigmaNP*Phi10NP-2*kappaNP*Phi11NP-kappabarNP*Phi02NP
    BI2rhs=(4*gammaNP-muNP)*Psi0NP - 2*(2*tauNP + betaNP)*Psi1NP + 3*sigmaNP*Psi2NP + (2*epsilonNP - 2*epsilonbarNP + rhobarNP)*Phi02NP + 2*(pibarNP - betaNP)*Phi01NP + 2*sigmaNP*Phi11NP - 2*kappaNP*Phi12NP - lambdabarNP*Phi00NP
    BI3rhs=(4*epsilonNP - rhoNP)*Psi4NP - 2*(2*piNP + alphaNP)*Psi3NP + 3*lambdaNP*Psi2NP + (2*gammaNP - 2*gammabarNP + mubarNP)*Phi20NP + 2*(taubarNP - alphaNP)*Phi21NP + 2*lambdaNP*Phi11NP - 2*nuNP*Phi10NP - sigmabarNP*Phi22NP
    BI4rhs=(4*betaNP - tauNP)*Psi4NP - 2*(2*muNP + gammaNP)*Psi3NP + 3*nuNP*Psi2NP + (taubarNP - 2*betabarNP - 2*alphaNP)*Phi22NP + 2*(gammaNP+mubarNP)*Phi21NP + 2*lambdaNP*Phi12NP - 2*nuNP*Phi11NP - nubarNP*Phi20NP
    BI5rhs=-lambdaNP*Psi0NP + 2*(piNP-alphaNP)*Psi1NP + 3*rhoNP*Psi2NP - 2*kappaNP*Psi3NP + (2*gammaNP + 2*gammabarNP - mubarNP)*Phi00NP - 2*(taubarNP + alphaNP)*Phi01NP - 2*tauNP*Phi10NP + 2*rhoNP*Phi11NP + sigmabarNP*Phi02NP
    BI6rhs=sigmaNP*Psi4NP + 2*(betaNP-tauNP)*Psi3NP - 3*muNP*Psi2NP + 2*nuNP*Psi1NP + (rhobarNP - 2*epsilonNP - 2*epsilonbarNP)*Phi22NP + 2*(pibarNP + betaNP)*Phi21NP + 2*piNP*Phi12NP - 2*muNP*Phi11NP - lambdabarNP*Phi20NP
    BI7rhs=-kappaNP*Psi4NP + 2*(rhoNP-epsilonNP)*Psi3NP + 3*piNP*Psi2NP - 2*lambdaNP*Psi1NP + (2*alphabarNP-2*betaNP - pibarNP)*Phi20NP - 2*(rhobarNP-epsilonNP)*Phi21NP - 2*piNP*Phi11NP + 2*muNP*Phi10NP + kappabarNP*Phi22NP
    BI8rhs=nuNP*Psi0NP + 2*(gammaNP-muNP)*Psi1NP - 3*tauNP*Psi2NP + 2*sigmaNP*Psi3NP + (taubarNP - 2*betabarNP + 2*alphaNP)*Phi02NP + 2*(mubarNP - gammaNP)*Phi01NP + 2*tauNP*Phi11NP-2*rhoNP*Phi12NP-nubarNP*Phi00NP
    BI9rhs=(2*gammaNP - muNP + 2*gammabarNP - mubarNP)*Phi00NP + (piNP - 2*alphaNP - 2*taubarNP)*Phi01NP + (pibarNP - 2*alphabarNP - 2*tauNP)*Phi10NP + 2*(rhoNP+rhobarNP)*Phi11NP + sigmabarNP*Phi02NP + sigmaNP*Phi20NP - kappabarNP*Phi12NP - kappaNP*Phi21NP
    BI10rhs=(-2*alphaNP + 2*betabarNP + piNP - taubarNP)*Phi02NP + (rhobarNP + 2*rhoNP - 2*epsilonbarNP)*Phi12NP + 2*(pibarNP - tauNP)*Phi11NP + (2*gammaNP - 2*mubarNP - muNP)*Phi01NP + nubarNP*Phi00NP - lambdabarNP*Phi10NP + sigmaNP*Phi21NP - kappaNP*Phi22NP
    BI11rhs=(rhoNP + rhobarNP - 2*epsilonNP - 2*epsilonbarNP)*Phi22NP + (2*betabarNP + 2*piNP - taubarNP)*Phi12NP + (2*betaNP + 2*pibarNP - tauNP)*Phi21NP - 2*(muNP+mubarNP)*Phi11NP + nuNP*Phi01NP + nubarNP*Phi10NP - lambdabarNP*Phi20NP - lambdaNP*Phi02NP

def show_Bianchi():
    show("Bianchi1=",BI1.expr())
    show("Bianchi2=",BI2.expr())
    show("Bianchi3=",BI3.expr())
    show("Bianchi4=",BI4.expr())
    show("Bianchi5=",BI5.expr())
    show("Bianchi6=",BI6.expr())
    show("Bianchi7=",BI7.expr())
    show("Bianchi8=",BI8.expr())
    show("Bianchi9=",BI9.expr())
    show("Bianchi10=",BI10.expr())
    show("Bianchi11=",BI11.expr())
    
def show_Bianchirhs():
    show("Bianchi1rhs=",BI1rhs.expr())
    show("Bianchi2rhs=",BI2rhs.expr())
    show("Bianchi3rhs=",BI3rhs.expr())
    show("Bianchi4rhs=",BI4rhs.expr())
    show("Bianchi5rhs=",BI5rhs.expr())
    show("Bianchi6rhs=",BI6rhs.expr())
    show("Bianchi7rhs=",BI7rhs.expr())
    show("Bianchi8rhs=",BI8rhs.expr())
    show("Bianchi9rhs=",BI9rhs.expr())
    show("Bianchi10rhs=",BI10rhs.expr())
    show("Bianchi11rhs=",BI11rhs.expr())
##########################################################################
# Ricci tensor in matrix notation (R{_a ^b})
def calculate_RiccimatrixNP():
    global RiccimatrixNP
    RiccimatrixNP=2*matrix([[-(Phi11NP-3*LambdaNP),-Phi00NP,Phi10NP,Phi01NP],
                          [-Phi22NP,-(Phi11NP-3*LambdaNP),Phi21NP,Phi12NP],
                          [-Phi12NP,-Phi01NP,(Phi11NP+3*LambdaNP),Phi02NP],
                          [-Phi21NP,-Phi10NP,Phi20NP,(Phi11NP+3*LambdaNP)]])
    RiccimatrixNP=matrix(4,4,lambda i,j: RiccimatrixNP[i,j].expr())
    #Note: Array slicing: RiccimatrixNP[0:2,0:2]
    #Note: Eigenvalues and eigenvectors: eigs=RiccimatrixNP.eigenvectors_right()
##########################################################################
# First Pontrjagin class:
def calculate_FirstPontrNP():
    global firstpontrNP
    firstpontrNP=(3*Psi2NP^2-4*Psi1NP*Psi3NP+Psi0NP*Psi4NP)+(-2*Phi11NP^2+6*LambdaNP^2+2*Phi01NP*Phi21NP+2*Phi10NP*Phi12NP-Phi02NP*Phi20NP-Phi00NP*Phi22NP)
##########################################################################
# Petrov Type
# Petrov invariant I and J (Kramer p.54, 4.19)
PetrovinvINPcalculated=0
def calculate_PetrovinvINP():
    global PetrovinvINP,PetrovinvINPcalculated
    PetrovinvINPcalculated=1
    PetrovinvINP=Psi0NP*Psi4NP-4*Psi1NP*Psi3NP+3*Psi2NP^2

PetrovinvJNPcalculated=0    
def calculate_PetrovinvJNP():
    global PetrovinvJNP,PetrovinvJNPcalculated
    PetrovinvJNPcalculated=1
    PetrovinvJmatrix=matrix([[Psi4NP,Psi3NP,Psi2NP],
                            [Psi3NP,Psi2NP,Psi1NP],
                            [Psi2NP,Psi1NP,Psi0NP]])
    PetrovinvJNP=PetrovinvJmatrix.determinant()

# Petrov invariant K, L, and N (Kramer p.121, 9.6)
# also check diagram Fig. 9.1 on p. 122
PetrovinvKNPcalculated=0
def calculate_PetrovinvKNP():
    global PetrovinvKNP,PetrovinvKNPcalculated
    PetrovinvKNPcalculated=1
    PetrovinvKNP=Psi1NP*Psi4NP^2-3*Psi4NP*Psi3NP*Psi2NP+2*Psi3NP^3

PetrovinvLNPcalculated=0
def calculate_PetrovinvLNP():
    global PetrovinvLNP,PetrovinvLNPcalculated
    PetrovinvLNPcalculated=1
    PetrovinvLNP=Psi2NP*Psi4NP-Psi3NP^2

PetrovinvNNPcalculated=0
def calculate_PetrovinvNNP():
    global PetrovinvNNP,PetrovinvNNPcalculated
    PetrovinvNNPcalculated=1
    if PetrovinvLNPcalculated==0:
        calculate_PetrovinvLNP()
    if PetrovinvINPcalculated==0:
        calculate_PetrovinvINP()
    PetrovinvNNP=12*PetrovinvLNP^2-(Psi4NP^2)*PetrovinvINP
    
# Petrov type from the invariants
def Petrov_frominvariants():
    if weylcalculated==0:
        calculate_Weyl()
        Petrov_frominvariants()
    else:
        if PetrovinvINPcalculated==0:
            calculate_PetrovinvINP()
        if PetrovinvJNPcalculated==0:
            calculate_PetrovinvJNP()
        if PetrovinvKNPcalculated==0:
            calculate_PetrovinvKNP()
        if PetrovinvLNPcalculated==0:
            calculate_PetrovinvLNP()
        if PetrovinvNNPcalculated==0:
            calculate_PetrovinvNNP()
        print("Calculating Petrov Type...")
        LL=simplify_fullfull(PetrovinvLNP.expr())
        JJ=simplify_fullfull(PetrovinvJNP.expr())
        II=simplify_fullfull(PetrovinvINP.expr())
        KK=simplify_fullfull(PetrovinvKNP.expr())
        NN=simplify_fullfull(PetrovinvNNP.expr())
        IIJJterm=simplify_fullfull((II^3-27*JJ^2))
        if IIJJterm==0:
            if II==0 and JJ==0:
                if KK==0 and LL==0:
                    print("Petrov Type N")
                else:
                    print("Petrov Type III")
            else:
                if KK==0 and NN==0:
                    print("Petrov Type D")
                else:
                    print("Petrov Type II")
        else:
            print("Petrov Type I")
        print("Attention: This procedure depends on the simplification of the structures. Therefore the Petrov type can be simpler.")

# Petrov type from the Weyl components
def Petrov_fromWeyl():
    if weylcalculated==0:
        calculate_Weyl()
        Petrov_fromWeyl()
    else:
        print("Calculating Petrov Type...")
        psi0=simplify_fullfull(Psi0NP.expr())
        psi1=simplify_fullfull(Psi1NP.expr())
        psi2=simplify_fullfull(Psi2NP.expr())
        psi3=simplify_fullfull(Psi3NP.expr())
        psi4=simplify_fullfull(Psi4NP.expr())
        if psi0==0 and psi1==0 and psi2==0 and psi3==0 and psi4==0:
            print("Petrov Type O")
        elif psi0==0 and psi1==0 and psi2==0 and psi3==0:
            print("Petrov Type N")
        elif psi0==0 and psi1==0 and psi2==0:
            print("Petrov Type III")
        elif psi0==0 and psi1==0 and psi3==0 and psi4==0:
            print("Petrov Type D")
        elif psi0==0 and psi1==0:
            print("Petrov Type II")
        elif psi0==0:
            print("Petrov Type I")
        print("Attention: This procedure depends on the simplification of the structures. Therefore the Petrov type can be simpler.") 
##########################################################################
# sigma2 invariant: 
def calculate_sigma2invNP():
    global sigma2invNP
    sigma2invNP=54*LambdaNP^2-2*Phi11NP^2+2*Phi10NP*Phi12NP+2*Phi01NP*Phi21NP-Phi02NP*Phi20NP-Phi00NP*Phi22NP
##########################################################################
# Calculate everything about NP
def calculate_allNP():
    calculate_spincoefficients()
    calculate_Weyl()
    calculate_Ricci()
    calculate_NPeq()
    calculate_Bianchi()
def show_allNP():
    show_spincoefficients()
    show_Weyl()
    show_Ricci()
    show_NPeq()
    show_Bianchi()
##########################################################################
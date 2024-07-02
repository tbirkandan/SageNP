##########################################################################
# Reissner-Nordstrom (with contravariant null vectors)
# Chandrasekhar p. 225
# YOU CAN CHANGE THIS PART ACCORDING TO YOUR METRIC
# Define the coordinates
###################################################
# {t=0, r=1, theta(=th)=2, phi=3} with ranges (CO = Coordinates)
CO.<t,r,th,ph> = Man.chart(r't r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
###################################################
# Define the metric functions
var('M,Q')
Delta=function('Delta',imag_part_func=0)(r)
Delta=r^2-2*M*r+Q^2
###################################################
# Enter null tetrad elements
lNPvecup=[(r^2)/Delta,1,0,0]
nNPvecup=[1/2,-Delta/(2*r^2),0,0]
mNPvecup=[0,0,1/(r*sqrt(2)),I*csc(th)/(r*sqrt(2))]
mbarNPvecup=[0,0,1/(r*sqrt(2)),-I*csc(th)/(r*sqrt(2))]
invert_tetradNP()
##########################################################################
##########################################################################
# Schwarzschild
# Chandrasekhar, p. 135
# YOU CAN CHANGE THIS PART ACCORDING TO YOUR METRIC
# Define the coordinates
###################################################
# {t=0, r=1, theta(=th)=2, phi=3} with ranges (CO = Coordinates)
CO.<t,r,th,ph> = Man.chart(r't r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
###################################################
# Define the metric functions
var('M')
Delta=function('Delta',imag_part_func=0)(r)
Delta=r^2-2*M*r
###################################################
# Enter null tetrad elements
lNPvec=[1,-(r^2)/Delta,0,0]
nNPvec=[Delta/(2*r^2),1/2,0,0]
mNPvec=[0,0,(-r/sqrt(2)),(-I*r/sqrt(2))*sin(th)]
mbarNPvec=[0,0,(-r/sqrt(2)),(I*r/sqrt(2))*sin(th)]
##########################################################################
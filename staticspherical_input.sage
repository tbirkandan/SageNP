##########################################################################
# Static spherical metric
# YOU CAN CHANGE THIS PART ACCORDING TO YOUR METRIC
# Define the coordinates
###################################################
# {t=0, r=1, theta(=th)=2, phi=3} with ranges (CO = Coordinates)
CO.<t,r,th,ph> = Man.chart(r't r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
###################################################
# Define the metric functions
A=function('A',imag_part_func=0)(r)
B=function('B',imag_part_func=0)(r)
###################################################
# Enter null tetrad elements
lNPvec=[-B,-sqrt(A*B),0,0]
nNPvec=[-1/2,(1/2)*sqrt(A/B),0,0]
mNPvec=[0,0,(r/sqrt(2)),(r/sqrt(2))*I*sin(th)]
mbarNPvec=[0,0,(r/sqrt(2)),-I*sin(th)*(r/sqrt(2))]
##########################################################################
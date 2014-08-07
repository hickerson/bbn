from matplotlib import pyplot as plt
import numpy as np
from fractions import Fraction
f13=Fraction('1/3')
f23=Fraction('2/3')
f43=Fraction('4/3')
f53=Fraction('5/3')
f12=Fraction('1/2')
f32=Fraction('3/2')
f56=Fraction('5/6')
fm1=Fraction('-1')
fm23=Fraction('-2/3')
fm32=Fraction('-3/2')
#Powers of t9 from original code

t9=np.arange(0.01,2,0.01)
t9a=t9/(1+0.1071*t9)
rate=(4.817e+6)*(t9**fm23)*np.exp(-14.964/(t9**float(f13)))*(1+0.0325*(t9**f13)-(1.04e-3)*(t9**f23)-(2.37e-4)*t9-(8.11e-5)*(t9**f43)-(4.69e-5)*(t9**f53))+(5.938e+6)*(t9a**f56)*(t9**fm32)*np.exp(-12.859/(t9a**float(f13))) #Without floats, shows Attribute error:exp
plt.plot(t9, rate, label='Old Data')
plt.plot([.01, .011, .012, .013, .014, .015, .016, .018, .02, .025, .03, .04, .05, .06, .07, .08, .09, .1, .11, .12, .13, .14, .15, .16, .18, .2, .25, .3, .35, .4, .45, .5, .6, .7, .8, .9, 1, 1.25, 1.5, 1.75, 2], [1.715e-18, 1.035e-17, 5.079e-17, 2.104e-16, 7.578e-16, 2.426e-15, 7.028e-15, 4.609e-14, 2.325e-13, 5.918e-12, 6.951e-11, 2.493e-9, 3.151e-8, 2.168e-7, 1.007e-6, 3.56e-6, 1.033e-5, 2.578e-5, 5.726e-5, .0001159, .0002173, .0003826, .0006392, .001021, .002333, .004739, .01945, .05655, .1317, .2632, .4705, .7731, 1.739, 3.296, 5.55, 8.582, 12.45, 25.92, 44.9, 69.19, 98.47], 'ro', label='New Data')
plt.xlabel('T9')
plt.ylabel('Reaction Rates (cm3/s/mol')
plt.legend()
plt.title('He3 (a,g) Be7')
plt.show()

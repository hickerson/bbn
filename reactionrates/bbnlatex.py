#example of using sympy to create latex equations
from sympy import *
t9=Symbol('t9')
from fractions import Fraction
f13=Fraction('1/3')
f23=Fraction('2/3')
f43=Fraction('4/3')
f53=Fraction('5/3')
f12=Fraction('1/2')
f32=Fraction('3/2')
f56=Fraction('5/6')
f52=Fraction('5/2')
fm1=Fraction('-1')
fm23=Fraction('-2/3')
fm32=Fraction('-3/2')
#Powers of t9 from original code-NOT A COMPLETE LIST
rate=(7.29e+2)+2.40*((10**3)*(t9**fm32)*exp(-0.223/t9))  #B11 (n,g) B12
print latex(rate, mode='equation')

t9a=t9/(1+0.1071*t9)
rate=(4.817e+6)*(t9**fm23)*exp(-14.964/(t9**f13))*(1+0.0325*(t9**f13)-(1.04e-3)*(t9**f23)-(2.37e-4)*t9-(8.11e-5)*(t9**f43)-(4.69e-5)*(t9**f53))+(5.938e+6)*(t9a**f56)*(t9**fm32)*exp(-12.859/(t9a**f13)) #He3 (a,g) Be7
print latex(rate, mode='equation')

rate=(4.61e+5)/(t9**f23)*exp(-12.062/(t9**float(f13))-(t9/4.402)**2)*(1.0+0.035*(t9**f13)+0.426*(t9**f23)+0.103*t9+0.281*(t9**f43)+0.173*(t9**f53))+(1.93e+5)/(t9**f32)*exp(-12.041/t9)+(1.14e+4)/(t9**f32)*exp(-16.164/t9) #B10 (p,g) C11
print latex(rate, mode='equation')

rate=(1.26e+11)/(t9**f23)*exp(-12.062/(t9**float(f13))-(t9/4.402)**2)*(1.0+0.035*(t9**f13)-0.498*(t9**f23)-0.121*t9+0.300*(t9**f43)+0.184*(t9**f53))+(2.59e+09)/t9*exp(-12.260/t9) #B10 (p,a) Be7
print latex(rate, mode='equation')

t9a=t9/(1+13.076*t9)
rate=(2.675e+9)*(1-0.560*(t9**f12)+0.179*t9-0.0283*(t9**f32)+(2.21e-3)*(t9**2)-(6.85e-5)*(t9**f52))+(9.391e+8)*((t9a/t9)**f32)+(4.467e+7)*(t9**fm32)*exp(-0.07486/t9)
print latex(rate, mode='equation')






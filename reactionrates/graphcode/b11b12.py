from matplotlib import pyplot as plt
import numpy as np
from fractions import Fraction
f13=Fraction('1/3')
f23=Fraction('2/3')
f43=Fraction('4/3')
f53=Fraction('5/3')
f12=Fraction('1/2')
f32=Fraction('3/2')
fm1=Fraction('-1')
fm23=Fraction('-2/3')
fm32=Fraction('-3/2')
#Powers of t9 from original code
t9=np.arange(0.01,2,0.01)
rate=(7.29e+2)+2.40*((10**3)*(t9**fm32)*np.exp(-0.223/t9)) #Reaction from original code
plt.plot(t9, rate, label='Old Data')
plt.plot([0.058, 0.073, 0.087, 0.102, 0.116, 0.131, 0.145, 0.16, 0.174, 0.232, 0.29, 0.348, 0.406, 0.464, 0.522, 0.58, 0.638, 0.696, 0.754, 0.812, 0.87, 0.928, 0.986, 1.044, 1.102], [5334, 8511, 11170, 13130, 14450, 15270, 15700, 15850, 15800, 14590, 12940, 11370, 10030, 8900, 7957, 7168, 6507, 5950, 5479, 5079, 4738, 4447, 4197, 3982, 3797], 'ro', label='New Data') #Retrieved from EXFOR
plt.xlabel('T9')
plt.ylabel('Reaction Rates (cm/s/mol)')
plt.title('B11 (n,g) B12')
plt.legend()
plt.show()
#Graph for the B11 (n,g) B12 reaction. 


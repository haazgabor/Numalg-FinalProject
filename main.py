import numpy as np
import matplotlib.pyplot as plt

'''
This function simulates the chemical decomposition of a star.
There's a connection between the fraction of He
and fraction of metals in a star and we assume the ratio of different metals found in
a star are similar to the Sun.
'''
def star_decomposition(metalicity,sun_decomp,sun_metal_content):
    #primal He content of the universe:
    Yp=0.247

    #ratio of He and metal creation:
    a=1.5

    #metal to H ratio of the star:
    r=(sun_metal_content/sun_decomp["H"])*10**(metalicity)

    #metal content of the star:
    Z=r*(1-Yp)/(1+r*(a+1))

    #He content of the star:
    Y=Yp+a*Z

    #H content of the star:
    X=1-Z-Y

    star_decomp=sun_decomp.copy()
    for element, fraction in sun_decomp.items():
        if element=="H":
            star_decomp["H"]=X
        elif element=="He":
            star_decomp["He"]=Y
        else:
            star_decomp[element]=(sun_decomp[element]/sun_metal_content)*Z
    return star_decomp

'''
This function generates the electromagnetic spectrum emitted by a star 
simulated as an ideal black body whose radiation is partially absorbed by the elements found in the stars
atmosphere.
'''
def stellar_spectrum(T,star_decomp, atomic_mass, optical_depth_coefs, absorption_lines):
    #defining constants
    h=6.62607015*10**(-34)
    c=299792458
    e=2.718281828459045
    kb=1.380649*10**(-23)
    lambda_max=900*10**-9
    lambda_min=10*10**-9
    steps=100000
    lambdas=np.linspace(lambda_min,lambda_max,steps)
    I=np.array([])
    #creating the Gaussian dips:
    for i in lambdas:
        I=np.append(I,(2*h*c**2/i**5)*(1/(e**(h*c/(kb*T*i))-1)))

        for element, fraction in star_decomp.items():
            Ai=(1-e**(-optical_depth_coefs[element]*star_decomp[element]))*1000

            for wavelength in absorption_lines[element]:
                width=(wavelength/c)*(kb*T/atomic_mass[element])**0.5
                width = 0.1e-6
                deltaI=Ai*e**(-(i-wavelength)**2/(2*width**2))
                I[-1]-=deltaI

    return lambdas,I


#photosphere composition of the sun (percentage by number of atoms): H, He, Na, Mg, Ca, Fe
sun_decomp={"H":0.92,"He":0.078,"Na":0.000002,"Mg":0.00003,"Ca":0.0000021,"Fe":0.00003}

atomic_mass={"H":1.673*10**-27, "He":6.64*10**-27, "Na":3.82*10**-26, "Mg":4.03*10**-26, "Ca":6.655*10**-26, "Fe":9.27*10**-26}

#For the simulation the absorption properties of elements are simplified to K*X
#where X is the abundance and K is the absorption coefficient which represents the
#the tendency of absorption by a given element. These are the values of the K coefficients: 
optical_depth_coefs={"H":1, "He":11.75, "Na":575000, "Mg":29499, "Ca":478469, "Fe":35461}
absorption_lines = {"H":[656.28*10**-9],"He":[587.56*10**-9],"Na":[588.995*10**-9, 589.592*10**-9],"Mg":[516.73*10**-9],"Ca":[393.37*10**-9, 396.85*10**-9],"Fe":[438.35*10**-9, 495.76*10**-9, 526.95*10**-9]}

#precalculated value from the suns decomposition
sun_metal_content=3.341*10**-4

#For now a star is characterised by its metalicity and and temperature.
#Later a galaxy can be modeled by the PDF of metalicity in a galaxy instead of a predefined value
metalicity=-0.2
T=6000

star_decomp=star_decomposition(metalicity,sun_decomp,sun_metal_content)
print(star_decomp)

lambdas,I=stellar_spectrum(T,star_decomp, atomic_mass, optical_depth_coefs, absorption_lines)
plt.plot(lambdas,I)
plt.show()
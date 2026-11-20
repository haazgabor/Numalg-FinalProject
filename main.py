import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

'''
This function adds wavelength dependent gaussian noise to the original signal.
The relative error depends on the original intensity at a given wavelength.
Since at lower intensities less photons are measured the relative error is higher
at lower inetsities while lower at higher intensities. This phenomena is represented by this 
function.
'''
def add_noise(I, noise_factor=0.03):
    #standard deviation at each wavelength
    sigma=noise_factor*np.sqrt(I/np.max(I))

    noise=np.random.normal(loc=0,scale=sigma,size=len(I))

    I_noisy=I/np.max(I)+noise

    return I_noisy*I

'''This function uses Wien's law to estimate the temperature 
of a star based on the wavelength of maximum emission.'''

def estimate_temperature(wavelengths, flux):
    b = 2.898e-3
    max_index = np.argmax(flux)
    lambda_max = wavelengths[max_index]

    lambda_max = lambda_max * 1e-9 #we have to use meters, so we have to convert it
    temperature= b / lambda_max

    return temperature

"""
This function identifes the elements present in a star based 
on the detected absorption lines and the 
reference absorption lines of the given elements.
"""

def identify_elements(detected_wavelengths,absorption_lines,tolerance=0.05):
    identified = {}

    for element, reference_lines in absorption_lines.items():
        for reference in reference_lines:
            for detected in detected_wavelengths:
                difference = abs(detected - reference)
                if difference <= tolerance:
                    if element not in identified:
                        identified[element] = []
                    identified[element].append({"observed": detected,"reference": reference,"difference": difference})

    return identified

"""This function detects the absorption lines in a given spectrum
 by finding the local minima in the flux data."""

def detect_spectral_lines(wavelengths, flux):
    #find absorption lines
    peaks, properties = find_peaks(-flux, prominence=0.01,distance=5)

    #get wavelength and intensity of detected lines
    detected_wavelengths = wavelengths[peaks]
    detected_intensities = flux[peaks]

    return detected_wavelengths, detected_intensities

"""This function simulates the spectrum of a star
 by randomly generating a metalicity and temperature for the star
 , calculating its chemical decomposition, and generating its spectrum 
 based on the black body radiation and absorption lines of the elements 
 present in the star's atmosphere."""

def sim_random_spectrum(sun_decomp, sun_metal_content, atomic_mass, optical_depth_coefs, absorption_lines):
    #randomly generates a metalicity and temperature for the star

    metalicity=np.random.uniform(-2,0.5)
    T=np.random.uniform(2000,10000)
    
    #calculates the chemical decomp of the random star
    star_decomp=star_decomposition(metalicity,sun_decomp,sun_metal_content)

    #generates the spectrum of the star
    lambdas, I = stellar_spectrum(T,star_decomp, atomic_mass, optical_depth_coefs, absorption_lines)
    return lambdas, I, T, metalicity, star_decomp



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
    sigma=1400000
    lambdas=np.linspace(lambda_min,lambda_max,steps)
    I=np.array([])
    #creating the Gaussian dips:
    for i in lambdas:
        I=np.append(I,(2*h*c**2/i**5)*(1/(e**(h*c/(kb*T*i))-1)))

        for element, fraction in star_decomp.items():
            Ai=(1-e**(-optical_depth_coefs[element]*star_decomp[element]))

            for wavelength in absorption_lines[element]:
                width=(wavelength/c)*(kb*T/atomic_mass[element]+sigma**2)**0.5
                deltaI=Ai*e**(-(i-wavelength)**2/(2*width**2))
                I[-1]=I[-1]*(1-deltaI)

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
lambdas, I, true_T, true_metalicity, true_decomp = sim_random_spectrum(sun_decomp,sun_metal_content,atomic_mass,optical_depth_coefs,absorption_lines)


print("True temperature:", true_T)
print("True metallicity:", true_metalicity)
print("True composition:", true_decomp)
star_decomp=star_decomposition(true_metalicity,sun_decomp,sun_metal_content)
print(star_decomp)

lambdas,I=stellar_spectrum(true_T,star_decomp, atomic_mass, optical_depth_coefs, absorption_lines)
I_noisy=add_noise(I, noise_factor=0.03)

plt.plot(lambdas,I_noisy)
plt.show()
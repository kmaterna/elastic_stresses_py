"""
This is a set of functions that show important scaling relationships between 
earthquake magnitude and other parameters. 
FAULT_TYPE = SS, R, N, All
a, b = parameters in equation
sa, sb = parameter uncertainties
s = estimate uncertainty
r = R correlation coefficient

I haven't yet coded the Average displacements, Max displacements, or the inverse relationships. 

"""

import numpy as np
import sys

def check_fault_types(fault_type):
	if fault_type == "SS" or fault_type == "N" or fault_type == "R" or fault_type == "ALL":
		return;
	else:
		sys.exit("Error! Fault Type should be SS, N, R, or ALL!  Exiting...");
		return;
def check_magnitude(M, mmin, mmax):
	if M >=mmin and M<=mmax:
		return;
	else:
		sys.exit("Error!  You're asking for a magnitude outside the valid range of %f to %f. Exiting..." % (mmin, mmax))
	return;


def rectangular_slip(length,width,magnitude):
# Rectangular slip : length in m, width in m, magnitude, and slip in m. 
	mu = 30*1e9;    # 30 GPa rigidity. 
	area = length * width; # meters^2
	moment = np.power(10,((1.5 * magnitude) + 16.05));  # dyne-cm
	moment_nm = moment * 1e-7;   # Newton-meters
	slip = moment_nm / (area*mu);      # Moment = Area * slip * mu. 
	return slip;

def get_magnitude(length, width, slip):
	mu=30*1e9;
	area=length*width; # meters^2
	moment_nm = slip * area * mu; # newton-meters
	moment_dcm=moment_nm*1e7; # dyne-cm
	magnitude = (2/3.0)*np.log10(moment_dcm) - 10.7;
	return magnitude;


def SLR_from_M(M, fault_type):
	"""
	SLR = surface rupture length (km);
	M = Magnitude; 
	Form: log(SLR)=a + b*M
	"""
	check_fault_types(fault_type);
	if fault_type=="SS":
		[a,sa,b,sb,s,r,mmin,mmax] = [-3.55, 0.37, 0.74, 0.05, 0.23, 0.91, 5.6, 8.1];
	if fault_type=="R":
		[a,sa,b,sb,s,r,mmin,mmax] = [-2.86, 0.55, 0.63, 0.08, 0.20, 0.88, 5.4, 7.4];
	if fault_type=="N":
		[a,sa,b,sb,s,r,mmin,mmax] = [-2.01, 0.65, 0.50, 0.10, 0.21, 0.81, 5.2, 7.3];
	if fault_type=="ALL":
		[a,sa,b,sb,s,r,mmin,mmax] = [-3.22, 0.27, 0.69, 0.04, 0.22, 0.89, 5.2, 8.1];
	check_magnitude(M, mmin, mmax);
	SLR = np.power(10, (a + b*M));
	return SLR;


def RW_from_M(M, fault_type):
	"""
	RW = downdip rupture width (km);
	M = Magnitude; 
	Form: log(RW)=a + b*M
	"""
	check_fault_types(fault_type);
	if fault_type=="SS":
		[a,sa,b,sb,s,r,mmin,mmax] = [-0.76, 0.12, 0.27, 0.02, 0.14, 0.84, 4.8, 8.1];
	if fault_type=="R":
		[a,sa,b,sb,s,r,mmin,mmax] = [-1.61, 0.20, 0.41, 0.03, 0.15, 0.90, 4.8, 7.6];
	if fault_type=="N":
		[a,sa,b,sb,s,r,mmin,mmax] = [-1.14, 0.28, 0.35, 0.05, 0.12, 0.86, 5.2, 7.3];
	if fault_type=="ALL":
		[a,sa,b,sb,s,r,mmin,mmax] = [-1.01, 0.10, 0.32, 0.02, 0.15, 0.84, 4.8, 8.1];
	check_magnitude(M, mmin, mmax);
	RW = np.power(10, (a + b*M));
	return RW;


def RLD_from_M(M, fault_type):
	"""
	RLD = subsurface rupture length (km);
	M = Magnitude; 
	Form: log(RLD)=a + b*M
	"""
	check_fault_types(fault_type);
	if fault_type=="SS":
		[a,sa,b,sb,s,r,mmin,mmax] = [-2.57, 0.12, 0.62, 0.02, 0.15, 0.96, 4.8, 8.1];
	if fault_type=="R":
		[a,sa,b,sb,s,r,mmin,mmax] = [-2.42, 0.21, 0.58, 0.03, 0.16, 0.93, 4.8, 7.6];
	if fault_type=="N":
		[a,sa,b,sb,s,r,mmin,mmax] = [-1.88, 0.37, 0.50, 0.06, 0.17, 0.88, 5.2, 7.3];
	if fault_type=="ALL":
		[a,sa,b,sb,s,r,mmin,mmax] = [-2.44, 0.11, 0.59, 0.02, 0.16, 0.94, 4.8, 8.1];
	check_magnitude(M, mmin, mmax);
	RLD = np.power(10, (a + b*M));
	return RLD;

def RA_from_M(M, fault_type):
	"""
	RA = rupture area (km^2);
	M = Magnitude; 
	Form: log(RA)=a + b*M
	"""
	check_fault_types(fault_type);
	if fault_type=="SS":
		[a,sa,b,sb,s,r,mmin,mmax] = [-3.42, 0.18, 0.90, 0.03, 0.22, 0.96, 4.8, 7.9];
	if fault_type=="R":
		[a,sa,b,sb,s,r,mmin,mmax] = [-3.99, 0.36, 0.98, 0.06, 0.26, 0.94, 4.8, 7.6];
	if fault_type=="N":
		[a,sa,b,sb,s,r,mmin,mmax] = [-2.87, 0.50, 0.82, 0.08, 0.22, 0.92, 5.2, 7.3];
	if fault_type=="ALL":
		[a,sa,b,sb,s,r,mmin,mmax] = [-3.49, 0.16, 0.91, 0.03, 0.24, 0.95, 4.8, 7.9];
	check_magnitude(M, mmin, mmax);
	RA = np.power(10, (a + b*M));
	return RA;



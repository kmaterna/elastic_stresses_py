# Conversion functions


import numpy as np 
import math


def get_strain_tensor(dUidUj):
	# Starts with displacement gradient tensor (3x3 2D array);
	# Returns a strain tensor (3x3 2d array). 
	rows,cols=np.shape(dUidUj);
	strain_tensor=np.zeros(np.shape(dUidUj));
	for i in range(rows):
		for j in range(cols):
			strain_tensor[i][j]=0.5* (dUidUj[i][j] + dUidUj[j][i]);
	return strain_tensor;

def get_stress_tensor(eij, lamda, mu):
	# Starts with displacement gradient tensor (3x3 2D array);
	# Returns a strain tensor (3x3 2d array). 
	# lamda and mu are Lame parameters
	rows,cols=np.shape(eij);
	stress_tensor=np.zeros(np.shape(eij));
	for i in range(rows):
		for j in range(cols):
			if i==j:
				stress_tensor[i][j]=lamda*(eij[0][0]+eij[1][1]+eij[2][2]) + 2*mu*eij[i][j];
			else:
				stress_tensor[i][j]=2*mu*eij[i][j];
	return stress_tensor;

def get_coulomb_stresses(tau, strike, dip, rake, friction):
	# THIS IS IMPORTANT. WILL DO TOMORROW. 
	normal=0;
	shear=0;
	coulomb=0;
	return normal, shear, coulomb;

def get_strike(deltax, deltay):
	# Returns the strike of a line (in cw degrees from north) given the deltax and deltay in km. 
	slope = math.atan2(deltay,deltax);
	strike= 90-np.rad2deg(slope);
	if strike<0:
		strike=strike+360;
	return strike;

def get_rtlat_dip_slip(slip, rake):
	strike_slip = slip * np.cos(np.deg2rad(rake));
	dip_slip = slip * np.sin(np.deg2rad(rake));
	return strike_slip, dip_slip;

def get_strike_length(x0,x1,y0,y1):
	# Just the pythagorean theorem
	length=np.sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
	return length;

def get_downdip_width(top,bottom,dip):
	W = abs(top-bottom)/np.sin(np.deg2rad(dip));  # guaranteed to be between 0 and 90
	return W;

def add_vector_to_point(x0,y0,vector_mag,vector_heading):
	# Vector heading defined as strike- CW from north.
	theta=np.deg2rad(90-vector_heading);
	x1 = x0 + vector_mag*np.cos(theta);
	y1 = y0 + vector_mag*np.sin(theta);
	return x1, y1;

def get_rake(strike_slip, dip_slip):
	# This is something that I'm too tired to think about right now. 
	# WILL DO THE TRIGONOMETRY TOMORROW. 

	return 0;


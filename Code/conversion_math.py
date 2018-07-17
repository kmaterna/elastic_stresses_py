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
	# Positive slip is right lateral, and reverse. 
	# Range is -180 to 180.
	rake = np.rad2deg(math.atan2(dip_slip,strike_slip));
	return rake;

def get_fault_center(fault_object, index):
	# Compute the x-y-z coordinates of the center of a fault patch. 
	# Index is the i'th fault patch in this fault_object
	W = get_downdip_width(fault_object.top[index],fault_object.bottom[index],fault_object.dipangle[index]);
	center_z = (fault_object.top[index]+fault_object.bottom[index])/2.0;
	updip_center_x=(fault_object.xstart[index]+fault_object.xfinish[index])/2.0;
	updip_center_y=(fault_object.ystart[index]+fault_object.yfinish[index])/2.0;
	vector_mag = W*np.cos(np.deg2rad(fault_object.dipangle[index]))/2.0;  # how far the middle is displaced downdip from map-view
	center_point = add_vector_to_point(updip_center_x,updip_center_y,vector_mag, fault_object.strike[index]+90);  # strike+90 = downdip direction. 
	center = [center_point[0],center_point[1],center_z]; 
	return center; 

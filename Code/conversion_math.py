# Stress and strain functions
# Conversion functions
# Fault plane geometric functions


import numpy as np 
import math
import haversine


def get_strain_tensor(dUidUj):
	# Starts with displacement gradient tensor (3x3 2D array);
	# Returns a strain tensor (3x3 2D array). 
	rows,cols=np.shape(dUidUj);
	strain_tensor=np.zeros(np.shape(dUidUj));
	for i in range(rows):
		for j in range(cols):
			strain_tensor[i][j]=0.5* (dUidUj[i][j] + dUidUj[j][i]);
	return strain_tensor;

def get_stress_tensor(eij, lamda, mu):
	# Starts with strain tensor (3x3 2D array);
	# Returns a stress tensor (3x3 2D array). 
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



def get_coulomb_stresses(tau, strike, rake, dip, friction):
	# Given a stress tensor, strike, rake, and dip
	# Resolve the stress changes on the fault plane. 
	# Return in KPa

	strike_unit_vector= get_strike_vector(strike);  # a 3d vector in the horizontal plane. 
	dip_unit_vector = get_dip_vector(strike, dip);  # a 3d vector.
	plane_normal = get_plane_normal(strike, dip);  # a 3d vector.
	traction_vector = np.dot(tau,plane_normal);

	# The stress that's normal to the receiver fault plane:
	normal_stress = np.dot(plane_normal, traction_vector);  # positive for unclamping stress (same as Coulomb software)

	# The shear stress causing strike slip (in the receiver fault plane).
	shear_rtlat = np.dot(strike_unit_vector, traction_vector);

	# The shear stress causing reverse slip (in the receiver fault plane).
	shear_reverse = np.dot(dip_unit_vector, traction_vector);

	# The shear that we want (in the rake direction). 
	rake_rad=np.deg2rad(rake);
	R = np.array([[np.cos(rake_rad),-np.sin(rake_rad)],[np.sin(rake_rad),np.cos(rake_rad)]]);
	shear_vector = [-shear_rtlat, shear_reverse];  # There is a critical minus sign here to define right lateral as positive. 
	rotated_shear = np.dot(R,shear_vector);
	shear_in_rake_dir = rotated_shear[0];

	# Finally, do unit conversion	
	normal_stress=normal_stress/1000.0;  # convert to KPa
	shear_stress=shear_in_rake_dir/1000.0;

	# The Coulomb Failure Hypothesis
	coulomb_stress = shear_stress + (friction*normal_stress);   # the sign here is important. 

	return normal_stress, shear_stress, coulomb_stress;



def get_plane_normal(strike, dip):
	# Given a strike and dip, find the orthogonal unit vectors aligned with strike and dip directions that sit within the plane. 
	# The plane normal is their cross product. 
	# Returns in x, y, z coordinates. 
	strike_vector=get_strike_vector(strike); # unit vector
	dip_vector = get_dip_vector(strike, dip); # unit vector
	plane_normal = np.cross(dip_vector, strike_vector);  # dip x strike for outward facing normal, by right hand rule. 
	return plane_normal;


def get_strike_vector(strike):
	# Returns a unit vector in x-y-z coordinates 
	theta=np.deg2rad(90-strike);
	strike_vector=[np.cos(theta),np.sin(theta),0];
	return strike_vector;

def get_dip_vector(strike,dip):
	# Returns a unit vector in x-y-z coordinates
	downdip_direction_theta = np.deg2rad(-strike);  # theta(strike+90)
	dip_unit_vector_z=np.sin(np.deg2rad(dip))  # the vertical component of the downdip unit vector
	dip_unit_vector_xy = np.sqrt(1-dip_unit_vector_z*dip_unit_vector_z);  # the horizontal component of the downdip unit vector
	dip_vector = [ dip_unit_vector_xy*np.cos(downdip_direction_theta), dip_unit_vector_xy*np.sin(downdip_direction_theta), -dip_unit_vector_z];
	return dip_vector;

def get_vector_magnitude(vector):
	total=0;
	for i in range(len(vector)):
		total=total+vector[i]*vector[i];
		magnitude=np.sqrt(total);
	return magnitude;


def get_strike(deltax, deltay):
	# Returns the strike of a line (in cw degrees from north) given the deltax and deltay in km. 
	slope = math.atan2(deltay,deltax);
	strike= 90-np.rad2deg(slope);
	if strike<0:
		strike=strike+360;
	return strike;

def get_rtlat_dip_slip(slip, rake):
	strike_slip = -slip * np.cos(np.deg2rad(rake));  # negative sign for convention of right lateral slip
	dip_slip = slip * np.sin(np.deg2rad(rake));
	return strike_slip, dip_slip;

def get_strike_length(x0,x1,y0,y1):
	# Just the pythagorean theorem
	length=np.sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
	return length;

def get_downdip_width(top,bottom,dip):
	W = abs(top-bottom)/np.sin(np.deg2rad(dip));  # guaranteed to be between 0 and 90
	return W;

def get_top_bottom(center_depth, width,dip):
	# Given a fault, where is the top and bottom? 
	# Width is total downdip width of the fault. 
	top=center_depth-(width/2.0*np.sin(np.deg2rad(dip)));
	bottom=center_depth+(width/2.0*np.sin(np.deg2rad(dip)));
	return top,bottom;

def get_top_bottom_from_top(top_depth, width,dip):
	bottom=top_depth+(width*np.sin(np.deg2rad(dip)));
	return top_depth,bottom;

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

def get_fault_four_corners(fault_object, i):
	# Get the four corners of the object, including updip and downdip. 
	W = get_downdip_width(fault_object.top[i],fault_object.bottom[i],fault_object.dipangle[i]);
	depth       = fault_object.top[i];
	strike      = fault_object.strike[i];
	dip         = fault_object.dipangle[i];

	updip_point0 = [fault_object.xstart[i],fault_object.ystart[i]];
	updip_point1 = [fault_object.xfinish[i],fault_object.yfinish[i]];
	vector_mag = W*np.cos(np.deg2rad(fault_object.dipangle[i]));  # how far the bottom edge is displaced downdip from map-view
	downdip_point0 = add_vector_to_point(fault_object.xstart[i],fault_object.ystart[i],vector_mag, strike+90);  # strike+90 = downdip direction. 
	downdip_point1 = add_vector_to_point(fault_object.xfinish[i],fault_object.yfinish[i], vector_mag, strike+90);

	x_total = [updip_point0[0], updip_point1[0], downdip_point1[0], downdip_point0[0],updip_point0[0]];
	y_total = [updip_point0[1], updip_point1[1], downdip_point1[1], downdip_point0[1],updip_point0[1]];
	x_updip = [updip_point0[0], updip_point1[0]];
	y_updip = [updip_point0[1], updip_point1[1]];
	return [x_total, y_total, x_updip, y_updip];

def xy2lonlat(xi,yi,reflon,reflat):
	lat=reflat+( yi*1/111.000 );
	lon=reflon+( xi*1/(111.000*abs(np.cos(np.deg2rad(reflat)))) );
	return lon, lat;


def latlon2xy(loni,lati,lon0,lat0):
	# returns the distance between a point and a reference in km. 
	radius = haversine.distance([lat0,lon0], [lati,loni]);
	bearing = haversine.calculate_initial_compass_bearing((lat0, lon0),(lati, loni))
	azimuth = 90 - bearing;
	x = radius * np.cos(np.deg2rad(azimuth));
	y = radius * np.sin(np.deg2rad(azimuth));
	return [x, y];




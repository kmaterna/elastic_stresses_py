# Running dc3d, given an input namedtuple.


import numpy as np 
import matplotlib.pyplot as plt 
import sys
from okada_wrapper import dc3dwrapper
import coulomb_collections
import conversion_math


def do_stress_computation(params, inputs, disp_points):
	# Step 0. Split receiver fault into many sub-faults if necessary
	# Step 1. Compute strains and displacements
	# Step 2. Resolve stresses on receiver faults

	print("Beginning stress calcultaion.");
	print("Number of sources: %d " % len(inputs.source_object.xstart));
	print("Number of receivers: %d " % len(inputs.receiver_object.xstart));
	subfaulted_inputs = split_subfaults(params, inputs);
	
	[x, y, x2d, y2d, u_displacements, v_displacements, w_displacements] = compute_surface_disp_grid(params, subfaulted_inputs);
	[u_ll, v_ll, w_ll] = compute_surface_disp_ll(params, subfaulted_inputs, disp_points);
	[source_object, receiver_object, receiver_normal, receiver_shear, receiver_coulomb] = compute_strains_stresses(params, subfaulted_inputs);

	MyOutObject = coulomb_collections.Out_object(x=x,y=y,x2d=x2d, y2d=y2d, u_disp=u_displacements, v_disp=v_displacements, w_disp=w_displacements, 
		u_ll=u_ll, v_ll=v_ll, w_ll=w_ll, 
		source_object=source_object, receiver_object=receiver_object, receiver_normal=receiver_normal, receiver_shear=receiver_shear, 
		receiver_coulomb=receiver_coulomb); 
	return MyOutObject;




def split_subfaults(params,inputs):
	receiver_object = inputs.receiver_object;
	strike_split = params.strike_num_receivers;
	dip_split = params.dip_num_receivers;

	if strike_split==1 and dip_split==1:
		# If we're not splitting the subfaults...
		subfaulted_receivers=inputs.receiver_object;
		print("Not subdividing receiver faults further.");
	else:
		print("Splitting %d receiver faults into %d subfaults each." % (len(receiver_object.xstart), strike_split*dip_split) );
		new_xstart=[];
		new_xfinish=[];
		new_ystart=[];
		new_yfinish=[];
		new_Kode=[];
		new_rtlat=[];
		new_reverse=[];
		new_strike=[];
		new_dipangle=[];
		new_rake=[];
		new_top=[];
		new_bottom=[];
		new_comment=[];

		for i in range(len(receiver_object.xstart)):
			# For each receiver... 

			# We find the depths corresponding to the tops and bottoms of our new sub-faults
			zsplit_array = get_split_z_array(inputs.receiver_object.top[i],inputs.receiver_object.bottom[i],dip_split);
			x0=inputs.receiver_object.xstart[i];
			y0=inputs.receiver_object.ystart[i];
			x1=inputs.receiver_object.xfinish[i];
			y1=inputs.receiver_object.yfinish[i];			

			for j in range(dip_split):
				# First we split it up by dip. 

				# Get the new coordinates of the top of the fault plane. 
				W = conversion_math.get_downdip_width(inputs.receiver_object.top[i],zsplit_array[j],inputs.receiver_object.dipangle[i]);
				vector_mag = W*np.cos(np.deg2rad(inputs.receiver_object.dipangle[i]));  # how far the bottom edge is displaced downdip from map-view

				# Get the starting points for the next row of fault subpatches. 
				[start_x_top, start_y_top] = conversion_math.add_vector_to_point(x0,y0,vector_mag,inputs.receiver_object.strike[i]+90);
				[finish_x_top, finish_y_top] = conversion_math.add_vector_to_point(x1,y1,vector_mag,inputs.receiver_object.strike[i]+90);

				[xsplit_array, ysplit_array] = get_split_x_y_arrays(start_x_top, finish_x_top, start_y_top, finish_y_top, strike_split);

				for k in range(strike_split):
					new_xstart.append(xsplit_array[k]);
					new_xfinish.append(xsplit_array[k+1]);
					new_ystart.append(ysplit_array[k]);
					new_yfinish.append(ysplit_array[k+1]);
					new_Kode.append(receiver_object.Kode[i]);
					new_rtlat.append(receiver_object.rtlat[i]);
					new_reverse.append(receiver_object.reverse[i]);
					new_strike.append(receiver_object.strike[i]);
					new_dipangle.append(receiver_object.dipangle[i]);
					new_rake.append(receiver_object.rake[i]);
					new_top.append(zsplit_array[j]);
					new_bottom.append(zsplit_array[j+1]);
					new_comment.append(receiver_object.comment[i]);

		subfaulted_receivers = coulomb_collections.Faults_object(xstart=new_xstart, xfinish=new_xfinish, ystart=new_ystart, yfinish=new_yfinish, Kode=new_Kode, rtlat=new_rtlat, 
			reverse=new_reverse, strike=new_strike, dipangle=new_dipangle, rake=new_rake, top=new_top, bottom=new_bottom, comment=new_comment);
	
	subfaulted_inputs = coulomb_collections.Input_object(PR1=inputs.PR1,FRIC=inputs.FRIC,depth=inputs.depth,start_gridx=inputs.start_gridx, finish_gridx=inputs.finish_gridx, 
		start_gridy=inputs.start_gridy, finish_gridy=inputs.finish_gridy, xinc=inputs.xinc, yinc=inputs.yinc, minlon=inputs.minlon,maxlon=inputs.maxlon,
		zerolon=inputs.zerolon,minlat=inputs.minlat,maxlat=inputs.maxlat,zerolat=inputs.zerolat,eqlon=inputs.eqlon, eqlat=inputs.eqlat,
		source_object=inputs.source_object,
		receiver_object=subfaulted_receivers);

	return subfaulted_inputs;



def get_split_x_y_arrays(start_x_top, finish_x_top, start_y_top, finish_y_top, strike_split):
	# Take the coordinates of the top of a receiver fault plane. 
	# Generate the list of coordinates that will help split it up along-strike
	if start_x_top==finish_x_top:
		xsplit_array = [start_x_top for j in range(strike_split+1)];
	else:
		xincrement = (finish_x_top-start_x_top)/strike_split;
		xsplit_array = np.arange(start_x_top,finish_x_top+xincrement,xincrement);  
		# length : xsplit+1. contains all the xlocations that could be used as start-stop points in the top row. 
	if start_y_top==finish_y_top:
		ysplit_array = [start_y_top for j in range(strike_split+1)];
	else:
		yincrement = (finish_y_top-start_y_top)/strike_split;
		ysplit_array = np.arange(start_y_top,finish_y_top+yincrement,yincrement); 
		# length : xsplit+1. contains all the ylocations that could be used as start-stop points in the top row. 
	return [xsplit_array, ysplit_array] ;


def get_split_z_array(top,bottom,dip_split):
	if top==bottom:
		zsplit_array = [top for j in range(dip_split+1)];
	else:
		zincrement = abs(top-bottom)/dip_split;
		zsplit_array = np.arange(top,bottom+zincrement,zincrement); 
	return zsplit_array;


def compute_surface_disp_grid(params, inputs):

	x=np.linspace(inputs.start_gridx,inputs.finish_gridx,(inputs.finish_gridx-inputs.start_gridx)/inputs.xinc);
	y=np.linspace(inputs.start_gridy,inputs.finish_gridy,(inputs.finish_gridy-inputs.start_gridy)/inputs.yinc);
	[x2d,y2d] = np.meshgrid(x,y);
	u_displacements = np.zeros((len(y), len(x)));
	v_displacements = np.zeros((len(y), len(x)));
	w_displacements = np.zeros((len(y), len(x)));
	numrows=np.shape(u_displacements)[0]
	numcols=np.shape(u_displacements)[1]

	# A major compute loop for each source object. 
	for i in range(len(inputs.source_object.xstart)):

		# Fault parameters
		L = conversion_math.get_strike_length(inputs.source_object.xstart[i],inputs.source_object.xfinish[i],inputs.source_object.ystart[i],inputs.source_object.yfinish[i]);
		W = conversion_math.get_downdip_width(inputs.source_object.top[i],inputs.source_object.bottom[i],inputs.source_object.dipangle[i]);
		depth       = inputs.source_object.top[i];
		strike      = inputs.source_object.strike[i];
		dip         = inputs.source_object.dipangle[i];
		strike_slip = inputs.source_object.rtlat[i]*-1;  # The dc3d coordinate system has left-lateral positive. 
		dip_slip    = inputs.source_object.reverse[i];		

		# Preparing to rotate to a fault-oriented coordinate system.
		theta=inputs.source_object.strike[i]-90;
		theta=np.deg2rad(theta);
		R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
		R2=np.array([[np.cos(-theta),-np.sin(-theta)],[np.sin(-theta),np.cos(-theta)]])
 
		for ky in range(numrows):
			for kx in range(numcols):

				# Compute the position relative to the translated, rotated fault. 
				# print("%d %d " %(kx, ky));
				translated_pos = np.array([[x2d[ky][kx]-inputs.source_object.xstart[i]],[y2d[ky][kx]-inputs.source_object.ystart[i]]]);
				xy=R.dot(translated_pos);
				success, u, grad_u = dc3dwrapper(params.alpha, [xy[0], xy[1], 0.0], depth, dip, [0, L], [-W, 0], [strike_slip, dip_slip, 0.0]);  # solve for displacements at the surface
				urot=R2.dot(np.array([[u[0]], [u[1]]]));

				# Update the displacements from all sources 
				u_displacements[ky][kx]=u_displacements[ky][kx] + urot[0];
				v_displacements[ky][kx]=v_displacements[ky][kx] + urot[1];
				w_displacements[ky][kx]=w_displacements[ky][kx] + u[2];  # vertical


	# OUTPUT GRIDS AND DISPLACEMENTS
	return [x, y, x2d, y2d, u_displacements, v_displacements, w_displacements];
	


def compute_surface_disp_ll(params, inputs, disp_points):
	x=[]; y=[];
	if disp_points==[]:
		return [ [], [], [] ];

	# convert here
	for i in range(len(disp_points[0])):
		[xi,yi]=conversion_math.latlon2xy(disp_points[0][i],disp_points[1][i],inputs.zerolon,inputs.zerolat);
		x.append(xi);
		y.append(yi);

	u_ll = np.zeros(len(x));
	v_ll = np.zeros(len(x));
	w_ll = np.zeros(len(x));

	# A major compute loop for each source object. 
	for i in range(len(inputs.source_object.xstart)):

		# Fault parameters
		L = conversion_math.get_strike_length(inputs.source_object.xstart[i],inputs.source_object.xfinish[i],inputs.source_object.ystart[i],inputs.source_object.yfinish[i]);
		W = conversion_math.get_downdip_width(inputs.source_object.top[i],inputs.source_object.bottom[i],inputs.source_object.dipangle[i]);
		depth       = inputs.source_object.top[i];
		strike      = inputs.source_object.strike[i];
		dip         = inputs.source_object.dipangle[i];
		strike_slip = inputs.source_object.rtlat[i]*-1;  # The dc3d coordinate system has left-lateral positive. 
		dip_slip    = inputs.source_object.reverse[i];		

		# Preparing to rotate to a fault-oriented coordinate system.
		theta=inputs.source_object.strike[i]-90;
		theta=np.deg2rad(theta);
		R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
		R2=np.array([[np.cos(-theta),-np.sin(-theta)],[np.sin(-theta),np.cos(-theta)]])
 
		for k in range(len(x)):

			# Compute the position relative to the translated, rotated fault. 
			# print("%d %d " %(kx, ky));
			translated_pos = np.array([[x[k]-inputs.source_object.xstart[i]],[y[k]-inputs.source_object.ystart[i]]]);
			xy=R.dot(translated_pos);
			success, u, grad_u = dc3dwrapper(params.alpha, [xy[0], xy[1], 0.0], depth, dip, [0, L], [-W, 0], [strike_slip, dip_slip, 0.0]);  # solve for displacements at the surface
			urot=R2.dot(np.array([[u[0]], [u[1]]]));

			# Update the displacements from all sources 
			u_ll[k]=u_ll[k] + urot[0];
			v_ll[k]=v_ll[k] + urot[1];
			w_ll[k]=w_ll[k] + u[2];  # vertical


	# OUTPUT GRIDS AND DISPLACEMENTS
	return [u_ll, v_ll, w_ll];





def compute_strains_stresses(params, inputs):

	# Pseudocode: 
	# For each receiver, at the center point, sum up the strain and stress for each source.
	# Return : source object, receiver object, shear stress, normal stress, and coulomb stress on each receiver. 

	number_of_receivers=len(inputs.receiver_object.xstart);
	
	# Where do we calculate the stress tensor? 
	receiver_center_x=[]; receiver_center_y=[]; receiver_center_z=[];

	# The values we're actually going to output. 
	receiver_shear=[];
	receiver_normal=[];
	receiver_coulomb=[];

	for m in range(number_of_receivers):
		centercoords = conversion_math.get_fault_center(inputs.receiver_object,m);
		receiver_center_x.append(centercoords[0]);
		receiver_center_y.append(centercoords[1]);
		receiver_center_z.append(centercoords[2]);
		receiver_strike = inputs.receiver_object.strike[m];
		receiver_dip    = inputs.receiver_object.dipangle[m];
		receiver_rake   = inputs.receiver_object.rake[m];
		normal_sum=0;
		shear_sum=0;
		coulomb_sum=0;

		for i in range(len(inputs.source_object.xstart)):
		# A major compute loop for each source object. 

			L = conversion_math.get_strike_length(inputs.source_object.xstart[i],inputs.source_object.xfinish[i],inputs.source_object.ystart[i],inputs.source_object.yfinish[i]);
			W = conversion_math.get_downdip_width(inputs.source_object.top[i],inputs.source_object.bottom[i],inputs.source_object.dipangle[i]);
			depth       = inputs.source_object.top[i];
			strike      = inputs.source_object.strike[i];
			dip         = inputs.source_object.dipangle[i];
			strike_slip = inputs.source_object.rtlat[i]*-1; # The dc3d coordinate system has left-lateral positive. 
			dip_slip    = inputs.source_object.reverse[i];

			# Preparing to rotate to a fault-oriented coordinate system.
			theta=inputs.source_object.strike[i]-90;
			theta=np.deg2rad(theta);
			R=np.array([[np.cos(theta),-np.sin(theta), 0],[np.sin(theta),np.cos(theta), 0],[0,0,1]]);  # horizontal rotation into strike-aligned coordinates.
			R2=np.array([[np.cos(-theta),-np.sin(-theta),0],[np.sin(-theta),np.cos(-theta),0],[0,0,1]]); 
			
			# Compute the position relative to the translated, rotated fault. 
			translated_pos = np.array([[centercoords[0]-inputs.source_object.xstart[i]],[centercoords[1]-inputs.source_object.ystart[i]],[-centercoords[2]] ]);
			xyz=R.dot(translated_pos);
			success, u, grad_u = dc3dwrapper(params.alpha, [xyz[0], xyz[1], xyz[2]], depth, dip, [0, L], [-W, 0], [strike_slip, dip_slip, 0.0]);  
			# Solve for displacement gradients at center of receiver fault

			# Rotate grad_u back into the unprimed coordinates.  Divide by 1000 because coordinate units (km) and slip units (m) are different by 1000. 
			desired_coords_grad_u = np.dot(R2, np.dot(grad_u, R2.T));
			desired_coords_grad_u = [k/1000.0 for k in desired_coords_grad_u];
			
			# Then rotate again into receiver coordinates. 
			strain_tensor=conversion_math.get_strain_tensor(desired_coords_grad_u);
			stress_tensor=conversion_math.get_stress_tensor(strain_tensor, params.lame1, params.mu);

			# Then compute shear, normal, and coulomb stresses. 
			[normal, shear, coulomb]=conversion_math.get_coulomb_stresses(stress_tensor,receiver_strike,receiver_rake,receiver_dip,inputs.FRIC);
			normal_sum=normal_sum+normal;
			shear_sum=shear_sum+shear;
			coulomb_sum=coulomb_sum+coulomb;

		receiver_normal.append(normal_sum);
		receiver_shear.append(shear_sum);
		receiver_coulomb.append(coulomb_sum);

		# Maybe return a source_object, receiver_object, and lists of normal, shear, coulomb values. 
	return [inputs.source_object, inputs.receiver_object, receiver_normal, receiver_shear, receiver_coulomb]; 
	


"""



"""


# The purpose of these functions is to read/write a convenient input file
# Since we're using an input format that involves lon/lat/depth/mag, we have to convert
# We will use Wells and Coppersmith 1994, in the same way that Coulomb does. 
# We have to do a few other conversions too. 
# This input file assumes a fixed rake specified initially in the configure file, same as the .inp file format

import numpy as np
import matplotlib.pyplot as plt 
import haversine
import conversion_math
import coulomb_collections
import wells_and_coppersmith
import sys


def read_intxt(input_file):
	print("Reading source and receiver fault information from file %s " % input_file);
	
	strike_src=[]; dipangle_src=[]; rake_src=[]; magnitude_src=[]; faulting_type_src=[]; fault_center_lon_src=[]; fault_center_lat_src=[]; fault_dep_src=[];
	strike_rec=[]; dipangle_rec=[]; rake_rec=[]; length_rec=[]; width_rec=[]; updip_corner_lon_rec=[]; updip_corner_lat_rec=[]; updip_corner_dep_rec=[];
	
	ifile=open(input_file,'r');
	for line in ifile:
		temp=line.split();
		if len(temp)==0:
			continue;
		if temp[0]=='S:':
			[strike,rake,dip,magnitude,faulting_type,fault_center_lon,fault_center_lat,fault_center_depth]=read_source_line(line);
			strike_src.append(strike);
			rake_src.append(rake);
			dipangle_src.append(dip);
			magnitude_src.append(magnitude);
			faulting_type_src.append(faulting_type);
			fault_center_lon_src.append(fault_center_lon);
			fault_center_lat_src.append(fault_center_lat);
			fault_dep_src.append(fault_center_depth);
		elif temp[0]=='R:':
			[strike,rake,dip,length,width,updip_corner_lon,updip_corner_lat,updip_corner_dep]=read_receiver_line(line);
			strike_rec.append(strike);
			rake_rec.append(rake);
			dipangle_rec.append(dip);
			length_rec.append(length);
			width_rec.append(width);
			updip_corner_lon_rec.append(updip_corner_lon);
			updip_corner_lat_rec.append(updip_corner_lat);
			updip_corner_dep_rec.append(updip_corner_dep);
		elif temp[0]=='G:':
			[PR1,FRIC,minlon,maxlon,zerolon,minlat,maxlat,zerolat]=read_general_line(line);
		else:
			continue;
	ifile.close();

	# The main computing functions. 
	[start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc] = compute_grid_parameters(minlon, maxlon, zerolon, minlat, maxlat, zerolat);
	[xstart_src, xfinish_src, ystart_src, yfinish_src, Kode_src, rtlat_src, reverse_src, top_src, bottom_src, comment_src] = compute_params_for_source(strike_src, dipangle_src, rake_src, magnitude_src, faulting_type_src, fault_center_lon_src, fault_center_lat_src, fault_dep_src, zerolon, zerolat);
	[xstart_rec, xfinish_rec, ystart_rec, yfinish_rec, Kode_rec, rtlat_rec, reverse_rec, top_rec, bottom_rec, comment_rec] = compute_params_for_receiver(strike_rec, dipangle_rec, rake_rec, length_rec, width_rec, updip_corner_lon_rec, updip_corner_lat_rec, updip_corner_dep_rec, zerolon, zerolat);

	receivers=coulomb_collections.Faults_object(xstart=xstart_rec, xfinish=xfinish_rec, ystart=ystart_rec, yfinish=yfinish_rec, Kode=Kode_rec, rtlat=rtlat_rec, reverse=reverse_rec, potency=[], strike=strike_rec, dipangle=dipangle_rec, rake=rake_rec, top=top_rec, bottom=bottom_rec, comment=comment_rec);
	sources=coulomb_collections.Faults_object(xstart=xstart_src, xfinish=xfinish_src, ystart=ystart_src, yfinish=yfinish_src, Kode=Kode_src, rtlat=rtlat_src, reverse=reverse_src, potency=[], strike=strike_src, dipangle=dipangle_src, rake=rake_src, top=top_src, bottom=bottom_src, comment=comment_src);

	input_obj=coulomb_collections.Input_object(PR1=PR1,FRIC=FRIC, depth=0, start_gridx=start_gridx, finish_gridx=finish_gridx, start_gridy=start_gridy, finish_gridy=finish_gridy, 
		xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon, minlat=minlat, maxlat=maxlat, zerolat=zerolat, eqlon=fault_center_lon_src[0], eqlat=fault_center_lat_src[0], 
		receiver_object=receivers, source_object=sources);

	return input_obj;





def read_source_line(line):
	strike=float(line.split()[1]);
	rake=float(line.split()[2]);
	dip=float(line.split()[3]);
	magnitude=float(line.split()[4]);
	faulting_type=line.split()[5];
	fault_center_lon=float(line.split()[6]);
	fault_center_lat=float(line.split()[7]);
	fault_center_dep=float(line.split()[8]);
	return [strike,rake,dip,magnitude,faulting_type,fault_center_lon,fault_center_lat,fault_center_dep];
def read_receiver_line(line):
	strike=float(line.split()[1]);
	rake=float(line.split()[2]);
	dip=float(line.split()[3]);
	length=float(line.split()[4]);
	width=float(line.split()[5]);
	updip_corner_lon=float(line.split()[6]);
	updip_corner_lat=float(line.split()[7]);
	updip_corner_dep=float(line.split()[8]);	
	return [strike,rake,dip,length,width,updip_corner_lon,updip_corner_lat,updip_corner_dep];
def read_general_line(line):
	PR1=float(line.split()[1]);
	FRIC=float(line.split()[2]);
	lon_min=float(line.split()[3]);
	lon_max=float(line.split()[4]);
	lon_zero=float(line.split()[5]);
	lat_min=float(line.split()[6]);
	lat_max=float(line.split()[7]);
	lat_zero=float(line.split()[8]);
	return [PR1,FRIC,lon_min,lon_max,lon_zero,lat_min,lat_max,lat_zero];





def compute_grid_parameters(minlon, maxlon, zerolon, minlat, maxlat, zerolat):
	start_gridx=0; start_gridy=0; finish_gridx=0; finish_gridy=0; xinc=0; yinc=0;
	# Compute the grid parameters that we'll be using, based on the size of the map. 
	# Assuming 100 increments in both directions. 

	deltalon=(maxlon-minlon)*111.00*np.cos(np.deg2rad(zerolat));  # in km. 
	deltalat=(maxlat-minlat)*111.00;  # in km. 
	start_gridx=-deltalon/2.0;
	finish_gridx=deltalon/2.0;
	start_gridy=-deltalat/2.0;
	finish_gridy=deltalat/2.0;
	xinc=deltalon/100.0;
	yinc=deltalat/100.0;

	# print("start_gridx: %f " % start_gridx);
	# print("finish_gridx: %f " % finish_gridx);
	# print("start_gridy: %f " % start_gridy);
	# print("finish_gridy: %f " % finish_gridy);
	# print("xinc: %f " % xinc);
	# print("yinc: %f " % yinc);
	return [start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc];


def compute_params_for_source(strike, dip, rake, mag, eq_type, eqlon, eqlat, eqdep, zerolon, zerolat):
	xstart=[]; xfinish=[]; ystart=[]; yfinish=[]; Kode=[]; rtlat=[]; reverse=[]; top=[]; bottom=[]; comment=[];
	# Useful when you have earthquake catalog information, and you want the source information. 
	# The depth/eqlon/eqlat parameters refer to the center of the fault plane by assumption. 
	# Takes lists of parameters. 
	for i in range(len(strike)):
		[xcenter,ycenter]=conversion_math.latlon2xy(eqlon[i],eqlat[i],zerolon,zerolat);
		# print("EQ location: %f %f " % (eqlon[i],eqlat[i]) );
		# print("Ref location: %f %f " % (zerolon, zerolat) );
		# print("Coordinates: %f %f " % (x,y) );
		L=wells_and_coppersmith.RLD_from_M(mag[i],eq_type[i]);  # rupture length
		W=wells_and_coppersmith.RW_from_M(mag[i],eq_type[i]);   # rupture width
		slip = wells_and_coppersmith.rectangular_slip(L*1000,W*1000,mag[i]);  # must input in meters
		print("Fault slip: %f" % slip);
		
		# xistart,yistart=conversion_math.add_vector_to_point(xcenter,ycenter,-L/2,strike[i]);  # if the hypocenter is really the center of the rupture
		xistart,yistart=conversion_math.add_vector_to_point(xcenter,ycenter,0,strike[i]); # if the hypocenter is on one side of the rupture
		xifinish,yifinish=conversion_math.add_vector_to_point(xcenter,ycenter,L,strike[i]);
		rtlati,reversei=conversion_math.get_rtlat_dip_slip(slip, rake[i]);
		topi, bottomi=conversion_math.get_top_bottom(eqdep[i],W,dip[i]);


		xstart.append(xistart);
		ystart.append(yistart);
		xfinish.append(xifinish);
		yfinish.append(yifinish);
		Kode.append(100);
		rtlat.append(rtlati);
		reverse.append(reversei);
		top.append(topi);
		bottom.append(bottomi);
		comment.append('');

	return [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom, comment];


def compute_params_for_receiver(strike, dip, rake, length, width, lon, lat, depth, zerolon, zerolat):
	xstart=[]; xfinish=[]; ystart=[]; yfinish=[]; Kode=[]; rtlat=[]; reverse=[]; top=[]; bottom=[]; comment=[];
	# Useful when you have a receiver and you want to compute stresses on it. 
	# The depth/lon/lat parameters refer to the updip fault corner where you're looking down the strike direction, and the dip is off to your right. 
	# Takes lists of parameters. 
	for i in range(len(strike)):

		[xistart,yistart]=conversion_math.latlon2xy(lon[i],lat[i],zerolon,zerolat);
		xifinish,yifinish=conversion_math.add_vector_to_point(xistart,yistart,length[i],strike[i]);
		topi, bottomi=conversion_math.get_top_bottom_from_top(depth[i],width[i],dip[i]);

		xstart.append(xistart);
		ystart.append(yistart);
		xfinish.append(xifinish);
		yfinish.append(yifinish);
		Kode.append(100);
		rtlat.append(0);
		reverse.append(0);
		top.append(topi);
		bottom.append(bottomi);
		comment.append('');

	return [xstart, xfinish, ystart, yfinish, Kode, rtlat, reverse, top, bottom, comment];









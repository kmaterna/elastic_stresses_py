
import numpy as np 
import haversine
import conversion_math
import coulomb_collections
import conversion_math
import io_intxt


def read_intxt(input_file):
	print("Reading source and receiver fault information from file %s " % input_file);
	strike_src=[]; dipangle_src=[]; rake_src=[]; magnitude_src=[]; lon_src=[]; lat_src=[]; depth_src=[];
	strike_rec=[]; dipangle_rec=[]; rake_rec=[]; length_rec=[]; width_rec=[]; updip_corner_lon_rec=[]; updip_corner_lat_rec=[]; updip_corner_dep_rec=[];

	ifile=open(input_file,'r');
	for line in ifile:
		temp=line.split();
		if len(temp)==0:
			continue;
		if temp[0]=='S:':
			[strike,rake,dip,lon,lat,depth,magnitude, mu, lame1]=read_point_source_line(line);
			strike_src.append(strike);
			rake_src.append(rake);
			dipangle_src.append(dip);
			magnitude_src.append(magnitude);
			lon_src.append(lon);
			lat_src.append(lat);
			depth_src.append(depth);
		elif temp[0]=='R:':
			[strike,rake,dip,length,width,updip_corner_lon,updip_corner_lat,updip_corner_dep]=io_intxt.read_receiver_line(line);
			strike_rec.append(strike);
			rake_rec.append(rake);
			dipangle_rec.append(dip);
			length_rec.append(length);
			width_rec.append(width);
			updip_corner_lon_rec.append(updip_corner_lon);
			updip_corner_lat_rec.append(updip_corner_lat);
			updip_corner_dep_rec.append(updip_corner_dep);
		elif temp[0]=='G:':
			[PR1,FRIC,minlon,maxlon,zerolon,minlat,maxlat,zerolat]=io_intxt.read_general_line(line);
		else:
			continue;
	ifile.close();


	# The main computing functions. 
	[start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc] = io_intxt.compute_grid_parameters(minlon, maxlon, zerolon, minlat, maxlat, zerolat);	
	[xstart_rec, xfinish_rec, ystart_rec, yfinish_rec, Kode_rec, rtlat_rec, reverse_rec, top_rec, bottom_rec, comment_rec] = io_intxt.compute_params_for_receiver(strike_rec, dipangle_rec, rake_rec, length_rec, width_rec, updip_corner_lon_rec, updip_corner_lat_rec, updip_corner_dep_rec, zerolon, zerolat);

	# The point source formatting
	[x_src, y_src, Kode_src, rtlat_src, reverse_src, potency_src, comment_src] = compute_params_for_point_source(strike_src, dipangle_src, rake_src, magnitude_src, lon_src, lat_src, zerolon, zerolat, mu, lame1);
	sources=coulomb_collections.Faults_object(xstart=x_src, xfinish=x_src, ystart=y_src, yfinish=y_src, Kode=Kode_src, rtlat=rtlat_src, reverse=reverse_src, potency=potency_src, strike=strike_src, dipangle=dipangle_src, rake=rake_src, top=depth_src, bottom=depth_src, comment=comment_src);

	# Packaging
	receivers=coulomb_collections.Faults_object(xstart=xstart_rec, xfinish=xfinish_rec, ystart=ystart_rec, yfinish=yfinish_rec, Kode=Kode_rec, rtlat=rtlat_rec, reverse=reverse_rec, potency=[], strike=strike_rec, dipangle=dipangle_rec, rake=rake_rec, top=top_rec, bottom=bottom_rec, comment=comment_rec);	
	input_obj=coulomb_collections.Input_object(PR1=PR1,FRIC=FRIC, depth=0, start_gridx=start_gridx, finish_gridx=finish_gridx, start_gridy=start_gridy, finish_gridy=finish_gridy, 
		xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon, minlat=minlat, maxlat=maxlat, zerolat=zerolat, eqlon=lon_src[0], eqlat=lat_src[0], 
		receiver_object=receivers, source_object=sources);

	return input_obj;


def read_point_source_line(line):
	# Format: strike rake dip lon lat depth magnitude mu lamdba
	strike=float(line.split()[1]);
	rake=float(line.split()[2]);
	dip=float(line.split()[3]);
	lon=float(line.split()[4]);
	lat=float(line.split()[5]);
	depth=float(line.split()[6]);
	magnitude=float(line.split()[7]);
	mu=float(line.split()[8]);
	lame1=float(line.split()[9]);
	return [strike,rake,dip,lon,lat,depth,magnitude,mu,lame1];


def compute_params_for_point_source(strike_src, dipangle_src, rake_src, magnitude_src, lon_src, lat_src, zerolon, zerolat, mu, lame1):
	# Given information about point sources from focal mechanisms,
	# Return the right components that get packaged into input_obj. 
	x_src=[]; y_src=[]; Kode_src=[]; rtlat_src=[]; reverse_src=[]; potency_src=[]; comment=[];
	for i in range(len(strike_src)):
		[xcenter,ycenter]=conversion_math.latlon2xy(lon_src[i],lat_src[i],zerolon,zerolat);
	
		potency = get_DC_potency(rake_src[i], magnitude_src[i], mu, lame1);
		# The point source
		x_src.append(xcenter);
		y_src.append(ycenter);
		potency_src.append(potency); 

		# Filler variables for the point source case
		Kode_src.append(100);
		rtlat_src.append(0);
		reverse_src.append(0); 
		comment.append('');

	return [x_src, y_src, Kode_src, rtlat_src, reverse_src, potency_src, comment];


def get_DC_potency(rake, momentmagnitude, mu, lame1):
	# Given the basic double couple parameters, 
	# Return the four-vector used in Okada DC3D0. 
	# Pot1 = strike-slip moment of DC / mu
	# Pot2 = dip-slip moment of DC / mu
	# Pot3 = inflation = M_ISO / lambda
	# Pot4 = tensile = M_lineardipole / mu 
	# In a more general case, we would use a different MT format to handle non-DC parts. 
	# Right now, it only handles DC focal mechanisms. 
	total_moment = moment(momentmagnitude);
	strike_slip_fraction, dip_slip_fraction = conversion_math.get_rtlat_dip_slip(1.0, rake);
	print("strike_slip: ",strike_slip_fraction);
	print("dip_slip: ",dip_slip_fraction);
	strike_slip_fraction=-1*strike_slip_fraction; # DC3D0 wants left lateral slip. 
	p1 = total_moment*strike_slip_fraction/ mu;
	p2 = total_moment*dip_slip_fraction/ mu;
	# In the double-couple case, this is zero. 
	p3 = 0; 
	p4 = 0; 
	return [p1, p2, p3, p4];


def moment(Mw):
	# Current definition of moment magnitude, returning moment in newton meters
	exponent = 1.5*Mw + 1.5*10.7;
	moment = np.power(10, exponent);
	moment_newton_meters = moment * 1e-7;
	return moment_newton_meters;

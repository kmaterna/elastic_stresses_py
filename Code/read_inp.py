# Read Coulomb input files in the .inp format
# Important parameters: 
# 1. Poisson's Ratio
# 2. Coefficient of Friction
# 3. Depth
# 4. Fault params (x start, x finish, y start, y finish, Kode, rt. lat, reverse, dip angle, top, bottom, comment)
# 5. Grid parameters (start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc)
# 6. Map info (min lon, max lon, zero lon, min lat, max lat, zero lat)


import collections
import conversion_math

Input_object = collections.namedtuple('Input_object',
	['PR1','FRIC','depth','start_gridx', 'finish_gridx', 'start_gridy', 'finish_gridy', 'xinc', 'yinc', 'minlon','maxlon','zerolon','minlat','maxlat','zerolat','source_object','receiver_object'])
Faults_object = collections.namedtuple('Faults_object',
	['xstart','xfinish','ystart','yfinish','Kode','rtlat','reverse','strike','dipangle','rake','top','bottom','comment']);
# These objects will be common between the inp and inr formats. 


def read_inp(input_file,fixed_rake):
	# inp files require a fixed rake for the receiver faults, since they don't provide a fault-specific one inside the input file. 
	read_faults=0;
	read_maps=0;
	read_grid=0;
	minlon=[]; maxlon=[]; zerolon=[]; minlat=[]; maxlat=[]; zerolat=[]; start_grix=[]; start_gridy=[]; finish_gridx=[]; finish_gridy=[]; xinc=[]; yinc=[];
	xstart_rec=[]; xfinish_rec=[]; ystart_rec=[]; yfinish_rec=[]; Kode_rec=[]; rtlat_rec=[]; reverse_rec=[]; strike_rec=[]; dipangle_rec=[]; rake_rec=[]; top_rec=[]; bottom_rec=[]; comment_rec=[];
	xstart_src=[]; xfinish_src=[]; ystart_src=[]; yfinish_src=[]; Kode_src=[]; rtlat_src=[]; reverse_src=[]; strike_src=[]; dipangle_src=[]; rake_src=[]; top_src=[]; bottom_src=[]; comment_src=[];


	ifile=open(input_file);
	for line in ifile:
		temp=line.split();
		if 'PR1=' in line:
			prspot=temp.index('PR1=');
			PR1=float(temp[prspot+1]);
			dpspot=temp.index('DEPTH=');
			depth=float(temp[dpspot+1]);
		if 'FRIC=' in line:
			fric_place = temp.index('FRIC=');
			FRIC=float(temp[fric_place+1]);
		if 'xxxxxxxxxx' in line: # Moving into the fault definitions
			read_faults=1;
			continue;
		if read_faults==1:
			if len(temp)==0:
				read_faults=0;
				continue;  # Getting out of fault definitions
			else:
				slip = abs(float(temp[6]))+abs(float(temp[7]));
				if slip>0.0000001:  # Here we have a source fault
					xstart_src.append(float(temp[1]));
					ystart_src.append(float(temp[2]));
					xfinish_src.append(float(temp[3]));
					yfinish_src.append(float(temp[4]));
					Kode_src.append(int(temp[5]));
					rtlat_src.append(float(temp[6]));
					reverse_src.append(float(temp[7]));
					dipangle_src.append(float(temp[8]));
					top_src.append(float(temp[9]));
					bottom_src.append(float(temp[10]));
					comment_src.append(" ".join(temp[11:-1]));
					deltax=xfinish_src[-1]-xstart_src[-1];
					deltay=yfinish_src[-1]-ystart_src[-1];
					strike = conversion_math.get_strike(deltax, deltay);
					strike_src.append(strike);
					rake = conversion_math.get_rake(rtlat_src[-1],reverse_src[-1]);
					rake_src.append(rake);
				else:  # here we have a receiver fault
					xstart_rec.append(float(temp[1]));
					ystart_rec.append(float(temp[2]));
					xfinish_rec.append(float(temp[3]));
					yfinish_rec.append(float(temp[4]));
					Kode_rec.append(int(temp[5]));
					rtlat_rec.append(float(temp[6]));
					reverse_rec.append(float(temp[7]));
					dipangle_rec.append(float(temp[8]));
					top_rec.append(float(temp[9]));
					bottom_rec.append(float(temp[10]));
					comment_rec.append(" ".join(temp[11:-1]));
					deltax=xfinish_rec[-1]-xstart_rec[-1];
					deltay=yfinish_rec[-1]-ystart_rec[-1];
					strike = conversion_math.get_strike(deltax, deltay);
					strike_rec.append(strike);
					rake = fixed_rake;
					rake_rec.append(rake);
		# Moving into Grid Parameters and Map Information
		if 'Grid Parameters' in line:
			read_grid=1;
			continue;
		if read_grid==1:
			if '  1  -' in line:
				start_gridx=float(temp[-1]);
			elif '  2  -' in line:
				start_gridy=float(temp[-1]);
			elif '  3  -' in line:
				finish_gridx=float(temp[-1]);
			elif '  4  -' in line:
				finish_gridy=float(temp[-1]);
			elif '  5  -' in line:
				xinc=float(temp[-1]);
			elif '  6  -' in line:
				yinc=float(temp[-1]);
			else:
				read_grid=0;
				continue;

		# Reading Map Information
		if 'Map info' in line:
			read_maps=1;
			continue;
		if read_maps==1:
			if '  1  -' in line:
				minlon=float(temp[-1]);
			elif '  2  -' in line:
				maxlon=float(temp[-1]);
			elif '  3  -' in line:
				zerolon=float(temp[-1]);
			elif '  4  -' in line:
				minlat=float(temp[-1]);
			elif '  5  -' in line:
				maxlat=float(temp[-1]);
			elif '  6  -' in line:
				zerolat=float(temp[-1]);
			else:
				read_maps=0;
				continue;

	ifile.close();

	receivers=Faults_object(xstart=xstart_rec, xfinish=xfinish_rec, ystart=ystart_rec, yfinish=yfinish_rec, Kode=Kode_rec, rtlat=rtlat_rec, reverse=reverse_rec, strike=strike_rec, dipangle=dipangle_rec, rake=rake_rec, top=top_rec, bottom=bottom_rec, comment=comment_rec);
	sources=Faults_object(xstart=xstart_src, xfinish=xfinish_src, ystart=ystart_src, yfinish=yfinish_src, Kode=Kode_src, rtlat=rtlat_src, reverse=reverse_src, strike=strike_src, dipangle=dipangle_src, rake=rake_src, top=top_src, bottom=bottom_src, comment=comment_src);

	input_obj=Input_object(PR1=PR1,FRIC=FRIC, depth=depth, start_gridx=start_gridx, finish_gridx=finish_gridx, start_gridy=start_gridy, finish_gridy=finish_gridy, 
		xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon, minlat=minlat, maxlat=maxlat, zerolat=zerolat, receiver_object=receivers, source_object=sources);

	return input_obj;



def write_inp(outfile, data_obj):
	sources=data_obj.source_object;
	receivers=data_obj.receiver_object;
	num_faults=len(sources.xstart)+len(receivers.xstart);
	ofile=open(outfile,'w');
	ofile.write('This is a test file for Coulomb 3.0\n');
	ofile.write('This file is prepared by Python\n');
	ofile.write('#reg1=  0  #reg2=  0   #fixed=  %d  sym=  1\n' % (num_faults) );
	ofile.write(' PR1=       %4.3f     PR2=       .250    DEPTH=        %.1f\n' % (data_obj.PR1, data_obj.depth) );
	ofile.write('  E1=   0.800000E+06   E2=   0.800000E+06\n'); # The young's modulus, in bar
	ofile.write('XSYM=       .000     YSYM=       .000\n'); # An old pattern- ignore
	ofile.write('FRIC=       %4.3f\n' % (data_obj.FRIC) );
	ofile.write('S1DR=    19.0001     S1DP=     -0.0001    S1IN=    100.000     S1GD=   .000000\n');
	ofile.write('S3DR=    89.9999     S3DP=      89.999    S3IN=     30.000     S3GD=   .000000\n');
	ofile.write('S2DR=   109.0001     S2DP=     -0.0001    S2IN=      0.000     S2GD=   .000000\n\n');
	ofile.write('  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot\n');
	ofile.write('xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n');
	for i in range(len(sources.xstart)):
		ofile.write('  1 % 10.4f % 10.4f % 10.4f % 10.4f %3d % 10.4f % 10.4f % 10.4f % 10.2f % 10.2f %s\n' % (sources.xstart[i], sources.ystart[i], sources.xfinish[i], sources.yfinish[i], 
			sources.Kode[i], sources.rtlat[i], sources.reverse[i], sources.dipangle[i], sources.top[i], sources.bottom[i], sources.comment[i]) );
	for i in range(len(receivers.xstart)):
		ofile.write('  1 % 10.4f % 10.4f % 10.4f % 10.4f %3d % 10.4f % 10.4f % 10.4f % 10.2f % 10.2f %s\n' % (receivers.xstart[i], receivers.ystart[i], receivers.xfinish[i], receivers.yfinish[i], 
			receivers.Kode[i], receivers.rtlat[i], receivers.reverse[i], receivers.dipangle[i], receivers.top[i], receivers.bottom[i], receivers.comment[i]) );
	ofile.write('\n');
	ofile.write('    Grid Parameters\n');
	ofile.write('  1  ----------------------------  Start-x =    % 10.5f\n' % data_obj.start_gridx);
	ofile.write('  2  ----------------------------  Start-y =    % 10.5f\n' % data_obj.start_gridy);
	ofile.write('  3  --------------------------   Finish-x =    % 10.5f\n' % data_obj.finish_gridx);
	ofile.write('  4  --------------------------   Finish-y =    % 10.5f\n' % data_obj.finish_gridy);
	ofile.write('  5  ------------------------  x-increment =    % 10.5f\n' % data_obj.xinc);
	ofile.write('  6  ------------------------  y-increment =    % 10.5f\n' % data_obj.yinc);
	ofile.write('     Size Parameters\n');
	ofile.write('  1  --------------------------  Plot size =     2.000000\n');
	ofile.write('  2  --------------  Shade/Color increment =     1.000000\n');
	ofile.write('  3  ------  Exaggeration for disp.& dist. =     10000.00\n\n');
	ofile.write('Cross section default\n');  # Because the python program doesn't do cross sections yet, I'm leaving this hard-coded. 
	ofile.write('  1  ----------------------------  Start-x =    -36.00000\n');
	ofile.write('  2  ----------------------------  Start-y =     36.00000\n');
	ofile.write('  3  --------------------------   Finish-x =     38.00000\n');
	ofile.write('  4  --------------------------   Finish-y =    -36.00000\n');
	ofile.write('  5  ------------------  Distant-increment =     1.000000\n');
	ofile.write('  6  ----------------------------  Z-depth =     30.00000\n');
	ofile.write('  7  ------------------------  Z-increment =     1.000000\n');
	ofile.write('     Map infomation\n');
	ofile.write('  1  ---------------------------- min. lon =    % 9.4f \n' % data_obj.minlon);
	ofile.write('  2  ---------------------------- max. lon =    % 9.4f \n' % data_obj.maxlon);
	ofile.write('  3  ---------------------------- zero lon =    % 9.4f \n' % data_obj.zerolon);
	ofile.write('  4  ---------------------------- min. lat =    % 9.4f \n' % data_obj.minlat);
	ofile.write('  5  ---------------------------- max. lat =    % 9.4f \n' % data_obj.maxlat);
	ofile.write('  6  ---------------------------- zero lat =    % 9.4f \n' % data_obj.zerolat);
	ofile.close();
	return;




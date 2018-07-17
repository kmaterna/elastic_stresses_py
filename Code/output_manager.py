# output_manager

import numpy as np 
import matplotlib.pyplot as plt 
from subprocess import call
import conversion_math


def produce_outputs(params, inputs, out_object):
	# Do the rest. 
	call(['mkdir','-p',params.outdir],shell=False);
	surface_def_plot(params,inputs,out_object);

	stress_plot();  # TOMORROW will write a plot with polygons colored by stress. 'Coulomb' vs. 'shear' vs. 'normal' will be an argument
	map_plot();
	return;



def get_plotting_traces(inputs):
	# Get the updip and total traces of the input file faults, for plotting on an initial map view. 

	src_total_x = [];
	src_total_y = [];
	src_updip_x = [];
	src_updip_y = [];
	rec_total_x = [];
	rec_total_y = [];
	rec_updip_x = [];
	rec_updip_y = [];

	# The surface traces of the source faults
	for i in range(len(inputs.source_object.xstart)):
		W = conversion_math.get_downdip_width(inputs.source_object.top[i],inputs.source_object.bottom[i],inputs.source_object.dipangle[i]);
		depth       = inputs.source_object.top[i];
		strike      = inputs.source_object.strike[i];
		dip         = inputs.source_object.dipangle[i];

		updip_point0 = [inputs.source_object.xstart[i],inputs.source_object.ystart[i]];
		updip_point1 = [inputs.source_object.xfinish[i],inputs.source_object.yfinish[i]];
		vector_mag = W*np.cos(np.deg2rad(inputs.source_object.dipangle[i]));  # how far the bottom edge is displaced downdip from map-view
		downdip_point0 = conversion_math.add_vector_to_point(inputs.source_object.xstart[i],inputs.source_object.ystart[i],vector_mag, strike+90);  # strike+90 = downdip direction. 
		downdip_point1 = conversion_math.add_vector_to_point(inputs.source_object.xfinish[i],inputs.source_object.yfinish[i], vector_mag, strike+90);

		x_total = [updip_point0[0], updip_point1[0], downdip_point1[0], downdip_point0[0],updip_point0[0]];
		y_total = [updip_point0[1], updip_point1[1], downdip_point1[1], downdip_point0[1],updip_point0[1]];
		x_updip = [updip_point0[0], updip_point1[0]];
		y_updip = [updip_point0[1], updip_point1[1]];

		src_total_x.append(x_total);
		src_total_y.append(y_total);
		src_updip_x.append(x_updip);
		src_updip_y.append(y_updip);

	# The surface traces of the receiver faults
	for i in range(len(inputs.receiver_object.xstart)):
		W = conversion_math.get_downdip_width(inputs.receiver_object.top[i],inputs.receiver_object.bottom[i],inputs.receiver_object.dipangle[i]);
		depth       = inputs.receiver_object.top[i];
		strike      = inputs.receiver_object.strike[i];
		dip         = inputs.receiver_object.dipangle[i];

		updip_point0 = [inputs.receiver_object.xstart[i],inputs.receiver_object.ystart[i]];
		updip_point1 = [inputs.receiver_object.xfinish[i],inputs.receiver_object.yfinish[i]];
		vector_mag = W*np.cos(np.deg2rad(inputs.receiver_object.dipangle[i]));  # how far the bottom edge is displaced downdip from map-view
		downdip_point0 = conversion_math.add_vector_to_point(inputs.receiver_object.xstart[i],inputs.receiver_object.ystart[i],vector_mag, strike+90);  # strike+90 = downdip direction. 
		downdip_point1 = conversion_math.add_vector_to_point(inputs.receiver_object.xfinish[i],inputs.receiver_object.yfinish[i], vector_mag, strike+90);

		x_total = [updip_point0[0], updip_point1[0], downdip_point1[0], downdip_point0[0],updip_point0[0]];
		y_total = [updip_point0[1], updip_point1[1], downdip_point1[1], downdip_point0[1],updip_point0[1]];
		x_updip = [updip_point0[0], updip_point1[0]];
		y_updip = [updip_point0[1], updip_point1[1]];

		rec_total_x.append(x_total);
		rec_total_y.append(y_total);
		rec_updip_x.append(x_updip);
		rec_updip_y.append(y_updip);

	return [src_total_x, src_total_y, src_updip_x, src_updip_y, rec_total_x, rec_total_y, rec_updip_x, rec_updip_y];



def surface_def_plot(params, inputs, out_object):

	# Get the updip and downdip traces of the faults for plotting purposes. 
	[src_total_x, src_total_y, src_updip_x, src_updip_y, rec_total_x, rec_total_y, rec_updip_x, rec_updip_y] = get_plotting_traces(inputs);

	# Write output file
	ofile=open(params.outdir+'disps.txt','w')
	for i in np.arange(0,len(out_object.y)):
		for j in np.arange(0,len(out_object.x)):
			ofile.write("%f %f %f %f %f\n" % (out_object.x2d[i][j], out_object.y2d[i][j], out_object.u_disp[i][j], out_object.v_disp[i][j], out_object.w_disp[i][j]));	
	ofile.close();

	print("Max vertical displacement is %f m" % (out_object.w_disp.max()) );
	print("Making displacement plot.")

	#Plot of elastic surface deformation from the given input model. 
	plt.figure(figsize=(16,16))
	plt.pcolormesh(out_object.x2d, out_object.y2d, out_object.w_disp);
	plt.colorbar()
	for i in np.arange(0,len(out_object.y),5):
		for j in np.arange(0,len(out_object.x),5):
			plt.quiver(out_object.x2d[i][j],out_object.y2d[i][j],out_object.u_disp[i][j],out_object.v_disp[i][j],units='width',scale=0.2)
	for i in range(len(src_total_x)):
		plt.plot(src_total_x[i], src_total_y[i],'k',linewidth=1);
		plt.plot(src_updip_x[i], src_updip_y[i],'g',linewidth=3);
	for i in range(len(rec_total_x)):
		plt.plot(rec_total_x[i], rec_total_y[i],'b',linewidth=1);
		plt.plot(rec_updip_x[i], rec_updip_y[i],'b',linewidth=3);
	for i in range(len(inputs.receiver_object.xstart)):
		center=conversion_math.get_fault_center(inputs.receiver_object,i);
		plt.plot(center[0],center[1],'.b',markersize=8);
	for i in range(len(inputs.source_object.xstart)):
		center=conversion_math.get_fault_center(inputs.source_object,i);
		plt.plot(center[0],center[1],'.g',markersize=8);		
	plt.xlim([out_object.x.min(),out_object.x.max()])
	plt.ylim([out_object.y.min(),out_object.y.max()])
	plt.grid()
	plt.title('Surface Dipslacement')
	plt.savefig(params.outdir+"Displacement_model.eps")
	plt.close();

	return;



def stress_plot():
	return;



def map_plot():
	return


# output_manager

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from subprocess import call
import collections
import conversion_math


# For reference: 
Out_object = collections.namedtuple('Out_object',
	['x','y','x2d','y2d','u_disp','v_disp','w_disp','source_object','receiver_object','receiver_normal','receiver_shear','receiver_coulomb']);


def produce_outputs(params, out_object):
	call(['mkdir','-p',params.outdir],shell=False);
	surface_def_plot(params,out_object);
	# stress_plot(params,out_object,'shear');  
	# stress_plot(params,out_object,'normal');
	stress_plot(params,out_object,'coulomb');
	map_plot();
	write_output_files(params,out_object);

	return;



def get_plotting_traces(fault_object):
	# Get the updip and total traces of the input file faults, for plotting on map view. 
	# Each element of total_x and total_y contain 4 coordinates. 
	total_x = []; total_y = []; updip_x = []; updip_y = [];
	for i in range(len(fault_object.xstart)):
		[x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(fault_object, i);
		total_x.append(x_total);  # a 1x4 list
		total_y.append(y_total);
		updip_x.append(x_updip);  # a 1x2 list
		updip_y.append(y_updip);

	return [total_x, total_y, updip_x, updip_y];



def surface_def_plot(params, out_object):

	# Get the updip and downdip traces of the faults for plotting purposes. 
	[src_total_x, src_total_y, src_updip_x, src_updip_y] = get_plotting_traces(out_object.source_object);
	[rec_total_x, rec_total_y, rec_updip_x, rec_updip_y] = get_plotting_traces(out_object.receiver_object);

	print("Max vertical displacement is %f m" % (out_object.w_disp.max()) );
	print("Making displacement plot.")

	#Plot of elastic surface deformation from the given input model. 
	plt.figure(figsize=(16,16))
	plt.pcolormesh(out_object.x2d, out_object.y2d, out_object.w_disp);
	cb = plt.colorbar();	
	cb.set_label('Vertical displacement (meters)',fontsize=22);
	for l in cb.ax.yaxis.get_ticklabels():
		l.set_size(18)	
	for i in np.arange(0,len(out_object.y),5):
		for j in np.arange(0,len(out_object.x),5):
			plt.quiver(out_object.x2d[i][j],out_object.y2d[i][j],out_object.u_disp[i][j],out_object.v_disp[i][j],units='width',scale=0.2)
	for i in range(len(src_total_x)):
		plt.plot(src_total_x[i], src_total_y[i],'k',linewidth=1);
		plt.plot(src_updip_x[i], src_updip_y[i],'g',linewidth=3);
	for i in range(len(rec_total_x)):
		plt.plot(rec_total_x[i], rec_total_y[i],'b',linewidth=1);
		plt.plot(rec_updip_x[i], rec_updip_y[i],'b',linewidth=3);
	for i in range(len(out_object.receiver_object.xstart)):
		center=conversion_math.get_fault_center(out_object.receiver_object,i);
		plt.plot(center[0],center[1],'.b',markersize=8);
	for i in range(len(out_object.source_object.xstart)):
		center=conversion_math.get_fault_center(out_object.source_object,i);
		plt.plot(center[0],center[1],'.g',markersize=8);
	plt.xlim([out_object.x.min(),out_object.x.max()])
	plt.ylim([out_object.y.min(),out_object.y.max()])
	for l in plt.gca().yaxis.get_ticklabels():
		l.set_size(18);
	for l in plt.gca().xaxis.get_ticklabels():
		l.set_size(18);	
	plt.grid();
	plt.title('Surface Dipslacement',fontsize=28)
	plt.savefig(params.outdir+"Displacement_model.eps")
	plt.close();

	return;



def stress_plot(params, out_object, stress_type):
	# Here we will put plots of fault patches, colored by the magnitude of the stress component. 

	# Get the updip and downdip traces of the faults for plotting purposes. 
	[src_total_x, src_total_y, src_updip_x, src_updip_y] = get_plotting_traces(out_object.source_object);
	[rec_total_x, rec_total_y, rec_updip_x, rec_updip_y] = get_plotting_traces(out_object.receiver_object);

	if stress_type=='shear':
		stress_component=out_object.receiver_shear;
	elif stress_type=='normal':
		stress_component=out_object.receiver_normal;
	elif stress_type=='coulomb':
		stress_component=out_object.receiver_coulomb;
	else:
		print("Error! Invalid stress type : %s " % stress_type );

	# Select boundaries of color map. 
	# smallest_stress = np.min(stress_component);
	# largest_stress = np.max(stress_component);
	smallest_stress = -100;  # units: KPa
	largest_stress = 100;  # units: KPa
	color_boundary_object=matplotlib.colors.Normalize(vmin=smallest_stress,vmax=largest_stress, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='jet_r');

	# Figure of stresses. 
	plt.figure(figsize=(12,10));

	patches=[];
	colors=[];
	for i in range(len(stress_component)):
		xcoords=rec_total_x[i][0:4];
		ycoords=rec_total_y[i][0:4];
		fault_vertices=np.column_stack((xcoords, ycoords));
		patch_color=custom_cmap.to_rgba(stress_component[i]);
		
		mypolygon = Polygon(fault_vertices,color=patch_color,alpha=0.5);
		plt.gca().add_patch(mypolygon);

	custom_cmap.set_array(range(smallest_stress,largest_stress));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Kilopascals',fontsize=22);
	for l in cb.ax.yaxis.get_ticklabels():
		l.set_size(18)

	for i in range(len(src_total_x)):
		plt.plot(src_total_x[i], src_total_y[i],'k',linewidth=1);
		plt.plot(src_updip_x[i], src_updip_y[i],'g',linewidth=3);
	for i in range(len(rec_total_x)):
		plt.plot(rec_total_x[i], rec_total_y[i],'b',linewidth=1);
		plt.plot(rec_updip_x[i], rec_updip_y[i],'b',linewidth=3);
	for i in range(len(out_object.receiver_object.xstart)):
		center=conversion_math.get_fault_center(out_object.receiver_object,i);
		plt.plot(center[0],center[1],'.b',markersize=8);
	for i in range(len(out_object.source_object.xstart)):
		center=conversion_math.get_fault_center(out_object.source_object,i);
		plt.plot(center[0],center[1],'.g',markersize=8);			

	plt.grid();
	plt.axis('equal');
	for l in plt.gca().yaxis.get_ticklabels():
		l.set_size(18);
	for l in plt.gca().xaxis.get_ticklabels():
		l.set_size(18);	
	plt.title(stress_type+' stress from source faults',fontsize=22)
	plt.xlim([out_object.x.min(),out_object.x.max()])
	plt.ylim([out_object.y.min(),out_object.y.max()])		
	plt.savefig(params.outdir+'Stresses_'+stress_type+'.eps');
	plt.close();

	return;



def map_plot():
	# FUTURE FEATURE: Make Map

	return



def write_output_files(params, out_object):

	# Write displacement output file
	ofile=open(params.outdir+'disps.txt','w');
	ofile.write("Format: x y udisp vdisp wdisp (m) \n");
	for i in np.arange(0,len(out_object.y)):
		for j in np.arange(0,len(out_object.x)):
			ofile.write("%f %f %f %f %f\n" % (out_object.x2d[i][j], out_object.y2d[i][j], out_object.u_disp[i][j], out_object.v_disp[i][j], out_object.w_disp[i][j]));	
	ofile.close();

	# Write output file for stresses. 
	ofile=open(params.outdir+'stresses.txt','w');
	ofile.write("Format: centerx centery centerz normal shear coulomb (kpa)\n");
	for i in range(len(out_object.receiver_object.xstart)):
		# [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(fault_object, i);
		center=conversion_math.get_fault_center(out_object.receiver_object,i);
		ofile.write("%f %f %f %f %f %f \n" % (center[0], center[1], center[2], out_object.receiver_normal[i], out_object.receiver_shear[i], out_object.receiver_coulomb[i]) );
	ofile.close();
	print("Outputs written to file.")

	return;


# output_manager

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import pygmt
from subprocess import call
import coulomb_collections
import conversion_math
import io_inp
import io_additionals


def produce_outputs(params, inputs, disp_points, out_object):
	call(['mkdir','-p',params.outdir],shell=False);
	call(['cp','config.txt',params.outdir],shell=False);
	call(['cp',params.input_file,params.outdir],shell=False);	
	subfaulted_inputs=coulomb_collections.Input_object(PR1=inputs.PR1,FRIC=inputs.FRIC,depth=inputs.depth,
		start_gridx=inputs.start_gridx,start_gridy=inputs.start_gridy,finish_gridx=inputs.finish_gridx,
		finish_gridy=inputs.finish_gridy,xinc=inputs.xinc,yinc=inputs.yinc,minlon=inputs.minlon,maxlon=inputs.maxlon,
		zerolon=inputs.zerolon,minlat=inputs.minlat,maxlat=inputs.maxlat,zerolat=inputs.zerolat,eqlon=inputs.eqlon, eqlat=inputs.eqlat,
		source_object=out_object.source_object,receiver_object=out_object.receiver_object); # make a new object of the subfaulted configuration.
	io_inp.write_inp(params.outdir+'subfaulted.inp',subfaulted_inputs);
	surface_def_plot(params,out_object);
	stress_plot(params,out_object,'shear');  # can give vmin, vmax here if desired. 
	stress_plot(params,out_object,'normal');
	stress_plot(params,out_object,'coulomb');
	map_plot(params, inputs, out_object, 'coulomb');
	map_plot(params, inputs, out_object, 'normal');
	map_plot(params, inputs, out_object, 'shear');
	write_output_files(params,inputs, disp_points, out_object);
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

	print("Making plot of predicted displacement throughout model domain.")
	print("Max vertical displacement is %f m" % (out_object.w_disp.max()) );	

	#Plot of elastic surface deformation from the given input model. 
	plt.figure(figsize=(16,16))
	plt.pcolormesh(out_object.x2d, out_object.y2d, out_object.w_disp,cmap='jet');
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
	plt.axis('equal');
	plt.title('Surface Dipslacement',fontsize=28)
	plt.savefig(params.outdir+"Displacement_model_on_grid.eps")
	plt.close();
	return;



def stress_plot(params, out_object, stress_type, vmin="", vmax=""):
	# default vmin,vmax are in KPa
	# Here we will put plots of fault patches, colored by the magnitude of the stress component. 

	print("Making plot of %s stress on receiver fault patches. " % stress_type);
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
	if vmin=="":
		vmin = np.min(stress_component);
	if vmax=="":
		vmax = np.max(stress_component);
	color_boundary_object=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='RdYlBu_r');

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

	custom_cmap.set_array(np.arange(vmin,vmax,100));
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


def map_plot(params, inputs, out_object, stress_component):

	# Some options: 
	if stress_component=='shear':
		plotting_stress = out_object.receiver_shear;
		label='Shear';
	elif stress_component=='normal':
		plotting_stress = out_object.receiver_normal;
		label='Normal';
	else:
		plotting_stress = out_object.receiver_coulomb; # The default option
		label='Coulomb';

	# Make stress bounds for map. 
	stress_bounds=[abs(np.min(plotting_stress)),abs(np.max(plotting_stress))];
	stress_bound = np.max(stress_bounds);  # setting the scale to symmetric about zero
	smallest_stress = -stress_bound; # units: KPa
	largest_stress = stress_bound; # units: KPa
	smallest_stress = -1; # units: KPa
	largest_stress = 1; # units: KPa	

	color_boundary_object=matplotlib.colors.Normalize(vmin=smallest_stress,vmax=largest_stress, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='RdYlBu_r');

	# Make cpt
	pygmt.makecpt(C="jet",T=str(smallest_stress-0.1)+"/"+str(largest_stress+0.1)+"/0.05",H="mycpt.cpt",D=True);

	# Make Map
	region = [inputs.minlon, inputs.maxlon, inputs.minlat, inputs.maxlat];
	proj = "M7i"
	fig = pygmt.Figure()
	title="+t\""+stress_component+" stress\"";  # must put escaped quotations around the title. 
	fig.basemap(region=region,projection=proj,B=title);
	fig.coast(shorelines="1.0p,black",region=region,projection=proj,B="1.0");  # the boundary. 
	fig.coast(region=region,projection=proj,N='1',W='0.5p,black',S='white',L="g-125.5/39.6+c1.5+w50");

	# # Draw each source
	for i in range(len(out_object.source_object.xstart)):
		[x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(out_object.source_object,i);
		lons=[];
		lats=[];
		for j in range(len(x_total)):
			mylon, mylat = conversion_math.xy2lonlat(x_total[j],y_total[j],inputs.zerolon,inputs.zerolat);
			lons.append(mylon);
			lats.append(mylat);
		fig.plot( x=lons, y=lats, pen="thick,black");

	# Draw each receiver outline. This will eventually be fixed in a later pygmt, I hope. 
	for i in range(len(out_object.receiver_object.xstart)):
		lons=[];
		lats=[];
		[x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(out_object.receiver_object,i);
		for j in range(len(x_total)):
			mylon, mylat = conversion_math.xy2lonlat(x_total[j],y_total[j],inputs.zerolon,inputs.zerolat);
			lons.append(mylon);
			lats.append(mylat);
		fig.plot(x=lons, y=lats, pen="thick,black");  # Outline only

	# Draw data into each receiver. This will eventually be fixed in a later pygmt, I hope. 
	lons=[]; lats=[]; colors=[]; sizes=[];
	for i in range(len(out_object.receiver_object.xstart)):
		[x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(out_object.receiver_object,i);
		mylon_top, mylat_top = conversion_math.xy2lonlat(x_total[0],y_total[0],inputs.zerolon,inputs.zerolat);
		mylon_bot, mylat_bot = conversion_math.xy2lonlat(x_total[2],y_total[2],inputs.zerolon,inputs.zerolat);
		lons.append((mylon_top+mylon_bot)/2);
		lats.append((mylat_top+mylat_bot)/2);
		colors.append(plotting_stress[i]);
		sizes.append(5.0);
	fig.plot(x=lons, y=lats, color=colors, style="c0.8c", pen="thick,black", C="mycpt.cpt"); 

	# Colorbar annotation
	fig.colorbar(D="jBr+w3.5i/0.2i+o2.5c/1.5c+h",C="mycpt.cpt",I="0.8",G=str(smallest_stress)+"/"+str(largest_stress-0.1),B=["x"+str(0.2),"y+L\"KPa\""]); 

	# Annotate with earthquake location.
	fig.plot(inputs.eqlon, inputs.eqlat, style='s0.3c',G="purple",pen="thin,black");

	# Annotate with aftershock locations
	if len(params.aftershocks)>0:
		[lon, lat, depth, magnitude, time]=io_additionals.read_aftershock_table(params.aftershocks);
		for i in range(len(lon)):
			fig.plot(lon,lat,style='c0.1c',G='black',pen="thin,black");	

	fig.savefig(params.outdir+label+'_map.eps');
	plt.close();
	return


def write_output_files(params, inputs, disp_points, out_object):

	# Write displacement output file
	ofile=open(params.outdir+'disps_model_grid.txt','w');
	ofile.write("# Format: x y udisp vdisp wdisp (m) \n");
	for i in np.arange(0,len(out_object.y)):
		for j in np.arange(0,len(out_object.x)):
			ofile.write("%f %f %f %f %f\n" % (out_object.x2d[i][j], out_object.y2d[i][j], out_object.u_disp[i][j], out_object.v_disp[i][j], out_object.w_disp[i][j]));	
	ofile.close();

	# Write output file for stresses. 
	ofile=open(params.outdir+'stresses.txt','w');
	ofile.write("# Format: centerx centery centerz rake normal shear coulomb (kpa)\n");
	for i in range(len(out_object.receiver_object.xstart)):
		# [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(fault_object, i);
		center=conversion_math.get_fault_center(out_object.receiver_object,i);
		ofile.write("%f %f %f %f %f %f %f \n" % (center[0], center[1], center[2], out_object.receiver_object.rake[i], out_object.receiver_normal[i], out_object.receiver_shear[i], out_object.receiver_coulomb[i]) );
	ofile.close();
	print("Outputs written to file.")

	if disp_points != []:
		ofile = open(params.outdir+'ll_disps.txt','w');
		ofile.write("# Format: lon lat u v w (m)\n");
		for i in range(len(out_object.u_ll)):
			ofile.write("%f %f %f %f %f\n" % (disp_points[0][i], disp_points[1][i], out_object.u_ll[i], out_object.v_ll[i], out_object.w_ll[i]) );
		ofile.close();
	return;


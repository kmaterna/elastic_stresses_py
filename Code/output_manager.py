# output_manager

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from subprocess import call
import coulomb_collections
import conversion_math
import io_inp
import io_aftershocks


def produce_outputs(params, inputs, out_object):
	call(['mkdir','-p',params.outdir],shell=False);
	subfaulted_inputs=coulomb_collections.Input_object(PR1=inputs.PR1,FRIC=inputs.FRIC,depth=inputs.depth,
		start_gridx=inputs.start_gridx,start_gridy=inputs.start_gridy,finish_gridx=inputs.finish_gridx,
		finish_gridy=inputs.finish_gridy,xinc=inputs.xinc,yinc=inputs.yinc,minlon=inputs.minlon,maxlon=inputs.maxlon,
		zerolon=inputs.zerolon,minlat=inputs.minlat,maxlat=inputs.maxlat,zerolat=inputs.zerolat,
		source_object=out_object.source_object,receiver_object=out_object.receiver_object); # make a new object of the subfaulted configuration.
	io_inp.write_inp(params.outdir+'subfaulted.inp',subfaulted_inputs);
	surface_def_plot(params,out_object);
	stress_plot(params,out_object,'shear');  
	stress_plot(params,out_object,'normal');
	stress_plot(params,out_object,'coulomb');
	map_plot(params, inputs, out_object);
	write_output_files(params,out_object);
	side_on_plot(params);
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
	plt.savefig(params.outdir+"Displacement_model.eps")
	plt.savefig(params.outdir+"Displacement_model.png")
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
	smallest_stress = -14;  # units: KPa
	largest_stress = 14;  # units: KPa
	color_boundary_object=matplotlib.colors.Normalize(vmin=smallest_stress,vmax=largest_stress, clip=True);
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


def side_on_plot(params):
	[x,y,z,rake,normal,shear,coulomb]=np.loadtxt(params.outdir+'stresses.txt',skiprows=1,unpack=True)
	plt.figure(figsize=(10,6));
	plt.scatter(x,z,c=coulomb,s=1450,marker='s',cmap='jet',edgecolor='black');
	# plt.scatter(x,z,c=normal,s=1450,marker='s');
	# plt.scatter(x,z,c=shear,s=1450,marker='s');
	plt.ylim(z.min()-2,z.max()+2);
	plt.xlim(x.min()-10,x.max()+10);
	plt.xlabel('X axis (km)');
	plt.ylabel('Depth (km)')
	plt.gca().invert_yaxis();
	cb = plt.colorbar();
	if len(set(rake))==1:
		plt.title('Coulomb stress change on fault planes, rake = %.1f (KPa)' % rake[0],fontsize=20);
	else:
		plt.title('Coulomb stress change for variable rake (KPa)',fontsize=20);
	cb.set_label('Kilopascals',fontsize=18);
	plt.savefig(params.outdir+'side_view.eps');
	plt.savefig(params.outdir+'side_view.png');
	plt.close();

	return;


def map_plot(params, inputs, out_object):

	# Make stress bounds for map. 
	stress_bound=14;  # units: KPa
	smallest_stress = -stress_bound;  # units: KPa
	largest_stress = stress_bound;  # units: KPa
	color_boundary_object=matplotlib.colors.Normalize(vmin=smallest_stress,vmax=largest_stress, clip=True);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='RdYlBu_r');

	#Basemap: Make Map
	plt.figure(figsize=(12,10));

	# use low resolution coastlines.
	mymap=Basemap(projection='merc',llcrnrlat=inputs.minlat,llcrnrlon=inputs.minlon, urcrnrlat=inputs.maxlat, urcrnrlon=inputs.maxlon,resolution='i');
	mymap.drawcoastlines(linewidth=1.0,color='black');
	# draw coastlines, country boundaries, fill continents.
	# mymap.fillcontinents(color='white',lake_color='white')
	# draw the edge of the map projection region (the projection limb)

	mymap.drawmapboundary(fill_color='white')
	# draw lat/lon grid lines every degree.
	mymap.drawmeridians(np.arange(inputs.minlon,inputs.maxlon,1),labels=[1,0,0,1],fontsize=20);
	mymap.drawparallels(np.arange(inputs.minlat,inputs.maxlat,1),labels=[1,0,0,1],fontsize=20);

	# Draw each source
	for i in range(len(out_object.source_object.xstart)):
		[x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(out_object.source_object,i);
		lons=[];
		lats=[];
		for j in range(len(x_total)):
			mylon, mylat = conversion_math.xy2lonlat(x_total[j],y_total[j],inputs.zerolon,inputs.zerolat);
			lons.append(mylon);
			lats.append(mylat);
		draw_screen_poly( lats, lons, mymap, [0,0,0] );

	# Draw each receiver
	for i in range(len(out_object.receiver_object.xstart)):
		[x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(out_object.receiver_object,i);
		patch_color=custom_cmap.to_rgba(out_object.receiver_coulomb[i]);  # coloring the map by coulomb stresses. 
		lons=[];
		lats=[];
		for j in range(len(x_total)):
			mylon, mylat = conversion_math.xy2lonlat(x_total[j],y_total[j],inputs.zerolon,inputs.zerolat);
			lons.append(mylon);
			lats.append(mylat);
		draw_screen_poly( lats, lons, mymap, patch_color );

	# Annotate with earthquake location.
	draw_earthquake(params.eqlon, params.eqlat, mymap);

	# Annotate with aftershock locations
	if len(params.aftershocks)>0:
		draw_aftershocks(params.aftershocks,mymap);

	# Map colorbar. 
	custom_cmap.set_array(range(smallest_stress,largest_stress));
	cb = plt.colorbar(custom_cmap);
	cb.set_label('Coulomb Stress (Kilopascals)',fontsize=22);
	cb.ax.tick_params(labelsize=20);

	plt.title('Coulomb Stresses',fontsize=22);
	plt.savefig(params.outdir+'coulomb_basemap.eps');
	plt.close();

	return


def draw_earthquake(eqlon, eqlat, m):
	x,y=m(eqlon, eqlat);
	m.plot(x,y,marker='D',color='m');
	return;


def draw_screen_poly( lats, lons, m, color ):
    x, y = m( lons, lats )
    xy = list(zip(x,y));
    poly = Polygon( xy, facecolor=color, edgecolor='black',alpha=0.7 )
    plt.gca().add_patch(poly)
    return;

def draw_aftershocks(aftershock_file, m):
	[lon, lat, depth, magnitude, time]=io_aftershocks.read_aftershock_table(aftershock_file);
	for i in range(len(lon)):
		x,y=m(lon[i], lat[i]);
		m.plot(x,y,marker='D',color='b',markersize=1);
	return;

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
	ofile.write("Format: centerx centery centerz rake normal shear coulomb (kpa)\n");
	for i in range(len(out_object.receiver_object.xstart)):
		# [x_total, y_total, x_updip, y_updip] = conversion_math.get_fault_four_corners(fault_object, i);
		center=conversion_math.get_fault_center(out_object.receiver_object,i);
		ofile.write("%f %f %f %f %f %f %f \n" % (center[0], center[1], center[2], out_object.receiver_object.rake[i], out_object.receiver_normal[i], out_object.receiver_shear[i], out_object.receiver_coulomb[i]) );
	ofile.close();
	print("Outputs written to file.")

	return;


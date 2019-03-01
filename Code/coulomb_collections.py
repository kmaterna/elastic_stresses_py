# Definitions of objects used in this project

import collections

Params = collections.namedtuple('Params',[
	'input_file','aftershocks',
	'strike_num_receivers',
	'dip_num_receivers',
	'fixed_rake',
	'mu','lame1',
	'alpha',
	'outdir', 'title']);


Input_object = collections.namedtuple('Input_object',[
	'PR1','FRIC','depth',
	'start_gridx','finish_gridx',
	'start_gridy','finish_gridy',
	'xinc','yinc',
	'minlon','maxlon',
	'zerolon',
	'minlat','maxlat',
	'zerolat',
	'eqlon','eqlat',
	'source_object',
	'receiver_object'])


Faults_object = collections.namedtuple('Faults_object',[
	'xstart','xfinish',
	'ystart','yfinish',
	'Kode',
	'rtlat','reverse',
	'strike','dipangle','rake',
	'top','bottom','comment']);


Out_object = collections.namedtuple('Out_object',[
	'x','y',
	'x2d','y2d',
	'u_disp','v_disp','w_disp',
	'source_object','receiver_object',
	'receiver_normal','receiver_shear','receiver_coulomb']);


# Definitions of objects used in this project

import collections

Params = collections.namedtuple('Params', [
    'config_file', 'input_file', 'aftershocks', 'disp_points_file', 'strain_file',
    'strike_num_receivers',
    'dip_num_receivers',
    'fixed_rake',
    'mu', 'lame1',
    'alpha', 'outdir']);

Faults_object = collections.namedtuple('Faults_object', [
    'xstart', 'xfinish',
    'ystart', 'yfinish',
    'Kode', 'zerolon', 'zerolat',
    'rtlat', 'reverse', 'tensile', 'potency',
    'strike', 'dipangle', 'rake',
    'top', 'bottom', 'comment']);

Input_object = collections.namedtuple('Input_object', [
    'PR1', 'FRIC', 'depth',
    'start_gridx', 'finish_gridx',
    'start_gridy', 'finish_gridy',
    'xinc', 'yinc',
    'minlon', 'maxlon',
    'zerolon',
    'minlat', 'maxlat',
    'zerolat',
    'source_object',
    'receiver_object'])

Out_object = collections.namedtuple('Out_object', [
    'x', 'y',
    'x2d', 'y2d',
    'u_disp', 'v_disp', 'w_disp', 'strains',
    'zerolon', 'zerolat',
    'u_ll', 'v_ll', 'w_ll',
    'source_object', 'receiver_object',
    'receiver_normal', 'receiver_shear', 'receiver_coulomb']);

Displacement_points = collections.namedtuple('Disp_Points', [
    'lon', 'lat',
    'dE_obs', 'dN_obs', 'dU_obs',
    'Se_obs', 'Sn_obs', 'Su_obs',
    'name']);
# Currently, disp_points sub-elements are lists, not organized as a large list of individual objects.

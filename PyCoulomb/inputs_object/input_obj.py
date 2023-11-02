
import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


class Input_object:
    # Input object for the calculation of displacements and stresses.
    def __init__(self, PR1, FRIC, depth, start_gridx, finish_gridx, start_gridy, finish_gridy, xinc, yinc,
                 minlon, maxlon, zerolon, minlat, maxlat, zerolat, source_object, receiver_object,
                 receiver_horiz_profile):
        self.PR1 = PR1;
        self.FRIC = FRIC;
        self.depth = depth;
        self.start_gridx = start_gridx;
        self.finish_gridx = finish_gridx;
        self.start_gridy = start_gridy;
        self.finish_gridy = finish_gridy;
        self.xinc = xinc;
        self.yinc = yinc;
        self.minlon = minlon;
        self.maxlon = maxlon;
        self.zerolon = zerolon;
        self.minlat = minlat;
        self.maxlat = maxlat;
        self.zerolat = zerolat;
        self.source_object = source_object;  # list of pycoulomb Faults, with same zerolon/zerolat as overall system
        self.receiver_object = receiver_object;  # list of pycoulomb Faults, with same zerolon/zerolat as overall system
        self.receiver_horiz_profile = receiver_horiz_profile;

    def define_map_region(self):
        """
        Define bounding box for map [W, E, S, N] based on sources and receivers, if bigger than coord system
        """
        region = [self.minlon, self.maxlon, self.minlat, self.maxlat];
        allfaults = self.receiver_object + self.source_object
        for item in allfaults:
            lon, lat = fault_vector_functions.xy2lonlat(item.xstart, item.ystart, self.zerolon, self.zerolat);
            if lon < region[0]:
                region[0] = lon;
            if lon > region[1]:
                region[1] = lon;
            if lat < region[2]:
                region[2] = lat;
            if lat > region[3]:
                region[3] = lat;
        return region;

    def check_each_fault_has_same_coord_system(self, tol=0.00001):
        """
        Check that all faults have the same coord system
        """
        for item in self.source_object + self.receiver_object:
            if np.abs(item.zerolon-self.zerolon) > tol:
                raise(ValueError, "input or receiver faults lack a good longitude coordinate system.");
            if np.abs(item.zerolat-self.zerolat) > tol:
                raise(ValueError, "input or receiver faults lack a good latitude coordinate system.");
        return;

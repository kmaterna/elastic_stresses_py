import numpy as np
from Tectonic_Utils.geodesy import fault_vector_functions


class Input_object:
    # Input object for the calculation of displacements and stresses.
    # source_object is required.
    def __init__(self, xinc, yinc, minlon, maxlon, zerolon, minlat, maxlat, zerolat, source_object,
                 PR1=0.25, FRIC=0.4, depth=0, start_gridx=-20, finish_gridx=20, start_gridy=-20, finish_gridy=20,
                 receiver_object=(), receiver_horiz_profile=None):
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
        if len(self.source_object) == 0:
            raise ValueError("Error! Valid Input_objects must have nonzero source_object.");

    def define_map_region(self):
        """
        Define bounding box for map [W, E, S, N] based on sources and receivers, if bigger than coord system
        """
        region = [self.minlon, self.maxlon, self.minlat, self.maxlat];
        allfaults = self.source_object;
        if len(self.receiver_object) > 0:
            allfaults = self.receiver_object + self.source_object
        if len(allfaults) == 0:
            raise ValueError("Error! No faults given, so automatic region cannot be determined.")
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
            if np.abs(item.zerolon - self.zerolon) > tol:
                raise (ValueError, "input or receiver faults lack a good longitude coordinate system.");
            if np.abs(item.zerolat - self.zerolat) > tol:
                raise (ValueError, "input or receiver faults lack a good latitude coordinate system.");
        return;

    def modify_inputs_object(self, PR1=None, FRIC=None, depth=None, start_gridx=None, finish_gridx=None,
                             start_gridy=None, finish_gridy=None, xinc=None, yinc=None, minlon=None, maxlon=None,
                             zerolon=None, minlat=None, maxlat=None, zerolat=None, source_object=None,
                             receiver_object=None, receiver_horiz_profile=None):
        """
        Modify the fields in a PyCoulomb.Input_object.
        """
        PR1 = self.PR1 if PR1 is None else PR1;
        FRIC = self.FRIC if FRIC is None else FRIC;
        depth = self.depth if depth is None else depth;
        start_gridx = self.start_gridx if start_gridx is None else start_gridx;
        finish_gridx = self.finish_gridx if finish_gridx is None else finish_gridx;
        start_gridy = self.start_gridy if start_gridy is None else start_gridy;
        finish_gridy = self.finish_gridy if finish_gridy is None else finish_gridy;
        xinc = self.xinc if xinc is None else xinc;
        yinc = self.yinc if yinc is None else yinc;
        minlon = self.minlon if minlon is None else minlon;
        maxlon = self.maxlon if maxlon is None else maxlon;
        zerolon = self.zerolon if zerolon is None else zerolon;
        minlat = self.minlat if minlat is None else minlat;
        maxlat = self.maxlat if maxlat is None else maxlat;
        zerolat = self.zerolat if zerolat is None else zerolat;
        source_object = self.source_object if source_object is None else source_object;
        receiver_object = self.receiver_object if receiver_object is None else receiver_object;
        rec_profile = self.receiver_horiz_profile if receiver_horiz_profile is None else receiver_horiz_profile;

        modified_inputs = Input_object(PR1=PR1, FRIC=FRIC, depth=depth, start_gridx=start_gridx,
                                       finish_gridx=finish_gridx, start_gridy=start_gridy, finish_gridy=finish_gridy,
                                       xinc=xinc, yinc=yinc, minlon=minlon, maxlon=maxlon, zerolon=zerolon,
                                       minlat=minlat, maxlat=maxlat, zerolat=zerolat,
                                       source_object=source_object, receiver_object=receiver_object,
                                       receiver_horiz_profile=rec_profile);
        return modified_inputs;


def configure_default_displacement_input(source_object, zerolon, zerolat, bbox, domainsize=20,
                                         num_points_x=20, num_points_y=20):
    """
    Build a default Input object for displacement-only calculations using a simpler-parameterized interface.

    :param source_object: a list of PyCoulomb faults
    :param zerolon: float
    :param zerolat: float
    :param bbox: [w, e, s, n]
    :param domainsize: in km. default is 20.
    :param num_points_x: default 20
    :param num_points_y: default 20
    """
    inputs = Input_object(start_gridx=-domainsize, finish_gridx=domainsize,
                          start_gridy=-domainsize, finish_gridy=domainsize,
                          xinc=domainsize / num_points_x, yinc=domainsize / num_points_y,
                          minlon=bbox[0], maxlon=bbox[1], zerolon=zerolon,
                          minlat=bbox[2], maxlat=bbox[3], zerolat=zerolat,
                          source_object=source_object);
    return inputs;

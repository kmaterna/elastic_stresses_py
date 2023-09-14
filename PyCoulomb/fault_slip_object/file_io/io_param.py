
from .. import fault_slip_object


def read_param_file(infile):
    """
    Read faults in the param file format
    """
    print("Reading Param distribution %s " % infile);
    fault_list = [];
    ifile = open(infile, 'r');
    start_read = 0;
    current_segment = 1;
    current_length, current_width = 0, 0
    for line in ifile:
        if start_read and len(line.split()) > 8 and line.split()[0][0] != '#':
            temp = line.split();
            lat = float(temp[0])
            lon = float(temp[1])
            depth = float(temp[2])
            slip = float(temp[3])/100;
            rake = float(temp[4]);
            strike = float(temp[5]);
            dip = float(temp[6]);
            one_fault = fault_slip_object.FaultSlipObject(strike=strike, dip=dip, length=current_length,
                                                          width=current_width,
                                                          depth=depth, rake=rake, slip=slip, tensile=0,
                                                          lon=lon, lat=lat, segment=current_segment);
            fault_list.append(one_fault);
        if '#Fault_segment = ' in line:
            current_segment = int(line.split()[2]);
            current_length = float(line.split()[4].split('km')[0]);
            current_width = float(line.split()[9].split('km')[0]);
            start_read = 0;
        if ' #Lat. Lon. depth slip rake strike' in line:
            start_read = 1;
    ifile.close();
    print("--> Returning %d fault patches" % len(fault_list));
    return fault_list;

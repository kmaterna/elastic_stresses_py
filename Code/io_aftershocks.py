# The purpose of this script is to read aftershock tables and return plotting-friendly arrays. 



def read_aftershock_table(infile):
	lon=[]; lat=[]; time=[]; depth=[]; magnitude=[];

	ifile=open(infile);
	for line in ifile:
		temp=line.split();
		if temp[0][0]=='#':
			continue;
		else:
			time.append(temp[0]);
			lon.append(float(temp[3]));
			lat.append(float(temp[2]));
			depth.append(float(temp[4]));
			magnitude.append(float(temp[5]));

	return [lon, lat, depth, magnitude, time];
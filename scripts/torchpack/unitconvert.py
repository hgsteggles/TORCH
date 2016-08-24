

def ra_to_deg(hrs, mins, secs):
	return (hrs * 15.0) + (mins * 15.0 / 60.0) + (secs * 15.0 / 3600.0)

def deg_to_ra(deg):
	hrs = int(deg / 15.0)
	mins = int((deg - 15.0 * hrs) / (15.0 / 60.0))
	secs = int((deg - 15.0 * hrs - (15.0 / 60.0) * mins) / (15.0 / 3600.0) + 0.5)
	return hrs, mins, secs

def dec_to_deg(deg, arcm, arcs):
	return deg + (arcm / 60.0) + (arcs / 3600.0)

def deg_to_dec(deg):
	hrs = int(deg)
	mins = int((deg - hrs) * 60.0)
	secs = int((deg - hrs - (mins / 60.0)) * 3600.0)
	return hrs, mins, secs
import matplotlib.colors as mcolors
import sys, getopt
import os

def get_marc_map():
	m0 =   0.0/255.0,  0.0/255.0,   0.0/255.0
	m1 =  90.0/255.0,  0.0/255.0, 180.0/255.0
	m2 = 127.0/255.0,  3.0/255.0, 255.0/255.0
	m3 = 156.0/255.0, 13.0/255.0, 180.0/255.0
	m4 = 179.0/255.0, 30.0/255.0,  12.0/255.0
	m5 = 202.0/255.0, 61.0/255.0,   0.0/255.0
	m6 = 221.0/255.0,108.0/255.0,   0.0/255.0
	m7 = 239.0/255.0,170.0/255.0,   0.0/255.0
	m8 = 255.0/255.0,250.0/255.0,   0.0/255.0
	return make_colormap([m0, m1, 0.125, m1, m2, 0.25, m2, m3, 0.375, m3, m4, 0.5, m4, m5, 0.625, m5, m6, 0.75, m6, m7, 0.875, m7, m8])

def get_green_map():
	green0 =   0.0/255.0,   0.0/255.0,   0.0/255.0
	green1 =   0.0/255.0,  28.0/255.0,  64.0/255.0
	green2 =   0.0/255.0, 128.0/255.0, 125.0/255.0
	green3 =   6.0/255.0, 191.0/255.0,   0.0/255.0
	green4 = 255.0/255.0, 221.0/255.0,   0.0/255.0
	return make_colormap([green0, green1, 0.25, green1, green2, 0.5, green2, green3, 0.75, green3, green4])

def get_bgyr_map():
	bgyr0 =   0.0/255.0,   0.0/255.0,   0.0/255.0
	bgyr1 =   0.0/255.0,   0.0/255.0,  51.0/255.0
	bgyr2 =   0.0/255.0, 102.0/255.0, 102.0/255.0
	bgyr3 =   0.0/255.0, 152.0/255.0,   0.0/255.0
	bgyr4 = 206.0/255.0, 206.0/255.0,   0.0/255.0
	bgyr5 = 255.0/255.0,   0.0/255.0,   0.0/255.0
	return make_colormap([bgyr0, bgyr1, 0.2, bgyr1, bgyr2, 0.4, bgyr2, bgyr3, 0.6, bgyr3, bgyr4, 0.8, bgyr4, bgyr5])

def get_water_map():
	fb0 =   0/255.0,  0/255.0,  0/255.0
	fb1 =  21/255.0,  4/255.0,168/255.0
	fb2 =   0/255.0,247/255.0,214/255.0
	return make_psychodelic_colormap([fb0, fb1, 0.5, fb1, fb2])

def get_grey_map():
	black = 0, 0, 0
	white = 1, 1, 1
	return make_colormap([white, black])

def make_colormap(seq):
	"""Return a LinearSegmentedColormap
	seq: a sequence of floats and RGB-tuples. The floats should be increasing
	and in the interval (0,1).
	"""
	seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
	cdict = {'red': [], 'green': [], 'blue': []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			r1, g1, b1 = seq[i - 1]
			r2, g2, b2 = seq[i + 1]
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
	return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def make_psychodelic_colormap(seq):
	"""Return a psychodelic LinearSegmentedColormap
	seq: a sequence of floats and RGB-tuples. The floats should be increasing
	and in the interval (0,1).
	"""
	seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
	cdict = {'red': [], 'green': [], 'blue': []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			r1, g1, b1 = tuple(255*x if x is not None else x for x in seq[i - 1])
			r2, g2, b2 = tuple(255*x if x is not None else x for x in seq[i + 1])
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
	return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def parse_plot_args(argv, figformat):
	try:
		opts, args = getopt.getopt(argv[1:],"i:phdt")
	except getopt.GetoptError:
		print 'Error:', argv[0], '-i <inputfile>'
		sys.exit(2)
	inputfile = ''
	outputfile = ''
	var_type = ''
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		if opt in ("-h"):
			var_type = "hii"
		if opt in ("-p"):
			var_type = "pre"
		if opt in ("-d"):
			var_type = "den"
		if opt in ("-t"):
			var_type = "tem"
	if inputfile == '':
		print argv[0], '-i <inputfile>'
		sys.exit(2)
	outputfile = os.path.splitext(inputfile)[0] + '.' + figformat
	return inputfile, outputfile, var_type



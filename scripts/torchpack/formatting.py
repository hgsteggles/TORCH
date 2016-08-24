import math

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("convert", "./unitconvert.py")
import convert

def latexify(str):
	return r'${}$'.format(str)

def fmt(x, pos):
	return latexify(fmt_nolatex(x, pos))

def fmt_nolatex(x, pos):
	if x >= -100 and x <= 100:
		a = '{:.1f}'.format(x)
		return '{}'.format(a)
	else:
		a, b = '{:.1e}'.format(x).split('e')
		b = int(b)
		return '{} \\times 10^{{{}}}'.format(a, b)

def fmt_power(x, fmt, power):
	a = x / math.pow(10.0, int(power))
	a = fmt.format(a)
	if int(power) == 0:
		return a
	else:
		return '{} \\times 10^{{{}}}'.format(a, int(power))

def fmt_mass(x):
	a = '{:.0f}'.format(x)
	return '{}'.format(a)

def fmt_ra(x, pos):
	hrs, mins, secs = convert.deg_to_ra(x)
	if hrs < 10:
		hrs = "0" + str(hrs)
	return '{0:02d}h{1:02d}m{2:02d}s'.format(hrs, abs(mins), abs(secs))

def fmt_dec(x, pos):
	hrs, mins, secs = convert.deg_to_dec(x)
	return r"{0:02d}$^\circ${1:02d}'{2:02d}''".format(hrs, abs(mins), abs(secs))
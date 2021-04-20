import concerto as conc
import numpy as np
from scipy import signal

def m(ds, sr): 
	""" 
	The function to be used as a macro (it can be renamed)
	
	:param ds: the first argument is always dataset as dictionary (it can have additional arguments)
	:return: transformed dataset as dictionary
	"""
	
	y=ds
	z,p,k = signal.ellip(8,0.1,60,1,'highpass',analog=False,output='zpk',fs=sr)
	sos = signal.zpk2sos(z,p,k)
	filtered = signal.sosfiltfilt(sos,y)
	
	return filtered

"""Strip a folder of all files that aren't DICOMs, and that don't look like VMAT images."""
from machinelearning.tools import strip

path = r'C:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files\VMATs\TrueBeam 3'
strip(path, 'vmat', [1, 2, 3], incorrect_names=['pf', 'vmat', 't1'], correct_names=['t2', 't3', 'open', 'gantry', 'speed'])

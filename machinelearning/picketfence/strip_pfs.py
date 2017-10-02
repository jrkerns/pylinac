"""Strip a folder of all files that aren't DICOMs, and that don't look like VMAT images."""
from machinelearning.tools import strip

good_names = ('pf', 'vmat', '90', '270', '180', 'ra', 't1', 'picket')
bad_names = ('open', 't2', 't3', 'speed', 'dose', 'rate', 'drgs', 'mlcs', 'jaw', 'coll', 'strip')
path = r'C:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files\Picket Fences'
strip(path, 'picketfence', correct_prediction=[1], incorrect_names=bad_names, correct_names=good_names)

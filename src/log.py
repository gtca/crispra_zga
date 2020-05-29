# 
# Collection of logging utility functions
# 

import yaml

def write_log(statlogfile: str, stats):
	with open(statlogfile, 'a') as slf:
		yaml.dump(stats, slf, default_flow_style=False)
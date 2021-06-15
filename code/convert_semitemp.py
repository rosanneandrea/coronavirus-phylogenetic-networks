from funcs_format_conversions import semitemp_extract_tcs
from analyse_results import list_files
import subprocess

dir1 = "../results/temporal/rerun/"

# EXTRACT TCS FROM SEMITEMPORAL OUTPUT FILES
"""
input_list = list_files(dir1, ".log")
for infilename in input_list:
	semitemp_extract_tcs(infilename)

"""
# TURN TCS INTO NETWORK
input_list = list_files(dir1, ".seq")

for infilename in input_list:
	
	cmd1 = f"python3 TCS_to_Nw.py {dir1}{infilename} -o {dir1}{infilename.replace('.seq','.net')}"
	subprocess.run(cmd1, shell=True)

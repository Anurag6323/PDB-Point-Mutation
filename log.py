import sys
import subprocess


def logger():
    output = subprocess.check_output([sys.executable, './point_mutator.py'])
    with open('logfile.txt', 'wb') as outfile:
        outfile.write(output)


logger()

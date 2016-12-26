#!/usr/bin/env python
""" a python script to create a script init template
"""

import sys

def main():
    pathout = sys.argv[1]

    with open(pathout, "w") as outf:
        outf.write("""#!/usr/bin/env python
\"\"\" 
\"\"\"

import os, sys, argparse  

def get_cmd():
    \"\"\" read input line arguments
    \"\"\"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfile", help="the input file")
    parser.add_argument("-o", action="store", dest="outputfile", help="the output file")
    params = parser.parse_args()
    return params

def main():
    params = get_cmd()

    sys.exit(0)  

if __name__ == "__main__":
    main()
""")
    sys.exit(0)

if __name__ == "__main__":
    main()

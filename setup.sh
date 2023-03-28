#!/bin/bash
# creatation: 1-25-2022
# Author: Ruitu Lyu (lvruitu@gmail.com)

# Stop on error
set -e
###
### setup.sh - This script is used to setup the 'KAS-Analyzer' pipeline.
###
### Usage: ./setup.sh	
###
### -h or --help Print the help.
###

# Help message for shell scripts

help() { 
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    help
fi

if [[ "$#" -eq 1 ]]; then
    exit 1
elif [[ "$#" -gt 1 ]]; then
    echo "please refer the usage: ./setup.sh"
    exit 1
elif [[ "$#" -eq 0 ]]; then 
    echo "setup KAS-Analyzer!"
fi    

# make these scripts executable and add their directory to your PATH
# get the path of shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

chmod u+x ${SH_SCRIPT_DIR}/KAS-Analyzer
chmod 755 ${SH_SCRIPT_DIR}/scripts/*sh
echo export PATH=\"${SH_SCRIPT_DIR}:"$"PATH\" >> $HOME/.bashrc
. $HOME/.bashrc

echo "All the shell scripts have been made to be executable and added to the \$path variable, please enjoy the KAS-Analyzer!"

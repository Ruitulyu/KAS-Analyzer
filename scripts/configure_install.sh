#!/bin/bash
# 'KAS-pipe2 install' was developed by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-pipe2 install [ -h\--help ] [ -conda ] [ -check ] [ -t tools ] [ -KAS-pipe2 ]"
exampleHelp="Example: KAS-pipe2 install or KAS-pipe2 install -check"
condaHelp="-conda: check the installation or install anaconda in your computer."
checkHelp="-check: list the installed and uninstalled tools."
KASpipe2Help="-KAS-pipe2: install and configure the KAS-pipe2 conda environment."
toolsHelp="-t [tools]: check the installation of specific tool, if not, will install automaticially."
helpHelp="-h\-help: print this help and exit.
Note: this subcommand is used to install conda environment and specific tool that needed in KAS-pipe2."


printHelpAndExit() {
    echo -e ""	
    echo -e "$usageHelp"
    echo -e "" 
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$condaHelp"
    echo -e ""
    echo -e "$checkHelp"
    echo -e ""
    echo -e "$toolsHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 install' will print the help.
if [[ $# == 0 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# test the installation of conda.
installanaconda() {
if ! type conda >/dev/null 2>&1 ;then
   echo ""
   echo "conda was not installed, anaconda3 will be installed..."
   echo ""
   wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
   chmod 755 ./Anaconda3-2021.11-Linux-x86_64.sh
   bash ./Anaconda3-2021.11-Linux-x86_64.sh
else
   echo "conda have been installed."
fi     	
}

# function to check the installed and uninstalled tools.
checkinstall() {
all_tools="bedtools bwa bowtie2 macs2 macs14 deeptools homer fastqc cutadapt trim_galore bedGraphToBigWig samtools picard R"
cat /dev/null > .installed_tools.txt
cat /dev/null > .uninstalled_tools.txt
for ((i=1; i<=14; i++))	
do
echo 	
tool_selected=$( echo $all_tools | awk -v x=$i '{print $x}' )	
if ! type $tool_selected > /dev/null 2>&1 ;then	
   echo "$tool_selected" >> .uninstalled_tools.txt
else 
   echo "$tool_selected" >> .installed_tools.txt
fi
done

echo "Installed_tools:"
cat .installed_tools.txt
echo ""
echo "Uninstalled_tools:"
cat .uninstalled_tools.txt

rm -f .installed_tools.txt
rm -f .uninstalled_tools.txt
}

# function to install tools.
installtools() {
if [[ $tools == "bedGraphToBigWig" ]]; then
   tools=ucsc-bedgraphtobigwig
elif [[ $tools == "trim_galore" ]]; then
   tools=trim-galore
fi
	
if ! type $tools >/dev/null 2>&1 ;then
   echo ""
   echo "$tools was not installed, install $tools with 'conda'"
   echo ""
   conda install -c bioconda $tools
   echo "done."
   echo ""
else
   echo "$tools have been installed."
fi   
}	

installKAS-pipe2() {
CONDA_ENV_PY3=KAS-pipe2
REQ_TXT_PY3=${SH_SCRIPT_DIR}/../requirements.txt

# conda --version  # check if conda exist.
echo "=== Installing KAS-pipe2 pipeline's conda environments ==="
conda create -n ${CONDA_ENV_PY3} --file ${REQ_TXT_PY3} -y -c defaults -c bioconda -c conda-forge -c r

echo "=== Configuring for KAS-pipe2 pipeline's conda environments ==="
CONDA_PREFIX_PY3=$(conda env list | grep -P "\b${CONDA_ENV_PY3}\s" | awk '{if (NF==3) print $3; else print $2}')

if [ ! "${CONDA_PREFIX_PY3}" ] ;then
   echo "Error: 'KAS-pipe2' conda environments not found."
   echo "Try to reinstall pipeline's Conda environments."
   echo
   echo "1) $ bash configure_uninstall.sh"
   echo "2) $ KAS-pipe2 install -KAS-pipe2"
   exit 1
fi

# make activate.d to init pipeline's Conda envs
CONDA_LIB="${CONDA_PREFIX_PY3}/lib"
CONDA_BIN="${CONDA_PREFIX_PY3}/bin"
CONDA_ACTIVATE_D="${CONDA_PREFIX_PY3}/etc/conda/activate.d"
CONDA_DEACTIVATE_D="${CONDA_PREFIX_PY3}/etc/conda/deactivate.d"
CONDA_ACTIVATE_SH="${CONDA_ACTIVATE_D}/env_vars.sh"
CONDA_DEACTIVATE_SH="${CONDA_DEACTIVATE_D}/env_vars.sh"
mkdir -p ${CONDA_ACTIVATE_D}
mkdir -p ${CONDA_DEACTIVATE_D}
touch ${CONDA_ACTIVATE_SH}
touch ${CONDA_DEACTIVATE_SH}

# disable multithreading for BLAS
echo "export OPENBLAS_NUM_THREADS=1" > ${CONDA_ACTIVATE_SH}
echo "export MKL_NUM_THREADS=1" >> ${CONDA_ACTIVATE_SH}
echo "unset OPENBLAS_NUM_THREADS MKL_NUM_THREADS" >> ${CONDA_DEACTIVATE_SH}

# to prevent conflict between Conda's python packages and user's local one
echo "export PYTHONNOUSERSITE=True" >> ${CONDA_ACTIVATE_SH}
echo "unset PYTHONNOUSERSITE" >> ${CONDA_DEACTIVATE_SH}

# LD_LIBRARY_PATH due to libgcc problem
echo "export OLD_LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}" >> ${CONDA_ACTIVATE_SH}
echo "export LD_LIBRARY_PATH=${CONDA_LIB}:\${LD_LIBRARY_PATH}" >> ${CONDA_ACTIVATE_SH}

echo "export LD_LIBRARY_PATH=\${OLD_LD_LIBRARY_PATH}" >> ${CONDA_DEACTIVATE_SH}
echo "unset OLD_LD_LIBRARY_PATH" >> ${CONDA_DEACTIVATE_SH}

# to prevent conflict between Conda's R and global(local) R
echo "export OLD_R_HOME=\${R_HOME}" >> ${CONDA_ACTIVATE_SH}
echo "export OLD_R_LIBS=\${R_LIBS}" >> ${CONDA_ACTIVATE_SH}
echo "export R_HOME=${CONDA_LIB}/R" >> ${CONDA_ACTIVATE_SH}
echo "export R_LIBS=${CONDA_LIB}/R/library" >> ${CONDA_ACTIVATE_SH}

echo "export R_HOME=\${OLD_R_HOME}" >> ${CONDA_DEACTIVATE_SH}
echo "export R_LIBS=\${OLD_R_LIBS}" >> ${CONDA_DEACTIVATE_SH}
echo "unset OLD_R_HOME" >> ${CONDA_DEACTIVATE_SH}
echo "unset OLD_R_LIBS" >> ${CONDA_DEACTIVATE_SH}

echo "Configure KAS-pipe2 conda environment successfully!"
}


# get the value of options.
ARGS=`getopt -a -o hackt: --long help,conda,check,KAS-pipe2,tools: -- "$@"`

if [ $? != 0 ];then
        echo "Terminating..."
        exit 1
fi

eval set -- "${ARGS}"

while :
do
    case $1 in
        -h | --help | -help)
            printHelpAndExit
            shift
            ;;
        -a | --conda | -conda)
            installanaconda
            shift
            ;;
        -c | --check | -check)
            checkinstall
            shift
            ;;
        -k | --KAS-pipe2 | -KAS-pipe2)
            installKAS-pipe2
            shift
            ;;	    
        -t | --tools | -tools)
	    tools=$2
            installtools
            shift
            ;;
        --)
            shift
            break
            ;;
    esac
shift
done

if [[ $# == 1 ]] ;then
    printHelpAndExit
fi

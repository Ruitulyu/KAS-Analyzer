#!/bin/bash
# 'KAS-Analyzer uninstall' was developed by Ruitu Lyu on 1-24-2022.

conda deactivate

CONDA_ENV_PY3=KAS-Analyzer
CONDA_ENV_OLD_PY3=KAS-Analyzer-python3

conda env remove -n ${CONDA_ENV_PY3} -y
conda env remove -n ${CONDA_ENV_OLD_PY3} -y

echo "'KAS-Analyzer' conda environment was successfully removed!"

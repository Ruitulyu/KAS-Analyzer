#!/bin/bash
# 'KAS-pipe2 uninstall' was developed by Ruitu Lyu on 1-24-2022.

conda deactivate

CONDA_ENV_PY3=KAS-pipe2
CONDA_ENV_OLD_PY3=KAS-pipe2-python3

conda env remove -n ${CONDA_ENV_PY3} -y
conda env remove -n ${CONDA_ENV_OLD_PY3} -y

echo "'KAS-pipe2' conda environment was successfully removed!"

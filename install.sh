#!/bin/bash

pip install joblib
pip install varname
pip install jupyter_contrib_nbextensions
pip install rise
pip install jupyter_nbextensions_configurator
jupyter nbextensions_configurator enable --user
pip install hide-code
jupyter nbextension install --py hide_code
jupyter nbextension enable --py hide_code
jupyter serverextension enable --py hide_code
pip install xlrd
pip install sklearn

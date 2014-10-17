# -*- coding: utf-8 -*-
"""
Created on Mon Oct 07 09:34:29 2013

@author: James

This script converts the Qt Designer file "PyQA_ui.ui" file to a PySide programmatic file "PyQA_ui.py".
It is a specific script for use with PyQA such that a simple run will
recompile the *.py file when changes are made to the *.ui file in Qt Designer.
"""
import os
import os.path as osp

import pysideuic as uic


uifile = osp.join(os.getcwd(), 'PyQA.ui')
pyfile = osp.join(os.getcwd(), 'PyQAui.py')

# compile .py file from .ui file
with open(pyfile, 'w') as file:
    uic.compileUi(uifile, file)

# fixer for matplotlibwidget class location
with open(pyfile) as file:
    pyfile_text = file.read()
new_text = pyfile_text.replace("from matplotlibwidget", "from pyqa_gui.matplotlibwidget")
with open(pyfile, 'w') as file:
    file.write(new_text)





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 09:54:00 2021

@author: ben
"""

import svgutils.transform as sg
import sys 

import os 



os.chdir('/home/ben/OneDrive/Padova/Redaction/Articles/1b_InversionUsingGravityPotMethod/figs/')
#create new SVG figure
fig = sg.SVGFigure("16cm", "6.5cm")

# load matpotlib-generated figures
fig1 = sg.fromfile('fig2a.svg')
fig2 = sg.fromfile('fig2b.svg')
fig3 = sg.fromfile('fig2c.svg')

# get the plot objects
plot1 = fig1.getroot()
plot2 = fig2.getroot()
plot3 = fig3.getroot()
plot2.moveto(480, 0)
plot3.moveto(880, 0)

# add text labels
txt1 = sg.TextElement(25,20, "A", size=12, weight="bold")
txt2 = sg.TextElement(305,20, "B", size=12, weight="bold")
txt3 = sg.TextElement(705,20, "C", size=12, weight="bold")

# append plots and labels to figure
fig.append([plot1, plot2, plot3])
fig.append([txt1, txt2, txt3])

# save generated SVG files
fig.save("fig_final.svg")
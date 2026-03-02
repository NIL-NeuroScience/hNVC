# hNVC

This repository contains code to analyze simultaneously
acquired human fMRI, EEG, and pupil diameter. Types of 
analysis include modeling of fMRI BOLD with EEG and 
pupil diameter time-series using convolutional models.

## Contents

Main.m - main script which loads data and creates 
convolutional model to predict fMRI BOLD from EEG
and pupil time-series.

analysis - contains all plotted analysis.

+utils - utility functions

datasets - contains code to load fMRI, EEG, and pupil
data.
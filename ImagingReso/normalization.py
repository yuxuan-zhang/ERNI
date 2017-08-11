import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
import NeuNorm as neunorm
from NeuNorm.roi import ROI
from NeuNorm.normalization import Normalization


data_path = '/Users/y9z/Documents/Neutron_data/IPTS_13639/15_Resonance_OBs_ForDirs11-14/sum_all'
ob_path = '/Users/y9z/Documents/Neutron_data/IPTS_13639/15_Resonance_OBs_ForDirs11-14/sum_ob'
o_norm = Normalization()
o_norm.load(folder=ob_path, data_type='ob')
o_norm.load(folder=data_path)
o_norm.normalization()
o_norm.data['normalized']

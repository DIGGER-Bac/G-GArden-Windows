import argparse
import glob
import os
import seaborn as sns
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision as torchvision
from PIL import Image, ImageFilter
import matplotlib.pyplot as plt


class OhOptimizer(nn.module):
    def __init__(self):
        super(OhOptimizer, self).__init__()
        self.fc1 = nn.Linear()
        self.fc2 = nn.Linear()
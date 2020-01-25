#!/usr/bin/env python3
#
#   Program:
#       Visualize the combination of 2D-PES, 3D-PES and molecular evolution of one dynamic trajectory .
#
#   Input:
#       $1 = PES
#       $2 = Important point on this surface
#       $3 = Coordinate of trajectory
#       $4 = xyz structure file of trajectory
#
#   Output:
#
#
#   Keywords: 
#   plotly, animation, xyz2graph, MDanalysis
# 
# Reference of tutorial blog or youtube
#   2. basic example of creating animation
#   https://learndataanalysis.org/a-basic-example-how-to-create-animation-matplotlib-tutorial/
#   3. combine two 2D animations in one fig
#   https://pythonmatplotlibtips.blogspot.com/2018/01/combine-two-2d-animations-in-one-figure-matplotlib-artistanimation.html
#   4. combine 3D and two 2D animations in one figure
#   https://pythonmatplotlibtips.blogspot.com/2018/01/combine-3d-two-2d-animations-in-one-figure-artistdanimation.html
#   5. convert xyz file into molecular structure
#   https://github.com/zotko/xyz2graph
#   6. display static image of plotly lib; download orca
#   https://plot.ly/python/static-image-export/
#   https://plot.ly/python/renderers/
#
#   History:
#       2020/01/19, Grace
#       20

import sys
import os
import math
import numpy as np
import pandas as pd
from copy import deepcopy
import plotly.graph_objects as go 
import plotly.offline as py
from xyz2graph import MolGraph, to_plotly_figure
# from ipywidgets import interactive, HBox, VBox
# from IPython.display import display

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
# import matplotlib.tri as tri
# from mpl_toolkits.mplot3d import Axes3D
# # 
# import plotly
# from plotly.subplots import make_subplots
# import plotly.offline as py
# import plotly.io as pio
# from plotly.offline import offline, iplot

import matplotlib.animation as animation

# workDir = os.getcwd()
# PATH_fig = workDir

# title = ''
# xaxis = 'Modified IRC of TSS1'
# yaxis = 'IRC of TSS2'
# zaxis = 'kcal/mol'
# axis_fontsize = 6

# static figure objects
# fig = plt.figure()
# threeD = fig.add_subplot(1, 2, 2,projection='3d') # top and bottom right 
# mol = fig.add_subplot(2,2,1)
# mol.axis('off') # remove borders, ticks and labels
# twoD = fig.add_subplot(2, 2, 3)  # bottom left

# interact objects
# fig = make_subplots(rows = 2, cols = 2)
# threeD = fig.add_trace(rows = 1, cols = 2)
# mol = fig.add_trace(rows = 2, cols = 1)
# twoD = fig.add_trace(rows = 2, cols = 2)

def main():
    # 1. Import rawdata
    X, Y, E, Pts = getPES( sys.argv[1], sys.argv[2] ) # E = relative energy in kcal/mol
    traj = getTraj(sys.argv[2], sys.argv[3]) # format: x,y,E
    
    # 2. Plot the background: 2D-PES and 3D-PES
    twoD = twoDwPts(X, Y, E, Pts)
    threeD = threeDwPts(X, Y, E, Pts)
    mol = molGeo(sys.argv[4])
    
    # 3. Combine above  FIXME:
    # fig = ani(twoD,threeD,traj)
    # fig = VBox( [ twoD,threeD ] )
    # display(fig)
    # py.iplot(fig)
    
    
    # display(fig)
    # display(fig)
    # outines and export them into one video
    frames = ani(twoD,threeD,traj)
    # print('test')
    # fig.add_trace(threeD,row=1,col=2)
    # fig.add_trace(twoD,row=2,col=1)
    # fig.add_trace(mol,2,2)
    # py.plot(fig)

def getPES(file1,file2):
    PES = pd.read_csv(file1,' ',header=None)
    dim = int(math.sqrt(PES.shape[0]))
    X = np.reshape(PES.values[:,0], (dim, dim))
    Y = np.reshape(PES.values[:,1], (dim, dim))

    # import name and coord of pts
    Pts = pd.read_csv(file2,' ',header=None)

    # relative energy in kcal/mol
    refE = Pts.values[0][3]
    relE = np.arange(dim*dim)
    for i in range(dim*dim):
        relE[i] = ( PES.values[i][2] - refE ) * 627.5095
    E = np.reshape(relE, (dim, dim))

    tmpE = deepcopy(Pts[3])
    for i in range(Pts.shape[0]):
        Pts[3][i] = ( tmpE[i] - refE ) * 627.5095

    return X,Y,E,Pts

def getTraj(file1, file2):
    tmp = pd.read_csv(file1,' ',header=None)
    traj = pd.read_csv(file2,' ',header=None)
    refE = tmp.values[0][3]

    # relative energy in kcal/mol
    tmpE = deepcopy(traj[2])
    for i in range(traj.shape[0]):
        traj[2][i] = ( tmpE[i] - refE ) * 627.5095

    return traj

def twoDwPts( X, Y, E, Pts):
    dim = X.shape[0]
    X = np.unique(sorted( np.reshape( X, (dim*dim) ) ) )
    Y = np.unique(sorted( np.reshape( Y, (dim*dim) ) ) )
    
    contours = go.Contour(
        x=X, y=Y, z=E,
        line_smoothing=0.85,
        contours_coloring='none',
        line_width=2
    ) 

    pts = go.Scatter(
        x=Pts[1], y=Pts[2], 
        text=Pts[0],
        textposition="bottom center",
        mode='markers+text',
        marker_line_width=2,
        marker_size=10,
        marker_color='rgba(37,116,169,1)' #jelly blue
    )

    layout = go.Layout(
        showlegend=False
    )

    twoD = go.FigureWidget(data=[contours,pts],layout=layout)
    return twoD
    
def threeDwPts(X, Y, E, Pts):
    pes = go.Surface(
        x=X, y=Y, z=E,
        opacity=0.7
    ) 

    pts = go.Scatter3d(
        x=Pts[1], y=Pts[2], z=Pts[3],
        text=Pts[0],
        textposition="bottom center",
        mode='markers+text',
        marker_line_width=2,
        marker_size=10,
        marker_color='rgba(37,116,169,1)' #jelly blue
    )

    layout = go.Layout(
        showlegend = False,
        scene = dict(
            zaxis = dict(
                range=[-20,120]
            )
        )
    )

    threeD = go.FigureWidget( data=[pes,pts] , layout=layout)

    return threeD

def molGeo(filename):
    # cutConvert(str(sys.argv[4]),order)

    # Create the MolGraph object
    mg = MolGraph()

    # Read the data from the .xyz file
    mg.read_xyz(filename)

    # Create the Plotly figure object
    mol = to_plotly_figure(mg)
    
    return mol 

def ani(twoD,threeD,traj):
    # test for printing full trajectory 
    # twoD.plot(traj[:,0],traj[:,1],'r-')
    # threeD.plot(traj[:,0],traj[:,1],traj[:,2],'r-',zorder=10)

    ncol,nrow = np.shape(traj)
    t = np.linspace(0, 80, ncol)
    x = traj[:,0]
    y = traj[:,1]
    z = traj[:,2]

    frames = [ go.Frame(
        data = [ go.Scatter(
            x = x, y = x,
            mode='markers',
            marker=dict(color="red", size=10) )
        ] )  ]
    # FIXME: figtest = plt.figure()
    # i = 1
    # molGeo(i,mol)
    # im1 = mpimg.imread('test.png')
    # mol.imshow(im1)
    
    # lines = []
    # for i in range(ncol):
    # # for i in range(1):
    #     # molecular structure
    #     # molGeo(i,mol)
    #     # im1 = mpimg.imread('test.png')
    #     # mol.imshow(im1)
    #     # trajectory
    #     head = i - 1
    #     head_slice = (t > t[i] - 1.0) & (t < t[i])
    #     #print(x[head],y[head],z[head])
    #     # 3-D
    #     line1,  = threeD.plot(x[:i], y[:i], z[:i],color='black')
    #     line1a, = threeD.plot(x[head_slice], y[head_slice], z[head_slice],color='red', linewidth=2)
    #     line1e, = threeD.plot([x[head]], [y[head]], [z[head]],color='red', marker='o', markeredgecolor='r')
    #     # 2-D
    #     line2,  = twoD.plot(x[:i], y[:i],color='black')
    #     line2a, = twoD.plot(x[head_slice], y[head_slice],color='red', linewidth=2)
    #     line2e, = twoD.plot(x[head], y[head],color='red', marker='o', markeredgecolor='r')
    #     lines.append([line1,line1a,line1e,line2,line2a,line2e])
    #     # default: view_init(30,-60) FIXME: how to rotate it
    #     # threeD.view_init(elev=10.,azim=i)

    # ani = animation.ArtistAnimation(fig, lines, interval=100, blit=True)
    return frames

main()

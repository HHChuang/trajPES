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
#   1. MDanalysis
#   https://www.youtube.com/watch?v=zVQGFysYDew
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
from copy import deepcopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
# 
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objects as go 
from xyz2graph import MolGraph, to_plotly_figure
# import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import offline, iplot

import matplotlib.animation as animation

workDir = os.getcwd()
PATH_fig = workDir

title = ''
xaxis = 'Modified IRC of TSS1'
yaxis = 'IRC of TSS2'
zaxis = 'kcal/mol'
axis_fontsize = 6

# static figure objects
fig = plt.figure()
threeD = fig.add_subplot(1, 2, 2,projection='3d') # top and bottom right 
mol = fig.add_subplot(2,2,1)
mol.axis('off') # remove borders, ticks and labels
twoD = fig.add_subplot(2, 2, 3)  # bottom left

# interact objects
# fig = make_subplots(rows = 2, cols = 2)
# threeD = 


def main():
    # 1. Import rawdata
    X, Y, E, pts_name, pts_coord = getPES()
    traj = getTraj(pts_coord) # format: x,y,E
    
    # 2. Plot the background: 2D-PES and 3D-PES
    twoDwPts(twoD, X, Y, E, pts_name, pts_coord)
    threeDwPts(threeD, X, Y, E, pts_name, pts_coord)
    
    # 3. Combine above subroutines and export them into one video
    ani(fig,twoD,threeD,traj)
    fig.show()

def totline(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1

def cutConvert(filename,order):
    with open(filename) as f:
        natom = int(f.readline())
    f.close()
    

def getTraj(pts_coord):
    # read the coordinate of trajectory
    tline=totline(str(sys.argv[3]))
    traj=[]
    with open(str(sys.argv[3])) as f:
        for i,line in enumerate(f):
            traj.append(line.split())
    f.close()

    traj = np.array(traj)
    traj = traj.astype(np.float)

    # Calculate relative energy in kcal/mol, and use reactant energy as the energy 
    ncol,nrow = np.shape(traj)
    E_R = pts_coord[0][2]
    relE = deepcopy(traj) #copy by values, not by reference
    for i in range(ncol-1):
        relE[i][2] = (traj[i][2] - E_R) * 627.5095

    return relE

def getPES():
    # import PES w/ x, y and z coord
    # TODO: change square to recentagular
    tline = totline(str(sys.argv[1]))
    dim = int(math.sqrt(tline))
    x, y, e = np.loadtxt(str(sys.argv[1]), delimiter=' ', unpack=True)
    X = np.reshape(x, (dim, dim))
    Y = np.reshape(y, (dim, dim))
    E = np.reshape(e, (dim, dim))

    # import name and coord of pts
    # FIXME: pts_coord may have datatype problem
    pts_name = list()
    pts_coord = list()
    with open(str(sys.argv[2])) as f:
        for i, line in enumerate(f):
            pts_name.append(line.split()[0])
            pts_coord.append(line.split()[1:4])
    f.close()
    pts_name = np.array(pts_name)
    pts_coord = np.array(pts_coord)
    pts_coord = pts_coord.astype(np.float)
    return X, Y, E, pts_name, pts_coord

def twoDwPts(twoD, X, Y, E, pts_name, pts_coord):
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    plt.gca().set_aspect('equal')  # aspect ratio of x and y is equal

    # Calculate relative energy in kcal/mol, and use reactant energy as the energy reference
    E_R = pts_coord[0][2]
    dim = int(math.sqrt(np.size(E)))
    relE = np.arange(dim*dim).reshape(dim,dim)
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - E_R) * 627.5095

    ct = twoD.contour(X, Y, relE, 10, colors='k')
    twoD.set_xlabel(xaxis)
    twoD.set_ylabel(yaxis)

    # plot important points
    npts = totline(str(sys.argv[2]))
    # NCH1
    shift_x = [0, -5, 1, 0, 0]
    shift_y = [5, 5, 5, 5, -5]
    for n in range(npts):
        twoD.plot(pts_coord[n][0], pts_coord[n][1], 'bo',markersize=4,markeredgecolor='k')
        twoD.text(pts_coord[n][0] + shift_x[n], pts_coord[n][1]+shift_y[n],
                 pts_name[n], weight='bold', backgroundcolor='white',
                 verticalalignment='top', multialignment='right', fontsize=6)
    
def threeDwPts(threeD, X, Y, E, pts_name, pts_coord):
    # change hartree to relative energy
    npts = totline(str(sys.argv[2]))
    #   energy reference : reactant
    Renergy = pts_coord[0][2]
    dim2 = np.size(E)
    dim = int(math.sqrt(dim2))
    relE = np.arange(dim*dim).reshape(dim,dim)
    for i in range(dim):
        for j in range(dim):
            relE[i][j] = (E[i][j] - Renergy) * 627.5095

    # ratio of 3D box FIXME:
    # threeD.pbaspect = [2.0, 0.6, 0.25]
    # threeD.set_aspect('equal')  # ratio of x and y axes, 'auto', 'equal'

    threeD.set_xlabel(xaxis)
    threeD.set_ylabel(yaxis)
    threeD.set_zlabel(zaxis)
    # range of z-axis; range of energy profile
    zmin = -20
    zmax = 100
    threeD.set_xlim(X.min(), X.max())
    threeD.set_ylim(Y.min(), Y.max())
    threeD.set_zlim(zmin, zmax) 

    # color of surface; cm.coolwarm
    # surf = threeD.plot_surface(X, Y, relE,  rstride=1, cstride=1,cmap=mpl.cm.Spectral_r, linewidth=1, antialiased=False,vmin=zmin,vmax=zmax,shade=False)
    # fig.colorbar(surf, shrink=0.5, aspect=10)
    threeD.plot_wireframe(X,Y,relE,alpha=0.5,rcount=100,ccount=100,linewidth=0.5)

    # ct = plt.contour(X, Y, E, 25, colors='k')

    # plot important points
    npts = totline(str(sys.argv[2]))
    # 
    # shift_x = [-0.5, -0.6, -0.6, -0.8, -1.0]
    shift_x = [0,0,0,0,0]
    shift_y = [5,5,5,5,5]
    for i in range(npts):
        pts_coord[i][2] = (pts_coord[i][2] - Renergy) * 627.5095

    for n in range(npts):
        threeD.scatter(pts_coord[n][0], pts_coord[n][1],
                   pts_coord[n][2], marker='o', c='b', edgecolors='k', zorder=10)
        threeD.text(pts_coord[n][0] + shift_x[n], pts_coord[n][1] + shift_y[n], pts_coord[n]
                [2] + 2, pts_name[n], weight='bold', backgroundcolor='white', zorder=10, fontsize=6)

def molGeo(order,mol):
    # cutConvert(str(sys.argv[4]),order)
    mg = MolGraph()

    mg.read_xyz(str(sys.argv[4]))
    # mg.read_xyz( str(sys.argv[4]) )#,order)
    # Create the Plotly figure object
    mol = to_plotly_figure(mg)
    # pio.renderers.default = "png"
    # mol.show()
    # pio.renderers.default = 'png'
    # mol.show()
    # test = go.Figure(data=mol) 
    # test.show(renderer='png')

    # Plot the figure
    # iplot(mol)
    # Image(pio.to_image(mol,format='png'))
    #  convert the figure to a PNG bytes object
    # img_bytes = mol.to_image(format='png') 
    # Image(pio.to_image(mol,format='png'))
    # mol.show(renderer='png')
    # mol.write_image('test.png') #work

def ani(fig,twoD,threeD,traj):
    # test for printing full trajectory 
    # twoD.plot(traj[:,0],traj[:,1],'r-')
    # threeD.plot(traj[:,0],traj[:,1],traj[:,2],'r-',zorder=10)

    ncol,nrow = np.shape(traj)
    t = np.linspace(0, 80, ncol)
    x = traj[:,0]
    y = traj[:,1]
    z = traj[:,2]

    # FIXME: figtest = plt.figure()
    i = 1
    molGeo(i,mol)
    im1 = mpimg.imread('test.png')
    mol.imshow(im1)
    
    lines = []
    for i in range(ncol):
    # for i in range(1):
        # molecular structure
        # molGeo(i,mol)
        # im1 = mpimg.imread('test.png')
        # mol.imshow(im1)
        # trajectory
        head = i - 1
        head_slice = (t > t[i] - 1.0) & (t < t[i])
        #print(x[head],y[head],z[head])
        # 3-D
        line1,  = threeD.plot(x[:i], y[:i], z[:i],color='black')
        line1a, = threeD.plot(x[head_slice], y[head_slice], z[head_slice],color='red', linewidth=2)
        line1e, = threeD.plot([x[head]], [y[head]], [z[head]],color='red', marker='o', markeredgecolor='r')
        # 2-D
        line2,  = twoD.plot(x[:i], y[:i],color='black')
        line2a, = twoD.plot(x[head_slice], y[head_slice],color='red', linewidth=2)
        line2e, = twoD.plot(x[head], y[head],color='red', marker='o', markeredgecolor='r')
        lines.append([line1,line1a,line1e,line2,line2a,line2e])
        # default: view_init(30,-60) FIXME: how to rotate it
        # threeD.view_init(elev=10.,azim=i)

    ani = animation.ArtistAnimation(fig, lines, interval=100, blit=True)

main()

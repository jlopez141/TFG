import plotly.offline as py
import plotly.graph_objs as go
import math
import numpy as np

magnitudes = ['young','compr']

for archivo in magnitudes:

    if (archivo == 'young'):
        titulo = 'Modulo de Young'
    else:
    	titulo = 'Coeficiente de compresibilidad lineal'


    f = open('3D_plot_'+archivo+'.txt', 'r')

    dx = [[0.0 for i in range(250)] for j in range(250)]
    dy = [[0.0 for i in range(250)] for j in range(250)]
    dz = [[0.0 for i in range(250)] for j in range(250)]
    dr = [['0' for i in range(250)] for j in range(250)]

    for line in f:
    	l = line.split()
    	teta = np.pi * (int(l[0])-1)/248
    	phi = 2*np.pi * (int(l[1])-1)/248
    	dx[int(l[0])-1][int(l[1])-1] = float(l[2])
    	dy[int(l[0])-1][int(l[1])-1] = float(l[3])
    	dz[int(l[0])-1][int(l[1])-1] = float(l[4])
    	dr[int(l[0])-1][int(l[1])-1] = '('+str(teta*180.0/np.pi)+', '+str(phi*180.0/np.pi)+') :: '+str(math.sqrt(float(l[2])**2+float(l[3])**2+float(l[4])**2))

    f.close()


    data = go.Surface (x=dx, y=dy, z=dz, text=dr)
    
    data['colorscale'] = [[0,'rgb(255,255,100)'],[0.2,'rgb(255,200,50)'],[0.95,'rgb(255,73,0)'],[1,'rgb(202,52,0)']]

    data['contours'] = contours=dict(x=dict(show=True, color='rgb(250,250,250)'), y=dict(show=True, color='rgb(250,250,250)'), z=dict(show=True, color='rgb(250,250,250)'))

    def dist_origin(x, y, z):
    	return math.sqrt((1.0 * x)**2 + (1.0 * y)**2 + (1.0 * z)**2)

    lx = len(dx)
    ly = len(dy)
    out=[]
    for i in xrange(lx):
    	temp = []
    	for j in xrange(ly):
       	    temp.append(
            	dist_origin(dx[i][j], dy[i][j], dz[i][j]))
    	out.append(temp)

    data['surfacecolor'] = out

    data = go.Data([data])

    layout = go.Layout (title = titulo)


    fig = go.Figure (data = data, layout = layout)
    py.plot (fig, filename=archivo+'.html')

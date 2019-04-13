
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
data = np.genfromtxt("/Users/anshumanacharya/Desktop/Figs/colden_Z1.005_Q14.txt", usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
hden=[]
z=np.arange(-5,-1.99,0.05)
for i in z:
    hden.append(i)
y0=go.Scatter(x=hden, y=data[:,0], mode = 'lines+markers', name='CII', line = dict(color = 'red',width = 2))
y1=go.Scatter(x=hden, y=data[:,1], mode = 'lines+markers', name='CIII', line = dict(color = 'cyan',width = 2))
y2=go.Scatter(x=hden, y=data[:,2], mode = 'lines+markers', name='CIV', line = dict(color = 'black',width = 2))
y3=go.Scatter(x=hden, y=data[:,3], mode = 'lines+markers', name='CV', line = dict(color = 'yellow',width = 2))
y4=go.Scatter(x=hden, y=data[:,4], mode = 'lines+markers', name='MgII', line = dict(color = 'darkorchid',width = 2))
y5=go.Scatter(x=hden, y=data[:,5], mode = 'lines+markers', name='NeVIII', line = dict(color = 'maroon',width = 2))
y6=go.Scatter(x=hden, y=data[:,6], mode = 'lines+markers', name='NII', line = dict(color = 'dodgerblue',width = 2))
y7=go.Scatter(x=hden, y=data[:,7], mode = 'lines+markers', name='NIII', line = dict(color = 'olive',width = 2))
y8=go.Scatter(x=hden, y=data[:,8], mode = 'lines+markers', name='NIV', line = dict(color = 'orange',width = 2))
y9=go.Scatter(x=hden, y=data[:,9], mode = 'lines+markers', name='OII', line = dict(color = 'hotpink',width = 2))
y10=go.Scatter(x=hden, y=data[:,10], mode = 'lines+markers', name='OIII', line = dict(color = 'peru',width = 2))
y11=go.Scatter(x=hden, y=data[:,11], mode = 'lines+markers', name='OIV', line = dict(color = 'blue',width = 2))
y12=go.Scatter(x=hden, y=data[:,12], mode = 'lines+markers', name='OV', line = dict(color = 'springgreen',width = 2))
y13=go.Scatter(x=hden, y=data[:,13], mode = 'lines+markers', name='OVI', line = dict(color = 'gold',width = 2))
y14=go.Scatter(x=hden, y=data[:,14], mode = 'lines+markers', name='SiII', line = dict(color = 'salmon',width = 2))
y15=go.Scatter(x=hden, y=data[:,15], mode = 'lines+markers', name='SiIII', line = dict(color = 'aquamarine',width = 2))
y16=go.Scatter(x=hden, y=data[:,16], mode = 'lines+markers', name='SiIV', line = dict(color = 'sienna',width = 2))
y17=go.Scatter(x=hden, y=data[:,17], mode = 'lines+markers', name='SIV', line = dict(color = 'deeppink',width = 2))
y18=go.Scatter(x=hden, y=data[:,18], mode = 'lines+markers', name='SV', line = dict(color = 'seagreen',width = 2))
y19=go.Scatter(x=hden, y=data[:,19], mode = 'lines+markers', name='SVI', line = dict(color = 'thistle',width = 2))
data =[y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19]
layout = {'xaxis': {'range': [-5, -2]},'yaxis': {'range': [0, 15]},'shapes':[
{'type': 'line','x0': -5,'y0': 10.84414403,'x1': -2,'y1': 10.84414403,'line': {'color': 'red','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 12.48897916,'x1': -2,'y1': 12.48897916,'line': {'color': 'cyan','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 12.35609281,'x1': -2,'y1': 12.35609281,'line': {'color': 'black','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 12.58180862,'x1': -2,'y1': 12.58180862,'line': {'color': 'yellow','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 7.735584529,'x1': -2,'y1': 7.735584529,'line': {'color': 'darkorchid','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 9.041017451,'x1': -2,'y1': 9.041017451,'line': {'color': 'maroon','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 9.94387028,'x1': -2,'y1': 9.94387028,'line': {'color': 'dodgerblue','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 11.91367189,'x1': -2,'y1': 11.91367189,'line': {'color': 'olive','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 12.229244,'x1': -2,'y1': 12.229244,'line': {'color': 'orange','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 10.26909947,'x1': -2,'y1': 10.26909947,'line': {'color': 'hotpink','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 12.34694142,'x1': -2,'y1': 12.34694142,'line': {'color': 'peru','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 13.04273692,'x1': -2,'y1': 13.04273692,'line': {'color': 'blue','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 12.65870571,'x1': -2,'y1': 12.65870571,'line': {'color': 'springgreen','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 11.90272774,'x1': -2,'y1': 11.90272774,'line': {'color': 'gold','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 9.182157766,'x1': -2,'y1': 9.182157766,'line': {'color': 'salmon','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 10.59815008,'x1': -2,'y1': 10.59815008,'line': {'color': 'aquamarine','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 10.68641728,'x1': -2,'y1': 10.68641728,'line': {'color': 'sienna','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 11.34396563,'x1': -2,'y1': 11.34396563,'line': {'color': 'deeppink','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 10.95958654,'x1': -2,'y1': 10.95958654,'line': {'color': 'seagreen','width': 2,'dash': 'dashdot',},},
{'type': 'line','x0': -5,'y0': 10.97589343,'x1': -2,'y1': 10.97589343,'line': {'color': 'thistle','width': 2,'dash': 'dashdot',},},
{
            'type': 'rect',
            'xref': 'x',
            'yref': 'y',
            'x0': '-3.987480028',
            'y0': 0,
            'x1': '-3.270872913',
            'y1': 15,
            'fillcolor': 'gray',
            'opacity': 0.4,
            'line': {'width': 0,}},
{
            'type': 'rect',
            'xref': 'x',
            'yref': 'y',
            'x0': '-4.345783586',
            'y0': 0,
            'x1': '-3.987480028',
            'y1': 15,
            'fillcolor': 'gray',
            'opacity': 0.2,
            'line': {'width': 0,}},
{
            'type': 'rect',
            'xref': 'x',
            'yref': 'y',
            'x0': '-3.270872913',
            'y0': 0,
            'x1': '-2.912569355',
            'y1': 15,
            'fillcolor': 'gray',
            'opacity': 0.2,
            'line': {'width': 0,}},
]}
fig = {'data': data,'layout': layout}
py.iplot(fig, filename='ForZ_1.005_Q14')
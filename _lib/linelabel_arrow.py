from math import atan2, degrees
import numpy as np

from matplotlib.dates import date2num
from matplotlib import patches
from datetime import datetime


# Label line with line2D label data
def labelLine(line, x, arrow=False, label=None, align=True,  **kwargs):
    '''Label a single matplotlib line at position x

    Parameters
    ----------
    line : matplotlib.lines.Line
       The line holding the label
    x : number
       The location in data unit of the label
    arrow : draw arrow under label
    label : string, optional
       The label to set. This is inferred from the line by default
    kwargs : dict, optional
       Optional arguments passed to ax.text
    '''
    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    # Convert datetime objects to floats
    if isinstance(x, datetime):
        x = date2num(x)

    if (x < xdata[0]) or (x > xdata[-1]):
        raise Exception('x label location is outside data range!')

    # Find corresponding y co-ordinate and angle of the
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1]) * \
        (x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        # Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy, dx))

        # Transform to screen co-ordinates
        pt = np.array([x, y]).reshape((1, 2))
        trans_angle = ax.transData.transform_angles(np.array((ang, )), pt)[0]

    else:
        trans_angle = 0

    # Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x, y, label, rotation=trans_angle, **kwargs)

    if arrow:

        # get axis position   
        if ax.yaxis.get_label_position() == 'right':
            x_dir = 1
            
        if ax.yaxis.get_label_position() == 'left':
            x_dir = -1       
        
        # get y axis scaling
        y_scaling = ax.get_yticks()[1]-ax.get_yticks()[0]
        
        # arrow direction
        x_arrow_start = x-(x_dir*0.25*x)
        x_arrow_end = x+(x_dir*0.25*x)
        y_arrow = y-0.3*y_scaling  
        
        # draw arrow
        linecolor = line.get_color()
        ax.annotate("", xytext=(x_arrow_start, y_arrow),
                    xy=(x_arrow_end, y_arrow),arrowprops=dict(lw=2, arrowstyle='->', color=linecolor), **kwargs)
        
        #draw ellipse
        #x_ellipse_width = 0.15*x
        #x_ellipse_height = 0.15*y
        #x_ellipse_center = x-(x_dir*0.25*x)-x_ellipse_width/2
        #e1 = patches.Ellipse((x_ellipse_center, y), x_ellipse_width, x_ellipse_height, fill=False, edgecolor=linecolor)
        #ax.add_patch(e1)

def labelLines(lines, align=True, xvals=None, draw_arrow=False, **kwargs):
    '''Label all lines with their respective legends.

    Parameters
    ----------
    lines : list of matplotlib lines
       The lines to label
    align : boolean, optional
       If True, the label will be aligned with the slope of the line
       at the location of the label. If False, they will be horizontal.
    xvals : (xfirst, xlast) or array of float, optional
       The location of the labels. If a tuple, the labels will be
       evenly spaced between xfirst and xlast (in the axis units).
    kwargs : dict, optional
       Optional arguments passed to ax.text
    '''
    ax = lines[0].axes
    labLines = []
    labels = []

    # Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xvals = ax.get_xlim() # set axis limits as annotation limits, xvals now a tuple
    if type(xvals) == tuple:
        xmin, xmax = xvals
        xscale = ax.get_xscale()
        if xscale == "log":
            xvals = np.logspace(np.log10(xmin), np.log10(xmax), len(labLines)+2)[1:-1]
        else:
            xvals = np.linspace(xmin, xmax, len(labLines)+2)[1:-1]

    for line, x, label in zip(labLines, xvals, labels):
        labelLine(line, x, draw_arrow, label, align, **kwargs)

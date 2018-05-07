import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def bottomaxisonly(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    
def leftaxisonly(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([]) 
        
def check_dim(*args):
    size = {len(args[i]) for i in range(len(args)) if args[i] != list()}
    assert len(size) == 1, "Different size between parameters"
    #for i in range(len(args)-1):
        #if args[i] != list() and args[i+1] != list():
            #assert len(args[i]) == len(args[i+1]), "Different size between parameters"
        
def plot_scatter_hist(data_x, data_y, alpha=0.8, labels=list(), xlabel="", ylabel="", savefig=None, show=True,
                      plot_corr=False):
    """ combine a scatter plot with two histograms, one for each of the axis
    """
    check_dim(data_x, data_y, labels)
    #colors = mpl.rcParams['axes.color_cycle']
    colors = [color['color'] for color in list(mpl.rcParams['axes.prop_cycle'])]
    nullfmt = NullFormatter()   
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    fig = plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    simpleaxis(axScatter)
    bottomaxisonly(axHistx)
    leftaxisonly(axHisty)
    
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    ps = list()
    xmax, xmin = -np.inf, np.inf
    ymax, ymin = -np.inf, np.inf
    k = 0
    if labels != list():
        for i in range(len(data_x)):
            p = axScatter.scatter(data_x[i], data_y[i], label=labels[i], color=colors[k], alpha=alpha)
            # now determine nice limits
            xmax = max(np.max(data_x[i]), xmax)
            xmin = min(np.min(data_x[i]), xmin)
            ymax = max(np.max(data_y[i]), ymax)
            ymin = min(np.min(data_y[i]), ymin)
            ps.append(p)
            k += 1
            if k >= len(colors):
                k = 0
    else:
        for i in range(len(data_x)):
            p = axScatter.scatter(data_x[i], data_y[i], color = colors[k], alpha=alpha)
            # now determine nice limits
            xmax = max(np.max(data_x[i]), xmax)
            xmin = min(np.min(data_x[i]), xmin)
            ymax = max(np.max(data_y[i]), ymax)
            ymin = min(np.min(data_y[i]), ymin)
            ps.append(p)
            k += 1
            if k >= len(colors):
                k = 0

    if plot_corr:
        k = 0
        for i in range(len(data_x)):
            x = np.linspace(min(data_x[i]), max(data_x[i]), 100)
            cor, pval = scst.pearsonr(data_x[i], data_y[i])
            slope, intercept, r_value, p_value, std_err = scst.linregress(data_x[i], data_y[i])
            r2 = r_value ** 2
            print(i, cor, r2)
            f = lambda x , a, b : x * a + b
            y = f(x, slope, intercept)
            axScatter.plot(x, y, linestyle="-.", color=colors[k])
        k += 1
    
    xsize = max(len(str(int(xmax))), len(str(int(xmin)))) - 1 
    ysize = max(len(str(int(ymax))), len(str(int(ymin)))) - 1
    xoffset = 0.1 * (10 ** xsize)
    yoffset = 0.1 * (10 ** ysize)
    
    axScatter.set_xlim((xmin-xoffset, xmax+xoffset))
    axScatter.set_ylim((ymin-yoffset, ymax+yoffset))

    xbins = np.linspace(xmin, xmax, 21)
    k = 0
    for i in range(len(data_x)):
        axHistx.hist(data_x[i], bins=xbins, alpha=alpha)
        axHistx.axvline(x=np.median(data_x[i]), ymin=0,  ymax=10, color=colors[k], linestyle="--")
        k += 1
        if k >= len(colors):
            k = 0
    
    ybins = np.linspace(ymin, ymax, 21)
    k = 0
    for i in range(len(data_y)):
        axHisty.hist(data_y[i], bins=ybins, orientation='horizontal', alpha=alpha)
        axHisty.axhline(y=np.median(data_y[i]), xmin=0,  xmax=10, color=colors[k], linestyle="--")
        k += 1
        if k >= len(colors):
            k = 0
            
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    
    if xlabel != list():
        axScatter.set_xlabel(xlabel)
    if ylabel != list():
        axScatter.set_ylabel(ylabel)
    if labels != list():
        plt.legend(handles=ps)
        
    if savefig:
        plt.savefig(savefig, dpi=600)
    if show:
        plt.show()
    return [axScatter, axHistx, axHisty]

def jitter_plot(data, labels=list(), xlabel=list(), ylabel="", alpha=0.5, savefig=None, show=True):
    """ create a jitter plot
    """
    fig, axJitter = plt.subplots()
    check_dim(data, labels, xlabel)
    simpleaxis(axJitter)
    
    indexes = np.arange(1, len(data)+1)
    for i in range(len(data)):
        x = np.random.rand(len(data[i])) * 0.5 + (i+1-0.25)
        labeli=labels[i] if i < len(labels) else ""
        axJitter.scatter(x, data[i], s=12, alpha=alpha, label=labeli)
        ymean = sum(data[i])/len(data[i])
        axJitter.plot([i+1-0.20, i+1.20], [ymean, ymean], linewidth=2, linestyle="--")
    
    axJitter.set_xticks(indexes)
    if xlabel != list():
        axJitter.set_xticklabels(xlabel)
    if ylabel:
        axBox.set_ylabel(ylabel)
        
    plt.xlim(0.5, len(indexes)+0.5)
   
    if savefig:
        plt.savefig(savefig, dpi=600)
    if show:
        plt.show()
    return axJitter

def fancy_box(data, labels=list(), xlabel=list(), ylabel="", alpha=0.5, savefig=None, show=True):
    """ make a fancy box plot -> jitter + box
    """
    fig, axBox = plt.subplots()
    check_dim(data, labels, xlabel)
    simpleaxis(axBox)
    
    indexes = np.arange(1, len(data)+1)
    axBox.boxplot([y_neu, y_hdel], labels=["neutral", "highly deleterious"], showfliers="",)
    for i in range(len(data)):
        x = np.random.rand(len(data[i])) * 0.5 + (i+1-0.25)
        labeli=labels[i] if i < len(labels) else ""
        axBox.scatter(x, data[i], s=12, alpha=alpha, label=labeli)
        #ymean = sum(data[i])/len(data[i])
        ymean = np.median(np.array(data[i]))
        axBox.plot([i+1-0.20, i+1.20], [ymean, ymean], linewidth=2, linestyle="--")
    
    axBox.set_xticks(indexes)
    if xlabel != list():
        axBox.set_xticklabels(xlabel)
    
    if ylabel:
        axBox.set_ylabel(ylabel)
        
    plt.xlim(0.5, len(indexes)+0.5)
    
    if savefig:
        plt.savefig(savefig, dpi=600)
    if show:
        plt.show()
    return axBox

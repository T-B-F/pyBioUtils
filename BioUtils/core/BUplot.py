import numpy as np
import scipy.stats as scst
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram


def simple_and_outward(ax, out=10):
    simpleaxis(ax)
    outward_spines(ax, out)

def outward_spines(ax, out=10):
    ax.spines["left"].set_position(("outward", out))
    ax.spines["bottom"].set_position(("outward", out))


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
            axScatter.plot(x, y, linestyle="-.", color="black")
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


def clustermap(x, draw_top=True, draw_left=True, colorbar_pad=0.5, cmap=plt.cm.viridis,
               col_labels = None, row_labels = None, xlabel_rotation = -45, 
               ylabel_rotation = 0, label_fontsize = 8, figsize=(12, 8)):
    ''' adapted from https://github.com/WarrenWeckesser/heatmapcluster/blob/master/heatmapcluster.py
    '''
    assert draw_top or draw_left, "Warning, must at least specify one hystogram (use standard heatmap otherwise)"

    if col_labels is None:
        col_labels = np.arange(x.shape[1])
    if row_labels is None:
        row_labels = np.arange(x.shape[0])
        
    fig, ax_heatmap = plt.subplots(figsize=figsize)
    ax_heatmap.yaxis.tick_right()
    divider = axes_grid1.make_axes_locatable(ax_heatmap)

    if draw_left:
        ax_dendleft = divider.append_axes("left", 1.2, pad=0.0, sharey=ax_heatmap)
        ax_dendleft.set_frame_on(False)
        left_threshold = -1
        side_orientation = 'left'
        lnk0 = linkage(pdist(x))
        dg0 = dendrogram(lnk0, ax=ax_dendleft, orientation=side_orientation, color_threshold=left_threshold, no_labels=True)
    if draw_top:
        ax_dendtop = divider.append_axes("top", 1.2, pad=0.0, sharex=ax_heatmap)
        ax_dendtop.set_frame_on(False)
        top_threshold = -1
        lnk1 = linkage(pdist(x.T))
        dg1 = dendrogram(lnk1, ax=ax_dendtop, color_threshold=top_threshold, no_labels=True)

    colorbar_width = 0.45
    ax_colorbar = divider.append_axes("right", colorbar_width, pad=colorbar_pad)

    # Reorder the values in x to match the order of the leaves of
    # the dendrograms.
    if draw_left:
        z = x[dg0['leaves'], :]
    else:
        z = x
    if draw_top:
        z = z[:, dg1['leaves']]
    
    if draw_top:
        ymax = ax_dendtop.get_xlim()[1]
    else:
        ymax = ax_dendleft.get_xlim()[1]

    if draw_left:
        im = ax_heatmap.imshow(z[::-1], aspect='auto', cmap=cmap, interpolation='nearest',
                              extent=(0, ymax, 0, ax_dendleft.get_ylim()[1]))
    else:
        im = ax_heatmap.imshow(z[::-1], aspect='auto', cmap=cmap, interpolation='nearest',
                               extent=(0, ymax, 0, ax_heatmap.get_ylim()[1]))

    
    xlim = ax_heatmap.get_xlim()[1]        
    ncols = len(col_labels)
    halfxw = 0.5*xlim/ncols
    ax_heatmap.xaxis.set_ticks(np.linspace(halfxw, xlim - halfxw, ncols))
    if draw_top:
        ax_heatmap.xaxis.set_ticklabels(np.array(col_labels)[dg1['leaves']])
    else:
        ax_heatmap.xaxis.set_ticklabels(col_labels)

    ylim = ax_heatmap.get_ylim()[1]
    nrows = len(row_labels)
    halfyw = 0.5*ylim/nrows
    ax_heatmap.yaxis.set_ticks(np.linspace(halfyw, ylim - halfyw, nrows))
    if draw_left:
        ax_heatmap.yaxis.set_ticklabels(np.array(row_labels)[dg0['leaves']])
    else:
        ax_heatmap.yaxis.set_ticklabels(row_labels)

    # Make the dendrogram labels invisible.
    if draw_left:
        plt.setp(ax_dendleft.get_yticklabels() + ax_dendleft.get_xticklabels(), visible=False)
    if draw_top:
        plt.setp(ax_dendtop.get_xticklabels() + ax_dendtop.get_yticklabels(), visible=False)

    # Hide all tick lines.
    lines = (ax_heatmap.xaxis.get_ticklines() +
             ax_heatmap.yaxis.get_ticklines())
    plt.setp(lines, visible=False)

    if draw_left:
        lines = (ax_dendleft.xaxis.get_ticklines() +
             ax_dendleft.yaxis.get_ticklines())
        plt.setp(lines, visible=False)
        
    if draw_top:
        lines = (ax_dendtop.xaxis.get_ticklines() +
                 ax_dendtop.yaxis.get_ticklines())
        plt.setp(lines, visible=False)

    xlbls = ax_heatmap.xaxis.get_ticklabels()
    plt.setp(xlbls, rotation=xlabel_rotation)
    plt.setp(xlbls, fontsize=label_fontsize)

    ylbls = ax_heatmap.yaxis.get_ticklabels()
    plt.setp(ylbls, rotation=ylabel_rotation)
    plt.setp(ylbls, fontsize=label_fontsize)

    cb = plt.colorbar(im, cax=ax_colorbar)
    # This code to draw the histogram in the colorbar can
    # probably be simplified.
    # Also, there are several values that someone, sometime,
    # will probably want to change, but for now, all the
    # details are hardcoded.
    nbins = min(80, max(int(x.size/10+0.5), 11))
    counts, edges = np.histogram(x.ravel(), bins=nbins)
    max_count = counts.max()
    counts = counts / max_count
    edges = (edges - edges[0])/(edges[-1] - edges[0])
    # cc and ee contain the values in counts and edges, repeated
    # as needed to draw the histogram curve similar to the 'steps-mid'
    # drawstyle of the plot function.
    cc = np.repeat(counts, 2)
    ee = np.r_[edges[0], np.repeat(edges[1:-1], 2), edges[-1]]

    ax_colorbar.plot(cc, ee, 'k', alpha=0.5)
    ax_colorbar.xaxis.set_ticks([0, 1])
    pctstr = '%.2g%%' % (100*max_count/x.size)
    ax_colorbar.xaxis.set_ticklabels(['0', pctstr])
    ax_colorbar.xaxis.set_label_text('Histogram\n(% count)')

    plt.show()


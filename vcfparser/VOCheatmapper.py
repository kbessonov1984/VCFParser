import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np


def heatmap(data, row_labels, col_labels, ax=None, title="",
            cbar_kw={}, axis_kw={}, cbarlabel="",  **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar legend
    cbar = ax.figure.colorbar(im, ax=ax, aspect=10, fraction=0.10, drawedges=True, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom",size=3)
    cbar.ax.tick_params(labelsize=3)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, **axis_kw)
    ax.set_yticklabels(row_labels, **axis_kw)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="grey", linestyle='-', linewidth=0.1)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title(title, fontsize='3', fontstyle='oblique', fontweight='bold')

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None,
                     covdata=[],
                     is_annotate = False,
                     **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if type(textcolors) is not list:
        raise Exception("textcolors should be a list and not string")
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)


    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []


    for i in range(data.shape[0]):
        for j in range(data.shape[1]):

            if threshold:
                kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)]) #return 0 or 1
            else:
                kw.update(color=textcolors[0])


            value =  valfmt(data[i, j], None)

            if is_annotate and covdata:
                if covdata[j][i] == 0:
                    value = "NC"
                elif float(value) == 0:
                    value = ""
            elif is_annotate:
                if float(value) == 0:
                    value = ""
            elif covdata:
                value = ""
                if covdata[j][i] == 0:
                    value = "NC"
            else:
                value=""


            text = im.axes.text(j, i, value, **kw)
            texts.append(text)





    return texts


def renderplot(data=pd.DataFrame(),debug=False, title="",
               is_plot_annotate = False,
               axis=None, read_coverages_2Darray=list()):

    if debug: #DEBUG read in input dataframe with x and y labels for testing only purposes
        data = pd.read_csv('test/SNVsamplesummary.tsv', sep="\t")
        SNVnames = data["NucName+AAName"].tolist()
        data = data.iloc[:, 1:]
    else:
        SNVnames = data["NucName+AAName"].tolist()
        column_number = data.columns.to_list().index("NucName+AAName")+1
        data = data.iloc[:, column_number:]




    #data = np.array(
    #    [[np.array(i, dtype=np.float64) for i in j] for j in summarydf.values],
    #    dtype=np.float16,
    #)
    #print(summarydf, summarydf.index, summarydf.columns);exit()



    #data = np.random.random_sample((6, 6))
    #data[0,0]=1; data[5,5]=0; data[0,1]=0.5


    #y = ["Position {}".format(i) for i in range(1, data.shape[0]+1, 1)]
    y = ["{}".format(i) for i in SNVnames]
    #x = ["Sample {}".format(i) for i in range(1, data.shape[1]+1)]
    x = ["{}".format(s) for s in data.columns]


    #qrates = list("ABCDEFGJKL")
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.00001,1.1,0.1), 11)


    #fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])
    # plt.figure(figsize=(2, 3.5), dpi=300)
    cmap=plt.get_cmap("YlGnBu", 11) #RdYlGn red green
    cmap.set_over('black')
    cmap.set_under('white')

    im, _ = heatmap(data, y, x,
                    ax=axis,
                    cmap=cmap,
                    norm=norm,
                    axis_kw=dict(size=3),
                    cbar_kw=dict(ticks=np.arange(0.00001,1.1,0.1),
                                 extend="both",
                                 orientation="vertical",
                                 spacing='proportional',
                                 shrink=1),#, format=fmt),
                    cbarlabel="ALT_FREQ",
                    title=title)


    #if is_text_annotate:
    annotate_heatmap(im, valfmt="{x:.3f}",
                         size=1.8,  textcolors=["red"],
                         covdata=read_coverages_2Darray,
                         is_annotate = is_plot_annotate
                    )



    #plt.title(title,loc='right', fontsize='5', fontstyle='oblique', fontweight='bold')
    #plt.tight_layout()
    #plt.savefig("heatmap.png")
    #return plt

if __name__ == "__main__":
    renderplot(debug=True)
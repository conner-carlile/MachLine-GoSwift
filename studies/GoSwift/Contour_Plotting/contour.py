""" Produces formatted contour plots for publication in .emf format.

This module is intended to allow a user to produce journal-quality contour
plots. It takes an input file containing many of the formatting options
in one location. The input file is meant to streamline the formatting of
the plot, and allows the inclusion of several data series. It is meant to
be run in conjunction with the bash script contour_save.sh, which allows
the user to save the figure in enhnaced windows metafile (.emf) format
and specify the desired name of the .emf file.


Parameters
------------------------------------------------------------------------
contour_plot_settings.json : input .json file
    Input file containing the following:
{
    "General_Format":{Contains general formatting parameters
        "N" : Integer, number of data points per axis, i.e. N=100 means 100x100
              grid.
        "Font_Size":, float, Size of the default figure font, in pt
        "Figure_Size":, array Size of the figure, in inches in the order
                        [width,height]
        {"Grid": Contains grid parameters
            "is_present":, {0,1} integer, 0 for grid off, any other integer
                           for grid on
            "Which":, string, either 'major' for gridlines on major ticks
                      only, or 'minor' for gridlines on minor ticks
            "Alpha":, float between 0 and 1 specifying transparency of
                      gridlines
            "Line_Color":, string, color of the gridlines, accepts any
                           valid python color string
            "Line_Style":, string, style of the gridlines, accepts any
                           valid python linestyle string
            "Dash_Style":, tuple, (Optional) accepts dash/space
                           specification, i.e. (0, (10.0, 5.0))
            "Line_Width": float, width of the gridlines, in pt.
        }
    },
    "Axis_Format":{ Contains axis formatting parameters
        "Minor_Ticks":, integer, 0 for no minor ticks, any other
                    integer to include minor ticks
        "Tick_Label_Size":, float, fontsize of the tick labels, in pt
        "X_Axis":{ Contains formatting parameters for the x axis.
                   Repeat for Y_Axis
            "Limits":, array, contains lower and upper bounds for the x
                       axis in the order [lower, upper]
            "Tick_Density":, array of values specifying tick values
            "Tick_Label_Format": string, format string specifying the
                                 output format of the tick labels. Accepts
                                 any valid python format string.
        }
    },
    "Data_Series_Format":{ Contains formatting parameters for the data series
        "Series_1":{ formatting parameters for the first series
                     duplicate for any subsequent series, i.e. Series_2, etc.
            "index":, order that the data series is listed in contour_data
            "Name":, string, name of the data series
            "Line_Width":, float, width of the data series, in pt
            "Line_Color":, string, color of the data series, accepts any
                           valid python color string
            "Line_Style": string, style of the gridlines, accepts any
                           valid python linestyle string
            "Dash_Style":, tuple, (Optional) accepts dash/space
                           specification, i.e. (0, (10.0, 5.0))
            "Contour_Density":, float specifying contour density or
                                array of floats specifying contour levels
            "Labels":{ contains formatting parameters for contour labels
                "is_present":, {0,1} integer, 0 for labels off, any other
                               integer for labels on
                "Fontsize":, float, fontsize of the labels, in pt
                "Label_Format":, string, format string specifying the
                                 output format of the tick labels. Accepts
                                 any valid python format string.
            }
        },
    }
}

Returns
------------------------------------------------------------------------
Figure.emf
    Returns a figure in .emf format, with name specified by the user.

Notes
------------------------------------------------------------------------
Data is imported into the conotour.py module via contour_data.py.
contour_data.py contains the function data_import, which returns the arrays
A, B, and C. A and B are 1-D arrays containing the X and Y values,
respectively. C is a 1-D array of 2-D arrays specifying each data series.

Example
------------------------------------------------------------------------
contour.run_plot(input.json)

"""
import sys
import numpy as np
import matplotlib

matplotlib.use('TKAgg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import contour_data as cdat
import json
from collections import OrderedDict
from matplotlib.ticker import FormatStrFormatter
from ast import literal_eval

rc('font', **{'family':'serif', 'serif':['Times']})

def run_plot(input_file):
    """ The run_plot function creates the contour plot and allows the user
    the option to save the file.

    Parameters:
    --------------------------------------------------------------------
    input_file: dinput file containing formatting data

    Returns:
    --------------------------------------------------------------------
    matplotlib figure in .svg format.

    """
    # Read in settings
    format_data = _settings_read_(input_file)

    # Create Plot with General Formatting
    fig, ax = _create_baseplot_(format_data)

    # Format Axis
    ax = _format_axis_(format_data, ax)

    # Plot Data
    ax = plot_data(format_data, ax)

    # Save Plot
    saveflag = _plot_save_(input_file)

    # End Terminal Operation if figure is not to be saved
    if(saveflag is False):
        sys.exit(0)


def _settings_read_(input_file):
    """ The settings_read function reads in the input file and stores
    the values in the data dict.

    Parameters:
    --------------------------------------------------------------------
    input_file: .json file containing the formatting data as described above.

    Returns:
    --------------------------------------------------------------------
    format_data: dictionary containing the formatting data from the input file.

    """
    with open(input_file+'.json') as file:
        format_data = json.load(file, object_pairs_hook=OrderedDict)

    return format_data


def _create_baseplot_(format_data):
    """ The create_baseplot function creates the baseline plot with general
    formatting.

    Parameters:
    --------------------------------------------------------------------
    format_data: dictionary containing the formatting data from the input file.

    Returns:
    --------------------------------------------------------------------
    fig: matplotlib figure object
    ax: matplotlib axis object

    """
    # Create Figure
    fig = plt.figure(figsize=(format_data["General_Format"]["Figure_Size"][0],
                              format_data["General_Format"]["Figure_Size"][1]))
    ax = fig.add_subplot(111)

    # Format Grid
    if(format_data["General_Format"]["Grid"]["is_present"] != 0):
        # Set Line Style
        line_style = format_data["General_Format"]["Grid"]["Line_Style"]
        if(line_style == 'dashed'):
            line_style = literal_eval(format_data["General_Format"]["Grid"]["Dash_Style"])

        ax.grid(which=format_data["General_Format"]["Grid"]["Which"],
                alpha=format_data["General_Format"]["Grid"]["Alpha"],
                color=format_data["General_Format"]["Grid"]["Line_Color"],
                linestyle=line_style,
                linewidth=format_data["General_Format"]["Grid"]["Line_Width"])

    # Set Border
    ax.spines['top'].set_linewidth(1.3)
    ax.spines['right'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)
    ax.spines['left'].set_linewidth(1.3)

    return fig, ax


def _format_axis_(format_data, ax):
    """ The format_axis function applies axis-specific formatting.

    Parameters:
    --------------------------------------------------------------------
    format_data: dictionary containing the formatting data from the input file.
    ax: matplotlib axis object

    Returns:
    --------------------------------------------------------------------
    ax: matplotlib axis object

    """
    # Set X and Y Limits
    ax.set_xlim(format_data["Axis_Format"]["X_Axis"]["Limits"])
    ax.set_ylim(format_data["Axis_Format"]["Y_Axis"]["Limits"])

    # Set Tick Density
    ax.set_xticks(format_data["Axis_Format"]["X_Axis"]["Tick_Density"])
    ax.set_yticks(format_data["Axis_Format"]["Y_Axis"]["Tick_Density"])

    # Toggle Minor Ticks
    if(format_data["Axis_Format"]["Minor_Ticks"] != 0):
        ax.minorticks_on()

    # Format Tick Marks
    ax.tick_params(which='major',
                   labelsize=format_data["Axis_Format"]["Tick_Label_Size"],
                   direction='in',
                   width=0.75,
                   length=3.0,
                   top=True,
                   right=True,
                   pad=7.25)
    ax.tick_params(which='minor',
                   labelsize=format_data["Axis_Format"]["Tick_Label_Size"],
                   direction='in',
                   width=0.25,
                   length=1.75,
                   top=True,
                   right=True,
                   pad=7.25)

    # Format Tick Labels
    ax.xaxis.set_major_formatter(FormatStrFormatter(format_data["Axis_Format"]["X_Axis"]["Tick_Label_Format"]))
    ax.yaxis.set_major_formatter(FormatStrFormatter(format_data["Axis_Format"]["Y_Axis"]["Tick_Label_Format"]))

    return ax


def plot_data(format_data, ax):
    """ The plot_data function plots and formats the data series

    Parameters:
    --------------------------------------------------------------------
    format_data: dictionary containing the formatting data from the input file.
    ax: matplotlib axis object

    Returns:
    --------------------------------------------------------------------
    ax: matplotlib axis object

    """
    # Import data from contour_data.py
    A, B, C = cdat.data_import(format_data["General_Format"]["N"],
                               format_data["Axis_Format"]["X_Axis"]["Limits"][0],
                               format_data["Axis_Format"]["X_Axis"]["Limits"][1],
                               format_data["Axis_Format"]["Y_Axis"]["Limits"][0],
                               format_data["Axis_Format"]["Y_Axis"]["Limits"][1])

    X, Y = np.meshgrid(A, B)

    CS = []
    # Data Series Formatting and Storage
    for series in format_data["Data_Series_Format"].keys():
        CS.append(ax.contour(X, Y, C[format_data["Data_Series_Format"][series]["index"]],
                             format_data["Data_Series_Format"][series]["Contour_Density"],
                             colors=format_data["Data_Series_Format"][series]["Line_Color"],
                             linewidths=format_data["Data_Series_Format"][series]["Line_Width"],
                             linestyles=format_data["Data_Series_Format"][series]["Line_Style"]))

        # Dash Formatting
        if(CS[-1].linestyles == 'dashed'):
            for c in CS[-1].collections:
                c.set_dashes([literal_eval(format_data["Data_Series_Format"][series]["Dash_Style"])])

        # Data Series Labeling
        if(format_data["Data_Series_Format"][series]["Labels"]["is_present"] != 0):
            ax.clabel(CS[-1], inline=1,
                      fontsize=format_data["Data_Series_Format"][series]["Labels"]["Fontsize"],
                      fmt=format_data["Data_Series_Format"][series]["Labels"]["Label_Format"])

    return ax


def _plot_save_(input_file):
    """ The plot_save function allows the user to save the figure as .svg
    It is meant to interface with the contour_save.sh bash script to save
    the figure as .emf and allow the user to specify the desired name of
    the figure.

    """
    plt.show(block=False)
    choiceflag = False
    saveflag = False
    while(choiceflag is False):
        print(' ')
        saving = input("Save Figure? (y/n) ")

        if(saving == 'y'):
            plt.savefig(input_file+'tempfile.svg', bbox_inches='tight', transparent=True)
            choiceflag = True
            saveflag = True
        elif(saving == 'n'):
            saveflag = False
            choiceflag = True
        else:
            print(saving+' is not a valid response. Try again.')

    return saveflag

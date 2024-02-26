"""
Dictionaries for plotly figure layouts.

Layout inspired by
https://github.com/materialsproject/crystaltoolkit/blob/main/crystal_toolkit/components/bandstructure.py

Contains dicts for
    1. general figure layout
    2. axes ((I)COHP / energy)
    3. traces (Spin.up / Spin.down)
    4. legend

"""

layout_dict = {
    "titlefont": {"size": 18},
    "showlegend": True,
    "title_x": 0.1,
    "title_y": 0.9,
    "hovermode": "closest",
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
}

"""
    General layout of Plotly figure.

    :param dict titlefont: Line style dictionary.
    :param bool showlegend: Legend hide or show.
    :param float title_x: x axis title size.
    :param float title_y: y axis title size.
    :param str hovermode: hover behaviour.
    :param str paper_bgcolor: background color.
    :param str plot_bgcolor: plot background color.
"""

cohp_axis_style_dict = {
    "titlefont": {"size": 20},
    "tickfont": {"size": 16},
    "ticks": "inside",
    "tickwidth": 2,
    "showgrid": False,
    "showline": True,
    "zeroline": True,
    "zerolinecolor": "black",
    "zerolinewidth": 3,
    "mirror": True,
    "linewidth": 2,
    "linecolor": "black",
}

"""
    COXX axis style.

    :param dict titlefont: axis title font style.
    :param dict tickfont: axis tick font style.
    :param str ticks: ticks style.
    :param float tickwidth: width of ticks.
    :param bool showgrid: show or hide axis grid.
    :param bool showline: show or hide axis.
    :param bool zeroline: show or hide zero line.
    :param str zerolinecolor: color of zero line.
    :param float zerolinewidth: width of zero line.
    :param bool mirror: mirror axis.
    :param float linewidth: width of axis line.
    :param str linecolor: color of axis line.
"""

energy_axis_style_dict = {
    "titlefont": {"size": 20},
    "tickfont": {"size": 16},
    "ticks": "inside",
    "tickwidth": 2,
    "showgrid": False,
    "showline": True,
    "zeroline": True,
    "zerolinecolor": "rgba(0,0,0,0.5)",
    "zerolinewidth": 3,
    "mirror": True,
    "linewidth": 2,
    "linecolor": "black",
}

"""
    Energy axis style.

    :param dict titlefont: axis title font style.
    :param dict tickfont: axis tick font style.
    :param str ticks: ticks style.
    :param float tickwidth: width of ticks.
    :param bool showgrid: show or hide axis grid.
    :param bool showline: show or hide axis.
    :param bool zeroline: show or hide zero line.
    :param str zerolinecolor: color of zero line.
    :param float zerolinewidth: width of zero line.
    :param bool mirror: mirror axis.
    :param float linewidth: width of axis line.
    :param str linecolor: color of axis line.
"""

spin_up_trace_style_dict = {
    "line": {"width": 3},
    "mode": "lines",
    "visible": "legendonly",
}

"""
    Line style for Spin.up traces.

    :param dict line: Line style dictionary.
    :param str mode: Plotly mode.
    :param str visible: Trace visibility.
"""

spin_down_trace_style_dict = {
    "line": {"width": 3, "dash": "dash"},
    "mode": "lines",
    "visible": "legendonly",
}

"""
    Line style for Spin.down traces.

    :param dict line: Line style dictionary.
    :param str mode: Plotly mode.
    :param str visible: Trace visibility.
"""

legend_style_dict = {
    "bordercolor": "black",
    "borderwidth": 2,
    "traceorder": "normal",
}

"""
    Legend style.

    :param str bordercolor: border color of legend.
    :param int borderwidth: width of border.
    :param str traceorder: order of trace.
"""

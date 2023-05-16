"""
Dictionaries for plotly figure layouts.

Layout inspired by
https://github.com/materialsproject/crystaltoolkit/blob/main/crystal_toolkit/components/bandstructure.py

Contains dicts for
    - general figure layout
    - axes ((I)COHP / energy)
    - traces (Spin.up / Spin.down)
    - legend
"""


layout_dict = {
    "titlefont": {"size": 18},
    "showlegend": True,
    "title_x": 0.1,
    "title_y": 0.9,
    # height=500,
    # width=1000,
    "hovermode": "closest",
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
    # margin=dict(l=60, b=50, t=50, pad=0, r=30),
}

cohp_axis_style_dict = {
    "titlefont": {"size": 20},  # dict(size=20),
    "tickfont": {"size": 16},  # dict(size=16),
    "ticks": "inside",
    "tickwidth": 2,
    "showgrid": False,
    "showline": True,
    "zeroline": True,
    "zerolinecolor": "black",
    "zerolinewidth": 3,
    "mirror": True,
    "linewidth": 2,
    # "linecolor": "rgb(71,71,71)", # maybe replace
    "linecolor": "black",
}

energy_axis_style_dict = {
    "titlefont": {"size": 20},  # dict(size=20),
    "tickfont": {"size": 16},  # dict(size=16),
    "ticks": "inside",
    "tickwidth": 2,
    "showgrid": False,
    "showline": True,
    "zeroline": True,
    "zerolinecolor": "rgba(0,0,0,0.5)",
    "zerolinewidth": 3,
    "mirror": True,
    "linewidth": 2,
    # linecolor="rgb(71,71,71)" # maybe replace
    "linecolor": "black",
}

spin_up_trace_style_dict = {
    "line": {"width": 3},
    "mode": "lines",
    "visible": "legendonly",  # dict(width=3),
}

spin_down_trace_style_dict = {
    "line": {"width": 3, "dash": "dash"},  # dict(width=3, dash="dash"),
    "mode": "lines",
    "visible": "legendonly",
}

legend_style_dict = {
    "bordercolor": "black",
    "borderwidth": 2,
    "traceorder": "normal",
}

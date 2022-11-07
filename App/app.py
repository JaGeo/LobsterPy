import os

import plotly.graph_objs as go

from helper import get_structure_plot, get_dummy_cohp_plot
from dash import Dash, dcc, html, Input, Output



# get path of structure to be visualized
dir = "mp-10143/"
#dir = "mp-510401"
#dir = "mp-2384"
path = os.path.join(os.path.expanduser("~/automationresults"), dir)

# get necessary file paths
path_to_poscar = os.path.join(path, "POSCAR")
path_to_charge = os.path.join(path, "CHARGE.lobster")
path_to_icobilist = os.path.join(path, "ICOBILIST.lobster")
path_to_icooplist = os.path.join(path, "ICOOPLIST.lobster")
path_to_icohplist = os.path.join(path, "ICOHPLIST.lobster")
path_to_cohpcar = os.path.join(path, "COHPCAR.lobster")
path_to_madelung = os.path.join(path, "MadelungEnergies.lobster")

# structure plot as Figure object
#vecs = [
#    [1, 1, 0, 0],
#    [1, 0, 1, 0]
#]
vecs = None
fig = get_structure_plot(
    path_to_poscar,
    path_to_charge,
    path_to_icobilist,
    path_to_icooplist,
    path_to_icohplist,
    path_to_cohpcar,
    path_to_madelung,
    vecs=vecs
)

app = Dash(__name__)

# Figure plot with empty plot for creating COHP plot
cohp_plot = get_dummy_cohp_plot()

app.layout = html.Div([
        # container for structure plot
        html.Div(
            # container for drawing line around section of structure plot
            html.Div(
                # structure plot as Figure object
                dcc.Graph(
                    id="structuregraph",
                    figure=fig,
                    clear_on_unhover=True,
                ),
                id="helper-container",
            ),
            className="container-class",
            id="structuregraph-container",
        ),
        # helper container
        html.Div(
            # dummy for updating camera position of structure plot, not visible on screen
            html.Pre(
                id="helper"
            ),
        ),
        # container for COHP plot
        html.Div(
            # COHP plot as Figure object
            dcc.Graph(
                id="plot",
                figure=cohp_plot,
            ),
            className="container-class",
            id="plot-container",
        ),
        # container for bond property data
        html.Div(
            # bond properties arranged as table
            html.Table([
                html.Tr([
                    html.Td(
                        "Bond Length",
                        className="property-name"
                    ),
                    html.Td(
                        id="bond_length",
                        className="property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOBI",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOBI",
                        className = "property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOOP",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOOP",
                        className="property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOHP",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOHP",
                        className="property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOHP Bonding Perc",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOHP_bonding_perc",
                        className="property-data"
                    )
                ]),
            ]),
            className="container-class",
            id="data-container",
        ),
    ],
    id="container",
)



@app.callback(
    Output("bond_length", "children"),
    Output("ICOBI", "children"),
    Output("ICOOP", "children"),
    Output("ICOHP", "children"),
    Output("ICOHP_bonding_perc", "children"),
    Output("structuregraph", "figure"),
    Output("plot", "figure"),
    Input("structuregraph", "hoverData")
)
def edge_hoverevent(hover_data):
    """
    function for extracting bond properties from an edge that the user hovers over

    :param hover_data: contains which is collected when hovering structure plot
    :return: bond_length: if edge is hovered bond length data of hovered edge, else "-"
             icobi: if edge is hovered icobi data of hovered edge, else "-"
             icoop: if edge is hovered icoop data of hovered edge, else "-"
             icohp: if edge is hovered icohp data of hovered edge, else "-"
             icohp_bonding_perc: if edge is hovered icoph bonding percent of hovered egde, else "-"
             fig: structure plot
             cohp: if edge is hovered COHP plot of hovered edge as Figure object, else empty plot
    """

    # current camera position
    global last_camera_position

    # default layout of COHP plot
    layout = go.Layout(
        showlegend=False,
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        height=440,
        width=600,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        # fig.update_layout(scene_camera=last_camera_position) # maybe prevents "jumping" of the plot !
    )

    cohp = go.Figure(layout=layout)

    # if the hover data contains edge/bond properties, extract them and include them in the data container and
    # COHP plot
    try:
        cohp_data, bond_length, icobi, icoop, icohp, icohp_bonding_perc = hover_data["points"][0]["customdata"]
        bond_length = "{:.5f}".format(bond_length) + " " + chr(8491)
        icobi = "{:.5f}".format(icobi)
        icoop = "{:.5f}".format(icoop) # Einheit "Anteil Elektronen"?
        icohp = "{:.5f}".format(icohp) + " eV"
        icohp_bonding_perc = "{:.1f}".format(icohp_bonding_perc*100) + " %"

        x_range = [min(cohp_data[0]) - 0.25, max(cohp_data[0]) + 0.25]

        cohp.add_trace(
            go.Scatter(
                x=cohp_data[0],
                y=cohp_data[1],
                line=dict(color="red")
            )
        )
        cohp.add_trace(
            go.Scatter(
                x=x_range,
                y=[0, 0],
                mode="lines",
                line=dict(
                    width=1,
                    color="black",
                    dash="dash"
                )
            )
        )

        cohp.update_xaxes(
            visible=True,
            title="-COHP",
            title_font_color="black",
            showline=True,
            linewidth=1.5,
            linecolor="black",
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor="black",
            tickmode="linear",
            ticks="outside",
            dtick=0.25,
            range=x_range,
            tickfont={"color": "black"}
        )
        cohp.update_yaxes(
            visible=True,
            title="E - Efermi [eV]",
            title_font_color="black",
            showline=True,
            linewidth=1.5,
            linecolor="black",
            tickmode="linear",
            ticks="outside",
            dtick=2.0,
            tickfont={"color": "black"}
        )
    # if hover data does not contain edge/bond properties leave data container and COHP plot empty
    except:
        cohp.add_trace(go.Scatter(x=[None], y=[None], line=dict(color="red")))
        bond_length = "-"
        icobi = "-"
        icoop = "-"
        icohp = "-"
        icohp_bonding_perc = "-"

    # set line width of edges to default
    for trace in fig.data:
        if "customdata" in trace:
            trace["line"]["width"] = 2
    # highlight the edge that the user hovers over (make it thicker)
    if hover_data:
        if "customdata" in hover_data["points"][0]:
            trace_index = hover_data["points"][0]["curveNumber"]
            fig.data[trace_index]["line"]["width"] = 5
            fig.data[trace_index]["opacity"] = 1

    # structure plot has to be reloaded everytime an edge is higlighted or de-highlighted, at every reloading the
    # camera position/perspective of the structure plot is set to a standard position, which would lead to "jumping"
    # of the plot, therefore set the camera position to last known position after reloading the plot
    fig.update_layout(scene_camera=last_camera_position)

    return bond_length, icobi, icoop, icohp, icohp_bonding_perc, fig, cohp



@app.callback(Output("helper", "children"), Input("structuregraph", "relayoutData"))
def get_current_camera_position(layout_data):
    """
    function for updating current camera position of structure plot

    :param layout_data: contains data of new perspective on structure plot if it has been moved or zoomed into,
           otherwise None
    :return: None
    """

    # current camera position
    global last_camera_position
    # initializatin of new camera position
    last_camera_position = dict()

    # check if there is new layout data
    if layout_data:
        # if 3D plot has been moved/rotated/etc. update current camera position
        try:
            camera_data = layout_data["scene.camera"]
            up_x = camera_data["up"]["x"]
            up_y = camera_data["up"]["y"]
            up_z = camera_data["up"]["z"]
            center_x = camera_data["center"]["x"]
            center_y = camera_data["center"]["y"]
            center_z = camera_data["center"]["z"]
            eye_x = camera_data["eye"]["x"]
            eye_y = camera_data["eye"]["y"]
            eye_z = camera_data["eye"]["z"]
            last_camera_position = dict(
                up=dict(x=up_x, y=up_y, z=up_z),
                center=dict(x=center_x, y=center_y, z=center_z),
                eye=dict(x=eye_x, y=eye_y, z=eye_z)
            )
        except:
            pass
    return None



if __name__ == '__main__':
    # declare current camera position
    global last_camera_position
    app.run_server(debug=True)

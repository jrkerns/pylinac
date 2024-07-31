from typing import Sequence

from plotly import graph_objects as go


def add_title(fig: go.Figure, title: str):
    """Set the title of a plotly figure at the center of the image"""
    fig.update_layout(title_text=title, title_x=0.5)


def set_axis_range(fig: go.Figure, x: Sequence[float], y: Sequence[float]):
    """Set the axis range of a plotly figure. There's some bug in ploty that won't
    correctly range the Y axis if autorange is already set. This works around that bug
    and in one spot vs trying to remember on each call."""
    fig.update_layout(xaxis_range=x, yaxis_range=y, yaxis_autorange=False)


def fixed_aspect_ratio(fig: go.Figure):
    """Set the aspect ratio of a plotly figure to be fixed."""
    fig.update_layout(
        yaxis_scaleanchor="x",
        yaxis_constrain="domain",
        xaxis_scaleanchor="y",
        xaxis_constrain="domain",
    )


def add_vertical_line(
    fig: go.Figure,
    x: float,
    color="black",
    width=1,
    opacity=1,
    name: str = "",
):
    """Add a vertical line to a plotly figure."""
    # get the current data limits
    # otherwise this can fall outside the image plot
    d = None
    for trace in fig.data:
        if trace.type == "heatmap":
            d = trace
            break
    if d:
        fig.add_scatter(
            x=[x, x],
            y=[0, d.z.shape[0]],
            mode="lines",
            line=dict(color=color, width=width),
            opacity=opacity,
            name=name,
        )


def add_horizontal_line(
    fig: go.Figure, y: float, color="black", width=1, opacity=1, name: str = ""
):
    """Add a horizontal line to a plotly figure."""
    d = None
    for trace in fig.data:
        if trace.type == "heatmap":
            d = trace
            break
    if d:
        fig.add_shape(
            dict(
                type="line",
                x0=0,
                x1=d.z.shape[1],
                y0=y,
                y1=y,
                xref="x",
                yref="y",
                line=dict(color=color, width=width),
                opacity=opacity,
                name=name,
            )
        )
    else:
        # it's a simple plot, just use paper reference
        fig.add_shape(
            dict(
                type="line",
                x0=0,
                x1=1,
                y0=y,
                y1=y,
                xref="paper",
                yref="y",
                line=dict(color=color, width=width),
                opacity=opacity,
                name=name,
            )
        )

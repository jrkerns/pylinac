from plotly import graph_objects as go


def add_title(fig: go.Figure, title: str):
    """Set the title of a plotly figure at the center of the image"""
    fig.update_layout(title_text=title, title_x=0.5)


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
    text: str = "",
    text_kwargs: dict | None = None,
    offset: int = -10,
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
        fig.add_shape(
            dict(
                type="line",
                x0=x,
                x1=x,
                y0=0,
                y1=d.z.shape[0],
                xref="x",
                yref="y",
                line=dict(color=color, width=width),
                opacity=opacity,
                name=name,
            )
        )
        text_kwargs = text_kwargs or {}
        if text:
            fig.add_annotation(
                x=x + offset, y=d.z.shape[0] / 2, text=text, **text_kwargs
            )
    else:
        # it's a simple plot, just use paper reference
        fig.add_shape(
            dict(
                type="line",
                x0=x,
                x1=x,
                y0=0,
                y1=1,
                xref="x",
                yref="paper",
                line=dict(color=color, width=width),
                opacity=opacity,
                name=name,
            )
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

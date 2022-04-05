"""
volcanic_island_model.py

Simulation model of an eroding volcanic island. Team-built by
CU Boulder GEOL5702 group, Spring semester 2022.
"""

import numpy as np
from landlab import imshow_grid, RasterModelGrid
from landlab.components import (
    PriorityFloodFlowRouter,
)


_DEFAULT_TIMING_PARAMS = {
    "run_duration": 2.0,  # duration of run, years
    "timestep_size": 1.0,  # duration of one global time step, years
}

_DEFAULT_GRID_PARAMS = {
    "num_rows": 10,  # number of node rows
    "num_cols": 10,  # number of node columns
    "spacing": 100.0,  # node spacing, m
}

_DEFAULT_CONE_PARAMS = {
    "relief": 3500,  # number of node rows
    "angle": 3,  # number of node columns

}

class VolcanicIslandSimulator:
    def __init__(self, params={}):
        """Initialize VolcanicIslandSimulator"""

        # Run timing parameters
        if "timing" in params:
            t_params = params["timing"]
        else:
            t_params = _DEFAULT_TIMING_PARAMS
        self.dt = t_params["timestep_size"]
        self.remaining_time = t_params["run_duration"]

        # Create and configure grid
        if "grid" in params:
            grid_params = params["grid"]
        else:
            grid_params = _DEFAULT_GRID_PARAMS

        self.grid = RasterModelGrid(
            (grid_params["num_rows"], grid_params["num_cols"]),
            xy_spacing=grid_params["spacing"],
        )

        if "cone" in params:
            cone_params = params["cone"]
        else:
            cone_params = _DEFAULT_CONE_PARAMS
        relief = cone_params['relief']
        angle = cone_params['angle']

        self.grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # Set up initial topography...
        self.topo = self.grid.add_zeros("topographic__elevation", at="node")
        # define initial cone
        self.topo[:] = make_volcano_topography(relief, angle,self.grid.x_of_node, self.grid.y_of_node)

        # ...and soil
        self.soil = self.grid.add_zeros("soil__depth", at="node")
        # For each process/phenomenon, parse parameters, create field(s),
        # instantiate components, and perform other initialization

        #   sea level and/or tectonics
        if "sea_level" in params:
            self.sea_level = params["sea_level"]
        else:
            self.sea_level = 0.0

        #   lithosphere flexure?

        #   hillslope weathering and transport

        #   precipitation

        #   flow routing
        self.flow_router = PriorityFloodFlowRouter(
            self.grid, flow_metric="D8", update_flow_depressions=True
        )

        #   fluvial erosion, transport, deposition

        #   submarine sediment transport

        #   submarine carbonate production

        pass

    def update(self, dt):
        """Update simulation for one global time step of duration dt"""

        # Update tectonics and/or sea level

        # Set boundaries for subaerial processes: all interior submarine nodes
        # flagged as FIXED_VALUE
        under_water = np.logical_and(
            self.topo < self.sea_level,
            self.grid.status_at_node == self.grid.BC_NODE_IS_CORE,
        )
        self.grid.status_at_node[under_water] = self.grid.BC_NODE_IS_FIXED_VALUE

        # Apply weathering and hillslope transport

        # Update precipitation

        # Update flow routing
        self.flow_router.run_one_step()

        # Apply fluvial erosion, transport, and deposition

        # Switch boundaries back to full grid

        # Apply submarine sediment transport

        # Produce marine carbonate

    def run(self):
        """Run simulation from start to finish"""
        while self.remaining_time > 0.0:
            self.update(min(self.dt, self.remaining_time))
            self.remaining_time -= self.dt
        pass
    def plot_elevation(self):
        imshow_grid(self.mg, self.z)
        plt.show()
        pass
def make_volcano_topography(relief, angle, x, y):
    """
    Parameters
    ----------
    relief: float
        height of initial volcano above datum

    angle: float
        average hillslope angle

    x: array
        array of x locations

    y: array
        array of y locations

    Returns
    -------
    topo: array
        array of elevation at x,y

    """

    slope = np.tan(np.pi*angle/180)
    midx = np.mean(x)
    midy = np.mean(y)
    dist = np.sqrt((midx - x)**2+(midy - y)**2)
    topo = relief - slope*dist

    return topo

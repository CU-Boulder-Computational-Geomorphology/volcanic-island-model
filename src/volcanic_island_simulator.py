"""
volcanic_island_model.py

Simulation model of an eroding volcanic island. Team-built by
CU Boulder GEOL5702 group, Spring semester 2022.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from landlab.plot.imshow import imshow_grid
class VolcanicIslandSimulator:

    def __init__(self, params={}):
        """Initialize VolcanicIslandSimulator"""

        # Parse general simulation parameters, create common fields, and
        # initialize
        initial_conditions = params["initial_conditions"]
        # Create grid
        self.relief = initial_conditions["relief"]
        self.angle = initial_conditions["angle"]
        self.spacing = initial_conditions["spacing"]

        self.mg, self.z = make_volcano_grid(self.relief, self.angle, self.spacing)
        # Set up initial conditions

        self.sea_level = self.mg.add_zeros('sea_level__elevation', at='node')
        self.soil_height = self.mg.add_zeros('soil__height', at='node')
        ## set boundary conditions
        self.mg.set_status_at_node_on_edges(
            right=self.mg.BC_NODE_IS_CLOSED,
            top=self.mg.BC_NODE_IS_CLOSED,
            left=self.mg.BC_NODE_IS_CLOSED,
            bottom=self.mg.BC_NODE_IS_CLOSED,
            )
        # For each process/phenomenon, parse parameters, create field(s),
        # instantiate components, and perform other initialization

        #   sea level and/or tectonics

        #   lithosphere flexure?

        #   hillslope weathering and transport

        #   precipitation

        #   flow routing

        #   fluvial erosion, transport, deposition

        #   submarine sediment transport

        #   submarine carbonate production

        pass

    def update(self):
        """Update simulation for one global time step"""

        # Update tectonics and/or sea level

        # Apply weathering and hillslope transport

        # Update precipitation

        # Update flow routing

        # Apply fluvial erosion, transport, and deposition

        # Apply submarine sediment transport

        # Produce marine carbonate

        pass

    def run(self):
        """Run simulation from start to finish"""
        pass
    def plot_elevation(self):
        imshow_grid(self.mg, self.z)
        plt.show()
        return
def make_volcano_grid(relief, angle, spacing):
    """
    Parameters
    ----------
    relief: float
        height of initial volcano above datum

    angle: float
        average hillslope angle

    spacing: float
        distance between nodes

    Returns
    -------
    mg: grid object
        landlab raster model grid object

    z: field object
        model grid elevation field

    """
    from landlab import RasterModelGrid

    slope = np.tan(np.deg2rad(angle))
    mid = 2.0*relief/slope

    num_rows = (mid*2.0+1.0)/spacing
    mg = RasterModelGrid((num_rows, num_rows), xy_spacing = spacing)
    dist = np.sqrt((mid - mg.x_of_node)**2+(mid - mg.y_of_node)**2)
    topo = relief - slope*dist
    z = mg.add_ones('topographic__elevation', at='node')*topo

    return mg, z

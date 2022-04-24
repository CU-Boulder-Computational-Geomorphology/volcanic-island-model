"""
volcanic_island_model.py

Simulation model of an eroding volcanic island. Team-built by
CU Boulder GEOL5702 group, Spring semester 2022.
"""

import numpy as np
from landlab import imshow_grid, RasterModelGrid
from landlab.components import PriorityFloodFlowRouter
from landlab.components import Space
from landlab.components import SimpleSubmarineDiffuser
from landlab.io.netcdf import write_netcdf


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
    "relief": 3500,  # maximum cone elevation, m
    "angle": 3,  # hillslope angle, degrees
    "noise": 1.0,  # amplitude of random noise, m
}

_DEFAULT_FLOW_PARAMS = {
    "surface": "topographic__elevation",
    "flow_metric": "D8",
    "update_flow_depressions": True,
}

_DEFAULT_SPACE_PARAMS = {
    "K_sed": 0.01,  # sediment erodibility
    "K_br": 0.0001,  # bedrock erodibility
    "F_f": 0.0,  # fraction of fines
    "phi": 0.3,  # sediment porosity
    "H_star": 0.1,  # characteristic sediment thickness (roughness height)
    "v_s": 1.0,  # settling velocity
    "m_sp": 0.5,  # area exponent in stream power equation
    "n_sp": 1.0,  # slope exponent in stream power equation
    "sp_crit_sed": 0.0,  # threshold to erode sediment?
    "sp_crit_br": 0.0,  # threshold to erode bedrock?
    "discharge_field": "surface_water__discharge",
    "solver": "basic",
    "dt_min": 0.001,
}

_DEFAULT_MARINE_PARAMS = {
    "sea_level": 0.0,  # water surface height
    "wave_base": 1.0,  # depth to wave base
    "shallow_water_diffusivity": 1.0,  # in m2/yr (this is very small)
}

_DEFAULT_SEA_LEVEL_PARAMS = {"mean": 0, "amplitude": 100, "period": 10000}

_DFAULT_TECTONIC_PARAMS = {"uplift": 0}  # m/yr of relative uplift (+) or subsidence (-)

_DEFAULT_OUTPUT_PARAMS = {
    "output_interval_fraction": 0.5,
    "output_file_basename": "volcanic_island",
}

_OUTPUT_FIELDS = [
    "topographic__elevation",
    "bedrock__elevation",
    "soil__depth",
    "water__unit_flux_in",
    "flow__link_to_receiver_node",
    "drainage_area",
    "flood_status_code",
    "flow__upstream_node_order",
    "flow__receiver_node",
    "surface_water__discharge",
    "topographic__steepest_slope",
    "flow__receiver_proportions",
    "depression_free_elevation",
    "sediment__influx",
    "sediment__outflux",
    "sediment__flux",
    "kd",
    "sediment_deposit__thickness",
    "water__depth",
]


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
        self.update_interval = self.remaining_time / 20.0
        self.next_update = self.remaining_time - self.update_interval
        self.final_time = t_params["run_duration"]

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
        relief = cone_params["relief"]
        angle = cone_params["angle"]
        noise = cone_params["noise"]

        self.grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # Set up initial topography...
        self.topo = self.grid.add_zeros("topographic__elevation", at="node")
        # define initial cone
        self.topo[:] = make_volcano_topography(
            relief, angle, self.grid.x_of_node, self.grid.y_of_node, noise
        )

        # ...and rock and soil
        self.rock = self.grid.add_zeros("bedrock__elevation", at="node")
        self.soil = self.grid.add_zeros("soil__depth", at="node") + 0.01
        self.rock[:] = self.topo - self.soil
        # For each process/phenomenon, parse parameters, create field(s),
        # instantiate components, and perform other initialization

        #   sea level and/or tectonics
        if "sea_level" in params:
            sea_level_params = params["sea_level"]
        else:
            sea_level_params = _DEFAULT_SEA_LEVEL_PARAMS

        self.sea_level_mean = sea_level_params["mean"]
        self.sea_level_amplitude = sea_level_params["amplitude"]
        self.sea_level_period = sea_level_params["period"]

        if "tectonics" in params:
            self.uplift = params["tectonics"]["uplift"]
        else:
            self.uplift = _DFAULT_TECTONIC_PARAMS["uplift"]

        #   lithosphere flexure?

        #   hillslope weathering and transport

        #   precipitation

        #   flow routing
        if "flow" in params:
            flow_params = params["flow"]
        else:
            flow_params = _DEFAULT_FLOW_PARAMS
        self.flow_router = PriorityFloodFlowRouter(self.grid, **flow_params)

        #   fluvial erosion, transport, deposition
        if "space" in params:
            space_params = params["space"]

        else:
            space_params = _DEFAULT_SPACE_PARAMS

        self.space = Space(self.grid, **space_params)
        self.fluvial_sed_influx = self.space.sediment_influx

        #   submarine sediment transport
        if "marine" in params:
            marine_params = params["marine"]
        else:
            marine_params = _DEFAULT_MARINE_PARAMS
        self.ssd = SimpleSubmarineDiffuser(self.grid, **marine_params)

        #   submarine carbonate production

        # Output parameters
        if "output" in params:
            out_params = params["output"]
        else:
            out_params = _DEFAULT_OUTPUT_PARAMS
        self.output_file_basename = out_params["output_file_basename"]
        self.output_interval = (
            self.remaining_time * out_params["output_interval_fraction"]
        )
        self.next_output = self.remaining_time - self.output_interval
        self.output_file_number = 0
        self.num_outfile_zeros = 1 + int(
            np.log10(self.remaining_time / max(self.output_interval, self.dt))
        )
        self.write_output()

    def write_output(self):
        """Write output to netcdf file."""
        filename = (
            self.output_file_basename
            + str(self.output_file_number).zfill(self.num_outfile_zeros)
            + ".nc"
        )
        write_netcdf(filename, self.grid, names=_OUTPUT_FIELDS)
        self.output_file_number += 1

    def update(self, dt):
        """Update simulation for one global time step of duration dt"""

        # Update tectonics and/or sea level
        self.change_sea_level()
        self.apply_tectonics(dt)
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
        self.space.run_one_step()

        # Deposit incoming sediment at underwater nodes
        sea_node_areas = self.grid.area_of_cell[self.grid.cell_at_node[under_water]]
        self.soil[under_water] += (
            dt * self.fluvial_sed_influx[under_water] / sea_node_areas
        )
        self.topo[:] = self.rock + self.soil
        self.fluvial_sed_influx[:] = 0.0

        # Switch boundaries back to full grid
        self.grid.status_at_node[under_water] = self.grid.BC_NODE_IS_CORE

        # Apply submarine sediment transport
        self.ssd.run_one_step(dt)

        # Produce marine carbonate

    def run(self):
        """Run simulation from start to finish"""
        while self.remaining_time > 0.0:
            self.update(min(self.dt, self.remaining_time))
            self.remaining_time -= self.dt
            if self.remaining_time <= self.next_update:
                print("Remaining time", self.remaining_time)
                self.next_update -= self.update_interval
            if self.remaining_time <= self.next_output:
                write_output()
                self.next_output -= self.output_interval

    def plot_elevation(self):
        imshow_grid(self.mg, self.z)
        plt.show()

    def change_sea_level(self):
        """update sea level based on sinuosoidal cycle"""
        time = self.final_time - self.remaining_time
        self.sea_level = self.sea_level_mean + (self.sea_level_amplitude / 2) * np.sin(
            2 * np.pi * time / (self.sea_level_period)
        )
        pass

    def apply_tectonics(self, dt):
        """update base level based on simple releative uplift/subsidence rate"""
        self.topo[:] += self.uplift * dt


def make_volcano_topography(relief, angle, x, y, noise=0.0):
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

    noise: float
        amplitude of random noise added to topography, m

    Returns
    -------
    topo: array
        array of elevation at x,y

    """

    slope = np.tan(np.pi * angle / 180)
    midx = np.mean(x)
    midy = np.mean(y)
    dist = np.sqrt((midx - x) ** 2 + (midy - y) ** 2)
    topo = relief - slope * dist
    topo += noise * np.random.rand(len(topo))

    return topo

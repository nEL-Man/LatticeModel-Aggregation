import os
import pytest
from Simulations.simulation import Simulation

class TestSimulation:
    def test_output_file_creation(self):
        # Cleanup previous test runs
        if os.path.exists("results/simulation_plot.png"):
            os.remove("results/simulation_plot.png")
        
        sim = Simulation()
        sim.run()
        
        assert os.path.exists("results/simulation_plot.png"), "Output file not created"

    def test_data_integrity(self):
        sim = Simulation()
        sim.run()
        
        assert len(sim.t) == len(sim.x), "Time and data arrays length mismatch"
        assert len(sim.t) > 0, "No data generated"
        assert not any(np.isnan(x) for x in sim.x), "NaN values in simulation data"
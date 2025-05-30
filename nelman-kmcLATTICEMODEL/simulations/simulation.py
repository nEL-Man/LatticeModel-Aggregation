import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

class Simulation:
    def __init__(self, m=1.0, c=0.1, k=1.0, x0=1.0, v0=0.0, 
                 t_start=0, t_end=20, dt=0.01):
        """
        Damped harmonic oscillator simulation
        Parameters:
        m : mass
        c : damping coefficient
        k : spring constant
        x0: initial position
        v0: initial velocity
        t_start: start time
        t_end: end time
        dt: time step
        """
        self.m = m
        self.c = c
        self.k = k
        self.x0 = x0
        self.v0 = v0
        self.t_start = t_start
        self.t_end = t_end
        self.dt = dt
        
        # Initialize results storage
        self.t = None
        self.x = None
        self.v = None

    def _ode_func(self, t, y):
        """
        System of ODEs for the damped harmonic oscillator:
        y[0] = position (x)
        y[1] = velocity (v)
        """
        x, v = y
        dxdt = v
        dvdt = (-self.c * v - self.k * x) / self.m
        return [dxdt, dvdt]

    def run(self):
        """Run the simulation and plot results"""
        # Time vector
        t_eval = np.arange(self.t_start, self.t_end, self.dt)
        
        # Initial conditions [x0, v0]
        y0 = [self.x0, self.v0]
        
        # Solve ODE
        sol = solve_ivp(
            self._ode_func, 
            [self.t_start, self.t_end], 
            y0, 
            t_eval=t_eval,
            method='RK45'
        )
        
        # Store results
        self.t = sol.t
        self.x = sol.y[0]
        self.v = sol.y[1]
        
        # Data validation
        if len(self.t) != len(self.x):
            raise ValueError(f"Data length mismatch: t({len(self.t)}) vs x({len(self.x)})")
        
        if np.isnan(self.x).any():
            raise ValueError("NaN values detected in position data")
        
        if np.isnan(self.v).any():
            raise ValueError("NaN values detected in velocity data")
        
        # Create plot
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.x, 'b-', label='Position')
        plt.plot(self.t, self.v, 'r-', label='Velocity')
        plt.title('Damped Harmonic Oscillator')
        plt.xlabel('Time (s)')
        plt.ylabel('State')
        plt.legend()
        plt.grid(True)
        
        # Ensure results directory exists
        os.makedirs("results", exist_ok=True)
        
        # Save plot
        plt.savefig("results/simulation_plot.png")
        plt.close()
        
        # Print success message
        print(f"Simulation completed successfully! Results saved to results/simulation_plot.png")

    def get_results(self):
        """Return simulation results as a tuple (t, x, v)"""
        if self.t is None or self.x is None or self.v is None:
            raise RuntimeError("Simulation not run yet. Call run() first.")
        return self.t, self.x, self.v

# Example usage
if __name__ == "__main__":
    # Create and run simulation
    sim = Simulation(
        m=1.0,
        c=0.1,
        k=1.0,
        x0=1.0,
        v0=0.0,
        t_end=30
    )
    sim.run()
    
    # Retrieve results
    t, x, v = sim.get_results()
    print(f"Simulation generated {len(t)} data points")
import abc
import numpy as np
from typing import Optional

class Generator(metaclass=abc.ABCMeta):
    """
    Abstract class to model a Trajectory Generator. Classes inheriting
    from this class should implement ``generate`` method.

    Parameters
    ----------
    T : float
        Total duration of each Trajectory.
    dim : int, optional
        Dimension of each Trajectory, by default 1.
    N : int, optional
        Number of trajectories, by default 1.
    dt : float, optional
        Time step of the Trajectory, by default 1.0.
    seed : int, optional
        Seed for the random number generator. If None, no seed is set.
        By default None.

    Attributes
    ----------
    T : float
        Total duration of each Trajectory.
    dim : int, optional
        Dimension of each Trajectory, by default 1.
    N : int, optional
        Number of trajectories, by default 1.
    dt : float, optional
        Time step of the Trajectory, by default 1.0.
    n  : int
        Number of samples on each Trajectory.
    rng : np.random.Generator
        Random number generator.
    """

    def __init__(
        self,
        T: float,
        dim: int = 1,
        N: int = 1,
        dt: float = 1.0,
        seed: Optional[int] = None,
    ):
        # Simulation parameters
        self.T = T  # Total time
        self.dim = dim  # Trajectory dimension
        self.N = N  # Number of trajectories
        self.dt = dt  # Time step of the simulation
        self.n = int(T / dt)  # Number of time steps
        self.rng = (
            np.random.default_rng(seed) if seed is not None else np.random.default_rng()
        )

    @abc.abstractmethod
    def generate(self):
        """
        Abstract method that is implemented on inheriting classes.
        It should compute a list of ``N`` Trajectory objects with the
        given parameters using a method specific to the inheriting class.
        """

class OverdampedHarmonicLangevin(Generator):
    """
    Random Walk class from a multidimensional Langevin Equation
    describing a Brownian particle in a harmonic potential.

    Parameters
    ----------
    tau : float
        Total duration of trajectories.
    dim : int, optional
        Trajectories dimension, by default 1.
    N : int, optional
        Number of simulated trajectories, by default 1.
    dt : float, optional
        Time step, by default 1.0.
    gamma : float, optional
        Friction coefficient, by default 1.
    kB : float, optional
        Boltzmann's constant, by default 1.
    T : np.ndarray, optional
        Array containing a time-series of temperatures, by default None.
    kappa : np.ndarray, optional
        Array containing a time-series of stiffness values, by default None.
    r0 : Optional[np.ndarray]
        Initial positions, by default None.
    alpha0 : Optional[float]
        Initial inverse variance, by default None.
    seed : int, optional
        Seed for random number generator, by default None.
    """

    def __init__(
        self,
        tau: float,
        dim: int = 1,
        N: int = 1,
        dt: float = 1.0,
        gamma: float = 1.0,
        kB: float = 1.0,
        T: Optional[np.ndarray] = None,
        kappa: Optional[np.ndarray] = None,
        r0: Optional[np.ndarray] = None,
        v0: Optional[np.ndarray] = None,
        alpha0: Optional[np.ndarray] = None,
        seed: Optional[int] = None,
    ):
        
        super().__init__(tau, dim, N, dt, seed)
        
        # Main ID of generated trajectories
        self.traj_id = "Overdamped Langevin with harmonic potential"
        
        if np.shape(T) == (self.n,):
            self.T = T          # temperature of the bath
        elif np.ndim(T) == 0:
            self.T = np.ones(self.n) * T
        else:
            raise ValueError("T is expected to be a float or an "
                            f"array of shape {(self.n,)}")
            
        if np.shape(kappa) == (self.n, self.N) or np.shape(kappa) == (self.n,):
            self.kappa = kappa  # trap stiffness
        elif np.ndim(kappa) == 0:
            self.kappa = np.ones(self.n) * kappa
        else:
            raise ValueError("kappa is expected to be a float or an "
                            f"array of shape {(self.n,)} or "
                            f"{(self.n, self.N)}")
        
        # Model parameters
        self.kB = kB        # Boltzmann's constant
        self.gamma = gamma  # friction coefficient
        
        # Intial conditions
        self.r0 = r0  # initial position
        self.v0 = v0  # initial velocity
        self.alpha0 = alpha0
        
        # Initialize variables before simulating
        self._set_scaling_params()
        self._set_simulation_vars()
        self._set_init_cond()
        self._set_noise()
        
    # Set characteristic values for each parameter
    def _set_scaling_params(self):
        self.T_scale = self.T[0]
        self.kappa_scale = np.mean(self.kappa[0,...])
        self.t_scale = self.gamma / self.kappa_scale
        self.r_scale = np.sqrt(2 * self.kB * self.T_scale / self.kappa_scale)
        self.v_scale = self.r_scale / self.t_scale
        self.alpha_scale = self.kappa_scale / (2 * self.kB * self.T_scale)
        
        # Simulation parameters and dynamic variables
    def _set_simulation_vars(self):
        # Simulation parameters
        self.dt = self.dt / self.t_scale          # Dimensionless timestep
        self.shape = (self.n, self.dim, self.N)   # Shape of dynamic variables
        self.kappa = self.kappa / self.kappa_scale  # Dimensionless stiffness
        self.T = self.T / self.T_scale            # Dimensionless temperature
        
        # Dynamic variables
        self.t = np.arange(self.n) * self.dt
        self.r = np.empty(self.shape)
        self.v = np.empty(self.shape)
        self.alpha = np.empty(self.kappa.shape)
        
        # Set initial conditions
    def _set_init_cond(self):
        # Initial inverse variance
        if self.alpha0 is None:
            self.alpha[0] = self.alpha_scale
        elif np.shape(self.alpha0) == (self.dim, self.N) or np.ndim(self.alpha0) == 0:
            self.alpha[0] = self.alpha0  # User input
        else:
            raise ValueError(
                "alpha0 is expected to be a float or an "
                f"array of shape {(self.dim, self.N)}"
            )
        self.alpha[0] /= self.alpha_scale
        self.alpha[1] = self.alpha[0] * self.kappa[0] - (self.alpha[0] ** 2) * self.T[0]
        
        # Initial positions
        if self.r0 is None:
            # Sample initial position from equilibrium distribution if none is given
            sigma_x0 = (2 * self.alpha_scale) ** (-1/2)
            self.r[0] = self.rng.normal(scale=sigma_x0, 
                                        size=(self.dim, self.N))   # Default
        elif np.shape(self.r0) == (self.dim, self.N) or np.ndim(self.r0) == 0:
            self.r[0] = self.r0  # User input
        else:
            raise ValueError(
                "r0 is expected to be a float or an "
                f"array of shape {(self.dim, self.N)}"
            )
        self.r[0] /= self.r_scale
        
        # Initial velocity
        if self.v0 is None:
            # Sample initial position from equilibrium distribution if none is given
            self.v[0] = self.rng.normal(scale=self.v_scale, 
                                        size=(self.dim, self.N))   # Default
        elif np.shape(self.v0) == (self.dim, self.N) or np.ndim(self.v0) == 0:
            self.v[0] = self.v0  # User input
        else:
            raise ValueError(
                "v0 is expected to be a float or an "
                f"array of shape {(self.dim, self.N)}"
            )
        self.v[0] /= self.v_scale
        
    # Fill noise array in advance
    def _set_noise(self):
        self.noise = self.rng.normal(size=self.shape)
    
    # Solve the dimensionless Langevin equation
    def _solve(self):
        sqrt_dt = np.sqrt(self.dt)
        
        for i in range(self.n - 1):
            # Position
            self.r[i + 1] = (self.r[i] - self.kappa[i] * self.r[i] * self.dt
                                       + self.noise[i] * sqrt_dt * np.sqrt(self.T[i]))
            
            # Velocity
            self.v[i + 1] = (self.r[i + 1] - self.r[i]) / self.dt
            
            # 4th-order Runge-Kutta for inverse variance
            # f(t_i, alpha[i]) = alpha[i] * kappa[i] - alpha[i]^2 * T[i]
            kappa_halfstep = (self.kappa[i + 1] + self.kappa[i]) / 2
            T_halfstep = (self.T[i + 1] + self.T[i]) / 2

            k1 = self.alpha[i] * self.kappa[i] - (self.alpha[i] ** 2) * self.T[i]

            m2 = self.alpha[i] + self.dt * k1 / 2
            k2 = m2 * kappa_halfstep - (m2 ** 2) * T_halfstep

            m3 = self.alpha[i] + self.dt * k2 / 2
            k3 = m3 * kappa_halfstep - (m3 ** 2) * T_halfstep

            m4 = self.alpha[i] + self.dt * k3
            k4 = m4 * self.kappa[i + 1] - (m4 ** 2) * self.T[i + 1]

            self.alpha[i + 1] = self.alpha[i] + self.dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    # Recover dimensions
    def _set_scale(self):
        self.r *= self.r_scale
        self.v *= self.v_scale
        self.alpha *= self.alpha_scale
        self.t *= self.t_scale
        self.dt *= self.t_scale
        self.kappa *= self.kappa_scale
        self.T *= self.T_scale
    
    # Simulate the process
    def generate(self):
        self._solve()      # solve the Langevin equation
        self._set_scale()  # recover dimensions

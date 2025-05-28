import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import animation
import meshing.tools as mesh_tools
from processing import readwritedatafiles

# Base path for file structure
BASE_PATH = "/Users/paxton/git"
SOURCE_DIR = f"{BASE_PATH}/quail_volcano/src"

# Set working directory to source
os.chdir(SOURCE_DIR)

def animate_conduit_pressure(
    solver_from_i,
    iterations=100,
    d_iterations=1,
    viscosity_index=2,
    wall_friction_index=1,
    p0=4.1,
    y_min=-1000,
    y_max=0,
    max_speed_of_sound=1000,
    max_pressure=12,
    max_velocity=1,
    slip_final=3.5,
    max_slip=20,
    max_tau=0.5,
    max_viscosity=10,
    max_water=1,
    max_density=2.6e3,
    max_fragmentation=600,
    max_crystal=600,
    show_p0_line=False,
):
    """
    Animates various state variables in the conduit over time.

    Parameters
    ----------
    solver_from_i : callable
        Function that returns the solver state for a given iteration.
    iterations : int, optional
        Number of iterations in the simulation. Default is 100.
    d_iterations : int, optional
        Step size for iterations. Default is 1.
    viscosity_index : int, optional
        Index of the viscosity source term. Default is 2.
    wall_friction_index : int, optional
        Index of the wall friction source term. Default is 1.
    p0 : float, optional
        Reference pressure in MPa. Default is 4.1.
    y_min : float, optional
        Minimum y-axis value (height in meters). Default is -1000.
    y_max : float, optional
        Maximum y-axis value (height in meters). Default is 0.
    max_speed_of_sound : float, optional
        Maximum speed of sound (m/s). Default is 1000.
    max_pressure : float, optional
        Maximum pressure (MPa). Default is 12.
    max_velocity : float, optional
        Maximum velocity (m/s). Default is 1.
    slip_final : float, optional
        Final slip value (m). Default is 3.5.
    max_slip : float, optional
        Maximum slip value (m). Default is 20.
    max_tau : float, optional
        Maximum tau value (MPa). Default is 0.5.
    max_viscosity : float, optional
        Maximum viscosity (MPa·s). Default is 10.
    max_water : float, optional
        Maximum water partial density (kg/m³). Default is 1.
    max_density : float, optional
        Maximum density (kg/m³). Default is 2600.
    max_fragmentation : float, optional
        Maximum fragmentation (kg/m³). Default is 600.
    max_crystal : float, optional
        Maximum crystal partial density (kg/m³). Default is 600.

    Returns
    -------
    matplotlib.animation.FuncAnimation
        Animation object for the conduit state variables.
    """
    # Create figure
    fig = plt.figure(figsize=(12, 8))

    # Define subplots with consistent axis limits
    ax1 = fig.add_subplot(251, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_pressure))
    ax2 = fig.add_subplot(252, autoscale_on=False, ylim=(y_max, y_min), xlim=(-3.0, max_velocity))
    ax3 = fig.add_subplot(253, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_speed_of_sound))
    ax4 = fig.add_subplot(254, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_viscosity))
    ax5 = fig.add_subplot(255, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_water))
    ax6 = fig.add_subplot(256, autoscale_on=False, ylim=(y_max, y_min), xlim=(-1, max_slip))
    ax7 = fig.add_subplot(257, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_density))
    ax8 = fig.add_subplot(258, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_crystal))
    ax9 = fig.add_subplot(259, autoscale_on=False, ylim=(y_max, y_min), xlim=(0, max_fragmentation))
    ax10 = fig.add_subplot(2, 5, 10, autoscale_on=False, ylim=(y_max, y_min), xlim=(-0.1, max_tau))

    # Invert y-axis for each subplot
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]:
        ax.invert_yaxis()
        ax.grid(True)

    # Initialize plot lines
    pressure_line, = ax1.plot([], [], color="blue", label="Pressure")
    velocity_line, = ax2.plot([], [], color="red", label="Velocity")
    sound_speed_line, = ax3.plot([], [], color="green", label="Speed of Sound")
    viscosity_line, = ax4.plot([], [], color="orange", label="Viscosity")
    total_water_line, = ax5.plot([], [], color="purple", label="Total Water")
    exsolved_water_line, = ax5.plot([], [], color="blue", label="Exsolved Water")
    new_state_line, = ax6.plot([], [], color="purple", label="Slip")
    rho_line, = ax7.plot([], [], label="Density")
    crystal_line, = ax8.plot([], [], label="Crystals")
    arhoF_line, = ax9.plot([], [], label="Melt")
    tau_line, = ax10.plot([], [], label="Tau Slip")
    tau_viscous_line, = ax10.plot([], [], label="Tau Viscous")
    
    p0_line = None
    if show_p0_line:
         p0_line = ax1.axvline(x=11.4, color='b', linestyle='--', linewidth=1, label="p0 (analytical)")
        
    plug_boundary_line = ax10.axhline(y=-500, color='b', linestyle='--', linewidth=1, label="Plug Boundary")

    # Set labels
    ax1.set_ylabel("Height [m]")
    ax6.set_ylabel("Height [m]")
    ax1.set_xlabel("Pressure [MPa]")
    ax2.set_xlabel("Velocity [m/s]")
    ax3.set_xlabel("Speed of Sound [m/s]")
    ax4.set_xlabel("Viscosity [MPa·s]")
    ax5.set_xlabel("Water Partial Density [kg/m³]")
    ax6.set_xlabel("Slip (slip/ρ_mix) [m]")
    ax7.set_xlabel("Density [kg/m³]")
    ax8.set_xlabel("Crystal Partial Density [kg/m³]")
    ax9.set_xlabel("Fragmentation [kg/m³]")
    ax10.set_xlabel("Tau [MPa]")

    # Add legends
    ax1.legend(loc="lower right")
    ax5.legend(loc="lower right")
    ax10.legend(loc="lower right")

    # Initialize text annotations
    time_text = ax1.text(0.4, 0.95, "", transform=ax1.transAxes)
    pl_text = ax1.text(0.4, 0.9, "", transform=ax1.transAxes)
    velocity_text = ax2.text(0.5, 0.9, "", transform=ax2.transAxes)
    boundary_text = ax10.text(0.4, 0.4, "", transform=ax10.transAxes)

    # Define text templates
    time_template = 'time = %.2f [s]'
    pl_template = 'P_L = %.2f [MPa]'
    velocity_template = 'V = %.2f [m/s]'
    boundary_template = 'Boundary\n = %.2f'

    def init():
        """Initialize the animation with empty data."""
        for line in [
            pressure_line, velocity_line, sound_speed_line, viscosity_line,
            total_water_line, exsolved_water_line, new_state_line, rho_line,
            crystal_line, arhoF_line, tau_line, tau_viscous_line
        ]:
            line.set_data([], [])
        for text in [time_text, pl_text, velocity_text, boundary_text]:
            text.set_text("")
        return (
            pressure_line, velocity_line, sound_speed_line, viscosity_line,
            total_water_line, exsolved_water_line, new_state_line, rho_line,
            crystal_line, plug_boundary_line, arhoF_line, tau_line, tau_viscous_line,
            time_text, pl_text, velocity_text, boundary_text, p0_line
        )

    def animate(i):
        """Update the animation for each frame."""
        solver = solver_from_i(i)
        flag_non_physical = True

        # Get nodal points
        nodal_pts = solver.basis.get_nodes(solver.order)
        x = np.empty((solver.mesh.num_elems,) + nodal_pts.shape)
        for elem_ID in range(solver.mesh.num_elems):
            x[elem_ID] = mesh_tools.ref_to_phys(solver.mesh, elem_ID, nodal_pts)
            
        # Compute variables
        p = solver.physics.compute_additional_variable("Pressure", solver.state_coeffs, flag_non_physical)
        sound_speed = solver.physics.compute_additional_variable("SoundSpeed", solver.state_coeffs, flag_non_physical)
        solver.physics.compute_additional_variable("Temperature", solver.state_coeffs, flag_non_physical)
        viscosity = solver.physics.source_terms[viscosity_index].compute_viscosity(solver.state_coeffs, solver.physics)
        tau_slip = solver.physics.source_terms[wall_friction_index].compute_tau(solver.state_coeffs, x, solver.physics)
        plug_boundary = solver.physics.source_terms[wall_friction_index].compute_plug_boundary(solver.state_coeffs, x, solver.physics)
        tau_viscous = solver.physics.source_terms[viscosity_index].compute_tau(solver.state_coeffs, x, solver.physics)

        # Extract state variables
        arhoA = solver.state_coeffs[:, :, solver.physics.get_state_index("pDensityA")]
        arhoWt = solver.state_coeffs[:, :, solver.physics.get_state_index("pDensityWt")]
        arhoWv = solver.state_coeffs[:, :, solver.physics.get_state_index("pDensityWv")]
        arhoC = solver.state_coeffs[:, :, solver.physics.get_state_index("pDensityC")]
        arhoM = solver.state_coeffs[:, :, solver.physics.get_state_index("pDensityM")]
        arhoF = solver.state_coeffs[:, :, solver.physics.get_state_index("pDensityFm")]
        rhoSlip = solver.state_coeffs[:, :, solver.physics.get_state_index("rhoSlip")]
        momentum = solver.state_coeffs[:, :, solver.physics.get_momentum_slice()]
        rho = np.sum(solver.state_coeffs[:, :, solver.physics.get_mass_slice()], axis=2, keepdims=True)

        # Compute derived quantities
        u = momentum.ravel() / rho.ravel()
        slip = rhoSlip.ravel() / rho.ravel()
        tau = 5e4 - (5e4 - 2e5) * np.exp(-slip / 10)

        # Update plot data (flipped x/y axes)
        pressure_line.set_data(p.ravel() / 1e6, x.ravel())
        velocity_line.set_data(u.ravel(), x.ravel())
        sound_speed_line.set_data(sound_speed.ravel(), x.ravel())
        viscosity_line.set_data(viscosity.ravel() / 1e6, x.ravel())
        total_water_line.set_data(arhoWt.ravel(), x.ravel())
        exsolved_water_line.set_data(arhoWv.ravel(), x.ravel())
        new_state_line.set_data(slip.ravel(), x.ravel())
        rho_line.set_data(rho.ravel(), x.ravel())
        crystal_line.set_data(arhoC.ravel(), x.ravel())
        arhoF_line.set_data(arhoF.ravel(), x.ravel())
        tau_line.set_data(tau_slip.ravel() / 1e6, x.ravel())
        tau_viscous_line.set_data(tau_viscous.ravel() / 1e6, x.ravel())
        plug_boundary_line.set_ydata([plug_boundary])

        # Update text annotations
        time_text.set_text(time_template % (solver.time))
        pl_text.set_text(pl_template % (p.ravel()[0] / 1e6))
        velocity_text.set_text(velocity_template % u.ravel()[0])
        boundary_text.set_text(boundary_template % plug_boundary)

        return (
            pressure_line, velocity_line, sound_speed_line, viscosity_line,
            total_water_line, exsolved_water_line, new_state_line, rho_line,
            crystal_line, plug_boundary_line, arhoF_line, tau_line, tau_viscous_line,
            time_text, pl_text, velocity_text, boundary_text, p0_line
        )

    # Adjust layout to minimize spacing
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    plt.grid(True)

    # Close the figure to prevent display during animation creation
    plt.close()

    # Create and return the animation
    return animation.FuncAnimation(
        fig,
        animate,
        frames=np.linspace(0, iterations, int(iterations / d_iterations), dtype=int),
        init_func=init,
        interval=100,
        blit=False
    )
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import numpy as np
import processing.mdtools as mdtools
from helper_code.helper_functions import get_local_solver_from_index_func

def atmosphere_2d_animate(iterations, d_iter, folder, file_name):
    solver2D_atm1 = get_local_solver_from_index_func(folder, file_name)
    
    # Set up figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Define ranges for colormaps
    pressure_clims = (-1e3, 1e3)  # Pressure range in Pa
    temp_clims = (0, 500)  # Typical atmospheric temperature range in K (adjust as needed)
    
    # Set up colorbars
    sm_pressure = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=pressure_clims[0], vmax=pressure_clims[1]),
        cmap=plt.get_cmap('viridis')
    )
    cb_pressure = fig.colorbar(sm_pressure, ax=ax1)
    cb_pressure.set_label("Pressure (Pa)")
    
    sm_temp = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=temp_clims[0], vmax=temp_clims[1]),
        cmap=plt.get_cmap('inferno')  # Different cmap for temperature
    )
    cb_temp = fig.colorbar(sm_temp, ax=ax2)
    cb_temp.set_label("Temperature (K)")
    
    # Get initial data for pressure baseline
    x1, p1_0 = mdtools.downsample(solver2D_atm1(0), plot_qty="Pressure")
    x1, t1_0 = mdtools.downsample(solver2D_atm1(0), plot_qty="Temperature")
    
    # Animation update function
    def update(frame):
        ax1.cla()  # Clear pressure axis
        ax2.cla()  # Clear temperature axis
        
        # Fetch data for current frame
        x1, p1 = mdtools.downsample(solver2D_atm1(frame), plot_qty="Pressure")
        x1, t1 = mdtools.downsample(solver2D_atm1(frame), plot_qty="Temperature")
        
        # Update plots
        mdtools.plot_mean(x1, p1 - p1_0, pressure_clims, cmap=plt.get_cmap('viridis'), ax=ax1)
        mdtools.plot_mean(x1, t1 - t1_0, temp_clims, cmap=plt.get_cmap('inferno'), ax=ax2)

        # Set titles with time rounded to 1 decimal place
        ax1.set_title(f"$\Delta$ Pressure Field at t= {round(solver2D_atm1(frame).time, 1)}")
        ax2.set_title(f"$\Delta$ Temperature Field at t= {round(solver2D_atm1(frame).time, 1)}")

        return ax1, ax2

    plt.subplots_adjust(bottom=0.2)

    plt.close()
    
    # Create animation
    return animation.FuncAnimation(
        fig,
        update,
        frames=range(0, iterations, d_iter),
        interval=100,
        blit=False
    )
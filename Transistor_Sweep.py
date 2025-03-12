# -*- coding: utf-8 -*-
"""
Transistor Characterization Program for Keysight B2912A SMU
------------------------------------------------------------
This program controls a Keysight B2912A SMU using PyVISA to measure transistor
characteristics (transfer and output) in alternating cycles.

Connection setup:
- SMU1 connected to the gate terminal
- SMU2 connected to the drain terminal
- Source terminal connected to ground

Features:
- Alternating transfer and output characteristic measurements
- User-configurable sweep parameters
- Field-effect mobility calculation
- CSV data storage with timestamps
- Visualization of measurement results and trends
"""

##########################################################################
######################## Import Packages #################################
import pyvisa as visa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import os
from datetime import datetime
from scipy import constants
from scipy import optimize
import signal
import sys
import matplotlib
matplotlib.rcParams['figure.figsize'] = (10, 8)
matplotlib.rcParams['font.size'] = 10

##########################################################################
####################### Signal Handler for Exit ##########################
def signal_handler(signal, frame):
    global interrupted
    interrupted = True
    print("\nMeasurement interrupted by user. Cleaning up...")

signal.signal(signal.SIGINT, signal_handler)
interrupted = False

##########################################################################
####################### User Configuration Parameters ####################
class Config:
    # File paths
    save_directory = "Transistor_Measurements"
    device_name = "TestTransistor"
    
    # Measurement parameters
    cycles = 3                     # Number of transfer-output measurement cycles
    
    # Transfer sweep parameters
    transfer_vds = 5.0             # Fixed Vds during transfer sweep [V]
    transfer_vgs_start = -5.0      # Starting Vgs for transfer sweep [V]
    transfer_vgs_stop = 10.0       # Ending Vgs for transfer sweep [V]
    transfer_vgs_step = 0.2        # Vgs step size [V]
    transfer_backward_sweep = True # Enable/disable backward sweep for hysteresis observation
    
    # Output sweep parameters
    output_vds_start = 0.0         # Starting Vds for output sweep [V]
    output_vds_stop = 10.0         # Ending Vds for output sweep [V]
    output_vds_step = 0.2          # Vds step size [V]
    output_vgs_values = [0, 2, 4, 6, 8, 10]  # Vgs values for output sweeps [V]
    
    # Timing settings
    sweep_delay = 0.01             # Delay between points in sweep loop [s]
    settling_time = 0.05           # Settling time between measurements [s]
    integration_time = 0.1         # Integration time for measurements [s]
    
    # Display settings
    live_plot = True               # Enable/disable live plotting during measurement
    plot_update_interval = 1       # Update live plot every N data points
    
    # Measurement settings
    compliance_current = 0.1       # Compliance current [A]
    
    # Device parameters for mobility calculation
    channel_length = 50e-6         # Channel length [m]
    channel_width = 1000e-6        # Channel width [m]
    oxide_capacitance = 15e-9      # Gate oxide capacitance per unit area [F/m²]
    
    # Mobility calculation parameters
    mobility_fit_low_vgs = 2.0     # Lower Vgs boundary for mobility extraction
    mobility_fit_high_vgs = 8.0    # Upper Vgs boundary for mobility extraction
    
    # SMU VISA Address
    smu_address = 'GPIB0::24::INSTR'  # Update with your SMU address

##########################################################################
####################### Utility Functions ################################
def create_save_directory():
    """Create timestamped directory for saving measurement data"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(Config.save_directory, f"{Config.device_name}_{timestamp}")
    os.makedirs(path, exist_ok=True)
    print(f"Data will be saved to: {path}")
    return path

def generate_sweep_values(start, stop, step):
    """Generate sweep values with given step size"""
    return np.arange(start, stop + step/2, step).tolist()

##########################################################################
####################### SMU Control Class ################################
class SMUControl:
    def __init__(self, visa_address):
        """Initialize connection to SMU"""
        self.rm = visa.ResourceManager()
        try:
            self.smu = self.rm.open_resource(visa_address)
            self.smu.timeout = 250000  # 250 seconds timeout
            self.reset()
            print(f"Connected to: {self.smu.query('*IDN?').strip()}")
        except Exception as e:
            print(f"Error connecting to SMU: {e}")
            sys.exit(1)
            
    def reset(self):
        """Reset the SMU and clear buffers"""
        self.smu.write('*RST')
        self.smu.write('TRIG:CLE')
        self.smu.write('TRACe:CLEAr')
        self.smu.write('*CLS')
        
    def configure_smu(self):
        """Configure basic SMU settings"""
        # Format settings
        self.smu.write(':FORMat:ELEMents:SENSe %s' % ('VOLT,CURRent,TIME'))
        self.smu.write(':FORMat:DATA %s,%d' % ('REAL', 64))
        self.smu.write(':FORMat:BORDer %s' % ('SWAPped'))
        
        # Channel settings
        for ch in [1, 2]:
            # Auto on/off settings
            self.smu.write(f':OUTPut{ch}:ON:AUTO %d' % (1))
            
            # Filter settings
            self.smu.write(f':OUTPut{ch}:FILTer:LPASs:STATe %d' % (1))
            self.smu.write(f':OUTPut{ch}:FILTer:LPASs:AUTO %d' % (1))
            
            # Low terminal setting to ground
            self.smu.write(f':OUTPut{ch}:LOW %s' % ('GROund'))
            
            # Disable high capacitance mode
            self.smu.write(f':OUTPut{ch}:HCAPacitance:STATe %d' % (0))
            
            # Current protection
            self.smu.write(f':SENSe{ch}:CURRent:DC:PROTection:LEVel %G' % (Config.compliance_current))
            
            # Auto range for current
            self.smu.write(f':SENSe{ch}:CURRent:DC:RANGe:AUTO ON')
            self.smu.write(f':SENSe{ch}:CURRent:DC:RANGe:AUTO:LLIMit %s' % ('MIN'))
            self.smu.write(f':SENSe{ch}:CURRent:DC:RANGe:AUTO:MODE %s' % ('RES'))
            self.smu.write(f':SENSe{ch}:CURRent:DC:RANGe:AUTO:THReshold %s' % ('MINimum'))
            
            # Integration time
            self.smu.write(f':SENSe{ch}:CURRent:DC:APERture %G' % (Config.integration_time))

        # Set both channels as voltage sources
        self.smu.write(':SOURce1:FUNCtion:MODE %s' % ('VOLTage'))
        self.smu.write(':SOURce2:FUNCtion:MODE %s' % ('VOLTage'))
        
        # Error check
        print("SMU Configuration Complete. Errors:", self.smu.query(':SYSTem:ERRor:ALL?'))
        self.smu.write('*CLS')

    def setup_transfer_sweep(self, vds):
        """Set up SMU for transfer sweep (fixed Vds, sweep Vgs)"""
        # Set fixed drain voltage
        self.smu.write(':SOURce2:VOLTage:LEVel:IMMediate:AMPLitude %G' % (vds))
        
        # Initialize gate to 0V
        self.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude %G' % (0))
        
        # Error check
        self.smu.query(':SYSTem:ERRor:ALL?')
        self.smu.write('*CLS')

    def setup_output_sweep(self, vgs):
        """Set up SMU for output sweep (fixed Vgs, sweep Vds)"""
        # Set fixed gate voltage
        self.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude %G' % (vgs))
        
        # Initialize drain to 0V
        self.smu.write(':SOURce2:VOLTage:LEVel:IMMediate:AMPLitude %G' % (0))
        
        # Error check
        self.smu.query(':SYSTem:ERRor:ALL?')
        self.smu.write('*CLS')
        
    def close(self):
        """Close connection to SMU"""
        self.reset()
        self.smu.close()
        self.rm.close()
        print("SMU connection closed")

##########################################################################
####################### Live Plot Setup Functions ########################
def setup_live_plots():
    """Set up figures for live data plotting"""
    plt.ion()  # Turn on interactive mode
    
    # Create figure with subplots for transfer curves
    fig_transfer = plt.figure(1)
    fig_transfer.clear()
    ax1 = fig_transfer.add_subplot(211)  # Linear scale
    ax2 = fig_transfer.add_subplot(212)  # Log scale
    
    ax1.set_xlabel('Gate Voltage, Vgs (V)')
    ax1.set_ylabel('Drain Current, Ids (A)')
    ax1.set_title('Transfer Characteristics - Linear Scale')
    ax1.grid(True)
    
    ax2.set_xlabel('Gate Voltage, Vgs (V)')
    ax2.set_ylabel('|Drain Current|, |Ids| (A)')
    ax2.set_title('Transfer Characteristics - Log Scale')
    ax2.set_yscale('log')
    ax2.grid(True)
    
    # Create figure for output curves
    fig_output = plt.figure(2)
    fig_output.clear()
    ax3 = fig_output.add_subplot(111)
    ax3.set_xlabel('Drain Voltage, Vds (V)')
    ax3.set_ylabel('Drain Current, Ids (A)')
    ax3.set_title('Output Characteristics')
    ax3.grid(True)
    
    plt.tight_layout()
    
    return (fig_transfer, ax1, ax2), (fig_output, ax3)

def update_transfer_plot(axes, data, draw=True):
    """Update the transfer plot with current data"""
    ax1, ax2 = axes
    
    ax1.clear()
    ax2.clear()
    
    # Split data into forward and backward sweeps if backward sweep is enabled
    if 'Sweep_Direction' in data and 'Backward' in data['Sweep_Direction']:
        # Forward sweep data
        forward_indices = [i for i, direction in enumerate(data['Sweep_Direction']) if direction == 'Forward']
        forward_vgs = [data['Vgs'][i] for i in forward_indices]
        forward_ids = [data['Ids'][i] for i in forward_indices]
        
        # Backward sweep data
        backward_indices = [i for i, direction in enumerate(data['Sweep_Direction']) if direction == 'Backward']
        backward_vgs = [data['Vgs'][i] for i in backward_indices]
        backward_ids = [data['Ids'][i] for i in backward_indices]
        
        # Plot linear scale
        ax1.plot(forward_vgs, forward_ids, 'b-', marker='o', label='Forward')
        ax1.plot(backward_vgs, backward_ids, 'r-', marker='s', label='Backward')
        ax1.set_xlabel('Gate Voltage, Vgs (V)')
        ax1.set_ylabel('Drain Current, Ids (A)')
        ax1.set_title('Transfer Characteristics - Linear Scale')
        ax1.grid(True)
        ax1.legend()
        
        # Plot log scale
        ax2.plot(forward_vgs, np.abs(forward_ids), 'b-', marker='o', label='Forward')
        ax2.plot(backward_vgs, np.abs(backward_ids), 'r-', marker='s', label='Backward')
        ax2.set_xlabel('Gate Voltage, Vgs (V)')
        ax2.set_ylabel('|Drain Current|, |Ids| (A)')
        ax2.set_title('Transfer Characteristics - Log Scale')
        ax2.set_yscale('log')
        ax2.grid(True)
        ax2.legend()
    else:
        # Just a single forward sweep
        # Plot linear scale
        ax1.plot(data['Vgs'], data['Ids'], 'b-', marker='o')
        ax1.set_xlabel('Gate Voltage, Vgs (V)')
        ax1.set_ylabel('Drain Current, Ids (A)')
        ax1.set_title('Transfer Characteristics - Linear Scale')
        ax1.grid(True)
        
        # Plot log scale
        ax2.plot(data['Vgs'], np.abs(data['Ids']), 'r-', marker='o')
        ax2.set_xlabel('Gate Voltage, Vgs (V)')
        ax2.set_ylabel('|Drain Current|, |Ids| (A)')
        ax2.set_title('Transfer Characteristics - Log Scale')
        ax2.set_yscale('log')
        ax2.grid(True)
    
    if draw:
        plt.figure(1)
        plt.tight_layout()
        plt.draw()
        plt.pause(0.001)

def update_output_plot(ax, data, draw=True):
    """Update the output plot with current data"""
    ax.clear()
    
    # Group data by Vgs values
    vgs_values = sorted(set(data['Vgs']))
    for vgs in vgs_values:
        indices = [i for i, v in enumerate(data['Vgs']) if v == vgs]
        vds = [data['Vds'][i] for i in indices]
        ids = [data['Ids'][i] for i in indices]
        
        if len(vds) > 0:  # Only plot if we have data
            ax.plot(vds, ids, marker='o', label=f'Vgs = {vgs} V')
    
    ax.set_xlabel('Drain Voltage, Vds (V)')
    ax.set_ylabel('Drain Current, Ids (A)')
    ax.set_title('Output Characteristics')
    ax.grid(True)
    ax.legend()
    
    if draw:
        plt.figure(2)
        plt.tight_layout()
        plt.draw()
        plt.pause(0.001)

##########################################################################
####################### Measurement Functions ############################
def perform_transfer_sweep(smu, cycle_num, save_dir, plot_axes=None):
    """Perform a transfer sweep (varying Vgs at fixed Vds) and save data"""
    print(f"\nPerforming Transfer Sweep (Cycle {cycle_num})...")
    
    # Setup fixed drain voltage
    smu.setup_transfer_sweep(Config.transfer_vds)
    
    # Generate Vgs values for sweep
    vgs_values = generate_sweep_values(Config.transfer_vgs_start, Config.transfer_vgs_stop, Config.transfer_vgs_step)
    
    # Data storage
    data = {
        'Vds': [],
        'Vgs': [],
        'Ids': [],
        'Igs': [],
        'Time_since_start': []
    }
    
    # Record start time for this sweep
    start_time = time.time()
    
    # Perform sweep
    for i, Vgs in enumerate(vgs_values):
        if interrupted:
            print("Transfer sweep interrupted!")
            break
            
        # Set gate voltage (SMU1)
        smu.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude %G' % (Vgs))
        
        # Initialize and trigger measurement
        smu.smu.write(':INIT (%s)' % ('@1,2'))
        smu.smu.write(':TRIGger:ALL (%s)' % ('@1,2'))
        
        # Read data
        Igs = smu.smu.query_binary_values(':SENSe1:DATA? %s,%d' % ('CURRent', 1), 'd', False)[0]
        Ids = smu.smu.query_binary_values(':SENSe2:DATA? %s,%d' % ('CURRent', 1), 'd', False)[0]
        
        # Error check
        smu.smu.query(':SYSTem:ERRor:ALL?')
        smu.smu.write('*CLS')
        
        # Record time since start
        elapsed_time = time.time() - start_time
        
        # Store data
        data['Vds'].append(Config.transfer_vds)
        data['Vgs'].append(Vgs)
        data['Ids'].append(Ids)
        data['Igs'].append(Igs)
        data['Time_since_start'].append(elapsed_time)
        
        # Print progress
        print(f"Vgs = {Vgs:.2f} V, Ids = {Ids:.8f} A, Elapsed time = {elapsed_time:.2f} s")
        
        # Update live plot if enabled
        if Config.live_plot and plot_axes and (i % Config.plot_update_interval == 0 or i == len(vgs_values) - 1):
            update_transfer_plot(plot_axes, data)
        
        # Wait for sweep delay
        time.sleep(Config.sweep_delay)
    
    # Save data to CSV
    filename = os.path.join(save_dir, f"Transfer_Sweep_{cycle_num}.csv")
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"Transfer data saved to {filename}")
    
    # Final plot update and save
    if Config.live_plot and plot_axes:
        update_transfer_plot(plot_axes, data)
        plt.figure(1)
        plt.savefig(os.path.join(save_dir, f"Transfer_Plot_{cycle_num}.png"))
    
    # Return data for analysis
    return data

def perform_output_sweep(smu, cycle_num, save_dir, plot_ax=None):
    """Perform output sweeps (varying Vds at different fixed Vgs) and save data"""
    print(f"\nPerforming Output Sweep (Cycle {cycle_num})...")
    
    # Generate Vds values for sweep
    vds_values = generate_sweep_values(Config.output_vds_start, Config.output_vds_stop, Config.output_vds_step)
    
    # Data storage
    data = {
        'Vds': [],
        'Vgs': [],
        'Ids': [],
        'Igs': [],
        'Time_since_start': []
    }
    
    # Record start time for this sweep
    start_time = time.time()
    
    # Perform sweep for each Vgs value
    for vgs_idx, Vgs in enumerate(Config.output_vgs_values):
        # Setup fixed gate voltage
        smu.setup_output_sweep(Vgs)
        
        for i, Vds in enumerate(vds_values):
            if interrupted:
                print("Output sweep interrupted!")
                break
                
            # Set drain voltage (SMU2)
            smu.smu.write(':SOURce2:VOLTage:LEVel:IMMediate:AMPLitude %G' % (Vds))
            
            # Initialize and trigger measurement
            smu.smu.write(':INIT (%s)' % ('@1,2'))
            smu.smu.write(':TRIGger:ALL (%s)' % ('@1,2'))
            
            # Read data
            Igs = smu.smu.query_binary_values(':SENSe1:DATA? %s,%d' % ('CURRent', 1), 'd', False)[0]
            Ids = smu.smu.query_binary_values(':SENSe2:DATA? %s,%d' % ('CURRent', 1), 'd', False)[0]
            
            # Error check
            smu.smu.query(':SYSTem:ERRor:ALL?')
            smu.smu.write('*CLS')
            
            # Record time since start
            elapsed_time = time.time() - start_time
            
            # Store data
            data['Vds'].append(Vds)
            data['Vgs'].append(Vgs)
            data['Ids'].append(Ids)
            data['Igs'].append(Igs)
            data['Time_since_start'].append(elapsed_time)
            
            # Print progress
            print(f"Vgs = {Vgs:.2f} V, Vds = {Vds:.2f} V, Ids = {Ids:.8f} A, Elapsed time = {elapsed_time:.2f} s")
            
            # Update live plot at intervals or when completing a Vgs sweep
            if Config.live_plot and plot_ax and (i % Config.plot_update_interval == 0 or i == len(vds_values) - 1):
                update_output_plot(plot_ax, data)
            
            # Wait for sweep delay
            time.sleep(Config.sweep_delay)
        
        if interrupted:
            break
            
        # Update live plot after each Vgs sweep
        if Config.live_plot and plot_ax:
            update_output_plot(plot_ax, data)
    
    # Save data to CSV
    filename = os.path.join(save_dir, f"Output_Sweep_{cycle_num}.csv")
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"Output data saved to {filename}")
    
    # Final plot update and save
    if Config.live_plot and plot_ax:
        update_output_plot(plot_ax, data)
        plt.figure(2)
        plt.savefig(os.path.join(save_dir, f"Output_Plot_{cycle_num}.png"))
    
    # Return data for analysis
    return data

##########################################################################
####################### Analysis Functions ###############################
def calculate_mobility(transfer_data):
    """Calculate field-effect mobility from transfer curve"""
    # Only use forward sweep data for mobility calculation if we have backward sweep
    if 'Sweep_Direction' in transfer_data and 'Backward' in transfer_data['Sweep_Direction']:
        # Get indices for forward sweep data
        forward_indices = [i for i, direction in enumerate(transfer_data['Sweep_Direction']) if direction == 'Forward']
        vgs = np.array([transfer_data['Vgs'][i] for i in forward_indices])
        ids = np.array([transfer_data['Ids'][i] for i in forward_indices])
        print("Using forward sweep data for mobility calculation")
    else:
        # Convert data to numpy arrays for easier manipulation
        vgs = np.array(transfer_data['Vgs'])
        ids = np.array(transfer_data['Ids'])
    
    # Select data within the fitting range for mobility calculation
    mask = (vgs >= Config.mobility_fit_low_vgs) & (vgs <= Config.mobility_fit_high_vgs)
    vgs_fit = vgs[mask]
    ids_fit = ids[mask]
    
    if len(vgs_fit) < 2:
        print("Not enough data points for mobility calculation!")
        return None, None
    
    # Define linear function for fitting
    def linear_func(x, a, b):
        return a * x + b
    
    # Perform linear fit to extract transconductance
    try:
        params, _ = optimize.curve_fit(linear_func, vgs_fit, ids_fit)
        transconductance = params[0]  # Slope of the Ids vs Vgs curve
        
        # Extract threshold voltage (x-intercept)
        threshold_voltage = -params[1] / params[0]
        
        # Calculate mobility using the standard equation
        mobility = (Config.channel_length / Config.channel_width) * (transconductance / (Config.oxide_capacitance * Config.transfer_vds))
        
        print(f"Calculated mobility: {mobility:.4e} cm²/(V·s)")
        print(f"Threshold voltage: {threshold_voltage:.4f} V")
        
        return mobility, threshold_voltage
        
    except Exception as e:
        print(f"Error during mobility calculation: {e}")
        return None, None

def plot_mobility_trend(mobilities, cycle_nums, save_dir):
    """Plot mobility vs measurement cycle"""
    plt.figure(figsize=(10, 6))
    plt.plot(cycle_nums, mobilities, 'bo-', linewidth=2)
    plt.xlabel('Measurement Cycle')
    plt.ylabel('Field-Effect Mobility (cm²/V·s)')
    plt.title('Field-Effect Mobility vs. Measurement Cycle')
    plt.grid(True)
    plt.xticks(cycle_nums)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "Mobility_vs_Cycle.png"))
    plt.close()
    
    # Save mobility data to CSV
    df = pd.DataFrame({
        'Cycle': cycle_nums,
        'Mobility': mobilities
    })
    df.to_csv(os.path.join(save_dir, "Mobility_vs_Cycle.csv"), index=False)

def plot_peak_current_trend(peak_currents, cycle_nums, save_dir):
    """Plot peak current vs measurement cycle"""
    plt.figure(figsize=(10, 6))
    plt.plot(cycle_nums, peak_currents, 'ro-', linewidth=2)
    plt.xlabel('Measurement Cycle')
    plt.ylabel('Peak Drain Current (A)')
    plt.title('Peak Drain Current vs. Measurement Cycle')
    plt.grid(True)
    plt.xticks(cycle_nums)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "PeakCurrent_vs_Cycle.png"))
    plt.close()
    
    # Save peak current data to CSV
    df = pd.DataFrame({
        'Cycle': cycle_nums,
        'PeakCurrent': peak_currents
    })
    df.to_csv(os.path.join(save_dir, "PeakCurrent_vs_Cycle.csv"), index=False)

##########################################################################
####################### Main Function ####################################
def main():
    """Main function to execute the measurement sequence"""
    # Create save directory
    save_dir = create_save_directory()
    
    # Save configuration to file
    with open(os.path.join(save_dir, "measurement_config.txt"), 'w') as f:
        for attr in dir(Config):
            if not attr.startswith('__'):
                f.write(f"{attr} = {getattr(Config, attr)}\n")
    
    # Initialize SMU
    print("Initializing SMU connection...")
    smu = SMUControl(Config.smu_address)
    
    # Configure SMU
    print("Configuring SMU...")
    smu.configure_smu()
    
    # Storage for mobility and peak current values
    mobilities = []
    peak_currents = []
    cycle_nums = []
    
    # Setup live plots if enabled
    if Config.live_plot:
        transfer_plots, output_plots = setup_live_plots()
        transfer_axes = (transfer_plots[1], transfer_plots[2])
        output_ax = output_plots[1]
    else:
        transfer_axes = None
        output_ax = None
    
    # Start time for the entire measurement
    overall_start_time = time.time()
    
    try:
        # Perform alternating transfer and output sweeps
        for cycle in range(1, Config.cycles + 1):
            if interrupted:
                break
                
            print(f"\n{'='*50}")
            print(f"Starting Measurement Cycle {cycle}/{Config.cycles}")
            print(f"{'='*50}")
            
            # Clear plots for this cycle
            if Config.live_plot:
                transfer_plots[0].clear()
                output_plots[0].clear()
                transfer_axes = (transfer_plots[0].add_subplot(211), transfer_plots[0].add_subplot(212))
                output_ax = output_plots[0].add_subplot(111)
                
                # Set log scale for the second transfer axis
                transfer_axes[1].set_yscale('log')
            
            # Perform transfer sweep with live plotting
            transfer_data = perform_transfer_sweep(smu, cycle, save_dir, transfer_axes)
            
            # Calculate mobility
            mobility, threshold = calculate_mobility(transfer_data)
            if mobility is not None:
                mobilities.append(mobility)
                cycle_nums.append(cycle)
                print(f"Cycle {cycle} Mobility: {mobility:.4e} cm²/(V·s), Vth: {threshold:.2f} V")
            
            # Get peak current
            peak_current = max(np.abs(transfer_data['Ids']))
            peak_currents.append(peak_current)
            print(f"Cycle {cycle} Peak Current: {peak_current:.6e} A")
            
            if interrupted:
                break
            
            # Perform output sweep with live plotting
            output_data = perform_output_sweep(smu, cycle, save_dir, output_ax)
            
            # Print time elapsed
            elapsed = time.time() - overall_start_time
            print(f"\nCycle {cycle} completed. Total elapsed time: {elapsed:.2f} seconds.")
            
        # Plot mobility and peak current trends
        if mobilities and not interrupted:
            # Close previous plots if they exist
            if Config.live_plot:
                plt.close(1)
                plt.close(2)
            
            # Create trend plots
            plot_mobility_trend(mobilities, cycle_nums, save_dir)
            plot_peak_current_trend(peak_currents, cycle_nums, save_dir)
            
            # Show final trend plots
            if Config.live_plot:
                plt.figure(3)  # Mobility trend
                plt.figure(4)  # Peak current trend
                plt.show()
            
    except Exception as e:
        print(f"Error during measurement: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Turn off interactive mode
        if Config.live_plot:
            plt.ioff()
        
        # Ensure SMU is properly closed
        print("\nClosing SMU connection...")
        smu.close()
        
        total_time = time.time() - overall_start_time
        print(f"\nMeasurement completed in {total_time:.2f} seconds.")
        
        if interrupted:
            print("Measurement was interrupted by user.")
        else:
            print("Measurement completed successfully.")

if __name__ == "__main__":
    print("\n" + "="*50)
    print("Transistor Characterization Program")
    print("="*50)
    
    # Create base directory if it doesn't exist
    os.makedirs(Config.save_directory, exist_ok=True)
    
    # Prompt user to start
    input("Press Enter to start the measurement...")
    
    # Run main function
    main()
    
    print("\n" + "="*50)
    print("Program execution completed")
    print("="*50)
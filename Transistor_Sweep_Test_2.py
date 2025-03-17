# -*- coding: utf-8 -*-
"""
Transistor Transfer and Output Characteristics Measurement Script
- Gate connected to Source1 (Channel 1) of the B2912A SMU
- Drain connected to Source2 (Channel 2) of the B2912A SMU
- Source grounded

This script can measure transfer curves, output curves, or cycle between both,
with automatic extraction of mobility, threshold voltage, and other parameters.
Works for both n-type and p-type transistors with automatic detection.
"""

import pyvisa as visa
import matplotlib.pyplot as plt
import numpy as np
import time
import os
from datetime import datetime
import pandas as pd
from scipy import stats
from scipy.signal import savgol_filter
import traceback
import signal

# ========================= MEASUREMENT SETTINGS =========================

# Measurement Mode - Choose one: "transfer", "output", or "transfer-output"
MEASUREMENT_MODE = "transfer"  # Set your desired measurement mode here

# Instrument configuration
SMU_ADDRESS = 'GPIB0::24::INSTR'  # Single SMU with two channels

# Transfer curve settings (Vg sweep with fixed Vd)
TRANSFER_VG_START = -50e-3       # Start gate voltage (V)
TRANSFER_VG_STOP = 1.0           # Stop gate voltage (V)
TRANSFER_VG_POINTS = 102         # Number of gate voltage points
TRANSFER_VD_FIXED = [1.0]        # Fixed drain voltages to test (V)

# Output curve settings (Vd sweep with fixed Vg)
OUTPUT_VD_START = -50e-3         # Start drain voltage (V)
OUTPUT_VD_STOP = 800e-3          # Stop drain voltage (V)
OUTPUT_VD_POINTS = 101           # Number of drain voltage points
OUTPUT_VG_FIXED = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # Fixed gate voltages to test (V)

# General settings
NPLC = 3.0                       # Integration time in power line cycles (0.1 to 10)
CURRENT_COMPLIANCE = 0.1         # Current compliance (A)
NUM_CYCLES = 5                   # Number of cycles
DELAY_BETWEEN_SWEEPS = 1.0       # Delay between sweeps (seconds)

# Device parameters
CHANNEL_LENGTH = 20e-6           # Channel length in meters
CHANNEL_WIDTH = 1000e-6          # Channel width in meters
DEVICE_AREA = CHANNEL_WIDTH * CHANNEL_LENGTH  # Device area in m²
# Oxide capacitance per unit area (F/m²) - adjust this based on your device
C_OX = 1.15e-3                  

# Mobility calculation parameters
SMOOTHING_WINDOW = 15            # Window size for Savitzky-Golay filter (must be odd)
SMOOTHING_POLY_ORDER = 3         # Polynomial order for Savitzky-Golay filter

# Setup interrupt handler
interrupted = False

def signal_handler(sig, frame):
    """Handle keyboard interrupt (Ctrl+C)"""
    global interrupted
    print("\nInterrupt detected! Finishing current measurement and exiting...")
    interrupted = True

# Register the signal handler for keyboard interrupt (CTRL+C)
signal.signal(signal.SIGINT, signal_handler)

# Output directory
SAVE_DIR = "Transistor_Measurements"
if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
BASE_FILENAME = f"{SAVE_DIR}/Transistor_{timestamp}"

# ========================= HELPER FUNCTIONS =========================

def setup_smu(smu):
    """Initialize SMU with basic settings"""
    smu.write('*RST')                       # Reset instrument
    smu.write('TRIG:CLE')                   # Clear trigger buffer
    smu.write('TRACe:CLEAr')                # Clear data buffer
    smu.write('*CLS')                       # Clear status registers
    
    # Setup data format
    smu.write(':FORMat:ELEMents:SENSe %s' % ('VOLT,CURRent,TIME'))
    smu.write(':FORMat:DATA %s,%d' % ('REAL', 64))
    smu.write(':FORMat:BORDer %s' % ('SWAPped'))
    
    return smu.query('*IDN?').strip()

def calculate_mobility_and_vth(vg_data, id_data, vd_fixed):
    """
    Calculate mobility and threshold voltage from transfer curve data.
    Works for both p-type and n-type transistors by automatically detecting
    device type from the data.
    
    Returns:
    - mobility: Field-effect mobility (cm²/Vs)
    - vth: Threshold voltage (V)
    - r_squared: R² of the linear fit
    - selected_range: (start_idx, end_idx) of the selected linear region
    - device_type: 'n-type' or 'p-type'
    """
    # Convert to numpy arrays
    vg = np.array(vg_data)
    id = np.array(id_data)
    
    # Determine device type by analyzing current behavior with gate voltage
    # For n-type, current increases as Vg increases
    # For p-type, current increases as Vg decreases (more negative)
    
    # Get indices of min and max gate voltages
    min_vg_idx = np.argmin(vg)
    max_vg_idx = np.argmax(vg)
    
    # If current at max Vg > current at min Vg, it's likely n-type
    # Otherwise, it's likely p-type
    if np.abs(id[max_vg_idx]) > np.abs(id[min_vg_idx]):
        device_type = 'n-type'
        # For n-type, conventional operation is in positive voltages
        working_id = id
    else:
        device_type = 'p-type'
        # For p-type, use absolute values and remember sign for later
        working_id = np.abs(id)
    
    print(f"Detected device type: {device_type}")
    
    # Calculate square root of drain current
    sqrt_id = np.sqrt(np.abs(working_id))
    
    # Smooth the data using Savitzky-Golay filter
    if len(sqrt_id) > SMOOTHING_WINDOW:
        sqrt_id_smooth = savgol_filter(sqrt_id, SMOOTHING_WINDOW, SMOOTHING_POLY_ORDER)
    else:
        sqrt_id_smooth = sqrt_id
    
    # Calculate numerical derivative (transconductance)
    d_sqrt_id_dvg = np.gradient(sqrt_id_smooth, vg)
    
    # Find the region of maximum slope
    if len(d_sqrt_id_dvg) > SMOOTHING_WINDOW:
        d_sqrt_id_dvg_smooth = savgol_filter(d_sqrt_id_dvg, SMOOTHING_WINDOW, SMOOTHING_POLY_ORDER)
    else:
        d_sqrt_id_dvg_smooth = d_sqrt_id_dvg
    
    # Find the region with maximum derivative (slope)
    # For n-type, find maximum (positive slope)
    # For p-type, find minimum (negative slope, steepest part of the curve)
    if device_type == 'n-type':
        max_deriv_idx = np.argmax(d_sqrt_id_dvg_smooth)
    else:  # p-type
        max_deriv_idx = np.argmin(d_sqrt_id_dvg_smooth)
    
    # Define a range around the maximum derivative
    # Use 30% of the data range centered around max_deriv_idx
    range_size = max(int(len(vg) * 0.3), 10)  # at least 10 points
    start_idx = max(0, max_deriv_idx - range_size//2)
    end_idx = min(len(vg), max_deriv_idx + range_size//2)
    
    # Linear fit in this range
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        vg[start_idx:end_idx], 
        sqrt_id_smooth[start_idx:end_idx]
    )
    
    # Calculate threshold voltage (x-intercept)
    vth = -intercept / slope
    
    # Calculate mobility using Equation (2) from the provided image
    # μ = (2*L / (W*Ci)) * (∂√Id / ∂Vgs)²
    # We use absolute value of slope squared so it works for both n and p-type
    mobility_factor = 2 * CHANNEL_LENGTH / (CHANNEL_WIDTH * C_OX)
    mobility = mobility_factor * (slope ** 2)
    
    # Convert mobility to cm²/Vs
    mobility_cm2_vs = mobility * 10000  # m²/Vs to cm²/Vs
    
    return mobility_cm2_vs, vth, r_value**2, (start_idx, end_idx), device_type

def perform_transfer_sweep(smu, vd_fixed_list, cycle_num):
    """
    Perform transfer sweep (Vg sweep with fixed Vd)
    Returns: List of DataFrames with measurement results for each Vd value
    """
    print(f"\nCycle {cycle_num}: Transfer Curve Measurement")
    
    # Reset SMU to default state
    smu.write('*RST')
    smu.write('*CLS')
    
    # Calculate total points
    total_points = TRANSFER_VG_POINTS * len(vd_fixed_list)
    
    # Setup gate voltage sweep (Source1 - primary sweep)
    smu.write(':SOURce1:FUNCtion:MODE %s' % ('VOLTage'))
    smu.write(':SOURce1:VOLTage:MODE %s' % ('SWEep'))
    smu.write(':SOURce1:VOLTage:STARt %G' % (TRANSFER_VG_START))
    smu.write(':SOURce1:VOLTage:STOP %G' % (TRANSFER_VG_STOP))
    smu.write(':SOURce1:VOLTage:POINts %d' % (TRANSFER_VG_POINTS))
    
    # Create a list of drain voltages for each gate voltage point
    vd_list_str = ""
    for vd in vd_fixed_list:
        # Repeat each Vd value for each Vg point
        vd_list_str += ','.join([str(vd)] * TRANSFER_VG_POINTS)
        if vd != vd_fixed_list[-1]:
            vd_list_str += ","
    
    # Setup drain voltage list (Source2 - secondary sweep)
    smu.write(':SOURce2:FUNCtion:MODE %s' % ('VOLTage'))
    smu.write(':SOURce2:VOLTage:MODE %s' % ('LIST'))
    smu.write(':SOURce2:LIST:VOLTage %s' % (vd_list_str))
    
    # Setup current measurement for both channels
    smu.write(':SENSe1:FUNCtion %s' % ('"CURRent"'))
    smu.write(':SENSe1:CURRent:NPLCycles %G' % (NPLC))
    smu.write(':SENSe1:CURRent:PROTection %G' % (CURRENT_COMPLIANCE))
    smu.write(':SENSe1:CURRent:RANGe:AUTO ON')
    
    smu.write(':SENSe2:FUNCtion %s' % ('"CURRent"'))
    smu.write(':SENSe2:CURRent:NPLCycles %G' % (NPLC))
    smu.write(':SENSe2:CURRent:PROTection %G' % (CURRENT_COMPLIANCE))
    smu.write(':SENSe2:CURRent:RANGe:AUTO ON')
    
    # Setup triggering
    smu.write(':TRIGger1:SOURce %s' % ('AINT'))
    smu.write(':TRIGger1:COUNt %d' % (total_points))
    smu.write(':TRIGger2:SOURce %s' % ('AINT'))
    smu.write(':TRIGger2:COUNt %d' % (total_points))
    
    # Turn on outputs
    smu.write(':OUTPut1 ON')
    smu.write(':OUTPut2 ON')
    
    # Start the sweep and measurement
    smu.write(':INITiate (@1,2)')
    
    # Wait for sweep to complete - adjust based on measurement time
    wait_time = NPLC * total_points * 0.02 + 2.0
    print(f"Waiting {wait_time:.1f} seconds for sweep to complete...")
    start_time = time.time()
    while time.time() - start_time < wait_time:
        time.sleep(0.5)  # Check for interrupts more frequently
        if interrupted:
            print("Interrupting current sweep!")
            break
    
    # Fetch measurements
    ig_data = smu.query_ascii_values(':FETCh:ARRay:CURRent? (@1)', container=np.array)
    id_data = smu.query_ascii_values(':FETCh:ARRay:CURRent? (@2)', container=np.array)
    
    # Fetch voltages (optional, but good for verification)
    vg_measured = smu.query_ascii_values(':FETCh:ARRay:VOLTage? (@1)', container=np.array)
    vd_measured = smu.query_ascii_values(':FETCh:ARRay:VOLTage? (@2)', container=np.array)
    
    # Turn off outputs
    smu.write(':OUTPut1 OFF')
    smu.write(':OUTPut2 OFF')
    
    # Split the data for each Vd value
    results = []
    vg_data = np.linspace(TRANSFER_VG_START, TRANSFER_VG_STOP, TRANSFER_VG_POINTS)
    
    # Create a unified dataframe for all VDs
    all_transfer_data = pd.DataFrame()
    
    for i, vd in enumerate(vd_fixed_list):
        start_idx = i * TRANSFER_VG_POINTS
        end_idx = (i + 1) * TRANSFER_VG_POINTS
        
        # Extract data for this Vd value
        vg_subset = vg_measured[start_idx:end_idx]
        id_subset = id_data[start_idx:end_idx]
        ig_subset = ig_data[start_idx:end_idx]
        
        # Calculate transconductance: gm = dId/dVg
        gm = np.gradient(id_subset, vg_subset)
        
        # Calculate mobility and threshold voltage with improved function
        mobility, vth, r_squared, selected_range, device_type = calculate_mobility_and_vth(vg_subset, id_subset, vd)
        
        # Calculate sqrt(Id) for plotting
        sqrt_id = np.sqrt(np.abs(id_subset))
        
        # Find maximum current
        max_current = np.max(np.abs(id_subset))
        
        # Calculate on/off ratio
        on_off_ratio = np.max(np.abs(id_subset)) / np.min(np.abs(id_subset[id_subset != 0]))
        
        # Calculate subthreshold swing in the subthreshold region
        # For n-type: look for region where current rises rapidly near threshold
        # For p-type: look for region where current rises rapidly as Vg becomes more negative
        
        # Define a general subthreshold region search
        # We'll define a wide current range to find the steepest part of the log(Id) vs Vg curve
        low_current_threshold = max_current * 1e-4  # 10⁻⁴ of max current
        high_current_threshold = max_current * 1e-1  # 10⁻¹ of max current
        
        subthreshold_mask = (np.abs(id_subset) > low_current_threshold) & (np.abs(id_subset) < high_current_threshold)
        
        if np.sum(subthreshold_mask) >= 3:  # Need at least 3 points for a good fit
            log_id_subth = np.log10(np.abs(id_subset[subthreshold_mask]))
            vg_subth = vg_subset[subthreshold_mask]
            slope_subth, _, _, _, _ = stats.linregress(vg_subth, log_id_subth)
            # Subthreshold swing in mV/decade
            subthreshold_swing = np.abs(1 / slope_subth) * 1000  # Convert to mV/decade, use abs for both types
        else:
            subthreshold_swing = np.nan
        
        # Create DataFrame
        data = pd.DataFrame({
            'Vg': vg_subset,
            'Vd': vd_measured[start_idx:end_idx],
            'Id': id_subset,
            'Ig': ig_subset,
            'Id_per_W': [i/CHANNEL_WIDTH for i in id_subset],  # Current per width (A/m)
            'sqrt_Id': sqrt_id,
            'gm': gm.tolist()
        })
        
        # Add to the unified dataframe
        # First, add a column to identify the Vd value
        data['Vd_fixed'] = vd
        all_transfer_data = pd.concat([all_transfer_data, data], ignore_index=True)
        
        # Print calculated parameters
        print(f"\nVd = {vd}V measurement results:")
        print(f"  Device Type: {device_type}")
        print(f"  Mobility: {mobility:.4f} cm²/Vs")
        print(f"  Threshold Voltage: {vth:.4f} V")
        print(f"  On/Off Ratio: {on_off_ratio:.2e}")
        print(f"  Subthreshold Swing: {subthreshold_swing:.2f} mV/decade")
        print(f"  R-squared of linear fit: {r_squared:.4f}")
        print(f"  Maximum current: {max_current:.10f} A")
        
        # Plot results
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Linear plot
        axes[0, 0].plot(vg_subset, np.array(id_subset))
        axes[0, 0].set_xlabel('Gate Voltage (V)')
        axes[0, 0].set_ylabel('Drain Current (A)')
        axes[0, 0].set_title(f'Transfer Curve (Vd = {vd}V)')
        axes[0, 0].grid(True)
        
        # Semilog plot
        axes[0, 1].semilogy(vg_subset, np.abs(np.array(id_subset)))
        axes[0, 1].set_xlabel('Gate Voltage (V)')
        axes[0, 1].set_ylabel('|Drain Current| (A)')
        axes[0, 1].set_title(f'Transfer Curve (log scale) - {device_type}')
        axes[0, 1].grid(True)
        
        # Sqrt(Id) plot with linear fit
        axes[1, 0].plot(vg_subset, sqrt_id, 'o', markersize=2, label='Data')
        start_range, end_range = selected_range
        
        # Fit line
        fit_vg = np.linspace(vg_subset[start_range], vg_subset[end_range-1], 100)
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            vg_subset[start_range:end_range], 
            sqrt_id[start_range:end_range]
        )
        fit_sqrt_id = slope * fit_vg + intercept
        
        # Highlight selected range
        axes[1, 0].plot(vg_subset[start_range:end_range], sqrt_id[start_range:end_range], 'ro', markersize=3)
        
        # Plot fit line
        axes[1, 0].plot(fit_vg, fit_sqrt_id, 'g-', linewidth=2, 
                   label=f'Fit: y = {slope:.4f}x + {intercept:.4f}')
        
        # Mark threshold voltage
        axes[1, 0].axvline(x=vth, color='r', linestyle='--', 
                      label=f'Vth = {vth:.4f} V')
        
        axes[1, 0].set_xlabel('Gate Voltage (V)')
        axes[1, 0].set_ylabel('√|Id| (√A)')
        axes[1, 0].set_title(f'Mobility: {mobility:.4f} cm²/Vs, Vth: {vth:.4f} V')
        axes[1, 0].legend()
        axes[1, 0].grid(True)
        
        # Gate leakage
        axes[1, 1].semilogy(vg_subset, np.abs(np.array(ig_subset))*1e9)
        axes[1, 1].set_xlabel('Gate Voltage (V)')
        axes[1, 1].set_ylabel('|Gate Current| (nA)')
        axes[1, 1].set_title('Gate Leakage Current')
        axes[1, 1].grid(True)
        
        # Add subthreshold swing to the semilog plot if available
        if not np.isnan(subthreshold_swing) and np.sum(subthreshold_mask) >= 3:
            axes[0, 1].plot(vg_subset[subthreshold_mask], np.abs(id_subset[subthreshold_mask]), 'ro', markersize=4)
            axes[0, 1].text(0.05, 0.05, f'SS = {subthreshold_swing:.1f} mV/dec\nOn/Off = {on_off_ratio:.1e}', 
                       transform=axes[0, 1].transAxes, bbox=dict(facecolor='white', alpha=0.7))
        
        fig.tight_layout()
        
        # Save individual plot only if explicitly requested
        if False:  # Change to True if you want these saved
            filename = f"{BASE_FILENAME}_Transfer_Vd{vd:.1f}V_Cycle{cycle_num}"
            fig.savefig(f"{filename}.png")
        plt.close(fig)
        
        # Add to results
        results.append({
            'vd_fixed': vd,
            'data': data,
            'mobility': mobility,
            'vth': vth,
            'r_squared': r_squared,
            'max_current': max_current,
            'on_off_ratio': on_off_ratio,
            'subthreshold_swing': subthreshold_swing,
            'device_type': device_type,
            'cycle': cycle_num
        })
        
        # Check for interrupt
        if interrupted:
            print(f"Interrupting transfer sweep at Vd = {vd}V")
            break
    
    # Save all data to a single file per cycle with only required columns
    all_transfer_data_slim = all_transfer_data[['Vg', 'Vd', 'Id', 'Ig', 'Vd_fixed']]
    all_transfer_data_slim.to_csv(f"{BASE_FILENAME}_Transfer_AllVd_Cycle{cycle_num}.csv", index=False)
    print(f"All transfer data saved to {BASE_FILENAME}_Transfer_AllVd_Cycle{cycle_num}.csv")
    
    return results

def perform_output_sweep(smu, vg_fixed_list, cycle_num):
    """
    Perform output sweep (Vd sweep with fixed Vg)
    Returns: List of DataFrames with measurement results for each Vg value
    """
    print(f"\nCycle {cycle_num}: Output Curve Measurement")
    
    # Reset SMU to default state
    smu.write('*RST')
    smu.write('*CLS')
    
    # Calculate total points
    total_points = OUTPUT_VD_POINTS * len(vg_fixed_list)
    
    # Setup drain voltage sweep (Source2 - primary sweep)
    smu.write(':SOURce2:FUNCtion:MODE %s' % ('VOLTage'))
    smu.write(':SOURce2:VOLTage:MODE %s' % ('SWEep'))
    smu.write(':SOURce2:VOLTage:STARt %G' % (OUTPUT_VD_START))
    smu.write(':SOURce2:VOLTage:STOP %G' % (OUTPUT_VD_STOP))
    smu.write(':SOURce2:VOLTage:POINts %d' % (OUTPUT_VD_POINTS))
    
    # Create a list of gate voltages for each drain voltage point
    vg_list_str = ""
    for vg in vg_fixed_list:
        # Repeat each Vg value for each Vd point
        vg_list_str += ','.join([str(vg)] * OUTPUT_VD_POINTS)
        if vg != vg_fixed_list[-1]:
            vg_list_str += ","
    
    # Setup gate voltage list (Source1 - secondary sweep)
    smu.write(':SOURce1:FUNCtion:MODE %s' % ('VOLTage'))
    smu.write(':SOURce1:VOLTage:MODE %s' % ('LIST'))
    smu.write(':SOURce1:LIST:VOLTage %s' % (vg_list_str))
    
    # Setup current measurement for both channels
    smu.write(':SENSe1:FUNCtion %s' % ('"CURRent"'))
    smu.write(':SENSe1:CURRent:NPLCycles %G' % (NPLC))
    smu.write(':SENSe1:CURRent:PROTection %G' % (CURRENT_COMPLIANCE))
    smu.write(':SENSe1:CURRent:RANGe:AUTO ON')
    
    smu.write(':SENSe2:FUNCtion %s' % ('"CURRent"'))
    smu.write(':SENSe2:CURRent:NPLCycles %G' % (NPLC))
    smu.write(':SENSe2:CURRent:PROTection %G' % (CURRENT_COMPLIANCE))
    smu.write(':SENSe2:CURRent:RANGe:AUTO ON')
    
    # Setup triggering
    smu.write(':TRIGger1:SOURce %s' % ('AINT'))
    smu.write(':TRIGger1:COUNt %d' % (total_points))
    smu.write(':TRIGger2:SOURce %s' % ('AINT'))
    smu.write(':TRIGger2:COUNt %d' % (total_points))
    
    # Turn on outputs
    smu.write(':OUTPut1 ON')
    smu.write(':OUTPut2 ON')
    
    # Start the sweep and measurement
    smu.write(':INITiate (@1,2)')
    
    # Wait for sweep to complete - adjust based on measurement time
    wait_time = NPLC * total_points * 0.02 + 2.0
    print(f"Waiting {wait_time:.1f} seconds for sweep to complete...")
    start_time = time.time()
    while time.time() - start_time < wait_time:
        time.sleep(0.5)  # Check for interrupts more frequently
        if interrupted:
            print("Interrupting current sweep!")
            break
    
    # Fetch measurements
    ig_data = smu.query_ascii_values(':FETCh:ARRay:CURRent? (@1)', container=np.array)
    id_data = smu.query_ascii_values(':FETCh:ARRay:CURRent? (@2)', container=np.array)
    
    # Fetch voltages (optional, but good for verification)
    vg_measured = smu.query_ascii_values(':FETCh:ARRay:VOLTage? (@1)', container=np.array)
    vd_measured = smu.query_ascii_values(':FETCh:ARRay:VOLTage? (@2)', container=np.array)
    
    # Turn off outputs
    smu.write(':OUTPut1 OFF')
    smu.write(':OUTPut2 OFF')
    
    # Split the data for each Vg value
    results = []
    vd_data = np.linspace(OUTPUT_VD_START, OUTPUT_VD_STOP, OUTPUT_VD_POINTS)
    
    # Create a unified dataframe for all VGs
    all_output_data = pd.DataFrame()
    
    # Plot a single figure with all Vg curves
    plt.figure(figsize=(10, 8))
    
    for i, vg in enumerate(vg_fixed_list):
        start_idx = i * OUTPUT_VD_POINTS
        end_idx = (i + 1) * OUTPUT_VD_POINTS
        
        # Extract data for this Vg value
        vd_subset = vd_measured[start_idx:end_idx]
        id_subset = id_data[start_idx:end_idx]
        ig_subset = ig_data[start_idx:end_idx]
        
        # Calculate output conductance: gd = dId/dVd
        gd = np.gradient(id_subset, vd_subset)
        
        # Create DataFrame for this Vg
        data = pd.DataFrame({
            'Vd': vd_subset,
            'Vg': vg_measured[start_idx:end_idx],
            'Id': id_subset,
            'Ig': ig_subset,
            'Id_per_W': [i/CHANNEL_WIDTH for i in id_subset],  # Current per width (A/m)
            'gd': gd.tolist()
        })
        
        # Add to the unified dataframe
        # First, add a column to identify the Vg value
        data['Vg_fixed'] = vg
        all_output_data = pd.concat([all_output_data, data], ignore_index=True)
        
        # Print some data points for verification
        print(f"\nVg = {vg}V measurement:")
        for j in range(0, len(vd_subset), max(1, len(vd_subset)//5)):
            print(f"  Vd = {vd_subset[j]:.2f}V, Id = {id_subset[j]:.8e}A, Ig = {ig_subset[j]:.8e}A")
        
        # Add to combined output plot
        plt.plot(vd_subset, np.array(id_subset), label=f'Vg = {vg:.1f}V')
        
        # Individual output curve (optional)
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        
        # Output curve
        axes[0].plot(vd_subset, np.array(id_subset))
        axes[0].set_xlabel('Drain Voltage (V)')
        axes[0].set_ylabel('Drain Current (A)')
        axes[0].set_title(f'Output Curve (Vg = {vg}V)')
        axes[0].grid(True)
        
        # Output conductance
        axes[1].plot(vd_subset, gd)
        axes[1].set_xlabel('Drain Voltage (V)')
        axes[1].set_ylabel('Output Conductance (A/V)')
        axes[1].set_title('Output Conductance (gd = dId/dVd)')
        axes[1].grid(True)
        
        fig.tight_layout()
        
        # Save individual plot only if explicitly requested
        if False:  # Change to True if you want these saved
            filename = f"{BASE_FILENAME}_Output_Vg{vg:.1f}V_Cycle{cycle_num}"
            fig.savefig(f"{filename}.png")
        plt.close(fig)
        
        # Add to results
        results.append({
            'vg_fixed': vg,
            'data': data,
            'cycle': cycle_num
        })
        
        # Check for interrupt
        if interrupted:
            print(f"Interrupting output sweep at Vg = {vg}V")
            break
    
    # Finalize and save the combined output plot
    plt.xlabel('Drain Voltage (V)')
    plt.ylabel('Drain Current (A)')
    plt.title(f'Output Curves - Cycle {cycle_num}')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    # Save combined plot only if explicitly requested
    if False:  # Change to True if you want these saved
        combined_filename = f"{BASE_FILENAME}_Output_Combined_Cycle{cycle_num}"
        plt.savefig(f"{combined_filename}.png")
    plt.close()
    
    # Save all data to a single file per cycle with only required columns
    all_output_data_slim = all_output_data[['Vd', 'Vg', 'Id', 'Ig', 'Vg_fixed']]
    all_output_data_slim.to_csv(f"{BASE_FILENAME}_Output_AllVg_Cycle{cycle_num}.csv", index=False)
    print(f"All output data saved to {BASE_FILENAME}_Output_AllVg_Cycle{cycle_num}.csv")
    
    return results

# Transistor Characterisation Suite

A Python-based measurement suite for automated characterisation of transistor devices using Keysight B2912A Source Measure Units (SMUs).

## Features

- **Complete Transistor Characterisfation**: Measure both transfer and output characteristics with a single script
- **Automated Parameter Extraction**: Calculate field-effect mobility, threshold voltage, and other key parameters
- **Advanced Plotting**: Visualize data with linear and semi-log plots for comprehensive analysis
- **Multiple Sweep Modes**:
  - Transfer sweeps (Vg vs Id with fixed Vd)
  - Output sweeps (Vd vs Id with fixed Vg)
  - Sequential transfer-output cycles

## Requirements

- Python 3.7 or higher
- Required Python packages:
  - pyvisa
  - matplotlib
  - numpy
  - pandas
  - scipy
- Keysight B2912A or compatible SMU with VISA connectivity

## Usage

1. Connect your transistor:
   - Gate to Source1 (Channel 1) of the SMU
   - Drain to Source2 (Channel 2) of the SMU
   - Source to ground

2. Configure measurement parameters at the top of the script:
   ```python
   # Measurement Mode - Choose one: "transfer", "output", or "transfer-output"
   MEASUREMENT_MODE = "transfer-output"
   
   # Voltage ranges and steps
   TRANSFER_VG_START = -50e-3
   TRANSFER_VG_STOP = 1.0
   TRANSFER_VG_POINTS = 102
   TRANSFER_VD_FIXED = [1.0]
   
   # Device parameters
   CHANNEL_LENGTH = 20e-6
   CHANNEL_WIDTH = 1000e-6
   ```

3. Run the script directly from your IDE:
   ```
   python transistor_measurement.py
   ```

4. Measurement data and parameters will be saved to the "Transistor_Measurements" directory.

## Measurement Overview

### Transfer Sweep
Measures drain current (Id) as a function of gate voltage (Vg) with fixed drain voltage (Vd):
- Extracts field-effect mobility using sqrt(Id) vs Vg
- Determines threshold voltage using linear extrapolation
- Calculates transconductance (gm)
- Provides linear and semi-log plots

### Output Sweep
Measures drain current (Id) as a function of drain voltage (Vd) with fixed gate voltage (Vg):
- Creates family of curves with multiple Vg values
- Calculates output conductance (gd)

### Cycle Testing
Performs multiple measurement cycles to analyze device stability or degradation:
- Tracks mobility, threshold voltage, and maximum current vs cycle number
- Creates comparison plots for all cycles

## Data Output

- **CSV Files**: Raw measurement data with Vg, Vd, Id, and Ig values for each sweep
- **Parameter File**: Extracted parameters (mobility, threshold voltage, max current) by cycle
- **Interactive Plots**: Displayed at the end of each measurement sequence

## Keyboard Control

Press Ctrl+C at any time to safely interrupt the measurements. The script will complete the current measurement step and perform a clean exit.

## Customization

- **Integration Time**: Adjust NPLC setting for measurement quality vs speed
- **Device Parameters**: Set device dimensions and oxide capacitance
- **Voltage Ranges**: Configure sweep parameters for your specific device

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Based on SCPI commands for Keysight B2900 SMU series
- Mobility calculation based on the standard FET model: μ = (2L/(W·Ci))·(∂√Id/∂Vgs)²

## Overview of Transfer/Output sweep 
This Python program provides automated characterization of transistor devices using a Keysight B2912A Source Measure Unit (SMU). It performs alternating transfer and output characteristic measurements in cycles, allowing for studies of device stability, degradation, and hysteresis effects.

## Key Features
- **Automated Measurement Cycles**: Alternates between transfer and output sweeps for a specified number of cycles
- **Transfer Sweep**: Measures Ids vs. Vgs at fixed Vds with optional forward/backward sweeps for hysteresis observation
- **Output Sweep**: Measures Ids vs. Vds at various fixed Vgs values
- **Field-Effect Mobility Calculation**: Automatically extracts mobility and threshold voltage
- **Real-Time Visualization**: Displays measurement data as it's acquired
- **Comprehensive Data Storage**: Saves well-labeled CSV files and images for each measurement

## Requirements
- Python 3.6 or higher
- PyVISA 1.11 or higher
- NumPy, Matplotlib, Pandas, SciPy
- Keysight B2912A Source Measure Unit (or compatible unit with SCPI interface)
- GPIB/USB interface for connecting to the SMU

## Installation

1. Clone this repository:
```bash
git clone https://github.com/adambickerdike/SMU-Transistor-Control-AB.git
cd transistor-characterization
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

Alternatively, install packages individually:
```bash
pip install pyvisa numpy matplotlib pandas scipy
```

3. Install the appropriate VISA backend. For NI-VISA:
```bash
pip install pyvisa-py
```

## Connection Setup
1. Connect the SMU to your computer via GPIB or USB
2. Connect the transistor to the SMU:
   - SMU1 (Channel 1) → Gate terminal
   - SMU2 (Channel 2) → Drain terminal
   - Ground → Source terminal
3. Ensure all connections are secure before starting measurements

## Usage

1. Update the `Config` class with your specific parameters:
```python
# Update SMU address
smu_address = 'GPIB0::24::INSTR'  # Change to match your setup

# Configure measurement parameters
cycles = 3                     # Number of measurement cycles
transfer_vds = 5.0             # Fixed Vds for transfer sweep
transfer_backward_sweep = True # Enable/disable hysteresis observation
```

2. Run the program:
```bash
python transistor_characterization.py
```

3. Press Enter when prompted to start the measurement

4. Data will be saved in a timestamped directory under "Transistor_Measurements/"

## Configuration Options

### General Settings
- `save_directory`: Base directory for saving data
- `device_name`: Name prefix for saving files
- `cycles`: Number of transfer-output measurement cycles
- `smu_address`: VISA address of the SMU

### Transfer Sweep Settings
- `transfer_vds`: Fixed drain-source voltage during transfer sweep [V]
- `transfer_vgs_start`: Starting gate voltage [V]
- `transfer_vgs_stop`: Ending gate voltage [V]
- `transfer_vgs_step`: Step size for gate voltage [V]
- `transfer_backward_sweep`: Enable/disable backward sweep for hysteresis observation

### Output Sweep Settings
- `output_vds_start`: Starting drain voltage [V]
- `output_vds_stop`: Ending drain voltage [V]
- `output_vds_step`: Step size for drain voltage [V]
- `output_vgs_values`: List of gate voltages for output sweeps [V]

### Timing Settings
- `sweep_delay`: Delay between points in sweep loop [s]
- `settling_time`: Settling time between measurements [s]
- `integration_time`: Integration time for measurements [s]

### Display Settings
- `live_plot`: Enable/disable live plotting during measurement
- `plot_update_interval`: Update live plot every N data points

### Device Parameters for Mobility Calculation
- `channel_length`: Channel length [m]
- `channel_width`: Channel width [m]
- `oxide_capacitance`: Gate oxide capacitance per unit area [F/m²]
- `mobility_fit_low_vgs`: Lower Vgs boundary for mobility extraction
- `mobility_fit_high_vgs`: Upper Vgs boundary for mobility extraction

## Output Files

For each measurement cycle, the program generates:

### Data Files
- `Transfer_Sweep_{cycle_num}.csv`: Transfer sweep data
- `Output_Sweep_{cycle_num}.csv`: Output sweep data
- `Mobility_vs_Cycle.csv`: Extracted mobility values for each cycle
- `PeakCurrent_vs_Cycle.csv`: Peak current values for each cycle
- `measurement_config.txt`: Configuration settings used for the measurement

### Plot Files
- `Transfer_Plot_{cycle_num}.png`: Transfer characteristics (linear and log scales)
- `Output_Plot_{cycle_num}.png`: Output characteristics
- `Mobility_vs_Cycle.png`: Mobility trend across cycles
- `PeakCurrent_vs_Cycle.png`: Peak current trend across cycles

## Example CSV Output

### Transfer_Sweep_1.csv (with backward sweep enabled):
```
Vds,Vgs,Ids,Igs,Time_since_start,Sweep_Direction
5.0,-5.0,2.14e-12,3.45e-10,0.127,Forward
5.0,-4.0,5.67e-12,3.47e-10,0.381,Forward
...
5.0,10.0,1.96e-02,3.74e-10,3.937,Forward
5.0,10.0,1.95e-02,3.75e-10,4.127,Backward
...
5.0,-5.0,8.96e-12,3.46e-10,7.937,Backward
```

### Output_Sweep_1.csv:
```
Vds,Vgs,Ids,Igs,Time_since_start
0.0,0.0,1.23e-11,2.34e-11,0.127
1.0,0.0,1.45e-10,2.34e-11,0.254
...
10.0,10.0,5.47e-05,2.19e-10,11.811
```

## Interrupting Measurements
Press Ctrl+C at any time to cleanly interrupt the measurement process. The program will save collected data before exiting.

## Troubleshooting

### Communication Issues
- Verify the correct GPIB/USB address in the configuration
- Check physical connections between computer and SMU
- Ensure only one application is accessing the SMU at a time

### Measurement Errors
- Check compliance settings if current/voltage is being limited
- Verify device connections (gate, drain, source)
- Increase settling time if measurements seem unstable

### Data Quality Issues
- Adjust integration time for better measurement precision
- Check for environmental noise sources
- Use shielded cables and proper grounding

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
- Based on PyVISA and the SCPI command set for Keysight instruments
- Developed for semiconductor device characterization research

## Overview of Transient Measurement (beta)
This Python program enables time-domain characterization of transistors by applying pulsed gate voltage signals while measuring the resulting drain current response. It uses a Keysight B2912A Source Measure Unit (SMU) to control both gate and drain voltages with precise timing, allowing for detailed analysis of transient behavior.

## Key Features
- **Configurable Pulse Train**: User-defined pulse period, duty cycle, amplitude, and count
- **High-Resolution Time-Domain Measurements**: Customizable sampling interval for capturing fast transient responses
- **Fixed Drain Bias**: Maintains constant drain-source voltage during the measurement
- **Real-Time Visualization**: Displays both gate voltage and drain current vs. time as the measurement progresses
- **Comprehensive Data Storage**: Saves time-series data with accurate timestamps for further analysis

## Usage

1. Update the `Config` class with your specific parameters:
```python
# Update SMU address
smu_address = 'GPIB0::24::INSTR'  # Change to match your setup

# Configure pulse parameters
gate_low_voltage = 0.0            # Gate voltage when pulse is low [V]
gate_high_voltage = 10.0          # Gate voltage when pulse is high [V]
pulse_period = 500e-3             # Complete pulse cycle period [s]
duty_cycle = 50                   # Duty cycle percentage [%]
num_pulses = 5                    # Number of pulses to measure

# Set drain bias
drain_bias_voltage = 5.0          # Fixed Vds during measurement [V]

# Set measurement timing
sampling_interval = 5e-3          # Time between measurements [s]
```

2. Run the program:
```bash
python transistor_transient_measurement.py
```

3. Press Enter when prompted to start the measurement

4. Data will be saved in a timestamped directory under "Transistor_Transient_Measurements/"

## Configuration Options

### General Settings
- `save_directory`: Base directory for saving data
- `device_name`: Name prefix for saving files
- `smu_address`: VISA address of the SMU

### Drain Bias Settings
- `drain_bias_voltage`: Fixed drain-source voltage during measurement [V]

### Gate Pulse Settings
- `gate_low_voltage`: Gate voltage when pulse is low [V]
- `gate_high_voltage`: Gate voltage when pulse is high [V]
- `pulse_period`: Complete pulse cycle period [s]
- `duty_cycle`: Duty cycle percentage [%]
- `num_pulses`: Number of pulses to measure

### Measurement Timing Settings
- `pre_pulse_delay`: Time before first pulse starts [s]
- `post_pulse_delay`: Time after last pulse ends [s]
- `sampling_interval`: Time between measurements [s]

### Measurement Settings
- `compliance_current`: Maximum allowed current [A]
- `integration_time`: Integration time for each measurement [s]
- `averaging_points`: Number of readings to average per point

### Plot Settings
- `live_plot`: Enable/disable live plotting
- `update_interval`: Update plot every N data points

## Output Files

The program generates the following files in the timestamped output directory:

### Data Files
- `Transient_Measurement.csv`: Primary data file containing all time-series measurements
- `Expected_Waveform.csv`: Reference file with the programmed gate voltage waveform
- `measurement_config.txt`: Configuration settings used for the measurement

### Plot Files
- `Transient_Measurement.png`: Dual-panel plot showing gate voltage and drain current vs. time

## Example CSV Output

### Transient_Measurement.csv:
```
Time,Vgs,Ids,Igs
0.000,0.00,3.25e-10,1.23e-11
0.005,0.00,3.24e-10,1.23e-11
...
0.995,0.00,3.23e-10,1.22e-11
1.000,10.00,2.14e-08,2.45e-11
1.005,10.00,7.89e-07,2.46e-11
...
1.245,10.00,8.75e-04,2.49e-11
1.250,0.00,8.72e-04,1.24e-11
1.255,0.00,8.54e-04,1.24e-11
...
```

### Expected_Waveform.csv:
```
Time,Expected_Vgs
0.000,0.00
0.005,0.00
...
0.995,0.00
1.000,10.00
1.005,10.00
...
1.245,10.00
1.250,0.00
...
```

## Analysis Tips

### Typical Transient Parameters to Extract
- **Rise Time**: Time for drain current to rise from 10% to 90% of its final value
- **Fall Time**: Time for drain current to fall from 90% to 10% of its peak value
- **Delay Time**: Time between gate voltage change and onset of drain current response
- **Overshoot/Undershoot**: Magnitude of current peaks beyond steady-state values
- **Settling Time**: Time required for current to stabilize within a certain percentage of final value

Example Python code for extracting these parameters:
```python
import pandas as pd
import numpy as np

# Load data
data = pd.read_csv('Transient_Measurement.csv')

# Find rising edge of first pulse
rising_edge = np.where(np.diff(data['Vgs']) > 5)[0][0]
falling_edge = np.where(np.diff(data['Vgs']) < -5)[0][0]

# Extract rise time
steady_current = data['Ids'][rising_edge+20:falling_edge].mean()
current_10pct = 0.1 * steady_current
current_90pct = 0.9 * steady_current

rise_10pct_idx = rising_edge + np.argmin(abs(data['Ids'][rising_edge:falling_edge] - current_10pct))
rise_90pct_idx = rising_edge + np.argmin(abs(data['Ids'][rising_edge:falling_edge] - current_90pct))

rise_time = data['Time'][rise_90pct_idx] - data['Time'][rise_10pct_idx]
print(f"Rise time: {rise_time*1000:.2f} ms")
```

## Interrupting Measurements
Press Ctrl+C at any time to cleanly interrupt the measurement process. The program will save collected data before exiting.

## Troubleshooting

### Timing/Sampling Issues
- If transients are too fast for your sampling rate, decrease `sampling_interval`
- For noisy measurements, increase `integration_time` (sacrifices timing resolution)
- If gate voltage changes aren't being captured properly, ensure `sampling_interval` is at least 10× smaller than pulse transitions

### Communication Issues
- Verify the correct GPIB/USB address in the configuration
- Check physical connections between computer and SMU
- Ensure only one application is accessing the SMU at a time

### Measurement Errors
- Check compliance settings if current/voltage is being limited
- Verify device connections (gate, drain, source)
- For very fast transients, consider the capacitance of your measurement setup

## Performance Optimization
- For faster measurements, decrease `integration_time` (may increase noise)
- To reduce system load, increase `update_interval` for live plotting
- For memory-constrained systems, increase `sampling_interval` to reduce data points

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
- Based on PyVISA and the SCPI command set for Keysight instruments


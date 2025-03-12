"""
Transistor Transient Measurement with Pulsed Gate Voltage
---------------------------------------------------------
This program controls a Keysight B2912A SMU to perform transient measurements
of transistor drain current in response to pulsed gate voltage.

Connection setup:
- SMU1 connected to gate terminal (pulsed voltage)
- SMU2 connected to drain terminal (fixed bias voltage)
- Source terminal connected to ground

Features:
- Configurable pulse train parameters (period, duty cycle, magnitude)
- Fixed drain bias voltage
- High-resolution time domain measurements
- Real-time visualization
- CSV data storage with timestamps
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
import signal
import sys

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
    save_directory = "Transistor_Transient_Measurements"
    device_name = "TestTransistor"
    
    # Drain bias settings
    drain_bias_voltage = 5.0          # Fixed Vds during measurement [V]
    
    # Gate pulse settings
    gate_low_voltage = 0.0            # Gate voltage when pulse is low [V]
    gate_high_voltage = 10.0          # Gate voltage when pulse is high [V]
    pulse_period = (500e-3)           # Complete pulse cycle period [s]
    duty_cycle = 50                   # Duty cycle percentage [%]
    num_pulses = 5                    # Number of pulses to measure
    
    # Measurement timing settings
    pre_pulse_delay = 1.0             # Time before first pulse starts [s]
    post_pulse_delay = 1.0            # Time after last pulse ends [s]
    sampling_interval = 5e-3          # Time between measurements [s]
    
    # Measurement settings
    compliance_current = 0.1          # Compliance current [A]
    integration_time = 0.5e-3         # Integration time for each measurement [s]
    averaging_points = 1              # Number of readings to average per point
    
    # Plot settings
    live_plot = True                  # Enable/disable live plotting
    update_interval = 10              # Update plot every N data points
    
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

def calculate_measurement_parameters():
    """Calculate timing and voltage parameters for the pulse train"""
    # Calculate pulse timing
    high_time = Config.pulse_period * Config.duty_cycle / 100
    low_time = Config.pulse_period - high_time
    
    # Calculate total measurement time
    total_time = (Config.pre_pulse_delay + 
                 (Config.pulse_period * Config.num_pulses) + 
                 Config.post_pulse_delay)
    
    # Calculate number of measurement points
    num_points = int(total_time / Config.sampling_interval) + 1
    
    # Create time points array
    time_points = np.linspace(0, total_time, num_points)
    
    # Create expected Vgs array for plotting comparison
    vgs_expected = np.ones(num_points) * Config.gate_low_voltage
    
    # Add pulses to the expected Vgs array
    for pulse in range(Config.num_pulses):
        pulse_start = Config.pre_pulse_delay + (pulse * Config.pulse_period)
        pulse_end = pulse_start + high_time
        
        # Set high voltage during pulse
        pulse_indices = np.where((time_points >= pulse_start) & 
                                (time_points < pulse_end))
        vgs_expected[pulse_indices] = Config.gate_high_voltage
    
    params = {
        'high_time': high_time,
        'low_time': low_time,
        'total_time': total_time,
        'num_points': num_points,
        'time_points': time_points,
        'vgs_expected': vgs_expected
    }
    
    return params

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
            
            # Integration time
            self.smu.write(f':SENSe{ch}:CURRent:DC:APERture %G' % (Config.integration_time))
            self.smu.write(f':SENSe{ch}:CURRent:DC:APERture:AUTO %d' % (0))  # Disable auto integration time

        # Set both channels as voltage sources
        self.smu.write(':SOURce1:FUNCtion:MODE %s' % ('VOLTage'))
        self.smu.write(':SOURce2:FUNCtion:MODE %s' % ('VOLTage'))
        
        # Error check
        print("SMU Configuration Complete. Errors:", self.smu.query(':SYSTem:ERRor:ALL?'))
        self.smu.write('*CLS')
    
    def setup_transient_measurement(self):
        """Set up SMU for transient measurement with gate pulse and drain bias"""
        # Set initial gate voltage to low level
        self.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude %G' % (Config.gate_low_voltage))
        
        # Set fixed drain bias voltage
        self.smu.write(':SOURce2:VOLTage:LEVel:IMMediate:AMPLitude %G' % (Config.drain_bias_voltage))
        
        # Error check
        self.smu.query(':SYSTem:ERRor:ALL?')
        self.smu.write('*CLS')
    
    def set_gate_voltage(self, voltage):
        """Set gate voltage to specified level"""
        self.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude %G' % (voltage))
    
    def measure_currents(self):
        """Measure gate and drain currents"""
        # Initialize and trigger measurement
        self.smu.write(':INIT (%s)' % ('@1,2'))
        self.smu.write(':TRIGger:ALL (%s)' % ('@1,2'))
        
        # Read data
        Igs = self.smu.query_binary_values(':SENSe1:DATA? %s,%d' % ('CURRent', 1), 'd', False)[0]
        Ids = self.smu.query_binary_values(':SENSe2:DATA? %s,%d' % ('CURRent', 1), 'd', False)[0]
        
        # Error check
        self.smu.query(':SYSTem:ERRor:ALL?')
        self.smu.write('*CLS')
        
        return Igs, Ids
        
    def close(self):
        """Close connection to SMU"""
        self.reset()
        self.smu.close()
        self.rm.close()
        print("SMU connection closed")

##########################################################################
####################### Plotting Functions ###############################
def setup_transient_plot():
    """Set up plot for transient measurements"""
    plt.ion()  # Turn on interactive mode
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Top plot for gate voltage
    ax1.set_ylabel('Gate Voltage, Vgs (V)')
    ax1.set_title('Pulsed Gate Voltage')
    ax1.grid(True)
    
    # Bottom plot for drain current
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Drain Current, Ids (A)')
    ax2.set_title('Drain Current Response')
    ax2.grid(True)
    
    plt.tight_layout()
    
    return fig, (ax1, ax2)

def update_transient_plot(axes, data, expected_vgs=None):
    """Update the transient measurement plot with latest data"""
    ax1, ax2 = axes
    
    # Clear previous data
    ax1.clear()
    ax2.clear()
    
    # Plot gate voltage
    ax1.plot(data['Time'], data['Vgs'], 'b-', label='Measured Vgs')
    if expected_vgs is not None:
        ax1.plot(data['Time'], expected_vgs[:len(data['Time'])], 'r--', label='Programmed Vgs')
    ax1.set_ylabel('Gate Voltage, Vgs (V)')
    ax1.set_title('Pulsed Gate Voltage')
    ax1.grid(True)
    ax1.legend(loc='best')
    
    # Plot drain current
    ax2.plot(data['Time'], data['Ids'], 'g-')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Drain Current, Ids (A)')
    ax2.set_title('Drain Current Response')
    ax2.grid(True)
    
    plt.tight_layout()
    plt.draw()
    plt.pause(0.001)

##########################################################################
####################### Measurement Function #############################
def perform_transient_measurement(smu, save_dir):
    """Perform transient measurement with pulsed gate voltage"""
    print("Preparing transient measurement...")
    
    # Setup SMU for transient measurement
    smu.setup_transient_measurement()
    
    # Calculate timing and expected voltage waveforms
    params = calculate_measurement_parameters()
    high_time = params['high_time']
    low_time = params['low_time']
    total_time = params['total_time']
    time_points = params['time_points']
    vgs_expected = params['vgs_expected']
    
    # Create data storage
    data = {
        'Time': [],
        'Vgs': [],
        'Ids': [],
        'Igs': []
    }
    
    # Setup plot if enabled
    if Config.live_plot:
        fig, axes = setup_transient_plot()
    
    # Start measurement
    print(f"Starting transient measurement...")
    print(f"Total measurement time: {total_time:.2f} s")
    print(f"Number of measurement points: {params['num_points']}")
    print(f"Gate voltage: Low = {Config.gate_low_voltage:.2f} V, High = {Config.gate_high_voltage:.2f} V")
    print(f"Pulse period: {Config.pulse_period*1000:.1f} ms, Duty cycle: {Config.duty_cycle}%")
    
    start_time = time.time()
    measurement_start = start_time
    
    # Initial gate voltage is low
    current_vgs = Config.gate_low_voltage
    smu.set_gate_voltage(current_vgs)
    
    try:
        # Wait for pre-pulse delay
        pre_pulse_end = start_time + Config.pre_pulse_delay
        print(f"Pre-pulse delay ({Config.pre_pulse_delay:.2f} s)...")
        
        # Main measurement loop
        point_count = 0
        next_measurement_time = start_time
        
        while time.time() < start_time + total_time:
            if interrupted:
                print("Measurement interrupted!")
                break
            
            current_time = time.time()
            elapsed = current_time - measurement_start
            
            # Only measure at the specified sampling intervals
            if current_time >= next_measurement_time:
                # Determine what the gate voltage should be at this time
                if elapsed < Config.pre_pulse_delay:
                    # Pre-pulse period - gate is low
                    if current_vgs != Config.gate_low_voltage:
                        current_vgs = Config.gate_low_voltage
                        smu.set_gate_voltage(current_vgs)
                elif elapsed > (Config.pre_pulse_delay + Config.num_pulses * Config.pulse_period):
                    # Post-pulse period - gate is low
                    if current_vgs != Config.gate_low_voltage:
                        current_vgs = Config.gate_low_voltage
                        smu.set_gate_voltage(current_vgs)
                else:
                    # During pulse train
                    pulse_time = elapsed - Config.pre_pulse_delay
                    pulse_number = int(pulse_time / Config.pulse_period)
                    phase = pulse_time % Config.pulse_period
                    
                    # Determine if we're in the high or low portion of the pulse
                    if phase < high_time:
                        # High portion of pulse
                        if current_vgs != Config.gate_high_voltage:
                            current_vgs = Config.gate_high_voltage
                            smu.set_gate_voltage(current_vgs)
                    else:
                        # Low portion of pulse
                        if current_vgs != Config.gate_low_voltage:
                            current_vgs = Config.gate_low_voltage
                            smu.set_gate_voltage(current_vgs)
                
                # Measure currents
                Igs, Ids = smu.measure_currents()
                
                # Store data
                data['Time'].append(elapsed)
                data['Vgs'].append(current_vgs)
                data['Ids'].append(Ids)
                data['Igs'].append(Igs)
                
                # Print progress occasionally
                point_count += 1
                if point_count % 20 == 0:
                    print(f"Time: {elapsed:.3f} s, Vgs: {current_vgs:.2f} V, Ids: {Ids:.6e} A")
                
                # Update live plot if enabled
                if Config.live_plot and point_count % Config.update_interval == 0:
                    update_transient_plot(axes, data, vgs_expected)
                
                # Schedule next measurement
                next_measurement_time = current_time + Config.sampling_interval
            
            # Small sleep to prevent CPU hogging
            time.sleep(min(Config.sampling_interval/10, 0.001))
    
    finally:
        # Set gate voltage back to low level
        smu.set_gate_voltage(Config.gate_low_voltage)
        
        # Final plot update
        if Config.live_plot and len(data['Time']) > 0:
            update_transient_plot(axes, data, vgs_expected)
            plt.savefig(os.path.join(save_dir, "Transient_Measurement.png"))
        
        # Save data to CSV
        if len(data['Time']) > 0:
            filename = os.path.join(save_dir, "Transient_Measurement.csv")
            df = pd.DataFrame(data)
            df.to_csv(filename, index=False)
            print(f"Transient data saved to {filename}")
            
            # Also save a separate file with expected voltage waveform for reference
            reference_data = {
                'Time': time_points,
                'Expected_Vgs': vgs_expected
            }
            ref_filename = os.path.join(save_dir, "Expected_Waveform.csv")
            pd.DataFrame(reference_data).to_csv(ref_filename, index=False)
        
        return data

##########################################################################
####################### Main Function ####################################
def main():
    """Main function to execute the transient measurement"""
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
    
    # Start time for the measurement
    start_time = time.time()
    
    try:
        # Perform transient measurement
        data = perform_transient_measurement(smu, save_dir)
        
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
        
        total_time = time.time() - start_time
        print(f"\nMeasurement completed in {total_time:.2f} seconds.")
        
        if interrupted:
            print("Measurement was interrupted by user.")
        else:
            print("Measurement completed successfully.")

if __name__ == "__main__":
    print("\n" + "="*50)
    print("Transistor Transient Measurement Program")
    print("="*50)
    
    # Create base directory if it doesn't exist
    os.makedirs(Config.save_directory, exist_ok=True)
    
    # Prompt user to start
    input("Press Enter to start the transient measurement...")
    
    # Run main function
    main()
    
    print("\n" + "="*50)
    print("Program execution completed")
    print("="*50)
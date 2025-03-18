clear all;
close all;
clc;

%% User Settings - Modify these variables
% Base filename pattern
baseFileName = 'Transistor_XXXX_20250317_124810_Transfer_AllVd';
% Parameters filename
parametersFileName = 'Transistor_XXXX_20250317_124810_Parameters_vs_Cycle.csv';
% Maximum cycle number to look for
maxCycle = 40;
% Set to true if you want to save the output figure
saveOutput = false;

%% Load Parameters vs Cycle Data
if exist(parametersFileName, 'file')
    try
        fprintf('Loading parameters file: %s\n', parametersFileName);
        paramsData = readtable(parametersFileName);
        hasParamsData = true;
        fprintf('Successfully loaded parameters data with %d cycles\n', height(paramsData));
    catch e
        fprintf('Error loading parameters file: %s\n', e.message);
        hasParamsData = false;
    end
else
    fprintf('Parameters file not found: %s\n', parametersFileName);
    hasParamsData = false;
end

%% Setup Figure for Transfer Curves
figure('Position', [100, 100, 850, 600], 'Name', 'Transfer Curves');
set(gcf, 'NumberTitle', 'off'); % Hide figure number

% Create the axes with proper logarithmic scale
ax = gca;
set(ax, 'YScale', 'log');  % Set y-axis to logarithmic scale
set(ax, 'Box', 'on');      % Turn on the box
set(ax, 'FontSize', 12);   % Set consistent font size
set(ax, 'YGrid', 'on');    % Turn on just the horizontal grid lines
hold on;

% Define a colormap for the curves (goes from blue to purple)
cmap = cool(maxCycle);

% Keep track of successful loads for the legend
successfulCycles = [];

% Create array to store all cycle curve handles for colorbar
allPlots = [];

%% Process each cycle file
for cycle = 1:maxCycle
    % Construct the filename for this cycle
    filename = sprintf('%s_Cycle%d.csv', baseFileName, cycle);
    
    % Check if the file exists
    if exist(filename, 'file')
        try
            % Load the data
            fprintf('Loading cycle %d: %s\n', cycle, filename);
            data = readtable(filename);
            
            % Extract the relevant columns (Vg for x-axis, Id for y-axis)
            Vg = data.Vg;
            Id = data.Id;
            
            % Get absolute value of current for log plotting
            Id_abs = abs(Id);
            
            % Plot the data with color based on cycle number
            colorIdx = round((cycle / maxCycle) * size(cmap, 1));
            if colorIdx < 1
                colorIdx = 1;
            end
            
            % Plot the data directly and save handle
            h = plot(Vg, Id_abs, 'LineWidth', 1.5, 'Color', cmap(colorIdx,:));
            
            % Store the handle for colorbar mapping
            allPlots(cycle) = h;
            
            % Add to the list of successful cycles
            successfulCycles = [successfulCycles, cycle];
            
            fprintf('  Successfully plotted cycle %d\n', cycle);
        catch e
            fprintf('  Error loading cycle %d: %s\n', cycle, e.message);
        end
    else
        fprintf('File not found for cycle %d: %s\n', cycle, filename);
    end
end

%% Finalize the transfer curves plot
% Add labels and title
xlabel('Gate Voltage (V_{gs}) (V)', 'FontSize', 14);
ylabel('Drain Current (|I_{ds}|) (A)', 'FontSize', 14);

% Set appropriate y-axis limits for logarithmic scale
% Find the min and max currents across all successfully loaded data
currMin = inf;
currMax = 0;
for cycle = successfulCycles
    filename = sprintf('%s_Cycle%d.csv', baseFileName, cycle);
    if exist(filename, 'file')
        data = readtable(filename);
        Id_abs = abs(data.Id);
        % Filter out zeros for log plotting
        Id_abs = Id_abs(Id_abs > 0);
        if ~isempty(Id_abs)
            currMin = min(currMin, min(Id_abs));
            currMax = max(currMax, max(Id_abs));
        end
    end
end

% Set y-axis limits with some padding
if ~isinf(currMin) && currMin > 0
    % Use one decade below the minimum as the lower limit
    yLowerLimit = 10^(floor(log10(currMin)));
    % Use one decade above the maximum as the upper limit
    yUpperLimit = 10^(ceil(log10(currMax)));
    ylim([yLowerLimit, yUpperLimit]);
end

% Customize grid appearance
grid on;
set(gca, 'YMinorGrid', 'on'); % Add minor grid lines for y-axis
set(gca, 'GridLineStyle', ':'); % Dotted grid lines
set(gca, 'MinorGridLineStyle', ':'); % Dotted minor grid lines

% Format the tick labels to use scientific notation
set(gca, 'YTickLabelMode', 'auto');

% Add colorbar
cb = colorbar;
colormap(cool);
caxis([1 maxCycle]);

% Create custom tick marks for the colorbar
if maxCycle > 10
    % If many cycles, use appropriate number of ticks
    tickStep = ceil(maxCycle/10); % Approximately 10 ticks
    cbTicks = 1:tickStep:maxCycle;
    if cbTicks(end) ~= maxCycle
        cbTicks = [cbTicks, maxCycle];
    end
else
    % For fewer cycles, show all ticks
    cbTicks = 1:maxCycle;
end
cb.Ticks = cbTicks;

% Set the label on the side and make it same size as axis labels
cb.Label.String = 'Cycle Number';
cb.Label.FontSize = 14;
cb.FontSize = 12;

% Position the colorbar optimally
cb.Position = [0.92 0.15 0.02 0.7];

hold off;

%% Plot Parameters vs Cycle if available
if hasParamsData
    % Create a figure with 3 subplots for mobility, threshold, and max current (vertical arrangement)
    figure('Position', [900, 100, 500, 800], 'Name', 'Parameters vs Cycle');
    set(gcf, 'NumberTitle', 'off'); % Hide figure number
    
    % 1. Mobility vs Cycle
    subplot(3, 1, 1);
    plot(paramsData.Cycle, paramsData.Mobility_Vd1_0V, '^-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerFaceColor', 'auto');
    xlabel('Cycle Number', 'FontSize', 12);
    ylabel('Mobility (cm^{2}/Vs)', 'FontSize', 12);
    grid on;
    ylim([min(paramsData.Mobility_Vd1_0V)*0.9, max(paramsData.Mobility_Vd1_0V)*1.1]);
    
    % 2. Threshold Voltage vs Cycle
    subplot(3, 1, 2);
    plot(paramsData.Cycle, paramsData.Vth_Vd1_0V, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerFaceColor', 'auto');
    xlabel('Cycle Number', 'FontSize', 12);
    ylabel('Threshold Voltage (V)', 'FontSize', 12);
    grid on;
    ylim([min(paramsData.Vth_Vd1_0V)*0.9, max(paramsData.Vth_Vd1_0V)*1.1]);
    
    % 3. Max Current vs Cycle
    subplot(3, 1, 3);
    plot(paramsData.Cycle, paramsData.MaxCurrent_Vd1_0V_A * 1000, 's-', 'LineWidth', 1.5, 'Color', 'g', 'MarkerFaceColor', 'auto');  % Convert to mA
    xlabel('Cycle Number', 'FontSize', 12);
    ylabel('Maximum Current (mA)', 'FontSize', 12);
    grid on;
    ylim([min(paramsData.MaxCurrent_Vd1_0V_A)*0.9*1000, max(paramsData.MaxCurrent_Vd1_0V_A)*1.1*1000]);
    
    % Save parameters figure if requested
    if saveOutput
        saveas(gcf, 'TransistorParameters.png');
        saveas(gcf, 'TransistorParameters.fig');
        fprintf('Parameters figures saved as TransistorParameters.png and TransistorParameters.fig\n');
    end
end

%% Save the transfer curves figure if requested
if saveOutput
    % Go back to the first figure
    figure(1);
    saveas(gcf, 'TransistorTransferCurves.png');
    saveas(gcf, 'TransistorTransferCurves.fig');
    fprintf('Transfer curves figures saved as TransistorTransferCurves.png and TransistorTransferCurves.fig\n');
end

% Display summary information
fprintf('\nScript completed.\n');
fprintf('Plotted %d cycles from %s\n', length(successfulCycles), baseFileName);
if hasParamsData
    fprintf('Plotted parameters data with %d cycles\n', height(paramsData));
end
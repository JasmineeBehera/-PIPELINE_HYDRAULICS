% =========================================================================
%      PIPELINE HYDRAULICS AND PRESSURE RELIEF SYSTEM ANALYSIS
% =========================================================================
clear;
clc;
close all;

fprintf('================================================================\n');
fprintf('  Pipeline Hydraulics & Pressure Relief System Analysis\n');
fprintf('================================================================\n\n');

%% 1. Define System Input Parameters
fprintf('--- [Step 1] Defining system parameters...\n');
% Create a structure to hold all parameters for easy passing to functions
% --- Fluid Properties (Water at 20 C) ---
P.rho = 998;         % Density (kg/m^3)
P.mu = 1.0e-3;       % Dynamic Viscosity (Pa.s)

% --- Pipe Properties (Commercial Steel) ---
P.L = 1200;          % Pipe Length (m)
P.D = 0.15;          % Pipe Inner Diameter (m) (approx. 6-inch pipe)
P.epsilon = 4.5e-5;  % Absolute Roughness (m)

% --- System Operating Conditions ---
P.Q_op = 0.04;       % Operating Volumetric Flow Rate (m^3/s) (approx. 634 GPM)
P.delta_z = 20;      % Elevation increase from pump to destination (m)
P.P_dest = 2.0e5;    % Required pressure at destination (Pa, gauge) (2 bar)
P.g = 9.81;          % Gravitational acceleration (m/s^2)

% --- Fittings and Minor Losses (K-values) ---
P.fittings_K = [
    0.75;  % 2x Long Radius 90-deg Elbow
    0.75;
    0.20;  % 1x Gate Valve (fully open)
    1.00   % 1x Check Valve
];

% --- Safety System Parameters for PRV Sizing ---
P.P_set = 12.0e5;    % PRV Set Pressure (Pa, gauge) (12 bar)
P.overpressure = 0.10; % Allowed overpressure (10% for single valve)
P.Kd = 0.975;        % Discharge coefficient for certified PRV (standard value)

%% 2. Perform Hydraulic Calculations for Normal Operation
fprintf('--- [Step 2] Calculating normal operation hydraulics...\n');
hydraulics_results = calculate_hydraulics(P.Q_op, P);

%% 3. Perform PRV Sizing for Blocked Outlet Scenario
fprintf('--- [Step 3] Sizing PRV for blocked outlet scenario...\n');
prv_results = size_pressure_relief_valve(P);

%% 4. Display Results in a Report
fprintf('\n================== SYSTEM HYDRAULICS REPORT ==================\n');
fprintf('OPERATING CONDITIONS:\n');
fprintf('  Flow Rate:                 %.3f m^3/s\n', P.Q_op);
fprintf('  Fluid Velocity:            %.2f m/s\n', hydraulics_results.velocity);
fprintf('  Reynolds Number:           %.0f (Turbulent)\n', hydraulics_results.Re);
fprintf('  Friction Factor (f):       %.4f\n', hydraulics_results.f);
fprintf('------------------------------------------------------------------\n');
fprintf('PRESSURE DROP ANALYSIS (Pa):\n');
fprintf('  Major Loss (Friction):     %.2e Pa (%.2f bar)\n', hydraulics_results.dP_major, hydraulics_results.dP_major/1e5);
fprintf('  Minor Losses (Fittings):   %.2e Pa (%.2f bar)\n', hydraulics_results.dP_minor, hydraulics_results.dP_minor/1e5);
fprintf('  Static Head Loss:          %.2e Pa (%.2f bar)\n', hydraulics_results.dP_static, hydraulics_results.dP_static/1e5);
fprintf('  --------------------------------------------------\n');
fprintf('  TOTAL PRESSURE TO OVERCOME:  %.2e Pa (%.2f bar)\n', hydraulics_results.dP_total, hydraulics_results.dP_total/1e5);
fprintf('------------------------------------------------------------------\n');
fprintf('PUMP REQUIREMENTS:\n');
fprintf('  Required Pump Head:        %.2f m\n', hydraulics_results.H_pump);
fprintf('  Required Pump Pressure:    %.2e Pa (%.2f bar)\n', hydraulics_results.P_pump, hydraulics_results.P_pump/1e5);
fprintf('==================================================================\n');

fprintf('\n================== SAFETY SYSTEM DESIGN REPORT ===================\n');
fprintf('SCENARIO: BLOCKED OUTLET\n');
fprintf('  PRV Set Pressure:          %.2f bar\n', P.P_set/1e5);
fprintf('  Relieving Pressure:        %.2f bar\n', prv_results.P_relieving/1e5);
fprintf('  Required Flow at Relief:   %.3f m^3/s\n', P.Q_op);
fprintf('------------------------------------------------------------------\n');
fprintf('PRV SIZING RESULTS (API 520):\n');
fprintf('  Calculated Orifice Area:   %.2f cm^2\n', prv_results.A_req_cm2);
fprintf('  Selected API Orifice Size:   "%s"\n', prv_results.API_orifice_letter);
fprintf('  Standard Area of "%s":     %.2f cm^2\n', prv_results.API_orifice_letter, prv_results.A_standard_cm2);
fprintf('==================================================================\n\n');

%% 5. Generate and Display the System Curve Plot
fprintf('--- [Step 4] Generating System Curve Plot...\n');
plot_system_curve(P, hydraulics_results);
fprintf('Plot generated successfully.\n');

% =========================================================================
%                       LOCAL FUNCTION DEFINITIONS
% =========================================================================

function results = calculate_hydraulics(Q, P)
    % Calculates all hydraulic parameters for a given flow rate Q.
    A_pipe = pi * (P.D/2)^2;
    results.velocity = Q / A_pipe;
    
    % Reynolds Number
    results.Re = (P.rho * results.velocity * P.D) / P.mu;
    
    % Friction factor (f) using the Haaland Equation (explicit approximation)
    if results.Re > 4000
        term1 = (P.epsilon / (3.7 * P.D));
        term2 = (6.9 / results.Re);
        inv_f_sqrt = -1.8 * log10(term1^1.11 + term2);
        results.f = (1 / inv_f_sqrt)^2;
    else
        results.f = 64 / results.Re; % Laminar flow
    end
    
    % Pressure Drops (in Pascals)
    % Major Loss (Darcy-Weisbach)
    results.dP_major = results.f * (P.L/P.D) * (P.rho * results.velocity^2 / 2);
    % Minor Losses
    K_total = sum(P.fittings_K);
    results.dP_minor = K_total * (P.rho * results.velocity^2 / 2);
    % Static Head
    results.dP_static = P.rho * P.g * P.delta_z;
    
    results.dP_total = results.dP_major + results.dP_minor + results.dP_static;
    
    % Pump Requirements
    results.P_pump = results.dP_total + P.P_dest; % Total pressure needed from pump
    results.H_pump = results.P_pump / (P.rho * P.g); % Convert to head (m)
end

function results = size_pressure_relief_valve(P)
    % Sizes a PRV for liquid service based on API 520 guidelines.
    
    % Relieving pressure is set pressure + overpressure
    results.P_relieving = P.P_set * (1 + P.overpressure);
    
    % Required relieving mass flow rate (kg/s)
    W_req = P.Q_op * P.rho;
    
    % Required orifice area (m^2) using fundamental liquid sizing equation
    % A = W / (Kd * sqrt(2 * rho * (P_relieving - P_backpressure)))
    % Assuming backpressure is atmospheric (0 Pa, gauge)
    P_back = 0;
    results.A_req_m2 = W_req / (P.Kd * sqrt(2 * P.rho * (results.P_relieving - P_back)));
    
    % Convert to cm^2 for comparison with API standards
    results.A_req_cm2 = results.A_req_m2 * 100^2;
    
    % --- Select next standard API 526 Orifice Size ---
    api_orifices.letters = {'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'T'};
    api_orifices.areas_cm2 = [0.710, 1.264, 1.978, 3.245, 5.065, 8.355, 11.87, 18.39, 23.23, 33.74, 41.29, 71.29, 103.2, 167.7];
    
    % Find the first standard orifice with an area larger than required
    idx = find(api_orifices.areas_cm2 >= results.A_req_cm2, 1, 'first');
    
    if isempty(idx)
        results.API_orifice_letter = 'T+ (Custom)';
        results.A_standard_cm2 = NaN;
    else
        results.API_orifice_letter = api_orifices.letters{idx};
        results.A_standard_cm2 = api_orifices.areas_cm2(idx);
    end
end

function plot_system_curve(P, op_results)
    % Generates and plots the system curve (Head vs. Flow Rate).
    
    % Create a vector of flow rates from 0 to 150% of operating flow
    Q_range = linspace(0, 1.5 * P.Q_op, 100);
    H_total = zeros(size(Q_range));
    
    % Calculate the total head required for each flow rate
    for i = 1:length(Q_range)
        if Q_range(i) > 0
            temp_results = calculate_hydraulics(Q_range(i), P);
            % System curve head = friction head + minor loss head + static head
            H_total(i) = (temp_results.dP_major + temp_results.dP_minor) / (P.rho * P.g) + P.delta_z;
        else
            H_total(i) = P.delta_z; % At zero flow, only static head remains
        end
    end
    
    % Create the plot
    figure('Name', 'Pipeline System Curve');
    plot(Q_range, H_total, 'b-', 'LineWidth', 2, 'DisplayName', 'System Curve');
    hold on;
    
    % Mark the operating point
    H_op = (op_results.dP_major + op_results.dP_minor) / (P.rho * P.g) + P.delta_z;
    plot(P.Q_op, H_op, 'r*', 'MarkerSize', 10, 'LineWidth', 2, ...
        'DisplayName', sprintf('Operating Point (%.2f m)', H_op));
    
    grid on;
    title('Pipeline System Head Curve');
    xlabel('Volumetric Flow Rate (m^3/s)');
    ylabel('Total System Head (m)');
    legend('show', 'Location', 'northwest');
    xlim([0, 1.5 * P.Q_op]);
    hold off;

end

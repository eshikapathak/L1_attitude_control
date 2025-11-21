%% Master Experiment Script: L1 Adaptive Control + LQR
% Supports: LQR Baseline, Reference Tracking, Ts Sweep, Delay Sweep
clear; clc; close all;

%% ================= 1. GLOBAL INITIALIZATION =================
model_name = 'l1_lqr'; 

% --- Geometry & Physics ---
P.m  = 0.76;      
P.L  = 0.14;      
P.Ix = 0.0045; P.Iy = 0.0045; P.Iz = 0.0113;

% --- Linear Model (FORCE BASED) ---
P.A_true = [zeros(3,3), eye(3); zeros(3,6)];
% Bm defined for inputs: [Force_Roll, Force_Pitch, Torque_Yaw]
P.Bm = [zeros(3,3); 
        P.L/P.Ix,  0,        0;      
        0,         P.L/P.Iy, 0;       
        0,         0,        1/P.Iz]; 
P.Bum    = [eye(3); zeros(3,3)];                         
P.Bbar   = [P.Bm P.Bum];  

% --- Default Simulation Settings ---
T_final     = 10;
P.loop_delay = 0;         % Default: No delay
P.ref_amp    = deg2rad(20); 
P.ref_start  = 0.0;       

% --- Disturbance (Standard) ---
P.dist_amp  = [1.0; 0.5; 0.2]; % Nm [Roll, Pitch, Yaw]
P.dist_freq = [2.0; 1.5; 0.5]; % rad/s

% --- LQR Baseline Tuning ---
% High gains enabled by Force-based Bm
Q = diag([100, 100, 80, 1, 1, 1])/2; 
R = diag([0.01, 0.01, 0.01])*10;    
P.K = lqr(P.A_true, P.Bm, Q, R);

% --- L1 Adaptive Defaults ---
P.Ts = 0.002;              % 2ms default
P.Ae = -10.0 * eye(6);      
P.wc = 40;                 
[P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3));

% Function to update L1 Matrices (Must call this if Ts changes)
update_L1 = @(p) struct('Phi_inv', inv((expm(p.Ae*p.Ts) - eye(6))/p.Ae), ...
                        'Exp_AeTs', expm(p.Ae*p.Ts));

% Initialize L1 matrices once
L1_Mats = update_L1(P);
P.Phi_inv = L1_Mats.Phi_inv;
P.Exp_AeTs = L1_Mats.Exp_AeTs;

fprintf('Standard Parameters Loaded.\n');
fprintf('--------------------------------------------------\n');
fprintf('SELECT EXPERIMENT:\n');
fprintf('  1. LQR Baseline vs L1+LQR (Disturbance Rejection)\n');
fprintf('  2. Reference Tracking Performance (Varying Amplitude)\n');
fprintf('  3. Effect of Sampling Time (Ts Sweep)\n');
fprintf('  4. Robustness to Time Delay (Delay Sweep)\n');
fprintf('--------------------------------------------------\n');
choice = input('Enter choice (1-4): ');

if isempty(choice), choice = 1; end

%% ================= 2. EXPERIMENT LOGIC =================

switch choice
    
    case 1 % --- LQR vs L1 ---
        fprintf('Running Exp 1: LQR vs L1+LQR...\n');
        
        % Run 1: LQR Only (Disable L1 by zeroing Filter Output)
        P_lqr = P;
        P_lqr.Clpf = zeros(3); % L1 control signal becomes 0
        assignin('base', 'P', P_lqr);
        try
            out1 = sim(model_name, 'StopTime', num2str(T_final));
            res_lqr = extract_sim_data(out1, P_lqr);
        catch ME
             fprintf('\nERROR during Run 1 (LQR Only).\n');
             rethrow(ME);
        end

        % Run 2: L1 + LQR (Enable L1)
        P_l1 = P;
        P_l1.Clpf = eye(3); % Normal L1
        assignin('base', 'P', P_l1);
        try
            out2 = sim(model_name, 'StopTime', num2str(T_final));
            res_l1 = extract_sim_data(out2, P_l1);
        catch ME
             fprintf('\nERROR during Run 2 (L1 + LQR).\n');
             rethrow(ME);
        end
        
        % Plot Comparison
        plot_compare_runs(res_lqr, res_l1, 'LQR Only', 'L1 + LQR');
        
    case 2 % --- Reference Tracking ---
        fprintf('Running Exp 2: Reference Tracking...\n');
        amps_deg = [10, 30, 60];
        results = {};
        
        % Setup Subplots for Tracking AND Torque
        fig_handle = figure('Color','w','Name','Reference Tracking', 'Position', [100, 100, 800, 600]);
        
        ax1 = subplot(2,1,1); hold on; grid on; 
        title('Roll Response', 'Interpreter', 'latex', 'FontSize', 14); 
        ylabel('Roll $\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
        
        ax2 = subplot(2,1,2); hold on; grid on; 
        title('Control Torque $\tau_\phi$ (Total)', 'Interpreter', 'latex', 'FontSize', 14); 
        ylabel('Torque (Nm)', 'Interpreter', 'latex', 'FontSize', 12); 
        xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
        
        colors = {'b', 'g', 'r'};
        
        for i = 1:length(amps_deg)
            P.ref_amp = deg2rad(amps_deg(i));
            assignin('base', 'P', P);
            
            out = sim(model_name, 'StopTime', num2str(T_final));
            res = extract_sim_data(out, P);
            
            % Plot Angle on Top
            plot(ax1, res.t, rad2deg(res.x(:,1)), 'LineWidth', 2, 'Color', colors{i}, ...
                'DisplayName', sprintf('Ref = %d deg', amps_deg(i)));
            plot(ax1, res.t, rad2deg(res.ref(:,1)), '--', 'Color', colors{i}, 'HandleVisibility','off');
            
            % Plot Torque on Bottom
            plot(ax2, res.t, res.tau_tot(:,1), 'LineWidth', 1.5, 'Color', colors{i}, ...
                'DisplayName', sprintf('Ref %d deg', amps_deg(i)));
        end
        legend(ax1, 'Location', 'best', 'Interpreter', 'latex');
        legend(ax2, 'Location', 'best', 'Interpreter', 'latex');
        
    case 3 % --- Ts Sweep ---
        fprintf('Running Exp 3: Ts Sweep...\n');
        Ts_vals = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05];
        rmse_log = [];
        
        for i = 1:length(Ts_vals)
            P.Ts = Ts_vals(i);
            
            % CRITICAL: Recompute L1 Matrices for new Ts
            mats = update_L1(P);
            P.Phi_inv = mats.Phi_inv;
            P.Exp_AeTs = mats.Exp_AeTs;
            
            assignin('base', 'P', P);
            fprintf('  Testing Ts = %.4f s ... ', P.Ts);
            
            try
                out = sim(model_name, 'StopTime', num2str(T_final));
                res = extract_sim_data(out, P);
                
                err = res.ref(:,1) - res.x(:,1);
                % RMSE in DEGREES
                rmse_rad = sqrt(mean(err.^2));
                rmse_deg = rad2deg(rmse_rad);
                
                if rmse_deg > 500, rmse_deg = NaN; end % Filter instability
                fprintf('RMSE: %.4f deg\n', rmse_deg);
            catch
                rmse_deg = NaN; fprintf('Failed.\n');
            end
            rmse_log(i) = rmse_deg;
        end
        
        figure('Color','w','Name','Ts Sweep');
        semilogx(Ts_vals, rmse_log, '-bo', 'LineWidth',2, 'MarkerFaceColor','b');
        set(gca, 'XDir', 'reverse'); % <--- FLIP X-AXIS HERE
        grid on; 
        xlabel('Sampling Time $T_s$ (s) [Decreasing $\rightarrow$]', 'Interpreter', 'latex', 'FontSize', 12); 
        ylabel('RMSE (deg)', 'Interpreter', 'latex', 'FontSize', 12);
        title('Effect of Sampling Rate on Tracking Accuracy', 'Interpreter', 'latex', 'FontSize', 14);
        
    case 4 % --- Delay Sweep ---
        fprintf('Running Exp 4: Time Delay Sweep...\n');
        
        % === MODIFIED DELAY SWEEP ===
        % Testing smaller increments to find exact failure point
        delays = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.010];
        rmse_log = [];
        
        for i = 1:length(delays)
            P.loop_delay = delays(i);
            assignin('base', 'P', P);
            fprintf('  Testing Delay = %.3f s ... ', P.loop_delay);
            
            try
                out = sim(model_name, 'StopTime', num2str(T_final));
                res = extract_sim_data(out, P);
                
                err = res.ref(:,1) - res.x(:,1);
                % RMSE in DEGREES
                rmse_rad = sqrt(mean(err.^2));
                rmse_deg = rad2deg(rmse_rad);
                
                if rmse_deg > 1000, rmse_deg = NaN; end % Unstable
                fprintf('RMSE: %.4f deg\n', rmse_deg);
            catch ME
                rmse_deg = NaN; 
                fprintf('Unstable (Crashed).\n');
                % We catch the crash so the script continues to plot valid points
            end
            rmse_log(i) = rmse_deg;
        end
        
        figure('Color','w','Name','Delay Sweep');
        plot(delays*1000, rmse_log, '-ro', 'LineWidth',2, 'MarkerFaceColor','r');
        grid on; 
        xlabel('Loop Delay (ms)', 'Interpreter', 'latex', 'FontSize', 12); 
        ylabel('RMSE (deg)', 'Interpreter', 'latex', 'FontSize', 12);
        title('Performance Degradation vs Time Delay', 'Interpreter', 'latex', 'FontSize', 14);
end

%% ================= 3. HELPER FUNCTIONS =================

function res = extract_sim_data(out, P)
    % Helper Function to Extract & Fix Data Orientation
    % This ensures data is always (Time x Channels)
    function data_fixed = fix_data(sim_data, t)
        d = squeeze(sim_data); % Remove singleton dims
        
        % If dimensions match Time length in the first dim, it's already correct
        if size(d, 1) == length(t)
            data_fixed = d;
        % If the second dim matches Time, we need to transpose
        elseif size(d, 2) == length(t)
            data_fixed = d.';
        else
            error('Data dimension does not match time vector length!');
        end
    end

    try
        % We first attempt to retrieve variables to trigger error if missing
        check_x = out.sim_x; 
        
        res.t = out.sim_x.Time;
        
        % Apply the fix function to each signal
        res.x      = fix_data(out.sim_x.Data, res.t);       
        res.ref    = fix_data(out.sim_ref.Data, res.t);     
        res.u_b    = fix_data(out.sim_u_base.Data, res.t);  
        res.u_l1   = fix_data(out.sim_u_L1.Data, res.t);    
        res.dist   = fix_data(out.sim_dist.Data, res.t);    
        
        % Optional: Sigma
        if isfield(out, 'sim_sigma')
             res.sigma = fix_data(out.sim_sigma.Data, res.t);
        end
        
        % Calculate Torques (Nm)
        res.tau_b  = [res.u_b(:,1)*P.L,  res.u_b(:,2)*P.L,  res.u_b(:,3)];
        res.tau_l1 = [res.u_l1(:,1)*P.L, res.u_l1(:,2)*P.L, res.u_l1(:,3)];
        res.tau_tot = res.tau_b + res.tau_l1;
        
    catch ME
        fprintf('\n!!! DATA EXTRACTION FAILED !!!\n');
        fprintf('The script could not find the required variables in "out".\n');
        fprintf('Current variables in "out" are:\n');
        disp(out);
        fprintf('Please rename your "To Workspace" blocks in Simulink to match: sim_x, sim_ref, sim_u_base, sim_u_L1, sim_dist\n');
        rethrow(ME);
    end
end

function plot_compare_runs(res1, res2, name1, name2)
    fig_size = [100, 100, 1000, 800];
    figure('Color','w','Name','LQR vs L1', 'Position', fig_size);
    
    % 1. Tracking
    subplot(3,1,1); hold on; grid on;
    plot(res1.t, rad2deg(res1.x(:,1)), 'g--', 'LineWidth', 2, 'DisplayName', name1);
    plot(res2.t, rad2deg(res2.x(:,1)), 'b-', 'LineWidth', 1.5, 'DisplayName', name2);
    plot(res1.t, rad2deg(res1.ref(:,1)), 'k--', 'DisplayName', 'Ref');
    ylabel('Roll $\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 
    title('Attitude Tracking Comparison', 'Interpreter', 'latex', 'FontSize', 14); 
    legend('Location','best', 'Interpreter', 'latex');
    
    % 2. Total Torque & Breakdown (The 4 Lines you requested)
    subplot(3,1,2); hold on; grid on;
    plot(res1.t, res1.tau_tot(:,1), 'g--', 'LineWidth', 2, 'DisplayName', [name1 ' Total']); % 1. LQR Total
    plot(res2.t, res2.tau_tot(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', [name2 ' Total']); % 2. L1 Total
    plot(res2.t, res2.tau_l1(:,1), 'r-.', 'LineWidth', 1.0, 'DisplayName', [name2 ' Adaptive']); % 3. L1 Adaptive Component
    plot(res1.t, res1.dist(:,1), 'm:', 'LineWidth', 2, 'DisplayName', 'Disturbance'); % 4. Disturbance
    
    ylabel('$\tau_\phi$ (Nm)', 'Interpreter', 'latex', 'FontSize', 12); 
    title('Control Effort Breakdown', 'Interpreter', 'latex', 'FontSize', 14); 
    legend('Location','best', 'Interpreter', 'latex');
    
    % 3. Error (Converted to Degrees)
    subplot(3,1,3); hold on; grid on;
    err1 = abs(res1.ref(:,1) - res1.x(:,1));
    err2 = abs(res2.ref(:,1) - res2.x(:,1));
    plot(res1.t, rad2deg(err1), 'g--', 'LineWidth', 1.5, 'DisplayName', [name1 ' Error']);
    plot(res2.t, rad2deg(err2), 'b-', 'LineWidth', 1.5, 'DisplayName', [name2 ' Error']);
    ylabel('$|$Error$|$ (deg)', 'Interpreter', 'latex', 'FontSize', 12); 
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); 
    title('Absolute Tracking Error', 'Interpreter', 'latex', 'FontSize', 14); 
    legend('Location','best', 'Interpreter', 'latex');
end
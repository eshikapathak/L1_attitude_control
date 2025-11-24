%% Master Experiment Script: L1 Adaptive Control + LQR + MRAC
clear; clc; close all;

%% ================= 1. GLOBAL INITIALIZATION =================
model_name = 'l1_lqr'; 

% --- Geometry & Physics ---
P.m  = 0.76; P.L  = 0.14;      
P.Ix = 0.0045; P.Iy = 0.0045; P.Iz = 0.0113;

% --- Linear Model (FORCE BASED) ---
P.A_true = [zeros(3,3), eye(3); zeros(3,6)];
P.Bm = [zeros(3,3); 
        P.L/P.Ix,  0,        0;      
        0,         P.L/P.Iy, 0;       
        0,         0,        1/P.Iz]; 
P.Bum    = [eye(3); zeros(3,3)];                         
P.Bbar   = [P.Bm P.Bum];  

% --- Default Simulation Settings ---
T_final     = 20;
P.loop_delay = 0;      
P.ref_amp    = deg2rad(20); 
P.ref_start  = 0.0;       
P.use_mrac   = 0;         

% --- Disturbance ---
P.dist_amp  = [0.5; 0.5; 0.2]; 
P.dist_freq = [1.0; 1.5; 1.0]; 

% --- LQR Baseline Tuning ---
Q = diag([100, 100, 80, 1, 1, 1]); 
R = diag([0.01, 0.01, 0.01]);
% Q = diag([10, 10, 10,  4, 4, 3]);   % [phi theta psi dphi dtheta dpsi]
% R = diag([2, 2, 2]);
P.K = lqr(P.A_true, P.Bm, Q, R);

% --- MRAC INITIALIZATION (Base) ---
P.MRAC.Am = P.A_true - P.Bm * P.K;
P.MRAC.Bm = P.Bm * P.K; 
Q_lyap = diag([100 100 100 175 175 125]);
% eye(6) * 50; 
P.MRAC.P_lyap = lyap(P.MRAC.Am', Q_lyap);

% Default Gains (Will be overwritten by Exp 5 logic)
P.MRAC.Gam_x = 100; P.MRAC.Gam_r = 100; P.MRAC.Gam_w = 100; 
P.MRAC.theta_x_max = 1.0; P.MRAC.theta_r_max = 1.0; P.MRAC.theta_w_max = 100.0;

% --- L1 Adaptive Defaults ---
P.Ts = 0.002;              
P.Ae = -1.0 * eye(6);      
P.wc = 40;                 
[P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3));
update_L1 = @(p) struct('Phi_inv', inv((expm(p.Ae*p.Ts) - eye(6))/p.Ae), 'Exp_AeTs', expm(p.Ae*p.Ts));
L1_Mats = update_L1(P);
P.Phi_inv = L1_Mats.Phi_inv;
P.Exp_AeTs = L1_Mats.Exp_AeTs;

fprintf('Parameters Loaded.\n');
fprintf('--------------------------------------------------\n');
fprintf('SELECT EXPERIMENT:\n');
fprintf('  1. LQR vs L1+LQR (Disturbance Rejection)\n');
fprintf('  2. Reference Tracking Performance\n');
fprintf('  3. Effect of Sampling Time (Ts Sweep)\n');
fprintf('  4. Robustness to Time Delay (Delay Sweep)\n');
fprintf('  5. MRAC Performance vs L1 (With Gain Tuning)\n');
fprintf('--------------------------------------------------\n');
choice = input('Enter choice (1-5): ');
if isempty(choice), choice = 1; end

%% ================= 2. EXPERIMENT LOGIC =================

switch choice
    case 1 % --- LQR vs L1 ---
        P.use_mrac = 0; 
        
        P_lqr = P; P_lqr.Clpf = zeros(3); 
        assignin('base', 'P', P_lqr);
        out1 = sim(model_name, 'StopTime', num2str(T_final));
        res_lqr = extract_sim_data(out1, P_lqr);
        
        P_l1 = P; P_l1.Clpf = eye(3); 
        assignin('base', 'P', P_l1);
        out2 = sim(model_name, 'StopTime', num2str(T_final));
        res_l1 = extract_sim_data(out2, P_l1);
        
        plot_compare_runs(res_lqr, res_l1, 'LQR Only', 'L1 + LQR');

    case 2 % --- Reference Tracking ---
        P.use_mrac = 0;
        amps_deg = [10, 30, 60];
        
        figure('Color','w','Name','Tracking');
        ax1 = subplot(2,1,1); hold on; grid on; title('Roll Response', 'FontSize', 14); ylabel('Deg');
        ax2 = subplot(2,1,2); hold on; grid on; title('Torque', 'FontSize', 14); ylabel('Nm');
        colors = {'b', [0.85 0.33 0.10], [0.49 0.18 0.56]}; 
        
        for i = 1:length(amps_deg)
            P.ref_amp = deg2rad(amps_deg(i));
            assignin('base', 'P', P);
            out = sim(model_name, 'StopTime', num2str(T_final));
            res = extract_sim_data(out, P);
            plot(ax1, res.t, rad2deg(res.x(:,1)), 'Color', colors{i}, 'LineWidth', 2);
            plot(ax1, res.t, rad2deg(res.ref(:,1)), '--', 'Color', colors{i}, 'LineWidth', 1.5);
            plot(ax2, res.t, res.tau_tot(:,1), 'Color', colors{i}, 'LineWidth', 2);
        end
        legend(ax1, 'Ref=10', 'Ref=10', 'Ref=30', 'Ref=30', 'Ref=60', 'Ref=60');

    case 3 % --- Ts Sweep ---
        fprintf('Running Ts Sweep...\n');
        Ts_vals = [0.0005, 0.001, 0.002, 0.005, 0.01];
        rmse_log = [];
        for i = 1:length(Ts_vals)
            P.Ts = Ts_vals(i);
            mats = update_L1(P); P.Phi_inv=mats.Phi_inv; P.Exp_AeTs=mats.Exp_AeTs;
            assignin('base', 'P', P);
            out = sim(model_name, 'StopTime', num2str(T_final));
            res = extract_sim_data(out, P);
            err = res.ref(:,1) - res.x(:,1);
            rmse_log(i) = rad2deg(sqrt(mean(err.^2)));
            fprintf('Ts=%.4f, RMSE=%.4f deg\n', P.Ts, rmse_log(i));
        end
        figure; semilogx(Ts_vals, rmse_log, '-o', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b'); 
        set(gca, 'XDir', 'reverse');
        xlabel('Ts (s)'); ylabel('RMSE (deg)'); title('Ts Sweep');

    case 4 % --- Delay Sweep ---
        fprintf('Running Delay Sweep...\n');
        delays = [0, 0.002, 0.004, 0.006, 0.008, 0.010];
        rmse_log = [];
        for i = 1:length(delays)
            P.loop_delay = delays(i);
            assignin('base', 'P', P);
            try
                out = sim(model_name, 'StopTime', num2str(T_final));
                res = extract_sim_data(out, P);
                err = res.ref(:,1) - res.x(:,1);
                val = rad2deg(sqrt(mean(err.^2)));
                if val > 1000, val=NaN; end
            catch
                val = NaN;
            end
            rmse_log(i) = val;
            fprintf('Delay=%.3f, RMSE=%.4f deg\n', P.loop_delay, val);
        end
        figure; plot(delays*1000, rmse_log, '-o', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r'); 
        xlabel('Delay (ms)'); ylabel('RMSE (deg)'); title('Delay Sweep');

    case 5 % --- MRAC vs L1 ---
        fprintf('Choose MRAC Learning Mode:\n');
        fprintf('  1: Conservative (Gamma = 2.0) -> Smooth Control\n');
        fprintf('  2: Aggressive (Gamma = 10) -> High Freq Oscillations (Slide 11)\n');
        mrac_mode = input('Choice (1/2): ');
        
        if mrac_mode == 2
            P.MRAC.Gam_x = 10.0; P.MRAC.Gam_r = 10.0; P.MRAC.Gam_w = 10.0;
            P.MRAC.theta_x_max = 5.0; % Relax bounds for high gain
            fprintf('>> Running AGGRESSIVE MRAC (Expect Oscillations)...\n');
        else
            % P.MRAC.Gam_x = 2.0; P.MRAC.Gam_r = 2.0; P.MRAC.Gam_w = 5.0;
            fprintf('>> Running CONSERVATIVE MRAC (Smooth)...\n');
        end
        
        % Run 1: L1 Adaptive (Always same reference)
        P_l1 = P; P_l1.use_mrac = 0; P_l1.Clpf = eye(3);
        assignin('base', 'P', P_l1);
        try
            out1 = sim(model_name, 'StopTime', num2str(T_final));
            res_l1 = extract_sim_data(out1, P_l1);
        catch ME
            disp('L1 Simulation Crashed!'); rethrow(ME);
        end
        
        % Run 2: MRAC (With chosen gains)
        P_mrac = P; P_mrac.use_mrac = 1; 
        assignin('base', 'P', P_mrac);
        try
            out2 = sim(model_name, 'StopTime', num2str(T_final));
            res_mrac = extract_sim_data(out2, P_mrac);
        catch ME
            disp('MRAC Simulation Crashed!'); rethrow(ME);
        end
        
        % Plot
        plot_compare_runs(res_l1, res_mrac, 'L1 Adaptive', 'MRAC');
        
    otherwise
        fprintf('Selected experiment logic is in previous versions.\n');
end

%% ================= 3. HELPER FUNCTIONS =================
function res = extract_sim_data(out, P)
    function data_fixed = fix_dim_robust(sim_data, t)
        d = squeeze(sim_data);
        if size(d,1) == length(t)
            data_fixed = d;
        elseif size(d,2) == length(t)
            data_fixed = d.';
        else
            error('Dimension mismatch. Length T=%d, Data size=%d', length(t), size(d,1));
        end
    end

    try
        res.t = out.sim_x.Time;
        res.x = fix_dim_robust(out.sim_x.Data, res.t);
        res.ref = fix_dim_robust(out.sim_ref.Data, res.t);
        res.dist = fix_dim_robust(out.sim_dist.Data, res.t);
        
        % LQR/L1 Signals
        try
            res.u_b = fix_dim_robust(out.sim_u_base.Data, res.t);
            res.tau_b  = [res.u_b(:,1)*P.L,  res.u_b(:,2)*P.L,  res.u_b(:,3)];
        catch
            res.tau_b = zeros(length(res.t), 3);
        end
        
        try
            res.u_l1 = fix_dim_robust(out.sim_u_L1.Data, res.t);
            res.tau_l1 = [res.u_l1(:,1)*P.L, res.u_l1(:,2)*P.L, res.u_l1(:,3)];
        catch
            res.tau_l1 = zeros(length(res.t), 3);
        end

        % MRAC Signals
        res.tau_mrac = zeros(length(res.t), 3);
        
        if isfield(P, 'use_mrac') && P.use_mrac == 1
            try 
                 u_mrac_raw = fix_dim_robust(out.sim_u_mrac.Data, res.t);
                 res.tau_mrac = [u_mrac_raw(:,1)*P.L, u_mrac_raw(:,2)*P.L, u_mrac_raw(:,3)];
                 
                 if max(abs(res.tau_mrac(:,1))) < 1e-6 && max(abs(res.tau_b(:,1))) > 1e-6
                     fprintf('WARNING: sim_u_mrac is 0, but sim_u_base has data. Using base as fallback.\n');
                     res.tau_mrac = res.tau_b;
                 end
            catch
                fprintf('\n!!! WARNING: sim_u_mrac NOT FOUND !!!\n');
                if max(abs(res.tau_b(:,1))) > 1e-6
                    res.tau_mrac = res.tau_b;
                end
            end
            res.tau_tot = res.tau_mrac;
            fprintf('MRAC Mode: Max Torque = %.4f Nm\n', max(abs(res.tau_tot(:,1))));
        else
            res.tau_tot = res.tau_b + res.tau_l1;
        end
        
    catch ME
        fprintf('\n!!! DATA EXTRACTION FAILED !!!\n');
        disp(out);
        rethrow(ME);
    end
end

function plot_compare_runs(res1, res2, name1, name2)
    fig_size = [100, 100, 1000, 800];
    figure('Color','w','Name',[name1 ' vs ' name2], 'Position', fig_size);
    
    c1 = [0.85 0.33 0.10]; % Dark Orange (L1)
    c2 = [0 0.45 0.74];    % Standard Blue (MRAC)
    
    % 1. TRACKING
    subplot(3,1,1); hold on; grid on;
    plot(res1.t, rad2deg(res1.x(:,1)), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', name1);
    plot(res2.t, rad2deg(res2.x(:,1)), '-', 'Color', c2, 'LineWidth', 2, 'DisplayName', name2);
    plot(res1.t, rad2deg(res1.ref(:,1)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ref');
    ylabel('Roll (deg)', 'FontSize', 12); title('Tracking Comparison', 'FontSize', 14); legend('Location','best');
    
    % 2. CONTROL EFFORT (Including L1 breakdown)
    subplot(3,1,2); hold on; grid on;
    % L1 Total
    plot(res1.t, res1.tau_tot(:,1), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', [name1 ' Total']);
    % L1 Base (Visual aid)
    if isfield(res1, 'tau_b')
        plot(res1.t, res1.tau_b(:,1), ':', 'Color', c1, 'LineWidth', 1.0, 'DisplayName', [name1 ' Base Only']);
    end
    % MRAC
    plot(res2.t, res2.tau_tot(:,1), '-', 'Color', c2, 'LineWidth', 2, 'DisplayName', name2);
    % Disturbance
    plot(res1.t, res1.dist(:,1), 'm:', 'LineWidth', 2, 'DisplayName', 'Disturbance');
    ylabel('Torque (Nm)', 'FontSize', 12); title('Control Effort Breakdown', 'FontSize', 14); legend('Location','best');
    
    % 3. ERROR
    subplot(3,1,3); hold on; grid on;
    plot(res1.t, rad2deg(abs(res1.ref(:,1)-res1.x(:,1))), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', name1);
    plot(res2.t, rad2deg(abs(res2.ref(:,1)-res2.x(:,1))), '-', 'Color', c2, 'LineWidth', 2, 'DisplayName', name2);
    ylabel('|Error| (deg)', 'FontSize', 12); title('Error', 'FontSize', 14); legend('Location','best');
end

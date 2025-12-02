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
T_final     = 10;
P.loop_delay = 1e-7;      
P.ref_amp    = deg2rad(30); 
P.ref_start  = 0.0; 
P.ref_type   = 0;         % 0=Step, 1=Sine
P.ref_freq   = 1.0;       % Frequency for Sine Reference (Used in Exp 10)
P.use_mrac   = 0;         
% --- Disturbance ---
P.dist_amp  = [1.0; 0.5; 0.2]; 
P.dist_freq = [2.0; 1.5; 0.5]; 
P.dist_type = 0;          % 0=Sine, 1=SumSines, 2=Triangle, 3=StateDep
% --- LQR Baseline Tuning ---
Q = diag([100, 100, 80, 1, 1, 1]); 
R = diag([0.01, 0.01, 0.01]);    
P.K = lqr(P.A_true, P.Bm, Q, R);
% --- MRAC INITIALIZATION ---
Q_att = diag([100, 100, 80, 1, 1, 1]);   % [phi theta psi dphi dtheta dpsi]
R_att = diag([0.01, 0.01, 0.01])*100;
P.K_att = lqr(P.A_true, P.Bm, Q_att, R_att);
P.MRAC.Am = P.A_true - P.Bm * P.K_att;
P.MRAC.Bm = P.Bm * P.K_att; 
Q_lyap = diag([100 100 100 175 175 125]);
P.MRAC.P_lyap = lyap(P.MRAC.Am', Q_lyap);
% Gains & Bounds
P.MRAC.Gam_x = 5.0; P.MRAC.Gam_r = 5.0; P.MRAC.Gam_w = 5.0; 
P.MRAC.theta_x_max = 100.0; P.MRAC.theta_r_max = 100.0; P.MRAC.theta_w_max = 100.0;
% --- L1 Adaptive Defaults ---
P.Ts = 0.002;              
P.Ae = -10.0 * eye(6);      
P.wc = 40;                 
[P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3));
update_L1 = @(p) struct('Phi_inv', inv((expm(p.Ae*p.Ts) - eye(6))/p.Ae), 'Exp_AeTs', expm(p.Ae*p.Ts));
L1_Mats = update_L1(P);
P.Phi_inv = L1_Mats.Phi_inv;
P.Exp_AeTs = L1_Mats.Exp_AeTs;
fprintf('Parameters Loaded.\n');
fprintf('--------------------------------------------------\n');
fprintf('SELECT EXPERIMENT:\n');
fprintf('  1. Disturbance Rejection (LQR vs L1 vs MRAC)\n');
fprintf('  2. Tracking Performance (Step/Sine: LQR vs L1 vs MRAC)\n');
fprintf('  3. Effect of Sampling Time (Ts Sweep)\n');
fprintf('  4. Robustness to Time Delay (Delay Sweep)\n');
fprintf('  5. MRAC Tuning Mode (Conservative vs Aggressive)\n');
fprintf('  6. MRAC Gain Sweep (RMSE vs Gamma)\n');
fprintf('  7. Filter BW Sweep (Dist Freq Variation) [Hz]\n');
fprintf('  8. Filter BW Sweep (2nd Order LPF) [Hz]\n');
fprintf('  9. Filter BW Sweep (Dist Amplitude Variation) [Hz]\n');
fprintf('  10. Filter BW Sweep (Reference Freq Variation) [Hz]\n');
fprintf('  11. Control Signal Analysis (Fixed Dist, Varying Filter BW) [Hz]\n');
fprintf('  12. Disturbance Type: Sum of High Freq Sines (BW Sweep & Signals)\n');
fprintf('  13. Disturbance Type: Triangle Wave (BW Sweep & Signals)\n');
fprintf('  14. LQR Failure Case: State-Dependent Disturbance (L1 ON)\n');
fprintf('  15. State-Dependent Disturbance: BW Sweep & Signal Breakdown\n');
% fprintf('  15. LQR Stress Test: High Amp/Resonant Disturbance (BW Sweep)\n');
fprintf('  16. Time Delay Margin vs Sampling Time Ts (Log Scale)\n');
fprintf('--------------------------------------------------\n');
choice = input('Enter choice (1-17): ');
if isempty(choice), choice = 1; end
%% ================= 2. EXPERIMENT LOGIC =================
switch choice
    % ... Cases 1-16 unchanged ...
    case 1 % --- EXP 1: DISTURBANCE REJECTION (3-Way) ---
        fprintf('Running Exp 1: Disturbance Rejection Comparison...\n');
        P.ref_type = 0; P.ref_start= 0;
        P_run = P; P_run.use_mrac = 0; P_run.Clpf = zeros(3); assignin('base', 'P', P_run);
        res_lqr = run_sim(model_name, T_final, P_run);
        P_run = P; P_run.use_mrac = 0; P_run.Clpf = eye(3); assignin('base', 'P', P_run);
        res_l1 = run_sim(model_name, T_final, P_run);
        P.MRAC.Gam_x = 50.0; P.MRAC.Gam_r = 50.0; P.MRAC.Gam_w = 50.0; 
        P_run = P; P_run.use_mrac = 1; assignin('base', 'P', P_run);
        res_mrac = run_sim(model_name, T_final, P_run);
        ref_str = sprintf('Reference: Step %.0f deg', rad2deg(P.ref_amp));
        dist_str = sprintf('Disturbance: $\\tau_\\phi = %.1f \\sin(%.1f t)$ Nm', P.dist_amp(1), P.dist_freq(1));
        plot_3way_comparison(res_lqr, res_l1, res_mrac, 'Disturbance Rejection Comparison', [ref_str ' | ' dist_str]);
    case 2 % --- EXP 2: TRACKING PERFORMANCE (3-Way x 3 Refs) ---
        fprintf('Running Exp 2: Tracking Performance (Single Figure)...\n');
        cases(1).type=0; cases(1).amp=deg2rad(10); cases(1).start=0; cases(1).name='Step 10 deg';
        cases(2).type=0; cases(2).amp=deg2rad(40); cases(2).start=0; cases(2).name='Step 40 deg';
        cases(3).type=1; cases(3).amp=deg2rad(30); cases(3).start=0; cases(3).name='Sine Wave'; 
        num_cases = length(cases);
        figure('Name', 'Experiment 2: Reference Tracking Comparison', 'Color', 'w', 'Position', [50, 50, 1400, 700]);
        tl = tiledlayout(2, num_cases, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, 'Tracking Performance Across Different References', 'FontSize', 16, 'Interpreter', 'latex');
        subtitle(tl, sprintf('Disturbance (Roll): $\\tau_\\phi = %.1f \\sin(%.1f t)$ Nm | (Applied to all cases)', P.dist_amp(1), P.dist_freq(1)), ...
            'Interpreter','latex', 'FontSize', 12);
        cLQR=[0.4 0.4 0.4]; cL1=[0.85 0.33 0.10]; cMRAC=[0 0.45 0.74];
        for i = 1:num_cases
            c = cases(i);
            fprintf('  -> Testing Reference: %s\n', c.name);
            P.ref_type = c.type; P.ref_amp  = c.amp; P.ref_start = c.start; 
            P_run = P; P_run.use_mrac=0; P_run.Clpf=zeros(3); assignin('base', 'P', P_run); r_lqr = run_sim(model_name, T_final, P_run);
            P_run = P; P_run.use_mrac=0; P_run.Clpf=eye(3); assignin('base', 'P', P_run); r_l1 = run_sim(model_name, T_final, P_run);
            P_run = P; P_run.use_mrac=1; assignin('base', 'P', P_run); r_mrac = run_sim(model_name, T_final, P_run);
            nexttile(i); hold on; grid on;
            if ~isempty(r_lqr), plot(r_lqr.t, rad2deg(r_lqr.x(:,1)), '--', 'Color', cLQR, 'LineWidth', 1.5); end
            if ~isempty(r_l1), plot(r_l1.t, rad2deg(r_l1.x(:,1)), '-.', 'Color', cL1, 'LineWidth', 2); end
            if ~isempty(r_mrac), plot(r_mrac.t, rad2deg(r_mrac.x(:,1)), '-', 'Color', cMRAC, 'LineWidth', 2); end
            if ~isempty(r_lqr), plot(r_lqr.t, rad2deg(r_lqr.ref(:,1)), 'k:', 'LineWidth', 1.5); end
            title(['\textbf{' c.name '}'], 'Interpreter', 'latex', 'FontSize', 12);
            if i==1, ylabel('$\phi$ (deg)', 'Interpreter', 'latex'); end
            if i==3, legend('LQR','L1+LQR','MRAC','Ref', 'Location','best','Interpreter','latex'); end
            nexttile(i + num_cases); hold on; grid on;
            if ~isempty(r_lqr), plot(r_lqr.t, r_lqr.tau_tot(:,1), '--', 'Color', cLQR, 'LineWidth', 1.5); end
            if ~isempty(r_l1), plot(r_l1.t, r_l1.tau_tot(:,1), '-.', 'Color', cL1, 'LineWidth', 2); end
            if ~isempty(r_mrac), plot(r_mrac.t, r_mrac.tau_tot(:,1), '-', 'Color', cMRAC, 'LineWidth', 2); end
            if ~isempty(r_lqr), plot(r_lqr.t, r_lqr.dist(:,1), 'm:', 'LineWidth', 1); end
            if i==1, ylabel('$\tau$ (Nm)', 'Interpreter', 'latex'); end
            xlabel('Time (s)', 'Interpreter', 'latex');
        end
    case 3 % --- Ts Sweep ---
        fprintf('Running Ts Sweep...\n');
        Ts_vals = [5e-5, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.1, 1.0];
        rmse_log = [];
        for i = 1:length(Ts_vals)
            P.Ts = Ts_vals(i);
            mats = update_L1(P); P.Phi_inv=mats.Phi_inv; P.Exp_AeTs=mats.Exp_AeTs;
            assignin('base', 'P', P);
            res = run_sim(model_name, T_final, P);
            if isempty(res), rmse_log(i)=NaN; continue; end
            err = res.ref(:,1) - res.x(:,1);
            rmse_log(i) = rad2deg(sqrt(mean(err.^2)));
            fprintf('Ts=%.4f, RMSE=%.4f deg\n', P.Ts, rmse_log(i));
        end
        figure('Color','w','Position',[100 100 800 600]);
        semilogx(Ts_vals, rmse_log, '-o', 'LineWidth', 2, 'Color', [0.49 0.18 0.56], 'MarkerFaceColor', [0.49 0.18 0.56]); 
        set(gca, 'XDir', 'reverse', 'TickLabelInterpreter', 'latex', 'FontSize', 12);
        xlabel('Sampling Time $T_s$ (s) (Log Scale)', 'Interpreter', 'latex'); 
        ylabel('RMSE (deg)', 'Interpreter', 'latex'); 
        subtitle('Tracking 30 deg Step', 'Interpreter', 'latex'); grid on;
    case 4 % --- Delay Sweep ---
        fprintf('Running Delay Sweep...\n');
        delays = [0, 0.002, 0.004, 0.006, 0.008, 0.010];
        rmse_log = [];
        for i = 1:length(delays)
            P.loop_delay = delays(i);
            assignin('base', 'P', P);
            res = run_sim(model_name, T_final, P);
            if isempty(res), val=NaN; else, val=rad2deg(sqrt(mean((res.ref(:,1)-res.x(:,1)).^2))); end
            if val > 1000, val=NaN; end
            rmse_log(i) = val;
            fprintf('Delay=%.3f, RMSE=%.4f deg\n', P.loop_delay, val);
        end
        figure('Color','w','Position',[100 100 800 600]);
        plot(delays*1000, rmse_log, '-o', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r'); 
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
        xlabel('Delay (ms)', 'Interpreter', 'latex'); ylabel('RMSE (deg)', 'Interpreter', 'latex'); 
        title('Time Delay Robustness', 'Interpreter', 'latex'); grid on;
    case 5 % --- MRAC Tuning ---
        fprintf('MRAC Tuning Mode:\n  1: Conservative\n  2: Aggressive\n');
        if input('Choice: ') == 2
            P.MRAC.Gam_x = 1000; P.MRAC.Gam_r = 1000; P.MRAC.Gam_w = 1000;
            P.MRAC.theta_x_max = 500; P.MRAC.theta_r_max = 500; P.MRAC.theta_w_max = 500;
            P.loop_delay = 0.004;
            fprintf('>> Aggressive Mode Active\n');
        else
            P.loop_delay = 1e-5;
            fprintf('>> Conservative Mode Active\n');
        end
        P_l1=P; P_l1.use_mrac=0; P_l1.Clpf=eye(3); assignin('base','P',P_l1); r_l1 = run_sim(model_name, T_final, P_l1);
        P_m=P; P_m.use_mrac=1; assignin('base','P',P_m); r_m = run_sim(model_name, T_final, P_m);
        ref_str = sprintf('Reference: Step %.0f deg', rad2deg(P.ref_amp));
        dist_str = sprintf('Disturbance: $\\tau_\\phi = %.1f \\sin(%.1f t)$ Nm', P.dist_amp(1), P.dist_freq(1));
        plot_3way_comparison([], r_l1, r_m, 'MRAC Tuning Check', [ref_str ' | ' dist_str]);
    case 6 % --- MRAC Gain Sweep (RMSE vs Gamma) ---
        fprintf('Running MRAC Gain Sweep (RMSE vs Gamma)...\n');
        gam_vals = [0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 2000, 10000];
        rmse_log = [];
        P.ref_type = 0; P.ref_amp  = deg2rad(30); P.ref_start = 0.0; P.use_mrac = 1; 
        P.MRAC.theta_x_max = 2000; P.MRAC.theta_r_max = 2000; P.MRAC.theta_w_max = 2000;
        for i = 1:length(gam_vals)
            g = gam_vals(i);
            P.MRAC.Gam_x = g; P.MRAC.Gam_r = g; P.MRAC.Gam_w = g;
            assignin('base', 'P', P);
            res = run_sim(model_name, T_final, P);
            if isempty(res), val = NaN; else, val = rad2deg(sqrt(mean((res.ref(:,1) - res.x(:,1)).^2))); end
            if val > 1000, val = NaN; end
            rmse_log(i) = val;
            fprintf('Gamma=%.1f, RMSE=%.4f deg\n', g, val);
        end
        figure('Color','w','Position',[100 100 800 600]);
        semilogx(gam_vals, rmse_log, '-sq', 'LineWidth', 2, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
        xlabel('Adaptation Gain $\Gamma$ (Log Scale)', 'Interpreter', 'latex');
        ylabel('RMSE (deg)', 'Interpreter', 'latex');
        subtitle('Tracking 30 deg Step', 'Interpreter', 'latex'); grid on;
        
    case 7 % --- Filter Bandwidth vs Disturbance Frequency (1st Order) ---
        fprintf('Running 1st Order Filter BW Sweep (Hz)...\n');
        dist_freqs_Hz = [0.1, 1, 2, 5, 10]; dist_freqs = dist_freqs_Hz * 2 * pi; 
        fc_vals = [0.1, 0.5, 1, 2, 5, 8, 12, 16, 20, 30, 40, 50]; 
        P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0; P.use_mrac = 0; 
        rmse_results = zeros(length(dist_freqs), length(fc_vals));
        colors = {'b', [0 0.5 0], [0.85 0.33 0.10], 'r', [0.49 0.18 0.56], [0.6 0.2 0], 'k', 'm'};
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title('Disturbance Rejection (1st Order LPF)', 'Interpreter','latex');
        subtitle(sprintf('Reference: Step 20 deg | Disturbance: Varying Frequency'), 'Interpreter','latex');
        xlabel('Filter Bandwidth $f_c$ (Hz)', 'Interpreter','latex'); ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        for i = 1:length(dist_freqs)
            d_freq = dist_freqs(i); P.dist_freq = [d_freq; d_freq; d_freq];
            for j = 1:length(fc_vals)
                f_c = fc_vals(j); w_c = f_c * 2 * pi; P.wc = w_c; 
                [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
                assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
                if isempty(res), val = NaN; else, val = rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2))); end
                rmse_results(i,j) = val;
            end
            plot(fc_vals, rmse_results(i,:), '-o', 'LineWidth', 2, 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$d_{freq} = %g$ Hz', dist_freqs_Hz(i)));
        end
        legend('Location','best', 'Interpreter','latex');
        
        FIXED_FREQ = 5 * 2 * pi; P.dist_freq = [FIXED_FREQ; FIXED_FREQ; FIXED_FREQ];
        figure('Color','w','Position',[100 50 1200 1000]); tl = tiledlayout(4,2, 'TileSpacing','compact','Padding','compact');
        title(tl, 'Control Signals: Disturbance Rejection (5 Hz)', 'Interpreter','latex', 'FontSize', 16);
        fc_plot = [1, 5, 10, 20];
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base','P',P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, 5);
        end

    case 8 % --- Filter Bandwidth vs Disturbance Frequency (2nd Order) ---
        fprintf('Running 2nd Order Filter BW Sweep (Hz)...\n');
        dist_freqs_Hz = [0.1, 1, 2, 5, 10]; dist_freqs = dist_freqs_Hz * 2 * pi;
        fc_vals = [0.1, 0.5, 1, 2, 5, 8, 12, 16, 20, 30, 40, 50]; 
        P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0; P.use_mrac = 0; 
        rmse_results = zeros(length(dist_freqs), length(fc_vals));
        colors = {'b', 'g', [0.85 0.33 0.10], 'r', 'm', 'k'};
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title('Disturbance Rejection (2nd Order LPF)', 'Interpreter','latex');
        subtitle(sprintf('Reference: Step 20 deg | Disturbance: Varying Frequency'), 'Interpreter','latex');
        xlabel('Filter Bandwidth $f_c$ (Hz)', 'Interpreter','latex'); ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        for i = 1:length(dist_freqs)
            d_freq = dist_freqs(i); P.dist_freq = [d_freq; d_freq; d_freq];
            for j = 1:length(fc_vals)
                f_c = fc_vals(j); w_c = f_c * 2 * pi;
                A_1ch = [0 1; -w_c^2 -2*w_c]; B_1ch = [0; w_c^2]; C_1ch = [1 0];
                P.Alpf = blkdiag(A_1ch, A_1ch, A_1ch); P.Blpf = blkdiag(B_1ch, B_1ch, B_1ch); 
                P.Clpf = blkdiag(C_1ch, C_1ch, C_1ch); P.Dlpf = zeros(3,3); P.wc = w_c; 
                assignin('base', 'P', P);
                try res = run_sim(model_name, T_final, P); if isempty(res), val=NaN; else, val=rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2))); end
                catch, val = NaN; end
                rmse_results(i,j) = val;
            end
            plot(fc_vals, rmse_results(i,:), '-^', 'LineWidth', 2, 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$d_{freq} = %g$ Hz', dist_freqs_Hz(i)));
        end
        legend('Location','best', 'Interpreter','latex');
        
        FIXED_FREQ = 5 * 2 * pi; P.dist_freq = [FIXED_FREQ; FIXED_FREQ; FIXED_FREQ];
        figure('Color','w','Position',[100 50 1200 1000]); tl = tiledlayout(4,2, 'TileSpacing','compact','Padding','compact');
        title(tl, 'Control Signals: 2nd Order LPF (Dist 5 Hz)', 'Interpreter','latex','FontSize',16);
        fc_plot = [1, 5, 10, 20];
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); w_c = f_c * 2 * pi;
            A_1ch = [0 1; -w_c^2 -2*w_c]; B_1ch = [0; w_c^2]; C_1ch = [1 0];
            P.Alpf = blkdiag(A_1ch, A_1ch, A_1ch); P.Blpf = blkdiag(B_1ch, B_1ch, B_1ch); 
            P.Clpf = blkdiag(C_1ch, C_1ch, C_1ch); P.Dlpf = zeros(3,3); P.wc = w_c; 
            assignin('base','P',P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, 5);
        end

    case 9 % --- Filter Bandwidth vs Disturbance Amplitude ---
        fprintf('Running Filter BW Sweep vs Dist Amplitude (Hz)...\n');
        dist_amps = [0.2, 0.5, 1.0, 2.0]; fc_vals = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50]; 
        P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0; P.use_mrac = 0; 
        FIXED_FREQ = 2.0 * 2 * pi; P.dist_freq = [FIXED_FREQ; FIXED_FREQ; FIXED_FREQ];
        rmse_results = zeros(length(dist_amps), length(fc_vals));
        colors = {'b', 'g', [0.85 0.33 0.10], 'r'}; 
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title(sprintf('Dist Rejection (Freq=2 Hz)'), 'Interpreter','latex');
        subtitle(sprintf('Reference: Step 20 deg | Disturbance: Varying Amplitude'), 'Interpreter','latex');
        xlabel('Filter Bandwidth $f_c$ (Hz)', 'Interpreter','latex'); ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        for i = 1:length(dist_amps)
            d_amp = dist_amps(i); P.dist_amp = [d_amp; d_amp; d_amp];
            for j = 1:length(fc_vals)
                f_c = fc_vals(j); w_c = f_c * 2 * pi; P.wc = w_c; 
                [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
                assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
                if isempty(res), val = NaN; else, val = rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2))); end
                rmse_results(i,j) = val;
            end
            plot(fc_vals, rmse_results(i,:), '-s', 'LineWidth', 2, 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$d_{amp} = %.1f$ Nm', d_amp));
        end
        legend('Location','best', 'Interpreter','latex');
        
        P.dist_amp = [1.0; 1.0; 1.0];
        figure('Color','w','Position',[100 50 1200 1000]); tl = tiledlayout(4,2, 'TileSpacing','compact','Padding','compact');
        title(tl, 'Control Signals: Dist Amp 1.0 Nm (Freq 2 Hz)', 'Interpreter','latex','FontSize',16);
        fc_plot = [1, 5, 10, 20];
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base','P',P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, 2.0);
        end

    case 10 % --- Filter Bandwidth vs Reference Frequency ---
        fprintf('Running Filter BW Sweep vs Reference Freq (Hz)...\n');
        ref_freqs_Hz = [0.1, 0.5, 1.0, 2.0]; ref_freqs = ref_freqs_Hz * 2 * pi;
        fc_vals = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50]; 
        P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 1; P.ref_amp = deg2rad(30); P.ref_start = 0; P.use_mrac = 0;
        P.dist_amp = [0.5; 0.5; 0.5]; P.dist_freq = [1.0; 1.0; 1.0];
        rmse_results = zeros(length(ref_freqs), length(fc_vals));
        colors = {'b', [0 0.5 0], [0.85 0.33 0.10], 'r', [0.49 0.18 0.56], 'k'};
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title('Tracking Performance (Sine Ref)', 'Interpreter','latex');
        subtitle(sprintf('Reference: Sine Amp 30 deg | Disturbance: Constant 0.5 Nm, 1 rad/s'), 'Interpreter','latex');
        xlabel('Filter Bandwidth $f_c$ (Hz)', 'Interpreter','latex'); ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        for i = 1:length(ref_freqs)
            r_freq = ref_freqs(i); P.ref_freq = r_freq;
            for j = 1:length(fc_vals)
                f_c = fc_vals(j); w_c = f_c * 2 * pi; P.wc = w_c; 
                [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
                assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
                if isempty(res), val = NaN; else, err = res.ref(:,1) - res.x(:,1); val = rad2deg(sqrt(mean(err.^2))); end
                rmse_results(i,j) = val;
            end
            plot(fc_vals, rmse_results(i,:), '-o', 'LineWidth', 2, 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$Ref_{freq} = %.1f$ Hz', ref_freqs_Hz(i)));
        end
        legend('Location','best', 'Interpreter','latex');
        
        P.ref_freq = 1.0 * 2 * pi;
        figure('Color','w','Position',[100 50 1200 1000]); tl = tiledlayout(4,2, 'TileSpacing','compact','Padding','compact');
        title(tl, 'Control Signals: Ref Freq 1 Hz', 'Interpreter','latex','FontSize',16);
        fc_plot = [1, 5, 10, 20];
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base','P',P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, nan);
        end

    case 11 % --- Control Signal Analysis ---
        fprintf('Plotting Control Signals for different Filter Bandwidths...\n');
        FIXED_D_AMP = 1.0; FIXED_D_FREQ_HZ = 5; FIXED_D_FREQ = FIXED_D_FREQ_HZ * 2 * pi;
        wc_vals_Hz = [0.5, 1, 5, 20]; 
        P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 0; P.ref_amp = 0; P.ref_start = 0; P.use_mrac = 0;
        P.dist_amp = [FIXED_D_AMP; FIXED_D_AMP; FIXED_D_AMP]; P.dist_freq = [FIXED_D_FREQ; FIXED_D_FREQ; FIXED_D_FREQ];
        figure('Color','w','Position',[100 50 1200 1000]); tl = tiledlayout(4,2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('Control Signals vs Bandwidth (Dist: Amp=%.1f, Freq=%d Hz)', FIXED_D_AMP, FIXED_D_FREQ_HZ), 'Interpreter','latex', 'FontSize', 16);
        for i = 1:length(wc_vals_Hz)
            f_c = wc_vals_Hz(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, FIXED_D_FREQ_HZ);
        end
        figure('Color','w','Position',[600 100 600 400]); hold on; grid on;
        title('Bode Magnitude Plot of L1 Low-Pass Filters', 'Interpreter','latex', 'FontSize', 14);
        xlabel('Frequency (Hz)', 'Interpreter','latex', 'FontSize', 12); ylabel('Magnitude (dB)', 'Interpreter','latex', 'FontSize', 12);
        f_bode = logspace(-1, 2, 1000); w_bode = f_bode * 2 * pi; colors_bode = {'r', 'g', 'b', 'k'}; 
        for i = 1:length(wc_vals_Hz)
            f_c = wc_vals_Hz(i); w_c = f_c * 2 * pi;
            mag = 20*log10( w_c ./ sqrt(w_bode.^2 + w_c^2) );
            semilogx(f_bode, mag, 'LineWidth', 2, 'Color', colors_bode{i}, 'DisplayName', sprintf('$f_c = %.1f$ Hz', f_c));
        end
        xline(FIXED_D_FREQ_HZ, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Dist Freq (5 Hz)');
        legend('Location','best', 'Interpreter','latex'); ylim([-40 5]);

    case 12 % --- Sum of High Freq Sines ---
        fprintf('Running Sum of Sinusoids Disturbance Sweep...\n');
        P.dist_type = 1; P.dist_amp = [1.0; 1.0; 1.0]; P.Ts = 0.002; P.Ae = -10.0 * eye(6);
        P.ref_type = 0; P.ref_amp = 0; P.ref_start = 0; P.use_mrac = 0;
        fc_vals = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50]; 
        rmse_log = [];
        for i = 1:length(fc_vals)
            f_c = fc_vals(i); w_c = f_c * 2 * pi; P.wc = w_c; 
            [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            if isempty(res), val=NaN; else, val=rad2deg(sqrt(mean(res.x(:,1).^2))); end
            rmse_log(i) = val;
        end
        fig1 = figure('Name', 'Exp 12: RMSE vs BW', 'Color', 'w'); 
        plot(fc_vals, rmse_log, '-o', 'LineWidth', 2);
        title('Disturbance: Sum of Sines'); 
        subtitle('Reference: Step 0 deg | Disturbance: Sum of Sines (20, 30, 50 rad/s)', 'Interpreter','latex');
        xlabel('Filter Bandwidth (Hz)'); ylabel('RMSE (deg)'); grid on;
        fc_plot = [1, 5, 10, 20];
        fig2 = figure('Name', 'Exp 12: Signals', 'Color','w','Position',[100 50 1200 1000]);
        tl = tiledlayout(fig2, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, 'Control Signals: Sum of Sines Disturbance', 'Interpreter','latex', 'FontSize', 16);
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, nan);
        end

    case 13 % --- Triangle Wave ---
        fprintf('Running Triangle Wave (BW in Hz)...\n');
        P.dist_type = 2; P.dist_amp = [1.0; 1.0; 1.0]; P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 0; P.ref_amp = 0; P.use_mrac = 0;
        fc_vals = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50]; 
        rmse_log = [];
        for i = 1:length(fc_vals)
            f_c = fc_vals(i); w_c = f_c * 2 * pi; P.wc = w_c; 
            [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            if isempty(res), val=NaN; else, val=rad2deg(sqrt(mean(res.x(:,1).^2))); end
            rmse_log(i) = val;
        end
        fig1 = figure('Name', 'Exp 13: RMSE vs BW', 'Color', 'w'); 
        plot(fc_vals, rmse_log, '-s', 'LineWidth', 2, 'Color', 'r');
        title('Disturbance: Triangle Wave (2 rad/s)'); 
        subtitle('Reference: Step 0 deg | Disturbance: Triangle Wave (2 rad/s, 1.0 Nm)', 'Interpreter','latex');
        xlabel('Filter Bandwidth (Hz)'); ylabel('RMSE (deg)'); grid on;
        fc_plot = [1, 5, 10, 20];
        fig2 = figure('Name', 'Exp 13: Signals', 'Color','w','Position',[100 50 1200 1000]);
        tl = tiledlayout(fig2, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, 'Control Signals: Triangle Wave Disturbance', 'Interpreter','latex', 'FontSize', 16);
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, nan);
        end

    case 14 % --- LQR Failure Case ---
        fprintf('Running LQR Failure Case...\n');
        amps = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        rmse_lqr = [];
        P.ref_type = 0; P.ref_amp = 0; P.ref_start = 0; P.dist_freq = [3.0; 3.0; 3.0]; P.dist_type = 0;
        P.use_mrac = 0; P.Clpf = eye(3); 
        for i = 1:length(amps)
            a = amps(i); P.dist_amp = [a; a; a]; assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            if isempty(res) || max(abs(res.x(:,1))) > 100, rmse_lqr(i) = NaN; else, rmse_lqr(i) = rad2deg(sqrt(mean(res.x(:,1).^2))); end
            fprintf('Amp=%.1f, RMSE=%.4f\n', a, rmse_lqr(i));
        end
        figure('Color','w'); plot(amps, rmse_lqr, '-rx', 'LineWidth', 2);
        xlabel('Disturbance Amplitude (Nm)'); ylabel('RMSE (deg)'); grid on;
        title('LQR Failure Case: High Amplitude Disturbance', 'FontSize', 14);
        subtitle('Reference: Regulator (0 deg) | Disturbance: Resonance (3 rad/s)', 'Interpreter','latex');
        
        fc_plot = [1, 5, 10, 20];
        figure('Color','w','Position',[100 50 1200 1000]); tl = tiledlayout(4,2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, 'Control Signals: LQR Failure (Amp=5.0 Nm)', 'Interpreter','latex', 'FontSize', 16);
        P.dist_amp = [5.0; 5.0; 5.0];
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, nan);
        end

    case 15 % --- State Dependent Disturbance ---
        fprintf('Running State-Dependent Analysis (BW in Hz)...\n');
        P.dist_type = 3; P.dist_amp = [2.0; 2.0; 2.0]; 
        P.Ts = 0.002; P.Ae = -10.0 * eye(6); P.ref_type = 0; P.ref_amp = deg2rad(30); P.ref_start = 0; P.use_mrac = 0;
        
        fc_vals = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 40, 50]; 
        rmse_log = [];
        for i = 1:length(fc_vals)
            f_c = fc_vals(i); w_c = f_c * 2 * pi;
            P.wc = w_c; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            if isempty(res), val=NaN; else, val=rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2))); end
            rmse_log(i) = val;
        end
        
        fig1 = figure('Name', 'Exp 15: RMSE vs BW', 'Color', 'w'); 
        plot(fc_vals, rmse_log, '-d', 'LineWidth', 2, 'Color', 'm');
        title('State-Dependent Disturbance'); 
        subtitle('Reference: Step 30 deg | Disturbance: d(x) = 2*sin(3\phi)*p', 'Interpreter','latex');
        xlabel('Filter Bandwidth (Hz)'); ylabel('RMSE (deg)'); grid on;
        
        fc_plot = [1, 5, 10, 20];
        fig2 = figure('Name', 'Exp 15: Signals', 'Color','w','Position',[100 50 1200 1000]);
        tl = tiledlayout(fig2, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, 'Control Signals: State-Dependent Disturbance', 'Interpreter','latex', 'FontSize', 16);
        for i = 1:length(fc_plot)
            f_c = fc_plot(i); P.wc = f_c*2*pi; [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3)); P.Clpf = eye(3); 
            assignin('base', 'P', P); res = run_sim(model_name, T_final, P);
            plot_control_row(res, tl, i, f_c, nan);
        end
    case 16 % --- Time Delay Margin vs Ts ---
        fprintf('Running Time Delay Margin vs Ts Sweep...\n');
        Ts_vals = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2];
        tdm_results = [];
        P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0; P.use_mrac = 0;
        
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title('Time Delay Margin vs Sampling Time', 'Interpreter','latex');
        xlabel('Sampling Time $T_s$ (s)', 'Interpreter','latex'); ylabel('Time Delay Margin (s)', 'Interpreter','latex');
        set(gca, 'XScale', 'log', 'TickLabelInterpreter','latex', 'FontSize',12);

        for i = 1:length(Ts_vals)
            ts = Ts_vals(i);
            P.Ts = ts; 
            update_L1 = @(p) struct('Phi_inv', inv((expm(p.Ae*p.Ts) - eye(6))/p.Ae), 'Exp_AeTs', expm(p.Ae*p.Ts));
            mats = update_L1(P); P.Phi_inv=mats.Phi_inv; P.Exp_AeTs=mats.Exp_AeTs;
            
            fprintf('Testing Ts = %.1e s... ', ts);
            
            % Find Max Delay
            d_min = 0; d_max = 0.1; % Search range 0 to 100ms
            stable_delay = 0;
            
            % Coarse search
            for delay = d_min:0.002:d_max
                P.loop_delay = delay;
                fprintf('Delay %.6e s', P.loop_delay);
                if P.loop_delay == 0, P.loop_delay = 1e-7; end
                assignin('base', 'P', P);
                
                try
                    out = sim(model_name, 'StopTime', '5'); % Short run is enough for stability check
                    % Check stability (norm of state)
                    if max(abs(out.sim_x.Data(:,1))) > 10 % 10 rad is definitely unstable
                        break;
                    end
                    stable_delay = delay;
                catch
                    break;
                end
            end
            tdm_results(i) = stable_delay;
            fprintf('TDM = %.3f s\n', stable_delay);
        end
        plot(Ts_vals, tdm_results, '-o', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b');
end
%% ================= 3. HELPER FUNCTIONS =================
function res = run_sim(model, T, P_struct)
    try out = sim(model, 'StopTime', num2str(T)); res = extract_sim_data(out, P_struct);
    catch ME, fprintf('Sim Failed: %s\n', ME.message); res = []; end
end
function res = extract_sim_data(out, P)
    function d_fix = fix_d(d, t)
        d = squeeze(d); if size(d,1)==length(t), d_fix=d; elseif size(d,2)==length(t), d_fix=d.'; else, error('Dim'); end
    end
    try
        res.t = out.sim_x.Time; res.x = fix_d(out.sim_x.Data, res.t); res.ref = fix_d(out.sim_ref.Data, res.t); res.dist = fix_d(out.sim_dist.Data, res.t);
        try res.u_b = fix_d(out.sim_u_base.Data, res.t); catch, res.u_b = zeros(length(res.t),3); end
        try res.u_l1 = fix_d(out.sim_u_L1.Data, res.t); catch, res.u_l1 = zeros(length(res.t),3); end
        res.tau_b = [res.u_b(:,1)*P.L, res.u_b(:,2)*P.L, res.u_b(:,3)]; res.tau_l1 = [res.u_l1(:,1)*P.L, res.u_l1(:,2)*P.L, res.u_l1(:,3)];
        if isfield(P, 'use_mrac') && P.use_mrac == 1
            try u_m = fix_d(out.sim_u_mrac.Data, res.t); res.tau_tot = [u_m(:,1)*P.L, u_m(:,2)*P.L, u_m(:,3)];
                if max(abs(res.tau_tot(:,1))) < 1e-6, res.tau_tot = res.tau_b; end 
            catch, res.tau_tot = res.tau_b; end
        else, res.tau_tot = res.tau_b + res.tau_l1; end
    catch, res = []; end
end
function plot_3way_comparison(r1, r2, r3, title_str, subtitle_str)
    fig_size = [100, 100, 1100, 450]; figure('Color','w', 'Position', fig_size); tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, title_str, 'Interpreter','latex', 'FontSize',16); 
    if nargin > 4 && ~isempty(subtitle_str), subtitle(tl, subtitle_str, 'Interpreter', 'latex', 'FontSize', 12); end
    cLQR=[0.4 0.4 0.4]; cL1=[0.85 0.33 0.10]; cMRAC=[0 0.45 0.74];
    nexttile(1); hold on; grid on;
    if ~isempty(r1), plot(r1.t, rad2deg(r1.x(:,1)), '--', 'Color', cLQR, 'LineWidth', 1.5, 'DisplayName', 'LQR'); end
    if ~isempty(r2), plot(r2.t, rad2deg(r2.x(:,1)), '-.', 'Color', cL1, 'LineWidth', 2, 'DisplayName', 'L1+LQR'); end
    if ~isempty(r3), plot(r3.t, rad2deg(r3.x(:,1)), '-', 'Color', cMRAC, 'LineWidth', 2, 'DisplayName', 'MRAC'); end
    if ~isempty(r1), plot(r1.t, rad2deg(r1.ref(:,1)), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Ref'); end
    ylabel('$\phi$ (deg)', 'Interpreter','latex', 'FontSize', 12); xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12); title('Tracking Angle', 'Interpreter','latex', 'FontSize', 14); legend('Location','best', 'Interpreter','latex');
    nexttile(2); hold on; grid on;
    if ~isempty(r1), plot(r1.t, r1.tau_tot(:,1), '--', 'Color', cLQR, 'LineWidth', 1.5, 'DisplayName', 'LQR'); end
    if ~isempty(r2), plot(r2.t, r2.tau_tot(:,1), '-.', 'Color', cL1, 'LineWidth', 2, 'DisplayName', 'L1+LQR'); end
    if ~isempty(r3), plot(r3.t, r3.tau_tot(:,1), '-', 'Color', cMRAC, 'LineWidth', 2, 'DisplayName', 'MRAC'); end
    if ~isempty(r1), plot(r1.t, r1.dist(:,1), 'm:', 'LineWidth', 1, 'DisplayName', 'Dist'); end
    ylabel('$\tau$ (Nm)', 'Interpreter','latex', 'FontSize', 12); xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12); title('Control Effort', 'Interpreter','latex', 'FontSize', 14); legend('Location','best', 'Interpreter','latex');
end
function plot_control_row(res, tl, row_idx, f_c, d_freq_hz)
    cTotal='k'; cBase='b'; cL1=[0 0.5 0]; cDist='m';
    nexttile(tl, 2*row_idx-1); hold on; grid on;
    if ~isempty(res)
        plot(res.t, res.tau_tot(:,1), 'Color', cTotal, 'LineWidth', 2, 'DisplayName', 'Total u');
        plot(res.t, res.tau_b(:,1), ':', 'Color', cBase, 'LineWidth', 1.5, 'DisplayName', 'Baseline u_{base}');
        plot(res.t, res.tau_l1(:,1), '--', 'Color', cL1, 'LineWidth', 1.5, 'DisplayName', 'Adaptive u_{l1}');
        plot(res.t, res.dist(:,1), ':', 'Color', cDist, 'LineWidth', 1.5, 'DisplayName', 'Disturbance');
    end
    title(sprintf('$f_c=%.1f$ Hz', f_c), 'Interpreter','latex'); xlim([0 5]);
    if row_idx==1, legend('Location','best', 'Interpreter','latex'); end
    nexttile(tl, 2*row_idx); hold on; grid on;
    if ~isempty(res)
        compute_fft = @(sig, Fs) abs(fft(sig)/length(sig));
        Ts_sim = mean(diff(res.t)); Fs = 1/Ts_sim; L=length(res.t); f=Fs*(0:(L/2))/L; 
        P_tot = compute_fft(res.tau_tot(:,1), Fs); P_tot=P_tot(1:floor(L/2)+1); P_tot(2:end-1)=2*P_tot(2:end-1);
        P_l1  = compute_fft(res.tau_l1(:,1), Fs); P_l1=P_l1(1:floor(L/2)+1); P_l1(2:end-1)=2*P_l1(2:end-1);
        P_base= compute_fft(res.tau_b(:,1), Fs); P_base=P_base(1:floor(L/2)+1); P_base(2:end-1)=2*P_base(2:end-1);
        plot(f, P_tot, '-', 'Color', cTotal); plot(f, P_l1, '--', 'Color', cL1); plot(f, P_base, ':', 'Color', cBase);
        if ~isnan(d_freq_hz), xline(d_freq_hz, 'r--', 'LineWidth', 1, 'Label', 'Dist Freq'); end
    end
    xlim([0 20]);
end

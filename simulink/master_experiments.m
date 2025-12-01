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
P.use_mrac   = 0;         
% --- Disturbance ---
P.dist_amp  = [1.0; 0.5; 0.2]; 
P.dist_freq = [2.0; 1.5; 0.5]; 
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
fprintf('  7. Filter BW Sweep (Dist Freq Variation)\n');
fprintf('  8. Filter BW Sweep (2nd Order LPF)\n');
fprintf('  9. Filter BW Sweep (Dist Amplitude Variation)\n');
fprintf('--------------------------------------------------\n');
choice = input('Enter choice (1-9): ');
if isempty(choice), choice = 1; end
%% ================= 2. EXPERIMENT LOGIC =================
switch choice
    case 1 % --- EXP 1: DISTURBANCE REJECTION (3-Way) ---
        fprintf('Running Exp 1: Disturbance Rejection Comparison...\n');
        P.ref_type = 0; % Step
        % P.ref_amp  = 0; % Regulate to 0
        P.ref_start= 0;
        
        % 1. Run LQR Only
        P_run = P; P_run.use_mrac = 0; P_run.Clpf = zeros(3); assignin('base', 'P', P_run);
        res_lqr = run_sim(model_name, T_final, P_run);
        
        % 2. Run L1 + LQR
        P_run = P; P_run.use_mrac = 0; P_run.Clpf = eye(3); assignin('base', 'P', P_run);
        res_l1 = run_sim(model_name, T_final, P_run);
        
        % 3. Run MRAC
        P.MRAC.Gam_x = 50.0; P.MRAC.Gam_r = 50.0; P.MRAC.Gam_w = 50.0; 
        P_run = P; P_run.use_mrac = 1; assignin('base', 'P', P_run);
        res_mrac = run_sim(model_name, T_final, P_run);
        
        plot_3way_comparison(res_lqr, res_l1, res_mrac, '');
    case 2 % --- EXP 2: TRACKING PERFORMANCE (3-Way x 3 Refs) ---
        fprintf('Running Exp 2: Tracking Performance (Single Figure)...\n');
        
        % Define Test Cases
        cases(1).type=0; cases(1).amp=deg2rad(10); cases(1).start=0; cases(1).name='Step 10 deg';
        cases(2).type=0; cases(2).amp=deg2rad(40); cases(2).start=0; cases(2).name='Step 40 deg';
        cases(3).type=1; cases(3).amp=deg2rad(30); cases(3).start=0; cases(3).name='Sine Wave'; 
        
        num_cases = length(cases);
        
        % Create ONE Big Figure (2 Rows x 3 Cols)
        figure('Name', 'Experiment 2: Reference Tracking Comparison', 'Color', 'w', 'Position', [50, 50, 1400, 700]);
        tl = tiledlayout(2, num_cases, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, 'Tracking Performance Across Different References', 'FontSize', 16, 'Interpreter', 'latex');
        
        % Colors
        cLQR  = [0.4 0.4 0.4];    % Grey
        cL1   = [0.85 0.33 0.10]; % Orange
        cMRAC = [0 0.45 0.74];    % Blue
        
        for i = 1:num_cases
            c = cases(i);
            fprintf('  -> Testing Reference: %s\n', c.name);
            P.ref_type = c.type; 
            P.ref_amp  = c.amp;
            P.ref_start = c.start; % Set Start Time
            
            % Run Sims
            P_run = P; P_run.use_mrac=0; P_run.Clpf=zeros(3); assignin('base', 'P', P_run);
            r_lqr = run_sim(model_name, T_final, P_run);
            
            P_run = P; P_run.use_mrac=0; P_run.Clpf=eye(3); assignin('base', 'P', P_run);
            r_l1 = run_sim(model_name, T_final, P_run);
            
            P_run = P; P_run.use_mrac=1; assignin('base', 'P', P_run);
            r_mrac = run_sim(model_name, T_final, P_run);
            
            % --- PLOT COLUMN 'i' ---
            
            % 1. Tracking (Row 1, Col i)
            nexttile(i); hold on; grid on;
            if ~isempty(r_lqr), plot(r_lqr.t, rad2deg(r_lqr.x(:,1)), '--', 'Color', cLQR, 'LineWidth', 1.5); end
            if ~isempty(r_l1), plot(r_l1.t, rad2deg(r_l1.x(:,1)), '-.', 'Color', cL1, 'LineWidth', 2); end
            if ~isempty(r_mrac), plot(r_mrac.t, rad2deg(r_mrac.x(:,1)), '-', 'Color', cMRAC, 'LineWidth', 2); end
            % Ref
            if ~isempty(r_lqr), plot(r_lqr.t, rad2deg(r_lqr.ref(:,1)), 'k:', 'LineWidth', 1.5); end
            
            title(['\textbf{' c.name '}'], 'Interpreter', 'latex', 'FontSize', 12);
            if i==1, ylabel('$\phi$ (deg)', 'Interpreter', 'latex'); end
            if i==3, legend('LQR','L1+LQR','MRAC','Ref', 'Location','best','Interpreter','latex'); end
            
            % 2. Torque (Row 2, Col i) -> Index is i + num_cases
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
        % title('Effect of $T_s$', 'Interpreter', 'latex');
        % ylabel('RMSE (deg)', 'Interpreter', 'latex');
        % title('MRAC Performance: RMSE vs Adaptation Gain', 'Interpreter', 'latex');
        subtitle('Tracking 30 deg Step', 'Interpreter', 'latex');
        grid on;
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
        xlabel('Delay (ms)', 'Interpreter', 'latex'); 
        ylabel('RMSE (deg)', 'Interpreter', 'latex'); 
        title('Time Delay Robustness', 'Interpreter', 'latex');
        grid on;
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
        P_l1=P; P_l1.use_mrac=0; P_l1.Clpf=eye(3); assignin('base','P',P_l1);
        r_l1 = run_sim(model_name, T_final, P_l1);
        P_m=P; P_m.use_mrac=1; assignin('base','P',P_m);
        r_m = run_sim(model_name, T_final, P_m);
        plot_3way_comparison([], r_l1, r_m, 'MRAC Tuning Check');
    case 6 % --- MRAC Gain Sweep (RMSE vs Gamma) ---
        fprintf('Running MRAC Gain Sweep (RMSE vs Gamma)...\n');
        
        % Sweep range for Gamma
        gam_vals = [0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 2000, 10000];
        rmse_log = [];
        
        % Experiment Setup: Step Reference 30 deg
        P.ref_type = 0; % Step
        P.ref_amp  = deg2rad(30);
        P.ref_start = 0.0;
        P.use_mrac = 1; % Force MRAC Mode
        
        % Relax bounds for high gain testing
        P.MRAC.theta_x_max = 2000; P.MRAC.theta_r_max = 2000; P.MRAC.theta_w_max = 2000;
        
        for i = 1:length(gam_vals)
            g = gam_vals(i);
            % Set all Gammas equal
            P.MRAC.Gam_x = g; 
            P.MRAC.Gam_r = g; 
            P.MRAC.Gam_w = g;
            
            assignin('base', 'P', P);
            res = run_sim(model_name, T_final, P);
            
            if isempty(res)
                val = NaN; % Crashed
            else
                % Calculate RMSE
                err = res.ref(:,1) - res.x(:,1);
                val = rad2deg(sqrt(mean(err.^2)));
            end
            
            % Filter out massive instabilities
            if val > 1000, val = NaN; end
            
            rmse_log(i) = val;
            fprintf('Gamma=%.1f, RMSE=%.4f deg\n', g, val);
        end
        
        % Plot Results
        figure('Color','w','Position',[100 100 800 600]);
        semilogx(gam_vals, rmse_log, '-sq', 'LineWidth', 2, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
        xlabel('Adaptation Gain $\Gamma$ (Log Scale)', 'Interpreter', 'latex');
        ylabel('RMSE (deg)', 'Interpreter', 'latex');
        % title('MRAC Performance: RMSE vs Adaptation Gain', 'Interpreter', 'latex');
        subtitle('Tracking 30 deg Step', 'Interpreter', 'latex');
        grid on;
        
    case 7 % --- Filter Bandwidth vs Disturbance Frequency (1st Order) ---
        fprintf('Running 1st Order Filter BW Sweep for Varying Disturbance Frequencies...\n');
        
        % Parameters to Sweep
        dist_freqs = [0.1, 1, 2, 5, 10, 20, 30, 60]; 
        wc_vals = [0.1, 1, 5, 10, 20, 30, 40, 50, 75, 100]; 
        
        % Setup: L1 Control, Step Reference = 20 deg
        P.Ts = 0.002;              
        P.Ae = -10.0 * eye(6);
        P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0;
        P.use_mrac = 0; 
        
        rmse_results = zeros(length(dist_freqs), length(wc_vals));
        % Standard distinct colors (avoiding yellow/cyan)
        % Blue, Green, Orange, Red, Purple, Brown, Black, Magenta
        colors = {'b', [0 0.5 0], [0.85 0.33 0.10], 'r', [0.49 0.18 0.56], [0.6 0.2 0], 'k', 'm'};
        
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title('Disturbance Rejection (Amp=1 Nm): RMSE vs Bandwidth', 'Interpreter','latex');
        xlabel('Filter Bandwidth $\omega_c$ (rad/s)', 'Interpreter','latex');
        ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        
        for i = 1:length(dist_freqs)
            d_freq = dist_freqs(i);
            P.dist_freq = [d_freq; d_freq; d_freq];
            fprintf('\n--- Testing Disturbance Freq: %g rad/s ---\n', d_freq);
            
            for j = 1:length(wc_vals)
                w_c = wc_vals(j);
                % 1st Order Filter
                P.wc = w_c; 
                [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3));
                P.Clpf = eye(3); 
                
                assignin('base', 'P', P);
                res = run_sim(model_name, T_final, P);
                
                if isempty(res)
                    val = NaN; 
                else
                    val = rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2)));
                end
                rmse_results(i,j) = val;
                fprintf('  wc=%d, RMSE=%.4f deg\n', w_c, val);
            end
            plot(wc_vals, rmse_results(i,:), '-o', 'LineWidth', 2, 'Color', colors{i}, ...
                'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$d_{freq} = %g$ rad/s', d_freq));
        end
        legend('Location','best', 'Interpreter','latex');

    case 8 % --- Filter Bandwidth vs Disturbance Frequency (2nd Order) ---
        fprintf('Running 2nd Order Filter BW Sweep for Varying Disturbance Frequencies...\n');
        fprintf('NOTE: This uses 6 filter states (2 per channel). Ensure L1_LPF block can handle dimensions!\n');
        
        % Parameters to Sweep
        dist_freqs = [2, 5, 10, 20, 30, 60]; 
        wc_vals = [0.1, 1, 5, 10, 20, 30, 40, 50, 75, 100]; 
        
        % Setup: L1 Control, Step Reference = 20 deg
        P.Ts = 0.002;              
        P.Ae = -10.0 * eye(6);
        P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0;
        P.use_mrac = 0; 
        
        rmse_results = zeros(length(dist_freqs), length(wc_vals));
        colors = {'b', 'g', [0.85 0.33 0.10], 'r', 'm', 'k'};
        
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title('Disturbance Rejection (2nd Order LPF): RMSE vs Bandwidth', 'Interpreter','latex');
        xlabel('Filter Bandwidth $\omega_c$ (rad/s)', 'Interpreter','latex');
        ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        
        for i = 1:length(dist_freqs)
            d_freq = dist_freqs(i);
            P.dist_freq = [d_freq; d_freq; d_freq];
            fprintf('\n--- Testing Disturbance Freq: %d rad/s ---\n', d_freq);
            
            for j = 1:length(wc_vals)
                w_c = wc_vals(j);
                
                % 2nd Order Critically Damped Filter Construction
                % C(s) = wc^2 / (s^2 + 2*wc*s + wc^2)
                % State Space for 1 channel: xdot = [0 1; -wc^2 -2wc]x + [0; wc^2]u, y = [1 0]x
                A_1ch = [0 1; -w_c^2 -2*w_c];
                B_1ch = [0; w_c^2];
                C_1ch = [1 0];
                D_1ch = 0;
                
                % Block Diagonal for 3 channels (Roll, Pitch, Yaw)
                P.Alpf = blkdiag(A_1ch, A_1ch, A_1ch); % 6x6
                P.Blpf = blkdiag(B_1ch, B_1ch, B_1ch); % 6x3
                P.Clpf = blkdiag(C_1ch, C_1ch, C_1ch); % 3x6
                P.Dlpf = zeros(3,3);                   % 3x3
                
                % Store w_c for checking
                P.wc = w_c; 
                
                assignin('base', 'P', P);
                try
                    res = run_sim(model_name, T_final, P);
                    if isempty(res)
                        val = NaN; 
                    else 
                        val = rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2))); 
                    end
                catch ME
                    fprintf('  Error: %s. Check L1_LPF block output dimensions.\n', ME.message);
                    val = NaN;
                end
                
                rmse_results(i,j) = val;
                fprintf('  wc=%d, RMSE=%.4f deg\n', w_c, val);
            end
            plot(wc_vals, rmse_results(i,:), '-^', 'LineWidth', 2, 'Color', colors{i}, ...
                'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$d_{freq} = %d$ rad/s', d_freq));
        end
        legend('Location','best', 'Interpreter','latex');

    case 9 % --- Filter Bandwidth vs Disturbance Amplitude ---
        fprintf('Running Filter BW Sweep for Varying Disturbance Amplitudes (Fixed Freq)...\n');
        
        % Parameters
        dist_amps = [0.2, 0.5, 1.0, 2.0, 4.0]; % Disturbance Amplitudes (Nm)
        wc_vals   = [0.1, 1, 5, 10, 20, 30, 40, 50, 75, 100]; 
        
        % Setup: L1 Control (1st Order), Step Ref = 20 deg
        P.Ts = 0.002;              
        P.Ae = -10.0 * eye(6);
        P.ref_type = 0; P.ref_amp = deg2rad(20); P.ref_start = 0;
        P.use_mrac = 0; 
        
        % Fixed Frequency for this experiment
        FIXED_FREQ = 2.0; % rad/s
        P.dist_freq = [FIXED_FREQ; FIXED_FREQ; FIXED_FREQ];
        
        rmse_results = zeros(length(dist_amps), length(wc_vals));
        colors = {[0.85 0.33 0.10], 'r', [0.49 0.18 0.56], [0.6 0.2 0], 'k', 'm'}; % Improved colors
        
        figure('Color','w','Position',[100 100 800 600]); hold on; grid on;
        title(sprintf('Disturbance Rejection (Freq=%d rad/s): RMSE vs Bandwidth', FIXED_FREQ), 'Interpreter','latex');
        xlabel('Filter Bandwidth $\omega_c$ (rad/s)', 'Interpreter','latex');
        ylabel('RMSE (deg)', 'Interpreter','latex');
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12);
        
        for i = 1:length(dist_amps)
            d_amp = dist_amps(i);
            P.dist_amp = [d_amp; d_amp; d_amp];
            fprintf('\n--- Testing Disturbance Amp: %.1f Nm ---\n', d_amp);
            
            for j = 1:length(wc_vals)
                w_c = wc_vals(j);
                % 1st Order Filter Logic
                P.wc = w_c; 
                [P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3));
                P.Clpf = eye(3); 
                
                assignin('base', 'P', P);
                res = run_sim(model_name, T_final, P);
                
                if isempty(res)
                    val = NaN; 
                else
                    val = rad2deg(sqrt(mean((res.x(:,1) - P.ref_amp).^2)));
                end
                rmse_results(i,j) = val;
            end
            plot(wc_vals, rmse_results(i,:), '-s', 'LineWidth', 2, 'Color', colors{i}, ...
                'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('$d_{amp} = %.1f$ Nm', d_amp));
        end
        legend('Location','best', 'Interpreter','latex');
end
%% ================= 3. HELPER FUNCTIONS =================
function res = run_sim(model, T, P_struct)
    try
        out = sim(model, 'StopTime', num2str(T));
        res = extract_sim_data(out, P_struct);
    catch ME
        fprintf('Simulation Failed: %s\n', ME.message);
        res = [];
    end
end
function res = extract_sim_data(out, P)
    function d_fix = fix_d(d, t)
        d = squeeze(d);
        if size(d,1)==length(t), d_fix=d; elseif size(d,2)==length(t), d_fix=d.'; else, error('Dim'); end
    end
    try
        res.t = out.sim_x.Time;
        res.x = fix_d(out.sim_x.Data, res.t);
        res.ref = fix_d(out.sim_ref.Data, res.t);
        res.dist = fix_d(out.sim_dist.Data, res.t);
        
        try res.u_b = fix_d(out.sim_u_base.Data, res.t); catch, res.u_b = zeros(length(res.t),3); end
        try res.u_l1 = fix_d(out.sim_u_L1.Data, res.t); catch, res.u_l1 = zeros(length(res.t),3); end
        
        res.tau_b = [res.u_b(:,1)*P.L, res.u_b(:,2)*P.L, res.u_b(:,3)];
        res.tau_l1 = [res.u_l1(:,1)*P.L, res.u_l1(:,2)*P.L, res.u_l1(:,3)];
        
        if isfield(P, 'use_mrac') && P.use_mrac == 1
            try 
                u_m = fix_d(out.sim_u_mrac.Data, res.t);
                res.tau_tot = [u_m(:,1)*P.L, u_m(:,2)*P.L, u_m(:,3)];
                if max(abs(res.tau_tot(:,1))) < 1e-6, res.tau_tot = res.tau_b; end 
            catch
                res.tau_tot = res.tau_b; 
            end
        else
            res.tau_tot = res.tau_b + res.tau_l1;
        end
    catch
        res = [];
    end
end
function plot_3way_comparison(r1, r2, r3, title_str)
    % r1=LQR, r2=L1, r3=MRAC
    % Define figure size: wider but balanced
    fig_size = [100, 100, 1100, 450];
    figure('Color','w', 'Position', fig_size);
    % Create a 1x2 tiled layout
    tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, title_str, 'Interpreter','latex', 'FontSize',16);
    % Colors for plots
    cLQR  = [0.4 0.4 0.4];    % Grey
    cL1   = [0.85 0.33 0.10]; % Orange
    cMRAC = [0 0.45 0.74];    % Blue
    % cMRAC = [0 0.7 0];        % Green
    % --- Left: Tracking (angles) ---
    nexttile(1); hold on; grid on;
    if ~isempty(r1), plot(r1.t, rad2deg(r1.x(:,1)), '--', 'Color', cLQR, 'LineWidth', 1.5, 'DisplayName', 'LQR'); end
    if ~isempty(r2), plot(r2.t, rad2deg(r2.x(:,1)), '-.', 'Color', cL1, 'LineWidth', 2, 'DisplayName', 'L1+LQR'); end
    if ~isempty(r3), plot(r3.t, rad2deg(r3.x(:,1)), '-', 'Color', cMRAC, 'LineWidth', 2, 'DisplayName', 'MRAC'); end
    % Reference line
    if ~isempty(r1), plot(r1.t, rad2deg(r1.ref(:,1)), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Ref'); end
    ylabel('$\phi$ (deg)', 'Interpreter','latex', 'FontSize', 12);
    xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
    title('Tracking Angle', 'Interpreter','latex', 'FontSize', 14);
    legend('Location','best', 'Interpreter','latex');
    % --- Right: Control Effort ---
    nexttile(2); hold on; grid on;
    if ~isempty(r1), plot(r1.t, r1.tau_tot(:,1), '--', 'Color', cLQR, 'LineWidth', 1.5, 'DisplayName', 'LQR'); end
    if ~isempty(r2), plot(r2.t, r2.tau_tot(:,1), '-.', 'Color', cL1, 'LineWidth', 2, 'DisplayName', 'L1+LQR'); end
    if ~isempty(r3), plot(r3.t, r3.tau_tot(:,1), '-', 'Color', cMRAC, 'LineWidth', 2, 'DisplayName', 'MRAC'); end
    if ~isempty(r1), plot(r1.t, r1.dist(:,1), 'm:', 'LineWidth', 1, 'DisplayName', 'Dist'); end
    ylabel('$\tau$ (Nm)', 'Interpreter','latex', 'FontSize', 12);
    xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12);
    title('Control Effort', 'Interpreter','latex', 'FontSize', 14);
    legend('Location','best', 'Interpreter','latex');
end

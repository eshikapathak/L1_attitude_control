%% nonlinear_plant_l1_LQR_v25_clean_solver.m
% Quad Attitude Control: LQR + L1
% Update: Modularized RK4 Solver for cleaner main loop.
clear; clc; close all;

%% ==================== USER CONFIGURATION ====================
% Set your reference amplitude here (Degrees).
USER_REF_AMP = 20; 
%% ============================================================

%% -------------------- Simulation Setup --------------------
P.dt_sim = 5e-5; % Physics
P.g = 9.81;

fprintf('Physics dt: %.6f s (%.0f Hz)\n', P.dt_sim, 1/P.dt_sim);
fprintf('Reference Amplitude defined as: %.1f deg\n', USER_REF_AMP);
fprintf('--------------------------------------------\n');
fprintf('Select Experiment:\n');
fprintf('  1: Single Run (Nominal)\n');
fprintf('  2: Exp 1 - Varying References (Step vs Sine)\n');
fprintf('  3: Exp 2 - Effect of Sampling Time (Ts)\n');
fprintf('  4: Exp 3 - LQR vs L1 (Disturbance Rejection)\n');
fprintf('  5: Exp 4 - Adding Time Delays\n');
choice = input('Enter choice (1-5): ');
if isempty(choice) || ~ismember(choice, 1:5), choice = 1; end

%% -------------------- Parameters --------------------
P.m  = 0.76; P.L  = 0.14;
P.Ix = 0.0045; P.Iy = 0.0045; P.Iz = 0.0113;

% Linear Model
P.A_true = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;
            0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
P.Bm   = [ zeros(3,3); diag([1/P.Ix, 1/P.Iy, 1/P.Iz]) ];      
P.Bum  = [ eye(3); zeros(3,3) ];                         
P.Bbar = [P.Bm P.Bum];                          

% --- LQR Tuning ---
Q_agg = diag([100, 100, 80, 1, 1, 1]); R_agg = diag([0.1, 0.1, 0.1]); 
P.K_agg = lqr(P.A_true, P.Bm, Q_agg, R_agg);

Q_con = diag([10, 10, 10, 4, 4, 3]); R_con = diag([2, 2, 2]);
P.K_con = lqr(P.A_true, P.Bm, Q_con, R_con);

% --- L1 Parameters ---
P.Ts_default = 0.0002;             
P.Ae = -10.0*eye(6);
P.wc = 50;           
[P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = lpf_ct_vec(P.wc);

% Disturbance
P.dist = @(t,x) [0.5*sin(3*t); 0; 0]; 

Tfinal = 10.0; 
fig_size = [100, 100, 1200, 800]; 

%% -------------------- Experiment Config --------------------
Configs = [];
switch choice
    case 1 
        Configs(1).type='step'; Configs(1).amp=USER_REF_AMP; Configs(1).Ts=P.Ts_default; 
        Configs(1).L1=true; Configs(1).K=P.K_con; Configs(1).delay=0; Configs(1).Name='Nominal';
    case 2 
        Configs(1).type='step'; Configs(1).amp=USER_REF_AMP; Configs(1).Ts=P.Ts_default; Configs(1).L1=true; Configs(1).K=P.K_con; Configs(1).delay=0; Configs(1).Name=sprintf('Step %.0f', USER_REF_AMP);
        Configs(2).type='sine'; Configs(2).amp=USER_REF_AMP; Configs(2).Ts=P.Ts_default; Configs(2).L1=true; Configs(2).K=P.K_con; Configs(2).delay=0; Configs(2).Name=sprintf('Sine %.0f', USER_REF_AMP);
    case 3 
        % Ts Sweep
        Ts_list = [0.003, 0.001, 0.0005, 0.0002, 0.0001, 1e-5];
        for i = 1:length(Ts_list)
            Configs(i).type='step'; Configs(i).amp=USER_REF_AMP; Configs(i).Ts=Ts_list(i); 
            Configs(i).L1=true; Configs(i).K=P.K_con; Configs(i).delay=0; Configs(i).Name=sprintf('Ts=%.3fs', Ts_list(i));
        end
    case 4 
        % P.dist = @(t,x) [0.5*sin(3*t); 0; 0]; 
        Configs(1).type='step'; Configs(1).amp=USER_REF_AMP; Configs(1).Ts=P.Ts_default; Configs(1).L1=false; Configs(1).K=P.K_con; Configs(1).delay=0; Configs(1).Name='LQR (Con)';
        Configs(2).type='step'; Configs(2).amp=USER_REF_AMP; Configs(2).Ts=P.Ts_default; Configs(2).L1=false; Configs(2).K=P.K_agg; Configs(2).delay=0; Configs(2).Name='LQR (Agg)';
        Configs(3).type='step'; Configs(3).amp=USER_REF_AMP; Configs(3).Ts=P.Ts_default; Configs(3).L1=true;  Configs(3).K=P.K_con; Configs(3).delay=0; Configs(3).Name='L1 (Con Base)';
    case 5 
        delays = [0, 0.002, 0.004, 0.006, 0.008];
        for i = 1:length(delays)
            Configs(i).type='step'; Configs(i).amp=USER_REF_AMP; Configs(i).Ts=P.Ts_default; 
            Configs(i).L1=true; Configs(i).K=P.K_con; Configs(i).delay=delays(i);
            Configs(i).Name=sprintf('Delay %.0fms', delays(i)*1000);
        end
end

%% -------------------- Simulation Loop --------------------
Results = {}; 
for i = 1:length(Configs)
    cfg = Configs(i);
    fprintf('Running %s (Ts=%.5fs)...', cfg.Name, cfg.Ts);
    
    SimSteps = round(Tfinal / P.dt_sim);
    ratio = round(cfg.Ts / P.dt_sim); if ratio < 1, ratio = 1; end 
    
    P.useL1 = cfg.L1; P.K_current = cfg.K;
    
    % Init
    x_plant = zeros(6,1); x_hat = zeros(6,1); chi = zeros(3,1); sigma = zeros(6,1);
    
    % Delay Buffer
    delay_ctrl_steps = round(cfg.delay / cfg.Ts);
    max_ctrl_steps = ceil(Tfinal/cfg.Ts) + 100;
    tau_buffer = zeros(max_ctrl_steps, 3);
    
    % Matrices
    Ad_pred = expm(P.A_true * cfg.Ts);
    Bd_pred = P.Bm * cfg.Ts + (P.A_true * P.Bm * cfg.Ts^2)/2; 
    Bud_pred = P.Bum * cfg.Ts;
    E_mx = expm(P.Ae * cfg.Ts);
    Phi_inv = inv(P.Ae \ (E_mx - eye(size(P.Ae))));
    
    H_t = []; H_x = []; H_ref = []; 
    H_tau = []; H_base = []; H_l1 = []; H_sigma = [];
    
    ctrl_step_count = 0;
    current_tau_cmd = [0;0;0]; current_tau_plant = [0;0;0];
    crashed = false;
    
    if strcmp(cfg.type, 'step'), ref_fun = @(t) deg2rad([cfg.amp; 10; 0]);
    else, ref_fun = @(t) deg2rad([cfg.amp*sin(2*pi*1*t); 10; 0]); end
    
    % Loop
    for k = 1:SimSteps
        t_now = (k-1) * P.dt_sim;
        
        % === CONTROLLER ===
        if mod(k-1, ratio) == 0
            ctrl_step_count = ctrl_step_count + 1;
            ref = ref_fun(t_now);
            
            tau_base = -P.K_current * ([x_plant(1:3); x_plant(4:6)] - [ref; 0;0;0]);
            
            if P.useL1, u_l1 = -(P.Clpf*chi + P.Dlpf*sigma(1:3)); 
            else, u_l1 = zeros(3,1); end
            
            current_tau_cmd = tau_base + u_l1; 
            
            tau_buffer(ctrl_step_count, :) = current_tau_cmd';
            read_idx = ctrl_step_count - delay_ctrl_steps;
            if read_idx < 1, current_tau_plant = zeros(3,1);
            else, current_tau_plant = tau_buffer(read_idx, :)'; end
            
            if P.useL1
                pred_input = current_tau_cmd + sigma(1:3);
                x_hat_next = Ad_pred * x_hat + Bd_pred * pred_input + Bud_pred * sigma(4:6) ...
                           + (E_mx - eye(6)) * (x_hat - x_plant); 
                chi_next = chi + cfg.Ts * (P.Alpf*chi + P.Blpf*sigma(1:3));
                mu = Phi_inv * (E_mx * (x_hat - x_plant));
                sigma = - (P.Bbar \ mu);
                x_hat = x_hat_next; chi = chi_next;
            end
            
            H_t(end+1,1) = t_now; 
            H_x(end+1,:) = x_plant';
            H_ref(end+1,:) = ref';
            H_tau(end+1,:) = current_tau_plant';
            H_base(end+1,:) = tau_base';
            H_l1(end+1,:) = u_l1';
            H_sigma(end+1,:) = sigma';
        end
        
        % === PHYSICS INTEGRATION (Using Helper) ===
        % Define the derivative function for this specific step
        derivs = @(t, x) plant_dynamics(x, current_tau_plant, t, P);
        
        % Perform Fixed Step Integration
        x_plant = RK4_step(derivs, t_now, x_plant, P.dt_sim);
        
        % === CRASH CHECK ===
        current_ref_check = ref_fun(t_now);
        current_error = x_plant(1:3) - current_ref_check;
        
        if any(isnan(x_plant)) || any(abs(current_error) > 5.0)
            crashed = true; 
            fprintf(' CRASHED (Error Explosion) at %.3fs\n', t_now);
            break; 
        end
    end
    if ~crashed, fprintf(' Done.\n'); end
    
    Res.t=H_t; Res.x=H_x; Res.tau=H_tau; Res.ref=H_ref; 
    Res.base=H_base; Res.l1=H_l1; Res.sigma=H_sigma;
    Res.name=cfg.Name; Res.hasL1 = cfg.L1; Res.Ts = cfg.Ts; Res.K = cfg.K;
    Results{i} = Res;
end

plot_comprehensive_results(Results, choice, fig_size, P);

%% ==================== SOLVERS & DYNAMICS ====================

% --- FIXED STEP RUNGE-KUTTA 4 SOLVER ---
function x_next = RK4_step(f, t, x, h)
    k1 = f(t, x);
    k2 = f(t + 0.5*h, x + 0.5*h*k1);
    k3 = f(t + 0.5*h, x + 0.5*h*k2);
    k4 = f(t + h, x + h*k3);
    x_next = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dx = plant_dynamics(x, tau, t, P)
    phi=x(1); theta=x(2); psi=x(3); p=x(4); q=x(5); r=x(6);
    d = P.dist(t,x);
    cth = cos(theta); if abs(cth)<1e-2, cth=sign(cth)*1e-2; end
    tt = tan(theta); cph = cos(phi); sph = sin(phi);
    phid = p + q*sph*tt + r*cph*tt;
    thtd = q*cph - r*sph;
    psid = q*sph/cth + r*cph/cth;
    pd = (P.Iy - P.Iz)*q*r/P.Ix + (tau(1) + d(1))/P.Ix;
    qd = (P.Iz - P.Ix)*p*r/P.Iy + (tau(2) + d(2))/P.Iy;
    rd = (P.Ix - P.Iy)*p*q/P.Iz + (tau(3) + d(3))/P.Iz;
    dx = [phid; thtd; psid; pd; qd; rd];
end

function [A,B,C,D] = lpf_ct_vec(wc)
    A = -wc*eye(3); B = wc*eye(3); C = eye(3); D = zeros(3);
end

function plot_comprehensive_results(Results, choice, fig_size, P)
    set(groot,'defaultTextInterpreter','latex');
    cols = lines(length(Results));
    
    % === FIGURE 1: PERFORMANCE SUMMARY ===
    figure('Color','w','Position',fig_size, 'Name', 'Performance Summary');
    subplot(2,2,[1 2]); hold on; grid on; title('\textbf{Attitude Tracking}'); ylabel('Angle (deg)');
    for i=1:length(Results)
        R = Results{i};
        if length(R.t) > 1
            plot(R.t, rad2deg(R.x(:,1)), 'Color', cols(i,:), 'LineWidth', 1.5, 'DisplayName', [R.name ' (\phi)']);
            if i==1, plot(R.t, rad2deg(R.ref(:,1)), 'k--', 'LineWidth', 1.0, 'DisplayName', 'Ref'); end
        end
    end
    legend('Location','best');
    subplot(2,2,3); hold on; grid on; title('Control (Full)'); ylabel('Nm');
    for i=1:length(Results)
        if length(Results{i}.t) > 1, plot(Results{i}.t, Results{i}.tau(:,1), 'Color', cols(i,:)); end
    end
    subplot(2,2,4); hold on; grid on; title('Control (Zoomed)'); ylabel('Nm');
    y_min=100; y_max=-100;
    for i=1:length(Results)
        R = Results{i};
        if length(R.t) > 5
            plot(R.t, R.tau(:,1), 'Color', cols(i,:));
            idx = max(2, round(length(R.t)*0.1)); 
            if idx < length(R.t), vals = R.tau(idx:end,1); y_min = min(y_min, min(vals)); y_max = max(y_max, max(vals)); end
        end
    end
    if y_max>y_min, ylim([y_min*1.2, y_max*1.2]); else, ylim([-5 5]); end
    
    % === FIGURE 2: SPECIALIZED ANALYSIS ===
    if choice == 3 
        % EXP 3: RMSE vs Ts & Margin vs Ts
        figure('Color','w','Position',[100, 100, 1000, 500], 'Name', 'Sampling Analysis');
        
        % -- Subplot 1: RMSE vs Ts (Log Scale) --
        subplot(1,2,1); hold on; grid on;
        Ts_vals = []; rmse_vals = [];
        for i=1:length(Results)
            R = Results{i};
            Ts_vals(end+1) = R.Ts;
            if length(R.t) > 50
                err = rad2deg(R.ref(:,1) - R.x(:,1));
                idx = find(R.t > 1.0, 1); if isempty(idx), idx=1; end
                rmse_vals(end+1) = sqrt(mean(err(idx:end).^2));
            else
                rmse_vals(end+1) = NaN; 
            end
        end
        semilogx(Ts_vals, rmse_vals, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b');
        title('\textbf{Tracking Performance vs Sampling Time}'); 
        xlabel('Sampling Time $T_s$ (s) [Log Scale]'); ylabel('RMSE (deg)');
        set(gca, 'XDir', 'reverse'); % Faster to the right
        
        % -- Subplot 2: Theoretical Delay Margin vs Ts --
        subplot(1,2,2); hold on; grid on;
        margin_vals = [];
        for i=1:length(Results)
            ts_curr = Results{i}.Ts;
            % Analytical Calculation for Linear LQR Baseline (Roll Axis)
            A_sub = [0 1; 0 0]; B_sub = [0; 1/P.Ix];
            sys_c = ss(A_sub, B_sub, [1 0], 0);
            sys_d = c2d(sys_c, ts_curr, 'zoh');
            K_sub = Results{i}.K(1, [1, 4]); % Roll gains
            
            % Loop Margin
            L = series(sys_d, tf(K_sub, 1, ts_curr));
            [Gm, Pm, Wcg, Wcp] = margin(L);
            
            if isinf(Gm) && isinf(Pm), time_margin = NaN;
            else, time_margin = (deg2rad(Pm) / Wcp); end
            margin_vals(end+1) = time_margin * 1000; % ms
        end
        
        semilogx(Ts_vals, margin_vals, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'r');
        title('\textbf{Theoretical LQR Delay Margin vs } $T_s$');
        xlabel('Sampling Time $T_s$ (s) [Log Scale]'); ylabel('Delay Margin (ms)');
        yline(0, 'k--');
        
    else
        % Standard Diagnostics
        figure('Color','w','Position',fig_size, 'Name', 'Diagnostics');
        subplot(2,1,1); hold on; grid on; title('\textbf{Tracking Error}'); ylabel('deg');
        for i=1:length(Results), R=Results{i}; if length(R.t)>1, plot(R.t, rad2deg(R.ref(:,1)-R.x(:,1)), 'Color', cols(i,:)); end; end
        subplot(2,1,2); hold on; grid on; title('\textbf{Phase Plane}'); ylabel('rate'); xlabel('angle');
        for i=1:length(Results), R=Results{i}; if length(R.t)>1, plot(R.x(:,1), R.x(:,4), 'Color', cols(i,:)); end; end
    end
    
    % === FIGURE 3: L1 DECOMPOSITION ===
    l1_indices = find(cellfun(@(x) x.hasL1, Results));
    if ~isempty(l1_indices)
        to_plot = [];
        if choice == 3 
             for k=1:length(l1_indices), if length(Results{l1_indices(k)}.t)>50, to_plot=[to_plot, l1_indices(k)]; end; end
             if length(to_plot)>2, to_plot=[to_plot(1), to_plot(end)]; end
             if isempty(to_plot), to_plot=l1_indices(end); end
        else, to_plot = l1_indices(1); end
        
        for idx = to_plot
            R = Results{idx}; if length(R.t)<5, continue; end
            figure('Color','w','Position',fig_size, 'Name', ['Decomposition: ' R.name]);
            subplot(2,1,1); hold on; grid on; title(['\textbf{Control Decomposition}: ' R.name]); ylabel('Nm');
            plot(R.t, R.base(:,1), 'b-', 'DisplayName', '$\tau_{LQR}$');
            plot(R.t, R.l1(:,1), 'r-', 'LineWidth', 1.5, 'DisplayName', '$u_{L1}$');
            plot(R.t, R.tau(:,1), 'k--', 'DisplayName', '$\tau_{Total}$');
            legend('Location','best');
            zoom_idx = max(2, round(length(R.t)*0.05)); vals=[R.base(zoom_idx:end,1); R.l1(zoom_idx:end,1)];
            if ~isempty(vals), mn=min(vals); mx=max(vals); if mx>mn, ylim([mn*1.2, mx*1.2]); end; end
            subplot(2,1,2); hold on; grid on; title('Adaptive Parameters'); ylabel('Mag');
            plot(R.t, R.sigma(:,1), 'r-', 'DisplayName', 'Matched');
            plot(R.t, R.sigma(:,4), 'b--', 'DisplayName', 'Unmatched');
            legend;
        end
    end
end
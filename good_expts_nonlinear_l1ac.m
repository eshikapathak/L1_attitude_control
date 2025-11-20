%% nonlinear_plant_l1_LQR_final_v8.m
% Quad attitude: LQR Baseline + L1 (PC closed-form).
% Update: Added Torque Difference Plot for Experiment 3.

clear; clc; close all;

%% -------------------- Experiment Selection --------------------
fprintf('Select Simulation Mode:\n');
fprintf('  1: Single Run (Detailed Breakdown)\n');
fprintf('  2: Exp 1 - Varying References (Step + Sine)\n');
fprintf('  3: Exp 2 - Effect of Sampling Time (Ts) on RMSE\n');
fprintf('  4: Exp 3 - Disturbance: LQR (Conservative) vs LQR (Aggressive) vs L1\n');
choice = input('Enter choice (1-4): ');

if isempty(choice) || ~ismember(choice, [1, 2, 3, 4])
    choice = 1; fprintf('Invalid choice. Defaulting to Mode 1.\n');
end

%% -------------------- Parameters --------------------
P.useL1 = true; 

% Geometry/inertias
P.m  = 0.76; P.L  = 0.14;
P.Ix = 0.0045; P.Iy = 0.0045; P.Iz = 0.0113;

% Linear Model (Predictor)
P.A_true = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;
            0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
P.Bm = [ zeros(3,3); diag([1/P.Ix, 1/P.Iy, 1/P.Iz]) ];      
P.Bum = [ eye(3); zeros(3,3) ];                         
P.Bbar = [P.Bm P.Bum];                          

P.Bd = [zeros(3,3); 1/P.Ix 0 0; 0 1/P.Iy 0; 0 0 1/P.Iz];

% --- LQR Designs (Multiple Tunings) ---
A_sys = P.A_true; B_sys = P.Bm;

% Tuning 1: Aggressive (High Q, Low R)
Q1 = diag([100, 100, 80, 1, 1, 1]); 
R1 = diag([0.01, 0.01, 0.01]);      
P.K_agg = lqr(A_sys, B_sys, Q1, R1);

% Tuning 2: Conservative (Low Q, High R)
% Q2 = diag([10, 10, 10, 1, 1, 1]); 
% R2 = diag([1, 1, 1]);      
Q2 = diag([10, 10, 10,  4, 4, 3]);   % [phi theta psi dphi dtheta dpsi]
R2 = diag([2, 2, 2]);
P.K_con = lqr(A_sys, B_sys, Q2, R2);

P.K_current = P.K_agg; 

% --- L1 Design ---
P.Ts_default = 0.002;              
P.Ae = -1.0*eye(6);                
P.A  = P.A_true;                   
P.wc = 25;                         
[P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = lpf_ct_vec(P.wc);

% --- Disturbance ---
P.dist = @(t,x) [1 * sin(3*t); 2 * sin(5*t); 0];

% Simulation settings
Tfinal = 8; 
fig_size = [100, 100, 900, 600]; 

%% -------------------- Experiment Configuration --------------------

switch choice
    case 1  % --- DEFAULT SINGLE RUN ---
        Configs(1).type = 'step'; Configs(1).amp = 10; Configs(1).freq = 0;
        Configs(1).Ts = P.Ts_default; Configs(1).L1 = true; Configs(1).K = P.K_agg;
        
    case 2  % --- EXP 1: VARYING REFERENCES ---
        Configs(1).type = 'step'; Configs(1).amp = 10; Configs(1).freq = 0; Configs(1).Ts = P.Ts_default; Configs(1).L1 = true; Configs(1).K = P.K_agg;
        Configs(2).type = 'step'; Configs(2).amp = 20; Configs(2).freq = 0; Configs(2).Ts = P.Ts_default; Configs(2).L1 = true; Configs(2).K = P.K_agg;
        Configs(3).type = 'step'; Configs(3).amp = 40; Configs(3).freq = 0; Configs(3).Ts = P.Ts_default; Configs(3).L1 = true; Configs(3).K = P.K_agg;
        Configs(4).type = 'sine'; Configs(4).amp = 15; Configs(4).freq = 1.0; Configs(4).Ts = P.Ts_default; Configs(4).L1 = true; Configs(4).K = P.K_agg;
        
    case 3  % --- EXP 2: VARYING Ts ---
        Ts_sweep = [0.1, 0.08, 0.05, 0.02, 0.01, 0.008, 0.005, 0.002, 0.001];
        for i = 1:length(Ts_sweep)
            Configs(i).type = 'step'; Configs(i).amp = 10; Configs(i).freq = 0; 
            Configs(i).Ts = Ts_sweep(i); Configs(i).L1 = true; Configs(i).K = P.K_agg;
        end
        
    case 4  % --- EXP 3: LQR Tunings vs L1 ---
        % Run 1: LQR Conservative (L1 OFF)
        Configs(1).type = 'step'; Configs(1).amp = 10; Configs(1).freq = 0; 
        Configs(1).Ts = P.Ts_default; Configs(1).L1 = false; Configs(1).K = P.K_con;
        Configs(1).Name = 'LQR (Conservative)';
        
        % Run 2: LQR Aggressive (L1 OFF)
        Configs(2).type = 'step'; Configs(2).amp = 10; Configs(2).freq = 0; 
        Configs(2).Ts = P.Ts_default; Configs(2).L1 = false; Configs(2).K = P.K_agg;
        Configs(2).Name = 'LQR (Aggressive)';
        
        % Run 3: L1 + Aggressive LQR
        Configs(3).type = 'step'; Configs(3).amp = 10; Configs(3).freq = 0; 
        Configs(3).Ts = P.Ts_default; Configs(3).L1 = true; Configs(3).K = P.K_agg;
        Configs(3).Name = 'L1 + LQR (Aggressive)';
end

%% -------------------- Plotting Prep --------------------
colors = {'r', 'g', 'b', 'k', 'm', 'c', [0.5 0 0.5], [0.9 0.5 0]};
Results = {}; 
rmse_results = [];

%% -------------------- Main Simulation Loop --------------------

for i = 1:length(Configs)
    
    % 1. Config
    cfg = Configs(i);
    current_Ts = cfg.Ts;
    P.useL1 = cfg.L1;
    P.K_current = cfg.K; 
    
    % 2. Define Reference Function
    if strcmp(cfg.type, 'step')
        ref_fun = @(t) deg2rad([cfg.amp; 15; 0]);
        ref_label = sprintf('Step %d^o', cfg.amp);
    else
        ref_fun = @(t) deg2rad([cfg.amp*sin(cfg.freq*t); 15; 0]);
        ref_label = sprintf('Sine %d^o', cfg.amp);
    end
    
    % 3. Time and Solver
    N = round(Tfinal/current_Ts);
    t = (0:N)'*current_Ts;
    opt = odeset('RelTol',1e-5,'AbsTol',1e-7, 'InitialStep', current_Ts, 'MaxStep', current_Ts);
    
    % 4. Initialize
    x = zeros(6,1); xhat = zeros(6,1); chi = zeros(3,1); sigma = zeros(6,1);
    x_hist = zeros(N+1,6); x_hist(1,:) = x.';
    tau_tot_hist = zeros(N+1,3);
    tau_b_hist = zeros(N+1,3);
    uL1_hist = zeros(N+1,3);
    
    [tau_b, uL1_tau, tau_tot, ~] = signals_at(0, x, chi, sigma, P, ref_fun);
    tau_b_hist(1,:) = tau_b.'; uL1_hist(1,:) = uL1_tau.'; tau_tot_hist(1,:) = tau_tot.';

    fprintf('Running Run %d/%d: Type=%s, L1=%d, Ts=%.4f...\n', ...
            i, length(Configs), cfg.type, cfg.L1, current_Ts);

    % 5. Integration Loop
    for k = 1:N
        t0 = t(k); t1 = t(k+1);
        
        [tau_b, uL1_tau, tau_tot, ~] = signals_at(t0, x, chi, sigma, P, ref_fun);
        tau_b_hist(k,:) = tau_b.';
        uL1_hist(k,:)   = uL1_tau.';
        tau_tot_hist(k,:) = tau_tot.';
        
        X0  = [x; xhat; chi];
        rhs = @(tt,XX) rhs_interval_slide(tt, XX, sigma, P, ref_fun);
        [~, Xseg] = ode45(rhs, [t0 t1], X0, opt); 
        X1  = Xseg(end,:).';
        
        x = X1(1:6); xhat = X1(7:12); chi = X1(13:15);
        
        if P.useL1
            sigma = pc_update_closedform(xhat, x, P.Ae, P.Bm, P.Bum, current_Ts);
        end
        x_hist(k+1,:) = x.';
    end
    [tau_b, uL1_tau, tau_tot, ~] = signals_at(t(end), x, chi, sigma, P, ref_fun);
    tau_b_hist(end,:) = tau_b.'; uL1_hist(end,:) = uL1_tau.'; tau_tot_hist(end,:) = tau_tot.';

    % 6. Robust Reference Generation
    REF_vals = zeros(length(t), 3);
    for k_ref = 1:length(t)
        r_val = ref_fun(t(k_ref));
        REF_vals(k_ref, :) = r_val(:).';
    end
    
    % Store results
    Res.t = t; Res.x = x_hist; Res.tau = tau_tot_hist; Res.tau_b = tau_b_hist; Res.uL1 = uL1_hist;
    Res.ref = REF_vals; Res.label = ref_label;
    if isfield(cfg, 'Name'), Res.name = cfg.Name; end
    Results{i} = Res;
    
    % RMSE
    if choice == 3
        s_idx = round(1.0/current_Ts); if s_idx<1, s_idx=1; end
        err = x_hist(s_idx:end, 1) - REF_vals(s_idx:end, 1);
        rmse = sqrt(mean(err.^2));
        rmse_results(i) = rmse;
    end
end

%% -------------------- Post-Simulation Plotting --------------------

if choice == 1
    plot_single_run_detailed(Results{1}, fig_size);
    
elseif choice == 2
    % EXP 1: Varying Reference
    f1 = figure('Color','w','Name','Exp 1: Angles', 'Position', fig_size); 
    ax1 = axes(f1); hold(ax1,'on'); grid(ax1,'on'); 
    title(ax1,'Roll ($\phi$) Tracking: Varying References'); xlabel(ax1,'Time (s)'); ylabel(ax1,'deg');
    
    f2 = figure('Color','w','Name','Exp 1: Torques', 'Position', fig_size);
    ax2 = axes(f2); hold(ax2,'on'); grid(ax2,'on');
    title(ax2,'Total Roll Torque ($\tau_\phi$)'); xlabel(ax2,'Time (s)'); ylabel(ax2,'Nm');
    
    for i = 1:length(Results)
        c = colors{mod(i-1, length(colors))+1};
        Res = Results{i};
        plot(ax1, Res.t, rad2deg(Res.x(:,1)), 'LineWidth', 2, 'Color', c, 'DisplayName', [Res.label ' (Actual)']);
        plot(ax1, Res.t, rad2deg(Res.ref(:,1)), 'LineStyle', '--', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(ax2, Res.t, Res.tau(:,1), 'LineWidth', 1.5, 'Color', c, 'DisplayName', Res.label);
    end
    pad_axes(ax1); legend(ax1, 'Location', 'best'); set(ax1, 'FontSize', 12);
    pad_axes(ax2); legend(ax2, 'Location', 'best'); set(ax2, 'FontSize', 12);
    
elseif choice == 3
    % EXP 2: RMSE
    figure('Color','w','Name','Effect of Ts on Performance', 'Position', fig_size);
    semilogx([Configs.Ts], rmse_results, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b');
    grid on; xlabel('Sampling Time $T_s$ (s) [Log Scale]'); ylabel('RMSE (rad)');
    title('Performance vs. Sampling Rate ($T_s$)');
    subtitle('Lower $T_s$ (Faster Adaptation) $\to$ Lower Error');
    set(gca, 'XDir','reverse', 'FontSize', 12); pad_axes(gca); set(groot,'defaultTextInterpreter','latex');
    
elseif choice == 4
    % EXP 3: Comparison
    plot_comparison_multi_LQR(Results, fig_size);
    
    % --- NEW PLOT: Torque Difference ---
    % We compare Run 3 (L1+Agg) vs Run 2 (LQR)
    plot_torque_difference(Results{3}, Results{1}, fig_size);
end

%% ==================== Helper functions ====================

function pad_axes(ax)
    axis(ax, 'tight'); yl = ylim(ax); range = yl(2) - yl(1);
    if range < 1e-6, range = 1; end
    ylim(ax, [yl(1)-0.15*range, yl(2)+0.15*range]);
end

function plot_torque_difference(ResL1, ResLQR, fig_size)
    set(groot,'defaultTextInterpreter','latex');
    t = ResL1.t;
    
    figure('Color','w','Name','Exp 3: Torque Difference', 'Position', fig_size);
    tiledlayout(3,1);
    lbl = {'$\tau_\phi$ (Nm)', '$\tau_\theta$ (Nm)', '$\tau_\psi$ (Nm)'};
    
    for i=1:3
        ax = nexttile; hold on; grid on;
        
        % Calculate Difference
        diff_tau = ResL1.tau(:,i) - ResLQR.tau(:,i);
        
        plot(t, diff_tau, 'k-', 'LineWidth', 1.5);
        ylabel(['\Delta ' lbl{i}]); 
        if i==1, title('Control Torque Difference (L1_{Total} - LQR_{Total})'); end
        if i==3, xlabel('Time (s)'); end
        
        yline(0, 'r:', 'LineWidth', 1);
        pad_axes(ax); set(gca,'FontSize',12);
    end
end

function plot_comparison_multi_LQR(Results, fig_size)
    set(groot,'defaultTextInterpreter','latex');
    
    % ANGLES
    figure('Color','w','Name','Exp 3: Angles Comparison', 'Position', fig_size); 
    tiledlayout(3,1);
    ylabs = {'$\phi$ (Roll)','$\theta$ (Pitch)','$\psi$ (Yaw)'};
    
    styles = {':', '--', '-'}; widths = {2, 2, 1.5}; cols = {'b', 'b', 'r'};
    
    for i=1:3
        ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
        plot(Results{3}.t, rad2deg(Results{3}.ref(:,i)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
        for r=1:3
            res = Results{r};
            plot(res.t, rad2deg(res.x(:,i)), 'LineStyle', styles{r}, 'LineWidth', widths{r}, ...
                 'Color', cols{r}, 'DisplayName', res.name);
        end
        ylabel([ylabs{i} ' (deg)']); 
        if i==1, title('Disturbance Rejection: LQR Tunings vs L1'); legend('Location','best'); end
        pad_axes(ax); set(gca,'FontSize',12);
    end
    
    % TORQUES
    figure('Color','w','Name','Exp 3: Torques Comparison', 'Position', fig_size); 
    tiledlayout(3,1);
    lbl = {'$\tau_\phi$ (Nm)', '$\tau_\theta$ (Nm)', '$\tau_\psi$ (Nm)'};
    
    for i=1:3
        ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
        for r=1:3
            res = Results{r};
            plot(res.t, res.tau(:,i), 'LineStyle', styles{r}, 'LineWidth', widths{r}, ...
                 'Color', cols{r}, 'DisplayName', res.name);
        end
        ylabel(lbl{i}); 
        if i==1, title('Control Effort Comparison'); legend('Location','best'); end
        pad_axes(ax); set(gca,'FontSize',12);
    end
end

function plot_single_run_detailed(Res, fig_size)
    t = Res.t; x_hist = Res.x; tau_b = Res.tau_b; uL1 = Res.uL1; tau_tot = Res.tau; REF = Res.ref;
    figure('Color','w','Name','Single Run Results', 'Position', fig_size); tiledlayout(3,1);
    set(groot,'defaultTextInterpreter','latex');
    
    ax1 = nexttile; plot(t,rad2deg(x_hist(:,1)),'LineWidth',2); hold on; grid on;
             plot(t,rad2deg(REF(:,1)),'k--','LineWidth',1.5); ylabel('$\phi$ (deg)'); title('Attitude Tracking'); legend('$\phi$','Ref'); pad_axes(ax1);
    ax2 = nexttile; plot(t,rad2deg(x_hist(:,2)),'LineWidth',2); hold on; grid on;
             plot(t,rad2deg(REF(:,2)),'k--','LineWidth',1.5); ylabel('$\theta$ (deg)'); legend('$\theta$','Ref'); pad_axes(ax2);
    ax3 = nexttile; plot(t,rad2deg(x_hist(:,3)),'LineWidth',2); hold on; grid on;
             plot(t,rad2deg(REF(:,3)),'k--','LineWidth',1.5); ylabel('$\psi$ (deg)'); xlabel('Time (s)'); legend('$\psi$','Ref'); pad_axes(ax3);
    figure('Color','w','Name','Control Torques', 'Position', fig_size); tiledlayout(3,1);
    lbl = {'$\tau_\phi$ (Nm)', '$\tau_\theta$ (Nm)', '$\tau_\psi$ (Nm)'};
    for i=1:3
        ax = nexttile; hold on; grid on;
        plot(t, tau_b(:,i), 'b--', 'LineWidth', 1.2); plot(t, uL1(:,i), 'r--', 'LineWidth', 1.2); plot(t, tau_tot(:,i), 'k-', 'LineWidth', 1.5);
        ylabel(lbl{i}); pad_axes(ax);
        if i==1, title('Base vs L1 Breakdown'); legend('$\tau_{Base}$', '$u_{L1}$', '$\tau_{Tot}$', 'Location','best'); end
    end
end

function [A,B,C,D] = lpf_ct_vec(wc)
    A = -wc*eye(3); B = wc*eye(3); C = eye(3); D = zeros(3);
end

function [tau_base, uL1_tau, tau_total, U_mix] = signals_at(tt, x, chi, sigma, P, ref_fun)
    ref     = ref_fun(tt);
    tau_base= baseline_controller_torque_LQR(x, ref, P); 
    
    if P.useL1
        sig_m   = sigma(1:3);                              
        uL1_tau = -(P.Clpf*chi + P.Dlpf*sig_m);            
    else
        uL1_tau = zeros(3,1);
    end
    tau_total = tau_base + uL1_tau;
    U_mix     = torque_to_mix(tau_total, P);
end

function dX = rhs_interval_slide(t, X, sigma_hold, P, ref_fun)
    x = X(1:6); xhat = X(7:12); chi  = X(13:15);
    ref     = ref_fun(t);
    tau_b   = baseline_controller_torque_LQR(x, ref, P);
    
    if P.useL1
        sig_m    = sigma_hold(1:3);
        chi_dot  = P.Alpf*chi + P.Blpf*sig_m;
        uL1_tau  = -(P.Clpf*chi + P.Dlpf*sig_m);
        xt       = xhat - x;
        sig_m    = sigma_hold(1:3); sig_um = sigma_hold(4:6);      
        xhat_dot = P.A*xhat + P.Bm*(tau_b + uL1_tau + sig_m) + P.Bum*sig_um + P.Ae*xt;
    else
        chi_dot = zeros(3,1); uL1_tau = zeros(3,1); xhat_dot = zeros(6,1);    
    end
    
    tau_tot  = tau_b + uL1_tau;
    U_mix    = torque_to_mix(tau_tot, P);
    xdot_nl  = nonlinear_plant_dynamics(x, U_mix, P);
    xdot     = xdot_nl + P.Bd*P.dist(t,x);           
    dX = [xdot; xhat_dot; chi_dot];
end

function sigma = pc_update_closedform(xhat, x, Ae, Bm, Bum, Ts)
    xt = xhat - x; E = expm(Ae*Ts); Phi = Ae \ (E - eye(size(Ae)));         
    z = Phi \ (E*xt); Bbar = [Bm Bum]; sigma= - (Bbar \ z);                     
end

function tau = baseline_controller_torque_LQR(x, ref, P)
    x_des = [ref(:); 0; 0; 0];
    e_full = x - x_des;
    tau = -P.K_current * e_full; 
end

function U = torque_to_mix(tau, P)
    U = [tau(1)/P.L; tau(2)/P.L; tau(3)];
end

function xdot = nonlinear_plant_dynamics(x, U_mix, P)
    Ix = P.Ix; Iy = P.Iy; Iz = P.Iz; L = P.L;
    p = x(4); q = x(5); r = x(6);
    U2 = U_mix(1); U3 = U_mix(2); U4 = U_mix(3);
    xdot_kin = [p; q; r];
    p_dot = q*r*(Iy-Iz)/Ix + (L/Ix)*U2;
    q_dot = p*r*(Iz-Ix)/Iy + (L/Iy)*U3;
    r_dot = p*q*(Ix-Iy)/Iz + (1/Iz)*U4;
    xdot = [xdot_kin; [p_dot; q_dot; r_dot]];
end
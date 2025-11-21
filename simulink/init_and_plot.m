%% Master Script: Init, Run, and Plot
% L1 Adaptive Control + LQR for Quadrotor Attitude
% ---------------------------------------------------------
clear; clc; close all;

%% ================= SECTION 1: INITIALIZATION =================

% --- 1. Simulation Settings ---
model_name = 'l1_lqr'; % <--- MAKE SURE THIS MATCHES YOUR FILE NAME
T_final    = 10;              % Simulation duration (s)

P.loop_delay = 0;

% --- 2. Geometry & Inertia ---
P.m  = 0.76;      % Mass (kg)
P.L  = 0.14;      % Arm length (m)
P.Ix = 0.0045;    % Inertia x
P.Iy = 0.0045;    % Inertia y
P.Iz = 0.0113;    % Inertia z

% % --- 3. Linear Model (Predictor Structure) ---
% % State: [phi theta psi p q r]'
% P.A_true = [zeros(3,3), eye(3); zeros(3,6)];
% P.Bm     = [zeros(3,3); diag([1/P.Ix, 1/P.Iy, 1/P.Iz])];      
% P.Bum    = [eye(3); zeros(3,3)];                         
% P.Bbar   = [P.Bm P.Bum];   
% % --- 3. Linear Model (Predictor Structure) ---
% % State: [phi theta psi p q r]'
% P.A_true = [zeros(3,3), eye(3); zeros(3,6)];
% 
% % CORRECTED Bm MATRIX
% % The input u is [RollForce, PitchForce, YawTorque].
% % We must include P.L in the B matrix for Roll/Pitch so the predictor
% % matches the plant physics (Acc = Force * L / Inertia).
% P.Bm = [zeros(3,3); 
%         P.L/P.Ix,  0,        0;       % Roll channel (Force->Acc)
%         0,         P.L/P.Iy, 0;       % Pitch channel (Force->Acc)
%         0,         0,        1/P.Iz]; % Yaw channel (Torque->Acc)
% 
% P.Bum    = [eye(3); zeros(3,3)];                         
% P.Bbar   = [P.Bm P.Bum];

% --- 3. Linear Model (Predictor Structure) ---
P.A_true = [zeros(3,3), eye(3); zeros(3,6)];

% DEFINING B MATRIX FOR FORCE INPUTS
% Roll/Pitch dynamics: Acc = (Force * L) / I
% Yaw dynamics:        Acc = Torque / I
P.Bm = [zeros(3,3); 
        P.L/P.Ix,  0,        0;      
        0,         P.L/P.Iy, 0;       
        0,         0,        1/P.Iz]; 

P.Bum    = [eye(3); zeros(3,3)];                         
P.Bbar   = [P.Bm P.Bum];  

% RE-CALCULATE GAINS WITH NEW Bm
% The solver will now see the L/I relationship and give you STRONGER gains.
% Q = diag([100, 100, 80, 1, 1, 1]); 
% R = diag([0.01, 0.01, 0.01]);     
% P.K = lqr(P.A_true, P.Bm, Q, R);

% --- 4. Reference & Disturbance Parameters ---
% Reference (Step)
P.ref_amp   = deg2rad(20); % 10 deg step
P.ref_start = 0.0;         % Step time

% Disturbance (Sinusoidal)
% Applied to Torque channels directly
P.dist_amp  = [0.01; 0.0; 0.0]; % Nm [Roll, Pitch, Yaw]
P.dist_freq = [2; 0.0; 0.0]; % rad/s

% --- 5. LQR Baseline Tuning ---
% Q = diag([100, 100, 80, 1, 1, 1]); 
% R = diag([0.1, 0.1, 0.1]);   
Q = diag([10, 10, 10,  4, 4, 3]);   % [phi theta psi dphi dtheta dpsi]
R = diag([2, 2, 2]);
P.K = lqr(P.A_true, P.Bm, Q, R);

% --- 6. L1 Adaptive Parameters ---
P.Ts = 2e-5;              % Sampling Time (Must match Adaptive Law Block)
P.Ae = -10.0 * eye(6);      % Error dynamics (A_sub_s)
P.wc = 40;                 % Filter cutoff (rad/s)

% LPF State Space
[P.Alpf, P.Blpf, P.Clpf, P.Dlpf] = deal(-P.wc*eye(3), P.wc*eye(3), eye(3), zeros(3));

% Pre-compute Closed-Form Adaptive Law Matrices
% This saves computation time inside the Simulink block
% Formula: Phi = (e^(Ae*Ts) - I) * inv(Ae)
Phi = (expm(P.Ae*P.Ts) - eye(6)) / P.Ae;
P.Phi_inv = inv(Phi);
P.Exp_AeTs = expm(P.Ae*P.Ts);

fprintf('Initialization Complete. Parameters loaded.\n');

%% ================= SECTION 2: RUN SIMULATION =================

% Check if model exists
if ~exist(model_name, 'file')
    error(['Model file "' model_name '" not found! Please save your Simulink model with this name.']);
end

fprintf('Running Simulink Model: %s ... ', model_name);

% Run Simulation
% 'out' will contain all the "To Workspace" data
out = sim(model_name, 'StopTime', num2str(T_final));

fprintf('Done.\n');

%% ================= SECTION 3: PLOTTING (FIXED) =================

% 1. Helper Function to Extract & Fix Data Orientation
%    This ensures data is always (Time x Channels)
function data_fixed = fix_data(sim_data, t)
    d = squeeze(sim_data); % Remove singleton dims (often becomes Channels x Time)
    
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

% 2. Extract Data
try
    t      = out.sim_x.Time;
    
    % Apply the fix function to each signal
    x      = fix_data(out.sim_x.Data, t);       
    ref    = fix_data(out.sim_ref.Data, t);     
    u_b    = fix_data(out.sim_u_base.Data, t);  
    u_l1   = fix_data(out.sim_u_L1.Data, t);    
    dist   = fix_data(out.sim_dist.Data, t);    
    sigma  = fix_data(out.sim_sigma.Data, t);   
    
catch ME
    % If it fails, show the REAL error message
    rethrow(ME); 
end

% 3. Unit Conversion (Mixing Inputs -> Torques)
% Plant inputs u are [Force/L, Force/L, Torque]. Convert to Nm.
tau_b_Nm  = [u_b(:,1)*P.L,  u_b(:,2)*P.L,  u_b(:,3)];
tau_l1_Nm = [u_l1(:,1)*P.L, u_l1(:,2)*P.L, u_l1(:,3)];
tau_tot   = tau_b_Nm + tau_l1_Nm;

% 4. Calculate RMSE (Roll)
err_phi = ref(:,1) - x(:,1);
rmse = sqrt(mean(err_phi.^2));
fprintf('--------------------------------------\n');
fprintf('Performance Metric (Roll RMSE): %.4f rad\n', rmse);
fprintf('--------------------------------------\n');

fig_size = [100, 100, 1000, 800];

%% --- PLOT 1: Attitude Tracking ---
figure('Color','w','Name','Attitude Tracking', 'Position', fig_size);
tiledlayout(3,1, 'Padding', 'compact');
ang_labels = {'\phi (Roll)', '\theta (Pitch)', '\psi (Yaw)'};

for i = 1:3
    nexttile; hold on; grid on;
    plot(t, rad2deg(x(:,i)), 'LineWidth', 2, 'Color', [0 0.447 0.741]);
    plot(t, rad2deg(ref(:,i)), 'k--', 'LineWidth', 1.5);
    ylabel([ang_labels{i} ' (deg)']);
    if i==1, title(['Attitude Tracking (RMSE_{\phi} = ' num2str(rmse, '%.4f') ')']); legend('Actual','Reference'); end
end
xlabel('Time (s)');

%% --- PLOT 2: Torque Breakdown ---
figure('Color','w','Name','Torque Breakdown', 'Position', fig_size);
tiledlayout(3,1, 'Padding', 'compact');
tau_labels = {'\tau_\phi', '\tau_\theta', '\tau_\psi'};

for i = 1:3
    nexttile; hold on; grid on;
    l1 = plot(t, tau_b_Nm(:,i), '--', 'Color', [0.466 0.674 0.188], 'LineWidth', 1.5); 
    l2 = plot(t, tau_l1_Nm(:,i), '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5); 
    l3 = plot(t, tau_tot(:,i), 'k-', 'LineWidth', 1.0);
    l4 = plot(t, dist(:,i), 'm:', 'LineWidth', 2.0);
    
    ylabel([tau_labels{i} ' (Nm)']);
    if i==1, title('Control Effort Breakdown'); legend([l1 l2 l3 l4], 'LQR', 'L1', 'Total', 'Disturbance', 'Location','best'); end
end
xlabel('Time (s)');

%% --- PLOT 3: Estimation Accuracy (CORRECTED SCALING) ---
figure('Color','w','Name','L1 Estimation', 'Position', fig_size);
tiledlayout(3,1, 'Padding', 'compact');

for i = 1:3
    nexttile; hold on; grid on;
    
    % SCALE SIGMA TO MATCH TORQUE
    % Sigma is estimated in "Input Space" (Force).
    % Disturbance is in "Torque Space" (Nm).
    % Force * L = Torque.
    
    if i == 3 
        scaler = 1.0; % Yaw input is already Torque
    else
        scaler = P.L; % Roll/Pitch inputs are Forces
    end
    
    estimate_Nm = sigma(:,i) * scaler; 
    
    plot(t, estimate_Nm, 'b', 'LineWidth', 2);     
    plot(t, dist(:,i), 'r--', 'LineWidth', 1.5);     
    
    ylabel(['Dist. ' tau_labels{i} ' (Nm)']);
    if i==1
        title('Uncertainty Estimation Accuracy'); 
        legend('L1 Estimate (Scaled)', 'True Disturbance', 'Location', 'best'); 
    end
end
xlabel('Time (s)');
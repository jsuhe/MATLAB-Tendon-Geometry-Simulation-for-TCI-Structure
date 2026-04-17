%% =========================================================================
%  TENDON-CONSTRAINED 3D RIGID PLATE — XZ BASE PLANE, PERIMETER RING
%  =========================================================================
%  Key fixes in this version:
%    - Orientation stored as a UNIT QUATERNION q = [w; x; y; z]
%      Integrated via  q_dot = 0.5 * quat_mult(q, [0; omega])
%      This avoids the drift error of the naive  rvec += omega*dt  approach
%      and keeps the plate normal physically correct at all times.
%    - Plate normal is re-derived from q at EVERY sub-step (no lag).
%    - Slack condition is explicit: tension = 0 whenever L <= L_slack(i),
%      where L_slack is a per-tendon threshold you set directly [m].
%      L_slack = L0 means zero-tension whenever the tendon is at or shorter
%      than its rest length (classical Hookean + unilateral constraint).
%      Set L_slack > L0 to introduce a "dead-band" where the tendon hangs
%      loose above its rest length but below L_slack.
%
%  Modular structure:
%    1. PARAMETERS          — plate, integration, pressure, damping
%    2. GEOMETRY            — XZ base anchors, plate local frame
%    3. TENDON PROPERTIES   — stiffness, absolute rest lengths, slack threshold
%    4. PERIMETER RING      — edge tether geometry and properties
%    5. INTEGRATION         — explicit Euler with quaternion orientation
%    6. DATA EXPORT         — JSON for tendon_visualizer.html
%    7. ANALYSIS PLOTS
%  =========================================================================

clear; clc; close all;

%% =========================================================================
%  SECTION 1 — PARAMETERS
%  =========================================================================

% --- Time integration ------------------------------------------------------
dt    = 0.002;   % Time step [s]  (reduce for higher accuracy)
t_end = 10.0;     % Total simulation time [s]

% --- Plate -----------------------------------------------------------------
plate_size = 0.25;     % Half-side of square plate [m]
plate_mass = 0.4;      % [kg]
I_plate    = plate_mass * plate_size^2 / 6;   % Isotropic MOI [kg·m²]

% --- Initial conditions ----------------------------------------------------
%  Y is "up". Plate starts 0.15 m above the XZ base plane, flat (no tilt).
pos0   = [0;  0.15;  0];
vel0   = [0;  0;     0];
q0     = [1;  0;  0;  0];   % Unit quaternion: w=1, x=y=z=0 → identity (no rotation)
omega0 = [0;  0;     0];

% --- Pressure --------------------------------------------------------------
%  Acts along the plate normal. Linear ramp: P(t) = P_rate * t [Pa].
P_rate     = 100.0;              % [Pa/s]
plate_area = (2*plate_size)^2;  % [m²]

% --- Gravity ---------------------------------------------------------------
%  Acts in -Y (down, away from the plate).
use_gravity = true;
g = 9.81;   % [m/s²]

% --- Damping ---------------------------------------------------------------
c_lin = 2.5;    % Linear damping  [N·s/m]
c_ang = 0.06;   % Angular damping [N·m·s/rad]

%% =========================================================================
%  SECTION 2 — THREE-TENDON GEOMETRY  (XZ base plane)
%  =========================================================================
%  Base anchors are fixed at Y = 0 (the XZ plane).
%  Plate attachment points are in the plate's LOCAL frame.
%  The plate face lies in local XZ (local Y = 0 for all attach points),
%  so the unrotated plate normal is local +Y = world +Y at t=0.
%
%  Triangular 120° layout — all three tendons geometrically equivalent.

angles_base  = deg2rad([0, 120, 240]);
angles_plate = deg2rad([0, 120, 240]);

r_base  = plate_size * 1.3;   % Radius of base anchor ring  [m]
r_plate = plate_size * 0.9;   % Radius of plate attach ring [m]

% Base anchors in XZ plane (Y = 0)
base_pts = [ r_base * cos(angles_base);   % world X
             zeros(1, 3);                 % world Y = 0
             r_base * sin(angles_base) ]; % world Z

% Plate attachment points in LOCAL plate frame (local XZ face, local Y = 0)
plate_local = [ r_plate * cos(angles_plate);
                zeros(1, 3);
                r_plate * sin(angles_plate) ];

n_tendons = 3;

%% =========================================================================
%  SECTION 3 — TENDON STIFFNESS, ABSOLUTE REST LENGTHS & SLACK THRESHOLD
%  =========================================================================
%
%  For each tendon i:
%    L         = current length
%    L0(i)     = natural rest length [m]  — Hookean zero-force point
%    L_slack(i)= slack threshold [m]      — tension is ZERO for L <= L_slack(i)
%
%  Force law:
%    if L <= L_slack(i)  → F = 0           (tendon is slack / unloaded)
%    if L >  L_slack(i)  → F = k(i)*(L - L0(i)) + k_nl(i)*(L-L0(i))^3
%
%  Setting L_slack(i) = L0(i) gives the classical unilateral Hookean spring
%  (slack when compressed, taut when stretched beyond rest length).
%
%  Setting L_slack(i) > L0(i) creates a dead-band: the tendon hangs loose
%  between L0 and L_slack before any tension develops.

% Compute initial tendon lengths at pos0 with identity rotation
init_L = zeros(1, n_tendons);
for i = 1:n_tendons
    init_L(i) = norm((pos0 + plate_local(:,i)) - base_pts(:,i));
end
fprintf('Initial tendon lengths:  %.4f  %.4f  %.4f m\n', init_L);

% --- Stiffness [N/m] -------------------------------------------------------
k = [60, 180, 350];

% --- Natural rest lengths [m] ----------------------------------------------
%  Set shorter than init_L so tendons are pre-stretched at t=0.
L0 = [0.10, 0.20, 0.30];

% --- Slack threshold [m] ---------------------------------------------------
%  Classical setting: L_slack = L0  (tension = 0 iff L <= L0).
%  To make tendon 2 have a dead-band above its rest length:
%    L_slack = [L0(1),  L0(2)+0.03,  L0(3)]
%  To make tendon 3 start fully slack at t=0 (since init_L(3) < L_slack(3)):
%    L_slack(3) = init_L(3) + 0.01
L_slack = L0;   % [m]  — set equal to L0 for classical behaviour

fprintf('Rest lengths (absolute): %.4f  %.4f  %.4f m\n', L0);
fprintf('Slack thresholds:        %.4f  %.4f  %.4f m\n', L_slack);
fprintf('Initial stretch:         %.4f  %.4f  %.4f m\n', init_L - L0);
fprintf('Taut at t=0:             %d  %d  %d\n\n', ...
    init_L(1)>L_slack(1), init_L(2)>L_slack(2), init_L(3)>L_slack(3));

% --- Nonlinear stiffness: F = k*(L-L0) + k_nl*(L-L0)^3 -------------------
k_nl = [0, 0, 0];

%% =========================================================================
%  SECTION 4 — PERIMETER RING
%  =========================================================================
%  N_perim elastic tethers connect points uniformly around the plate edge
%  to anchor points on the XZ base plane.
%  Each perimeter tether also obeys the slack threshold L_slack_perim:
%  zero tension when L_j <= L_slack_perim.

N_perim        = 12;            % Number of perimeter tethers
k_perim        = 30;            % Stiffness [N/m]
L0_perim       = 0.08;          % Natural rest length [m]
L_slack_perim  = L0_perim;      % Slack threshold [m] (= L0 for classical)
r_perim_plate  = plate_size * 1.0;   % Radius on plate local frame [m]
r_perim_base   = plate_size * 1.5;   % Radius on XZ base plane [m]

perim_angles = linspace(0, 2*pi*(1 - 1/N_perim), N_perim);

perim_local = [ r_perim_plate * cos(perim_angles);
                zeros(1, N_perim);
                r_perim_plate * sin(perim_angles) ];

perim_base = [ r_perim_base * cos(perim_angles);
               zeros(1, N_perim);
               r_perim_base * sin(perim_angles) ];

fprintf('Perimeter: %d tethers, k=%.0f N/m, L0=%.3f m, L_slack=%.3f m\n\n', ...
    N_perim, k_perim, L0_perim, L_slack_perim);

%% =========================================================================
%  SECTION 5 — TIME INTEGRATION  (Explicit Euler + quaternion orientation)
%  =========================================================================
%
%  State: pos [3], vel [3], q [4] (unit quaternion), omega [3]  → 13 values
%
%  Quaternion kinematics (exact, no drift):
%    q_dot = 0.5 * quat_mult(q, [0; omega_world])
%  where omega_world is the angular velocity in the WORLD frame.
%  After each step: re-normalise q to prevent floating-point drift.
%
%  Rotation matrix from quaternion (used to rotate local→world):
%    R = quat_to_R(q)   (see local functions below)
%
%  Plate normal in world frame:
%    normal = R * [0; 1; 0]   (local +Y = plate outward normal)

n_steps      = round(t_end / dt);
t_hist       = zeros(1, n_steps);
pos_hist     = zeros(3, n_steps);
vel_hist     = zeros(3, n_steps);
q_hist       = zeros(4, n_steps);
omega_hist   = zeros(3, n_steps);
tension_hist = zeros(3, n_steps);
perim_hist   = zeros(N_perim, n_steps);
normal_hist  = zeros(3, n_steps);
slack_hist   = zeros(3, n_steps);   % 1 = taut, 0 = slack (per tendon)

pos = pos0;  vel = vel0;  q = q0;  omega = omega0;
pos_hist(:,1)    = pos;
q_hist(:,1)      = q;
normal_hist(:,1) = [0; 1; 0];

fprintf('Integrating %d steps...\n', n_steps);

for step = 2:n_steps
    t = (step-1) * dt;

    % ---- Rotation matrix from current quaternion --------------------------
    R = quat_to_R(q);

    % ---- Plate normal: local +Y rotated to world -------------------------
    %  This is always correct because R is derived from the normalised
    %  quaternion — no approximation error accumulates here.
    normal = R * [0; 1; 0];

    % ---- Pressure (linear ramp, acts along plate normal) -----------------
    P = P_rate * t;
    F_pressure = P * plate_area * normal;

    % ---- Gravity (world -Y) ----------------------------------------------
    F_grav = [0; -plate_mass * g * use_gravity; 0];

    % ---- Main tendon forces ----------------------------------------------
    F_tend_sum = zeros(3,1);
    T_tend_sum = zeros(3,1);

    for i = 1:n_tendons
        p_att = pos + R * plate_local(:,i);   % Attachment point, world frame
        tvec  = p_att - base_pts(:,i);        % Vector from anchor to attachment
        L     = norm(tvec);

        % Slack condition: zero tension if L is at or below slack threshold
        if L > L_slack(i)
            ext   = L - L0(i);               % Extension beyond rest length
            F_mag = k(i)*ext + k_nl(i)*ext^3;
            tension_hist(i, step) = F_mag;
            slack_hist(i, step)   = 1;        % Taut
            F_i = -F_mag * (tvec / L);        % Pulls attachment toward anchor
            r_w = R * plate_local(:,i);       % Moment arm in world frame
            F_tend_sum = F_tend_sum + F_i;
            T_tend_sum = T_tend_sum + cross(r_w, F_i);
        end
        % else: L <= L_slack → tendon is slack, F = 0, no contribution
    end

    % ---- Perimeter tether forces -----------------------------------------
    F_perim_sum = zeros(3,1);
    T_perim_sum = zeros(3,1);

    for j = 1:N_perim
        p_att_j = pos + R * perim_local(:,j);
        tvec_j  = p_att_j - perim_base(:,j);
        L_j     = norm(tvec_j);

        if L_j > L_slack_perim
            ext_j   = L_j - L0_perim;
            F_mag_j = k_perim * ext_j;
            perim_hist(j, step) = F_mag_j;
            F_j   = -F_mag_j * (tvec_j / L_j);
            r_w_j = R * perim_local(:,j);
            F_perim_sum = F_perim_sum + F_j;
            T_perim_sum = T_perim_sum + cross(r_w_j, F_j);
        end
    end

    % ---- Damping ---------------------------------------------------------
    F_damp = -c_lin * vel;
    T_damp = -c_ang * omega;

    % ---- Net force & torque ----------------------------------------------
    F_net = F_pressure + F_grav + F_tend_sum + F_perim_sum + F_damp;
    T_net = T_tend_sum + T_perim_sum + T_damp;

    % ---- Translational update (Newton) -----------------------------------
    acc = F_net / plate_mass;
    vel = vel + acc * dt;
    pos = pos + vel * dt;

    % ---- Rotational update (quaternion kinematics) -----------------------
    %  Equation: q_dot = 0.5 * q ⊗ omega_quat
    %  where omega_quat = [0; omega]  (pure quaternion from angular velocity)
    %
    %  This is the exact kinematic equation for unit quaternions. It keeps
    %  the normal vector correct regardless of rotation magnitude.
    alpha     = T_net / I_plate;
    omega     = omega + alpha * dt;
    omega_q   = [0; omega];                   % Pure quaternion
    q_dot     = 0.5 * quat_mult(q, omega_q);
    q         = q + q_dot * dt;
    q         = q / norm(q);                  % Re-normalise every step

    % ---- Store -----------------------------------------------------------
    t_hist(step)        = t;
    pos_hist(:,step)    = pos;
    vel_hist(:,step)    = vel;
    q_hist(:,step)      = q;
    omega_hist(:,step)  = omega;
    normal_hist(:,step) = normal;
end

fprintf('Integration complete.\n\n');

% Convert quaternion history to Euler angles (for plots/export)
euler_hist = zeros(3, n_steps);   % [roll; pitch; yaw] in degrees
for s = 1:n_steps
    euler_hist(:,s) = quat_to_euler_deg(q_hist(:,s));
end

%% =========================================================================
%  SECTION 6 — DATA EXPORT  (JSON for tendon_visualizer.html)
%  =========================================================================

export_every = 4;
json_frames  = '';
frame_count  = 0;

for s = 1:export_every:n_steps
    pos_s  = pos_hist(:, s);
    q_s    = q_hist(:, s);
    norm_s = normal_hist(:, s);
    tens_s = tension_hist(:, s);
    slk_s  = slack_hist(:, s);
    t_s    = t_hist(s);
    perim_mean = mean(perim_hist(:, s));

    % Export quaternion directly — visualizer uses it to rebuild R
    frame_str = sprintf( ...
        ['{"pos":[%.6f,%.6f,%.6f],' ...
         '"q":[%.8f,%.8f,%.8f,%.8f],' ...
         '"normal":[%.6f,%.6f,%.6f],' ...
         '"tens":[%.4f,%.4f,%.4f],' ...
         '"slack":[%d,%d,%d],' ...
         '"perim_mean":%.4f,"t":%.4f}'], ...
        pos_s(1), pos_s(2), pos_s(3), ...
        q_s(1),   q_s(2),   q_s(3),   q_s(4), ...
        norm_s(1),norm_s(2),norm_s(3), ...
        tens_s(1),tens_s(2),tens_s(3), ...
        slk_s(1), slk_s(2), slk_s(3), ...
        perim_mean, t_s);

    if frame_count == 0; json_frames = frame_str;
    else;                json_frames = [json_frames, ',', frame_str];
    end
    frame_count = frame_count + 1;
end

bp_str       = mat_to_json_array(base_pts);
pl_str       = mat_to_json_array(plate_local);
pb_str       = mat_to_json_array(perim_base);
pl_perim_str = mat_to_json_array(perim_local);

json_out = sprintf( ...
    ['{\n' ...
     '  "plate_size": %.6f,\n' ...
     '  "plate_normal_axis": "y",\n' ...
     '  "k": [%s],\n' ...
     '  "L0": [%s],\n' ...
     '  "L_slack": [%s],\n' ...
     '  "k_perim": %.2f,\n' ...
     '  "L0_perim": %.4f,\n' ...
     '  "N_perim": %d,\n' ...
     '  "base_pts": %s,\n' ...
     '  "plate_loc": %s,\n' ...
     '  "perim_base": %s,\n' ...
     '  "perim_local": %s,\n' ...
     '  "frames": [%s]\n' ...
     '}'], ...
    plate_size, ...
    num2str(k,       '%.2f,'), ...
    num2str(L0,      '%.6f,'), ...
    num2str(L_slack, '%.6f,'), ...
    k_perim, L0_perim, N_perim, ...
    bp_str, pl_str, pb_str, pl_perim_str, ...
    json_frames);

json_out = regexprep(json_out, ',\s*}', '}');
json_out = regexprep(json_out, ',\s*]', ']');

fid = fopen('tendon_sim.json', 'w');
fprintf(fid, '%s', json_out);
fclose(fid);
fprintf('Exported %d frames to tendon_sim.json\n', frame_count);

%% =========================================================================
%  SECTION 7 — ANALYSIS PLOTS
%  =========================================================================

figure('Name','Tendon Analysis', 'Color','w', 'Position',[50 50 1200 800]);

subplot(2,3,1);
plot(t_hist, pos_hist(1,:),'r', t_hist, pos_hist(2,:),'g', ...
     t_hist, pos_hist(3,:),'b','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Position [m]');
legend('X','Y (up)','Z'); title('Plate centroid position'); grid on;

subplot(2,3,2);
plot(t_hist, euler_hist(1,:),'r', ...
     t_hist, euler_hist(2,:),'g', ...
     t_hist, euler_hist(3,:),'b','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Angle [deg]');
legend('Roll','Pitch','Yaw'); title('Plate orientation (Euler angles)'); grid on;

subplot(2,3,3);
cols = {'#E03030','#2060E0','#20A050'};
hold on;
for i = 1:n_tendons
    h = plot(t_hist, tension_hist(i,:), 'Color', cols{i}, 'LineWidth', 1.8, ...
         'DisplayName', sprintf('T%d  k=%d  L0=%.2f  Ls=%.2f', ...
         i, k(i), L0(i), L_slack(i)));
    % Shade slack regions
    slack_on = slack_hist(i,:) == 0;
    if any(slack_on)
        patch([t_hist(slack_on), fliplr(t_hist(slack_on))], ...
              [zeros(1,sum(slack_on)), zeros(1,sum(slack_on))], ...
              cols{i}, 'FaceAlpha', 0.08, 'EdgeColor','none', ...
              'HandleVisibility','off');
    end
end
xlabel('Time [s]'); ylabel('Tension [N]');
legend('Location','northwest'); title('Main tendon tensions  (0 = slack)'); grid on;
hold off;

subplot(2,3,4);
perim_total = sum(perim_hist, 1);
plot(t_hist, perim_total,    'k',  'LineWidth',1.5,'DisplayName','Total'); hold on;
plot(t_hist, mean(perim_hist,1),'--','Color',[0.5 0.5 0.5],'LineWidth',1.2,...
     'DisplayName','Mean/tether');
xlabel('Time [s]'); ylabel('Force [N]');
legend; title(sprintf('Perimeter ring (%d tethers)', N_perim)); grid on;

subplot(2,3,5);
plot(t_hist, normal_hist(1,:),'r', ...
     t_hist, normal_hist(2,:),'g', ...
     t_hist, normal_hist(3,:),'b','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Component');
legend('nx','ny (normal)','nz'); title('Plate normal vector'); grid on;
ylim([-1.2 1.2]);

subplot(2,3,6);
speed = sqrt(sum(vel_hist.^2, 1));
plot(t_hist, speed, 'k', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Plate centroid speed'); grid on;

fprintf('--- Equilibrium Summary ---\n');
fprintf('Final position:    [%.4f  %.4f  %.4f] m\n', pos_hist(:,end));
fprintf('Final orientation: [%.2f  %.2f  %.2f] deg (roll/pitch/yaw)\n', euler_hist(:,end));
fprintf('Final speed:       %.5f m/s\n', norm(vel_hist(:,end)));
fprintf('Final tensions:    [%.2f  %.2f  %.2f] N\n', tension_hist(:,end));
fprintf('Final slack state: [%d  %d  %d]  (1=taut, 0=slack)\n', slack_hist(:,end));
fprintf('Perim mean/total:  %.2f / %.2f N\n', ...
    mean(perim_hist(:,end)), sum(perim_hist(:,end)));

%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

function R = quat_to_R(q)
% QUAT_TO_R  Rotation matrix from unit quaternion q = [w; x; y; z].
%   Equivalent to Rodrigues but exact and singularity-free.
    w=q(1); x=q(2); y=q(3); z=q(4);
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),   2*(x*z+w*y);
           2*(x*y+w*z), 1-2*(x^2+z^2),   2*(y*z-w*x);
           2*(x*z-w*y),   2*(y*z+w*x), 1-2*(x^2+y^2)];
end

function qout = quat_mult(p, q)
% QUAT_MULT  Hamilton product of two quaternions p, q = [w; x; y; z].
    pw=p(1); px=p(2); py=p(3); pz=p(4);
    qw=q(1); qx=q(2); qy=q(3); qz=q(4);
    qout = [ pw*qw - px*qx - py*qy - pz*qz;
             pw*qx + px*qw + py*qz - pz*qy;
             pw*qy - px*qz + py*qw + pz*qx;
             pw*qz + px*qy - py*qx + pz*qw ];
end

function eul = quat_to_euler_deg(q)
% QUAT_TO_EULER_DEG  Convert unit quaternion to [roll; pitch; yaw] in degrees.
%   Convention: ZYX (yaw-pitch-roll), intrinsic rotations.
    w=q(1); x=q(2); y=q(3); z=q(4);
    roll  = atan2(2*(w*x+y*z), 1-2*(x^2+y^2));
    pitch = asin( max(-1, min(1, 2*(w*y-z*x))) );
    yaw   = atan2(2*(w*z+x*y), 1-2*(y^2+z^2));
    eul   = [roll; pitch; yaw] * 180/pi;
end

function S = skew(v)
    S = [ 0    -v(3)  v(2);
          v(3)  0    -v(1);
         -v(2)  v(1)  0   ];
end

function s = mat_to_json_array(M)
% Converts a 3×N matrix into a JSON array of [x,y,z] triplets.
    entries = '';
    for col = 1:size(M, 2)
        entry = sprintf('[%.6f,%.6f,%.6f]', M(1,col), M(2,col), M(3,col));
        if col == 1; entries = entry;
        else;        entries = [entries, ',', entry];
        end
    end
    s = ['[', entries, ']'];
end

%% =========================================================================
%  EXPERIMENT SUGGESTIONS
%  =========================================================================
%
%  1. SLACK DEAD-BAND
%     L_slack = [L0(1), L0(2)+0.03, L0(3)]
%     Tendon 2 hangs loose for 3 cm above its rest length before tensioning.
%     Watch the tension plot for the delayed onset.
%
%  2. TENDON STARTS FULLY SLACK
%     Set L_slack(3) = init_L(3) + 0.02  (computed after Section 2 runs).
%     Tendon 3 produces zero tension until the plate has risen/tilted enough.
%
%  3. ASYMMETRIC REST LENGTHS
%     L0 = [0.05, 0.10, 0.14] — T1 very tight, T3 lightly pre-loaded.
%     Combined with different k values → preferred tilt direction.
%
%  4. PERIMETER DENSITY
%     N_perim = 24  → smooth boundary torque
%     N_perim = 3   → triangular perimeter, pronounced angular asymmetry
%
%  5. NONLINEAR PERIMETER STIFFENING
%     Replace in the perimeter loop:
%       F_mag_j = k_perim * ext_j + 500 * ext_j^3;
%
%  6. PRESSURE RATE
%     P_rate = 20  → slow quasi-static rise, near-equilibrium at each step
%     P_rate = 200 → impulsive, ring-down transients visible in speed plot
%
%  7. VERIFY QUATERNION NORMALITY
%     After simulation: max(abs(vecnorm(q_hist) - 1)) should be < 1e-6.
%     Add this check: fprintf('Max quat norm error: %.2e\n', max(abs(vecnorm(q_hist)-1)));
%
%  =========================================================================
% tendon_3D_asymmetric.m
% Tendon-Constrained 3D Rigid Plate Simulation
% Exports tendon_sim.json for use with tendon_visualizer.html
%
% Physics:  Explicit Euler, 12-DOF (pos, vel, rvec, omega)
% Rotation: Rodrigues rotation vector (NOT quaternions)
% Axes:     XY base plane, Z up, gravity in -Z
% Plate:    square rigid body, normal along local +Z
%
% Perimeter walls:
%   - Contact distance measured in XY plane only (walls are vertical)
%   - Force magnitude scaled inversely with contact point height:
%     f = k_perim * (L_contact - L_xy) * (z0 / max(att_z, 0.01))
%     so a lower contact point receives a larger outward push
%
% Run this script, then open tendon_visualizer.html and load the JSON.

clear; clc;

%% PARAMETERS

% Plate
plate_size = 0.25;          % half-side length [m]
mass       = 0.4;           % [kg]
I_plate    = mass * plate_size^2 / 6;   % isotropic MOI

% Integration
dt    = 0.002;              % time step [s]
t_end = 6.0;                % end time [s]
skip  = 4;                  % export every Nth frame (~750 frames total)

% Environment
g      = 9.81;              % gravity [m/s^2]
P0     = 0;                 % initial pressure at t=0 [Pa]
P_rate = 50;                % pressure ramp rate [Pa/s]

% Damping
c_lin = 2.5;                % linear damping [N.s/m]
c_ang = 0.06;               % angular damping [N.m.s]

%% INITIAL STATE

z0    = 0.15;               % initial plate centroid height [m]
pos   = [0; 0; z0];         % plate centroid [m]
vel   = [0; 0; 0];          % linear velocity [m/s]
rvec  = [0; 0; 0];          % Rodrigues rotation vector
omega = [0; 0; 0];          % angular velocity [rad/s]

%% TENDON CONFIG

base_coords = [
     0.325,  0.000;
    -0.163,  0.282;
    -0.163, -0.282
];

plate_coords = [
     0.225,  0.000;
    -0.113,  0.195;
    -0.113, -0.195
];

n_tendons = size(base_coords, 1);

% Tendon stiffness [N/m]
k_ten = [60, 180, 350];

% Rest lengths [m] -- set directly as absolute values.
% Leave empty [] to auto-compute from L0_frac at the initial geometry.
L0_abs  = [];
L0_frac = [0.50, 0.60, 0.90]; % fractions of initial length (used only if L0_abs is empty)

% Nonlinear cubic term [N/m^3] -- set to 0 to disable
k_nl = [0, 0, 0];

%% PERIMETER COMPRESSION WALLS

N_perim         = 12;       % number of evenly-spaced contact points
k_perim         = 120;      % wall stiffness [N/m]
L_contact_perim = 0.12;     % XY-plane contact activation distance [m]

% Perimeter geometry
perim_angles = linspace(0, 2*pi, N_perim+1); perim_angles(end) = [];
perim_base  = [plate_size*1.5*cos(perim_angles)', plate_size*1.5*sin(perim_angles)', zeros(N_perim,1)];
perim_local = [plate_size*cos(perim_angles)',     plate_size*sin(perim_angles)',     zeros(N_perim,1)];

%% BUILD ANCHOR POINT ARRAYS

base_pts  = [base_coords,  zeros(n_tendons, 1)];  % [X, Y, 0]
plate_loc = [plate_coords, zeros(n_tendons, 1)];  % local [X, Y, 0]

%% COMPUTE TENDON REST LENGTHS

R0 = rot_rodrigues(rvec);
if isempty(L0_abs)
    L0 = zeros(1, n_tendons);
    for i = 1:n_tendons
        att   = pos + R0 * plate_loc(i,:)';
        L0(i) = norm(att - base_pts(i,:)') * L0_frac(i);
    end
else
    L0 = L0_abs(:)';
end
L_slack = L0;

fprintf('Rest lengths L0: '); fprintf('%.4f  ', L0); fprintf('\n');

%% TIME INTEGRATION

N_steps  = round(t_end / dt);
N_export = floor(N_steps / skip) + 1;

out_pos    = zeros(N_export, 3);
out_rv     = zeros(N_export, 3);
out_normal = zeros(N_export, 3);
out_tens   = zeros(N_export, n_tendons);
out_slack  = zeros(N_export, n_tendons);
out_perim  = zeros(N_export, 1);
out_t      = zeros(N_export, 1);
ei = 0;

fprintf('Running simulation: %d steps (%.1f s) ...\n', N_steps, t_end);

for step = 0:N_steps-1
    t = step * dt;
    R = rot_rodrigues(rvec);
    normal     = R * [0; 0; 1];
    plate_area = 4 * plate_size^2;
    P          = P0 + P_rate * t;

    % Forces
    F_pressure = P * plate_area * normal;
    F_gravity  = [0; 0; -mass * g];
    F_damp     = -c_lin * vel;
    T_damp     = -c_ang * omega;

    F_net = F_pressure + F_gravity + F_damp;
    T_net = T_damp;

    tens_vals  = zeros(1, n_tendons);
    slack_vals = zeros(1, n_tendons);

    % Tendon forces
    for i = 1:n_tendons
        r_w  = R * plate_loc(i,:)';
        att  = pos + r_w;
        tvec = att - base_pts(i,:)';
        L    = norm(tvec);

        if L <= L_slack(i)
            slack_vals(i) = 1;
        else
            dL    = L - L0(i);
            f_mag = k_ten(i)*dL + k_nl(i)*dL^3;
            F_i   = -f_mag * tvec / L;
            F_net = F_net + F_i;
            T_net = T_net + cross(r_w, F_i);
            tens_vals(i) = f_mag;
        end
    end

    % Perimeter wall forces
    % -- Contact distance is XY-only (walls are vertical).
    % -- Force is inversely proportional to contact point height:
    %    height_scale = z0 / max(att_z, 0.01)
    perim_fs = zeros(N_perim, 1);
    for j = 1:N_perim
        r_w = R * perim_local(j,:)';
        att = pos + r_w;

        % XY-only distance to base anchor
        dx   = att(1) - perim_base(j,1);
        dy   = att(2) - perim_base(j,2);
        L_xy = sqrt(dx*dx + dy*dy);

        if L_xy < L_contact_perim && L_xy > 1e-9
            % Inverse-height scaling: z0 height gives scale = 1
            att_z        = max(att(3), 0.01);
            height_scale = z0 / att_z;

            f_mag = k_perim * (L_contact_perim - L_xy) * height_scale;
            F_i   = f_mag * [dx; dy; 0] / L_xy;  % outward in XY, zero Z
            F_net = F_net + F_i;
            T_net = T_net + cross(r_w, F_i);
            perim_fs(j) = f_mag;
        end
    end
    perim_mean = mean(perim_fs);

    % Explicit Euler integration
    acc   = F_net / mass;
    vel   = vel + acc * dt;
    pos   = pos + vel * dt;
    alpha = T_net / I_plate;
    omega = omega + alpha * dt;
    rvec  = rvec + omega * dt;
    rvec  = wrap_rvec(rvec);

    % Export every skip-th frame
    if mod(step, skip) == 0
        ei = ei + 1;
        out_pos(ei,:)    = pos';
        out_rv(ei,:)     = rvec';
        out_normal(ei,:) = normal';
        out_tens(ei,:)   = tens_vals;
        out_slack(ei,:)  = slack_vals;
        out_perim(ei)    = perim_mean;
        out_t(ei)        = t;
    end
end

% Trim
out_pos    = out_pos(1:ei,:);
out_rv     = out_rv(1:ei,:);
out_normal = out_normal(1:ei,:);
out_tens   = out_tens(1:ei,:);
out_slack  = out_slack(1:ei,:);
out_perim  = out_perim(1:ei);
out_t      = out_t(1:ei);

fprintf('Done. Exported %d frames.\n', ei);

%% BUILD JSON

fprintf('Writing tendon_sim.json ...\n');

fid = fopen('tendon_sim.json', 'w');
fprintf(fid, '{\n');
fprintf(fid, '  "plate_size": %g,\n',         plate_size);
fprintf(fid, '  "n_tendons": %d,\n',           n_tendons);
fprintf(fid, '  "k": %s,\n',                   json_vec(k_ten));
fprintf(fid, '  "L0": %s,\n',                  json_vec(L0));
fprintf(fid, '  "L_slack": %s,\n',             json_vec(L_slack));
fprintf(fid, '  "k_perim": %g,\n',             k_perim);
fprintf(fid, '  "L_contact_perim": %g,\n',     L_contact_perim);
fprintf(fid, '  "N_perim": %d,\n',             N_perim);
fprintf(fid, '  "P0": %g,\n',                  P0);
fprintf(fid, '  "P_rate": %g,\n',              P_rate);
fprintf(fid, '  "z0": %g,\n',                  z0);
fprintf(fid, '  "base_pts": %s,\n',            json_mat(base_pts));
fprintf(fid, '  "plate_loc": %s,\n',           json_mat(plate_loc));
fprintf(fid, '  "perim_base": %s,\n',          json_mat(perim_base));
fprintf(fid, '  "perim_local": %s,\n',         json_mat(perim_local));
fprintf(fid, '  "frames": [\n');
for i = 1:ei
    fprintf(fid, '    {\n');
    fprintf(fid, '      "t": %.4f,\n',         out_t(i));
    fprintf(fid, '      "pos": %s,\n',         json_vec(out_pos(i,:)));
    fprintf(fid, '      "rv": %s,\n',          json_vec(out_rv(i,:)));
    fprintf(fid, '      "normal": %s,\n',      json_vec(out_normal(i,:)));
    fprintf(fid, '      "tens": %s,\n',        json_vec(out_tens(i,:)));
    fprintf(fid, '      "slack": %s,\n',       json_vec(out_slack(i,:)));
    fprintf(fid, '      "perim_mean": %.6f\n', out_perim(i));
    if i < ei
        fprintf(fid, '    },\n');
    else
        fprintf(fid, '    }\n');
    end
end
fprintf(fid, '  ]\n}\n');
fclose(fid);

fprintf('Saved tendon_sim.json  (%d frames, %d tendons)\n', ei, n_tendons);
fprintf('Open tendon_visualizer.html and click "Load tendon_sim.json"\n');

%% QUICK PLOT

figure('Name','Tendon Simulation Results','NumberTitle','off');

subplot(3,1,1);
plot(out_t, out_pos(:,3), 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Z [m]');
title('Plate height (centroid Z)'); grid on;

subplot(3,1,2);
hold on;
cols = lines(n_tendons);
for i = 1:n_tendons
    plot(out_t, out_tens(:,i), 'Color', cols(i,:), 'LineWidth', 1.2);
end
plot(out_t, out_perim, '--', 'Color', [0.4 0.7 1], 'LineWidth', 1);
lbl = arrayfun(@(i) sprintf('T%d', i), 1:n_tendons, 'UniformOutput', false);
legend([lbl, {'Wall mean'}], 'Location', 'northwest');
xlabel('Time [s]'); ylabel('Force [N]');
title('Tendon tensions + perimeter wall (height-scaled, XY contact)'); grid on; hold off;

subplot(3,1,3);
rpy = zeros(ei, 3);
for i = 1:ei
    Ri       = rot_rodrigues(out_rv(i,:)');
    rpy(i,1) = atan2(Ri(3,2), Ri(3,3))                       * 180/pi;
    rpy(i,2) = atan2(-Ri(3,1), sqrt(Ri(3,2)^2+Ri(3,3)^2))   * 180/pi;
    rpy(i,3) = atan2(Ri(2,1), Ri(1,1))                       * 180/pi;
end
plot(out_t, rpy, 'LineWidth', 1.2);
legend('Roll', 'Pitch', 'Yaw', 'Location', 'northwest');
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Plate orientation (roll/pitch/yaw)'); grid on;

sgtitle('Tendon-Constrained 3D Rigid Plate — Simulation Results');

%% LOCAL FUNCTIONS

function R = rot_rodrigues(rv)
    th = norm(rv);
    if th < 1e-12; R = eye(3); return; end
    k = rv / th;
    s = sin(th); c = cos(th); C = 1 - c;
    kx = k(1); ky = k(2); kz = k(3);
    R = [c+kx*kx*C,     kx*ky*C-kz*s, kx*kz*C+ky*s;
         ky*kx*C+kz*s,  c+ky*ky*C,    ky*kz*C-kx*s;
         kz*kx*C-ky*s,  kz*ky*C+kx*s, c+kz*kz*C  ];
end

function rv_out = wrap_rvec(rv)
    th = norm(rv);
    if th < 1e-12; rv_out = rv; return; end
    w = mod(th + pi, 2*pi) - pi;
    if abs(w) < 1e-12; w = pi; end
    rv_out = rv / th * w;
end

function s = json_vec(v)
    v = v(:)';
    s = ['[', strjoin(arrayfun(@(x) num2str(x,'%.8g'), v, 'UniformOutput', false), ', '), ']'];
end

function s = json_mat(M)
    rows = cell(1, size(M,1));
    for i = 1:size(M,1); rows{i} = json_vec(M(i,:)); end
    s = ['[', strjoin(rows, ', '), ']'];
end
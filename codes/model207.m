%% ================================================================
%  SHWS RI measurement:
%  Self-consistent iterative estimation of RI center value,
%  RI gradients and rotation angle
%
%  Iterative logic:
%
%  Initial n_center
%       ↓
%  Calculate gx, gy
%       ↓
%  Calculate rotation angle theta
%       ↓
%  Rotate wavefront gradients
%       ↓
%  Calculate nx' and ny'
%       ↓
%  Update n_center
%       ↓
%  Repeat until convergence
%
% ================================================================

clear;
clc;
close all;

rng(2026);


%% ================================================================
% 1. BASIC PARAMETERS
% ================================================================

lambda = 632.8e-6;       % mm

n_air = 1.0000;
n_ref = 1.3325;          % pure water
n_sample = 1.3595;       % true nominal sample RI

d0 = 0.5;                % mm


%% SHWS sampling

Nx = 40;
Ny = 40;

Lx = 6;
Ly = 6;

x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);

[X,Y] = meshgrid(x,y);

dx = x(2)-x(1);
dy = y(2)-y(1);


%% ================================================================
% 2. TRUE RI DISTRIBUTION
% ================================================================

gx_true = 4.0e-4;
gy_true = 1.0e-4;

qxx = 0;
qxy = 0;
qyy = 0;

n_true = n_sample ...
       + gx_true.*X ...
       + gy_true.*Y ...
       + qxx.*X.^2 ...
       + qxy.*X.*Y ...
       + qyy.*Y.^2;


%% True direction

theta_true = atan2(gy_true,gx_true)*180/pi;

fprintf('\n');
fprintf('================================================\n');
fprintf(' TRUE PARAMETERS\n');
fprintf('================================================\n');

fprintf('True center RI       = %.8f\n',n_sample);
fprintf('True gx              = %.4e RIU/mm\n',gx_true);
fprintf('True gy              = %.4e RIU/mm\n',gy_true);
fprintf('True angle           = %.6f deg\n',theta_true);


%% ================================================================
% 3. OPTICAL WEDGE
% ================================================================

alpha_x = 1.0e-4;
alpha_y = 2.0e-4;

d_wedge = d0 ...
        + alpha_x.*X ...
        + alpha_y.*Y;


%% ================================================================
% 4. GENERATE THREE WAVEFRONTS
% ================================================================

[W_empty_true,~,~] = ...
    generate_wavefront(d_wedge,n_air,lambda);

[W_ref_true,~,~] = ...
    generate_wavefront(d_wedge,n_ref,lambda);

[W_sample_true,~,~] = ...
    generate_wavefront(d_wedge,n_true,lambda);


%% ================================================================
% 5. SHWS PARAMETERS
% ================================================================

f_ml = 5;                    % mm
pixel_size = 0.005;          % mm
sigma_centroid_pixel = 0.03; % pixel


%% ================================================================
% 6. SIMULATE SHWS
% ================================================================

[sx_empty,sy_empty] = ...
    simulate_SHWS( ...
    W_empty_true,X,Y,dx,dy,...
    f_ml,pixel_size,...
    sigma_centroid_pixel);

[sx_ref,sy_ref] = ...
    simulate_SHWS( ...
    W_ref_true,X,Y,dx,dy,...
    f_ml,pixel_size,...
    sigma_centroid_pixel);

[sx_sample,sy_sample] = ...
    simulate_SHWS( ...
    W_sample_true,X,Y,dx,dy,...
    f_ml,pixel_size,...
    sigma_centroid_pixel);


%% ================================================================
% 7. RECONSTRUCT SHWS WAVEFRONTS
% ================================================================

W_empty = reconstruct_quadratic_wavefront( ...
    X,Y,sx_empty,sy_empty);

W_ref = reconstruct_quadratic_wavefront( ...
    X,Y,sx_ref,sy_ref);

W_sample = reconstruct_quadratic_wavefront( ...
    X,Y,sx_sample,sy_sample);


%% ================================================================
% 8. REMOVE CENTER PISTON
% ================================================================

ix0 = round((Nx+1)/2);
iy0 = round((Ny+1)/2);

W_empty = W_empty - W_empty(iy0,ix0);
W_ref = W_ref - W_ref(iy0,ix0);
W_sample = W_sample - W_sample(iy0,ix0);


%% ================================================================
% 9. WAVEFRONT GRADIENTS
% ================================================================

[dW_empty_x,dW_empty_y] = ...
    gradient(W_empty_true,dx,dy);

[dW_ref_x,dW_ref_y] = ...
    gradient(W_ref_true,dx,dy);

[dW_sample_x,dW_sample_y] = ...
    gradient(W_sample_true,dx,dy);


%% ================================================================
% 10. THREE-STATE DIFFERENCE
% ================================================================

D_x = dW_sample_x - dW_empty_x;
D_y = dW_sample_y - dW_empty_y;

R_x = dW_ref_x - dW_empty_x;
R_y = dW_ref_y - dW_empty_y;


%% ================================================================
% 11. REMOVE EDGE REGION
% ================================================================

edge_margin = 3;

mask_inner = true(size(X));

mask_inner(1:edge_margin,:) = false;
mask_inner(end-edge_margin+1:end,:) = false;

mask_inner(:,1:edge_margin) = false;
mask_inner(:,end-edge_margin+1:end) = false;


%% ================================================================
% 12. ITERATION PARAMETERS
% ================================================================

% Initial predicted center RI
n_center = 1.3600;

% Convergence tolerance
tol_n = 1e-4;

% Maximum iteration number
max_iter = 1000;

% Relaxation factor
%
% beta = 1:
% direct update
%
% beta < 1:
% more stable
%
beta = 0.5;


%% ================================================================
% 13. STORAGE FOR ITERATION
% ================================================================

iteration_history = [];

n_history = nan(max_iter,1);

gx_history = nan(max_iter,1);
gy_history = nan(max_iter,1);

theta_history = nan(max_iter,1);

nx_center_history = nan(max_iter,1);
ny_center_history = nan(max_iter,1);

error_history = nan(max_iter,1);


%% ================================================================
% 14. ITERATIVE SELF-CONSISTENT CALCULATION
% ================================================================

fprintf('\n');
fprintf('================================================\n');
fprintf(' ITERATIVE CALCULATION\n');
fprintf('================================================\n');

fprintf('\n');
fprintf('%5s %14s %14s %14s %14s\n',...
    'Iter','n_center','gx','gy','theta');

fprintf('%5s %14s %14s %14s %14s\n',...
    '-----','------------','------------',...
    '------------','------------');


for iter = 1:max_iter


    %% ------------------------------------------------------------
    % STEP 1
    % Use current center RI to remove wedge contribution
    % ------------------------------------------------------------

    wedge_ratio = ...
        (n_center-1)/(n_ref-1)


    % RI gradient
    Gx = ...
        (D_x - wedge_ratio.*R_x)./d0;

    Gy = ...
        (D_y - wedge_ratio.*R_y)./d0;


    %% ------------------------------------------------------------
    % STEP 2
    % Calculate gx and gy
    % ------------------------------------------------------------

    Gx_valid = Gx(mask_inner);
    Gy_valid = Gy(mask_inner);

    Gx_valid = Gx_valid(isfinite(Gx_valid));
    Gy_valid = Gy_valid(isfinite(Gy_valid));


    gx_est = median(Gx_valid);
    gy_est = median(Gy_valid);


    %% ------------------------------------------------------------
    % STEP 3
    % Calculate rotation angle
    %
    % New y' axis follows RI gradient
    % ------------------------------------------------------------

    theta = atan2(gy_est,gx_est);

    theta_deg = theta*180/pi;


    %% ------------------------------------------------------------
    % STEP 4
    % Rotate differential wavefront gradients
    % ------------------------------------------------------------

    c = cos(theta);
    s = sin(theta);


    % New x' direction
    D_x_rot = ...
        -s.*D_x + c.*D_y;

    R_x_rot = ...
        -s.*R_x + c.*R_y;


    % New y' direction
    D_y_rot = ...
        c.*D_x + s.*D_y;

    R_y_rot = ...
        c.*R_x + s.*R_y;


    %% ------------------------------------------------------------
    % STEP 5
    % Calculate directional RI
    % ------------------------------------------------------------
    %
    % nx' = 1 + (n_ref-1)*Dx'/Rx'
    %
    % ny' = 1 + (n_ref-1)*Dy'/Ry'
    %
    % A threshold is necessary because the denominator may be
    % close to zero.
    % ------------------------------------------------------------

    ref_threshold_x = ...
        0.20*max(abs(R_x_rot(mask_inner)));

    ref_threshold_y = ...
        0.20*max(abs(R_y_rot(mask_inner)));


    mask_x = ...
        mask_inner & ...
        isfinite(R_x_rot) & ...
        abs(R_x_rot)>ref_threshold_x;


    mask_y = ...
        mask_inner & ...
        isfinite(R_y_rot) & ...
        abs(R_y_rot)>ref_threshold_y;


    nx_rot = nan(size(X));

    ny_rot = nan(size(X));


    nx_rot(mask_x) = ...
        1 + ...
        (n_ref-1).* ...
        D_x_rot(mask_x)./R_x_rot(mask_x);


    ny_rot(mask_y) = ...
        1 + ...
        (n_ref-1).* ...
        D_y_rot(mask_y)./R_y_rot(mask_y);

    %% ------------------------------------------------------------
    % STEP 6
    % Calculate central RI values from nx' and ny'
    %
    % Priority:
    % 1. Use the exact central point
    % 2. If the central point is invalid, use the
    %    median value within the central ROI
    % ------------------------------------------------------------

    roi_radius = 2;     % SHWS sampling points

    roi_mask = ...
        (X-X(iy0,ix0)).^2 + ...
        (Y-Y(iy0,ix0)).^2 ...
        <= (roi_radius*max(dx,dy))^2;


    % ------------------------------------------------------------
    % Exact center point
    % ------------------------------------------------------------

    nx_center_point = nx_rot(iy0,ix0);
    ny_center_point = ny_rot(iy0,ix0);


    % ------------------------------------------------------------
    % If exact center point is invalid,
    % use central ROI median
    % ------------------------------------------------------------

    if ~isfinite(nx_center_point)

        nx_center_values = nx_rot(roi_mask & mask_x);

        nx_center_values = ...
            nx_center_values(isfinite(nx_center_values));

        if ~isempty(nx_center_values)
            nx_center = median(nx_center_values);
        else
            nx_center = NaN;
        end

    else

        nx_center = nx_center_point;

    end


    if ~isfinite(ny_center_point)

        ny_center_values = ny_rot(roi_mask & mask_y);

        ny_center_values = ...
            ny_center_values(isfinite(ny_center_values));

        if ~isempty(ny_center_values)
            ny_center = median(ny_center_values);
        else
            ny_center = NaN;
        end

    else

        ny_center = ny_center_point;

    end


    %% ------------------------------------------------------------
    % STEP 7
    % Select the directional RI closest to current predicted RI
    % ------------------------------------------------------------

    err_nx = abs(n_center - nx_center);
    err_ny = abs(n_center - ny_center);


    if isfinite(err_nx) && isfinite(err_ny)

        if err_nx <= err_ny

            n_target = nx_center;
            selected_direction = 'nx';

        else

            n_target = ny_center;
            selected_direction = 'ny';

        end

    elseif isfinite(err_nx)

        n_target = nx_center;
        selected_direction = 'nx';

    elseif isfinite(err_ny)

        n_target = ny_center;
        selected_direction = 'ny';

    else

        warning('No valid central RI estimate at iteration %d.',iter);

        break;

    end


    %% ------------------------------------------------------------
    % STEP 8
    % Update predicted center RI
    % ------------------------------------------------------------

    n_new_raw = n_target;

    n_new = ...
        (1-beta).*n_center ...
        + beta.*n_new_raw;


    %% ------------------------------------------------------------
    % STEP 9
    % Convergence error
    %
    % Convergence is determined by the difference between
    % the current predicted RI and the closer directional
    % central RI.
    % ------------------------------------------------------------

    dn = abs(n_center - n_target);


    %% ------------------------------------------------------------
    % STEP 10
    % Save iteration history
    % ------------------------------------------------------------

    n_history(iter) = n_center;

    gx_history(iter) = gx_est;
    gy_history(iter) = gy_est;

    theta_history(iter) = theta_deg;

    nx_center_history(iter) = nx_center;
    ny_center_history(iter) = ny_center;

    error_history(iter) = dn;


    %% ------------------------------------------------------------
    % STEP 11
    % Print current iteration
    % ------------------------------------------------------------

    fprintf('%5d %14.8f %14.5e %14.5e %14.5f ',...
        iter,...
        n_center,...
        gx_est,...
        gy_est,...
        theta_deg);

    fprintf('nx=%14.8f ny=%14.8f ',...
        nx_center,...
        ny_center);

    fprintf('target=%s dn=%10.3e\n',...
        selected_direction,...
        dn);


    %% ------------------------------------------------------------
    % STEP 12
    % Check convergence
    % ------------------------------------------------------------

    if dn < tol_n

        fprintf('\n');
        fprintf('Converged at iteration %d.\n',iter);

        fprintf('Selected direction = %s\n',selected_direction);

        fprintf('Predicted n_center  = %.10f\n',n_center);

        fprintf('nx'' center         = %.10f\n',nx_center);

        fprintf('ny'' center         = %.10f\n',ny_center);

        fprintf('Convergence error   = %.4e RIU\n',dn);

        % Use the selected central RI as final value
        n_center = n_target;

        break;

    end


    %% ------------------------------------------------------------
    % Update
    % ------------------------------------------------------------

    n_center = n_new;

end
   


%% ================================================================
% 15. FINAL RESULTS
% ================================================================

n_final = n_center;

gx_final = gx_est;
gy_final = gy_est;

theta_final = theta_deg;


fprintf('\n');
fprintf('================================================\n');
fprintf(' FINAL RESULTS\n');
fprintf('================================================\n');

fprintf('True center RI       = %.10f\n',n_sample);

fprintf('Estimated center RI  = %.10f\n',n_final);

fprintf('RI error             = %.4e RIU\n',...
    n_final-n_sample);


fprintf('\n');

fprintf('True gx              = %.6e RIU/mm\n',gx_true);

fprintf('Estimated gx         = %.6e RIU/mm\n',gx_final);


fprintf('\n');

fprintf('True gy              = %.6e RIU/mm\n',gy_true);

fprintf('Estimated gy         = %.6e RIU/mm\n',gy_final);


fprintf('\n');

fprintf('True angle           = %.6f deg\n',theta_true);

fprintf('Final angle          = %.6f deg\n',theta_final);

fprintf('Angle error          = %.6f deg\n',...
    wrapTo180(theta_final-theta_true));


%% ================================================================
% 16. RECOMPUTE FINAL RI MAPS
% ================================================================

wedge_ratio_final = ...
    (n_final-1)/(n_ref-1);


Gx_final = ...
    (D_x - wedge_ratio_final.*R_x)./d0;

Gy_final = ...
    (D_y - wedge_ratio_final.*R_y)./d0;


%% ================================================================
% 17. FINAL ROTATION
% ================================================================

theta_final_rad = theta_final*pi/180;

c = cos(theta_final_rad);
s = sin(theta_final_rad);


D_x_final = ...
    -s.*D_x + c.*D_y;

D_y_final = ...
    c.*D_x + s.*D_y;


R_x_final = ...
    -s.*R_x + c.*R_y;

R_y_final = ...
    c.*R_x + s.*R_y;


%% ================================================================
% 18. FINAL DIRECTIONAL RI
% ================================================================

threshold_x = ...
    0.20*max(abs(R_x_final(mask_inner)));

threshold_y = ...
    0.20*max(abs(R_y_final(mask_inner)));


mask_x_final = ...
    mask_inner & ...
    abs(R_x_final)>threshold_x;


mask_y_final = ...
    mask_inner & ...
    abs(R_y_final)>threshold_y;


nx_final = nan(size(X));
ny_final = nan(size(X));


nx_final(mask_x_final) = ...
    1 + ...
    (n_ref-1).* ...
    D_x_final(mask_x_final)./ ...
    R_x_final(mask_x_final);


ny_final(mask_y_final) = ...
    1 + ...
    (n_ref-1).* ...
    D_y_final(mask_y_final)./ ...
    R_y_final(mask_y_final);


%% ================================================================
% 19. FINAL CENTRAL RI
% ================================================================

roi_radius = 2;

roi_mask = ...
    (X-X(iy0,ix0)).^2 + ...
    (Y-Y(iy0,ix0)).^2 ...
    <= (roi_radius*max(dx,dy))^2;


center_nx = nx_final(roi_mask & mask_x_final);
center_ny = ny_final(roi_mask & mask_y_final);


center_nx = center_nx(isfinite(center_nx));
center_ny = center_ny(isfinite(center_ny));


fprintf('\n');
fprintf('================================================\n');
fprintf(' FINAL DIRECTIONAL RI\n');
fprintf('================================================\n');

fprintf('nx'' center = %.10f\n',...
    median(center_nx));

fprintf('ny'' center = %.10f\n',...
    median(center_ny));


%% ================================================================
% 20. ITERATION HISTORY
% ================================================================

valid_iter = find(isfinite(n_history));

n_history = n_history(valid_iter);

gx_history = gx_history(valid_iter);
gy_history = gy_history(valid_iter);

theta_history = theta_history(valid_iter);

nx_center_history = ...
    nx_center_history(valid_iter);

ny_center_history = ...
    ny_center_history(valid_iter);

error_history = ...
    error_history(valid_iter);


%% ================================================================
% 21. FIGURE: n CENTER CONVERGENCE
% ================================================================

figure('Name','Center RI convergence');

plot(1:length(n_history),...
    n_history,...
    'o-','LineWidth',1.5);

hold on;

yline(n_sample,'--');

hold off;

grid on;

xlabel('Iteration');
ylabel('Center RI');

title('Self-consistent convergence of center RI');

legend( ...
    'Estimated n_{center}',...
    'True n_{center}',...
    'Location','best');


%% ================================================================
% 22. FIGURE: gx gy CONVERGENCE
% ================================================================

figure('Name','Gradient convergence');

plot(1:length(gx_history),...
    gx_history,...
    'o-','LineWidth',1.5);

hold on;

plot(1:length(gy_history),...
    gy_history,...
    's-','LineWidth',1.5);

yline(gx_true,'--');
yline(gy_true,'--');

hold off;

grid on;

xlabel('Iteration');
ylabel('RI gradient (RIU/mm)');

legend( ...
    'g_x',...
    'g_y',...
    'True g_x',...
    'True g_y',...
    'Location','best');

title('Convergence of RI gradients');


%% ================================================================
% 23. FIGURE: ANGLE CONVERGENCE
% ================================================================

figure('Name','Rotation angle convergence');

plot(1:length(theta_history),...
    theta_history,...
    'o-','LineWidth',1.5);

hold on;

yline(theta_true,'--');

hold off;

grid on;

xlabel('Iteration');

ylabel('\theta (deg)');

title('Convergence of RI-gradient direction');

legend( ...
    'Estimated angle',...
    'True angle',...
    'Location','best');


%% ================================================================
% 24. FIGURE: RI MAP nx'
% ================================================================

figure('Name','RI x prime');

imagesc(x,y,nx_final);

axis image;
axis xy;

colorbar;

xlabel('x (mm)');
ylabel('y (mm)');

title(sprintf( ...
    'Reconstructed n_{x''}, \\theta = %.3f deg',...
    theta_final));


%% ================================================================
% 25. FIGURE: RI MAP ny'
% ================================================================

figure('Name','RI y prime');

imagesc(x,y,ny_final);

axis image;
axis xy;

colorbar;

xlabel('x (mm)');
ylabel('y (mm)');

title(sprintf( ...
    'Reconstructed n_{y''}, \\theta = %.3f deg',...
    theta_final));


%% ================================================================
% 26. FIGURE: TRUE RI
% ================================================================

figure('Name','True RI');

imagesc(x,y,n_true);

axis image;
axis xy;

colorbar;

xlabel('x (mm)');
ylabel('y (mm)');

title('True refractive-index distribution');


%% ================================================================
% 27. FIGURE: RI GRADIENT
% ================================================================

figure('Name','RI gradient');

quiver(X,Y,...
    Gx_final,...
    Gy_final,...
    'AutoScaleFactor',1);

axis image;

xlabel('x (mm)');
ylabel('y (mm)');

title('Reconstructed RI-gradient field');

hold on;

quiver(0,0,...
    gx_final,...
    gy_final,...
    0,...
    'LineWidth',2);

hold off;


%% ================================================================
% 28. SAVE
% ================================================================

save( ...
    'SHWS_RI_iterative_rotation.mat',...
    'X','Y','x','y',...
    'n_true',...
    'W_empty','W_ref','W_sample',...
    'D_x','D_y',...
    'R_x','R_y',...
    'Gx_final','Gy_final',...
    'nx_final','ny_final',...
    'n_history',...
    'gx_history','gy_history',...
    'theta_history',...
    'theta_true',...
    'theta_final',...
    'n_final',...
    'error_history');


fprintf('\n');
fprintf('================================================\n');
fprintf(' PROGRAM FINISHED\n');
fprintf(' Results saved as:\n');
fprintf(' SHWS_RI_iterative_rotation.mat\n');
fprintf('================================================\n');


%% ================================================================
% FUNCTION 1
% Generate wavefront
% ================================================================

function [W,phi,OPD] = ...
    generate_wavefront(d,n,lambda)

    OPD = (n-1).*d;

    [Ny,Nx] = size(OPD);

    ix0 = round((Nx+1)/2);
    iy0 = round((Ny+1)/2);

    W = OPD - OPD(iy0,ix0);

    phi = 2*pi/lambda.*W;

end


%% ================================================================
% FUNCTION 2
% SHWS simulation
% ================================================================

function [slope_x,slope_y] = ...
    simulate_SHWS( ...
    W,X,Y,dx,dy,...
    f_ml,pixel_size,...
    sigma_centroid_pixel)

    [dWdx,dWdy] = gradient(W,dx,dy);

    spot_x = f_ml.*dWdx;
    spot_y = f_ml.*dWdy;

    spot_x_pixel = spot_x./pixel_size;
    spot_y_pixel = spot_y./pixel_size;

    spot_x_pixel = ...
        spot_x_pixel + ...
        sigma_centroid_pixel.*randn(size(X));

    spot_y_pixel = ...
        spot_y_pixel + ...
        sigma_centroid_pixel.*randn(size(X));

    spot_x_measured = ...
        spot_x_pixel.*pixel_size;

    spot_y_measured = ...
        spot_y_pixel.*pixel_size;

    slope_x = spot_x_measured./f_ml;
    slope_y = spot_y_measured./f_ml;

end


%% ================================================================
% FUNCTION 3
% Quadratic wavefront reconstruction
% ================================================================

function W_fit = ...
    reconstruct_quadratic_wavefront( ...
    X,Y,slope_x,slope_y)

    N = numel(X);

    Ax = [ ...
        zeros(N,1),...
        ones(N,1),...
        zeros(N,1),...
        2*X(:),...
        Y(:),...
        zeros(N,1)];

    Ay = [ ...
        zeros(N,1),...
        zeros(N,1),...
        ones(N,1),...
        zeros(N,1),...
        X(:),...
        2*Y(:)];

    A = [Ax;Ay];

    B = [slope_x(:);slope_y(:)];

    coef = A\B;

    a0 = coef(1);
    a1 = coef(2);
    a2 = coef(3);
    a3 = coef(4);
    a4 = coef(5);
    a5 = coef(6);

    W_fit = ...
        a0 ...
        + a1.*X ...
        + a2.*Y ...
        + a3.*X.^2 ...
        + a4.*X.*Y ...
        + a5.*Y.^2;

    ix0 = round((size(X,2)+1)/2);
    iy0 = round((size(X,1)+1)/2);

    W_fit = W_fit - W_fit(iy0,ix0);

end
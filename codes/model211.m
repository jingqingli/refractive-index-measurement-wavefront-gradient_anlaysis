%% ============================================================
%  SHWS-based RI measurement simulation
%
%  Three-state differential measurement:
%       1. Empty cuvette
%       2. Pure water
%       3. Sample liquid
%
%  Main functions:
%       True wavefront generation
%       SHWS slope measurement
%       Centroid noise
%       Quadratic wavefront reconstruction
%       Rotation from 0 to 90 deg
%       nx / ny RI reconstruction
%       Automatic optimal-angle selection
%
%  ============================================================

clear;
clc;
close all;
rng(2026);
%% ============================================================
% 1. Optical parameters
% =============================================================
lambda = 632.8e-6;       % wavelength [mm]
n_air   = 1.0000;
n_water = 1.3317;

% Liquid optical path length
d0 = 0.5;                % [mm]
%% ============================================================
% 2. Spatial sampling
% =============================================================

dx = 0.15;dy =dx;

% for k1 = 1:25
%     for k2=1:25
       % tic
    k1=1;
    k2=5;
    k3=5;
    % filename = sprintf('1002\10\1\01\\%d.csv',k2);
    filename = fullfile('1002','10','1','01',sprintf('%d.csv',k2));
    II10 = double(readmatrix(filename));
    filename = fullfile('1002','10','1','10',sprintf('%d.csv',k3));
    II30= double(readmatrix(filename));
    filename = fullfile('1002','10','1','空1',sprintf('%d.csv',k1));
    II00= double(readmatrix(filename));

II3=II30(82:120,2:40);
II1=II10(82:120,2:40);
II0=II00(82:120,2:40);

[M,N]=size(II3(:,:,1));
%dnx=1500; dny=dnx;
dnx=38; dny=dnx;
nx1=round(M/2)-dnx/2+1+0; nx2=nx1+dnx-1;
ny1=round(N/2)-dny/2+1-0; ny2=ny1+dny-1;

Nx=dnx/2;Ny=Nx;
x = (-Nx:Nx-1)*dx;
y = (-Ny:Ny-1)*dx;
[X,Y] = meshgrid(x,y);

I1 = II1(nx1:nx2,ny1:ny2);
I3 = II3(nx1:nx2,ny1:ny2);
I0 = II0(nx1:nx2,ny1:ny2);
% 
% [p23, phi3] = extract_phase_fit1(I3);
% [p21, phi1] = extract_phase_fit0(I1);
% [p20, phi0] = extract_phase_fit0(I0);
% [p201, phi01] = extract_phase_fit1(I0);

W_sample=I3;
W_water=I1;
W_empty=I0;

% W_empty=W_empty-W_empty;
% W_water=W_water-W_empty;
% W_sample=W_sample-W_empty;

%%


[p_empty,W_empty_fit] = phase_fit1(W_empty,dx);
[p_water,W_water_fit] = phase_fit1(W_water,dx);
[p_sample,W_sample_fit] = phase_fit(W_sample,dx);

[dphix0,dphiy0]=gradient(W_empty_fit,dx);
% [dphix01,dphiy01]=gradient(phi01,dx);
[dphix1,dphiy1]=gradient(W_water_fit,dx);
[dphix3,dphiy3]=gradient(W_sample_fit,dx);
nref=n_water;
nxxx=(dphix3-dphix0)./(dphix1-dphix0)*(nref-1)+1;
nyyy=(dphiy3-dphiy0)./(dphiy1-dphiy0)*(nref-1)+1;
imagesc(nxxx)
mean(nxxx(:))
mean(nyyy(:))
mean(nxxx(:)+nyyy(:))/2

% p_sample.p20/p_sample.p02
ab=p_sample.p20/p_sample.p02/p_water.p10*p_water.p01;

theta_rad = -ab;
nref = n_water;
theta = atan(theta_rad);


%% --------------------------------------------------------
% Rotate wavefront coefficients
% --------------------------------------------------------
coef_empty = ...
    rotate_phase_coef1(theta,p_empty);
coef_water = ...
    rotate_phase_coef1(theta,p_water);
coef_sample = ...
    rotate_phase_coef(theta,p_sample);

%% --------------------------------------------------------
% Generate rotated wavefronts
% --------------------------------------------------------
W0_rot = ...
    polynomial_wavefront(coef_empty,X,Y);
W1_rot = ...
    polynomial_wavefront(coef_water,X,Y);
W3_rot = ...
    polynomial_wavefront(coef_sample,X,Y);

%% --------------------------------------------------------
% Calculate wavefront gradients
% --------------------------------------------------------
[dW0x,dW0y] = gradient(W0_rot,dx,dy);
[dW1x,dW1y] = gradient(W1_rot,dx,dy);
[dW3x,dW3y] = gradient(W3_rot,dx,dy);
%% --------------------------------------------------------
% Differential wavefront slopes
% --------------------------------------------------------
D_x = dW3x - dW0x;
D_y = dW3y - dW0y;
R_x = dW1x - dW0x;
R_y = dW1y - dW0y;

%% --------------------------------------------------------
% Avoid unstable division
% --------------------------------------------------------
threshold_x = ...
    0.05 * max(abs(R_x(:)));
threshold_y = ...
    0.05 * max(abs(R_y(:)));
mask_x = abs(R_x) > threshold_x;
mask_y = abs(R_y) > threshold_y;
%% --------------------------------------------------------
% RI reconstruction
% --------------------------------------------------------
nx_map = nan(size(X));
ny_map = nan(size(X));
nx_map(mask_x) = ...
    1 + (nref-1).* ...
    D_x(mask_x)./R_x(mask_x);
ny_map(mask_y) = ...
    1 + (nref-1).* ...
    D_y(mask_y)./R_y(mask_y);
%% --------------------------------------------------------
% Common valid region
% --------------------------------------------------------
mask_common = ...
    mask_x & mask_y & ...
    isfinite(nx_map) & ...
    isfinite(ny_map);
nx_all = nx_map;
ny_all = ny_map;

%% ============================================================
% 25. Plot optimal RI maps
% =============================================================
figure('Name','RI maps at optimal angle');
subplot(1,3,1);
imagesc(x,y,nx_map);
axis image;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar;
subplot(1,3,2);
imagesc(x,y,ny_map);
axis image;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar;
subplot(1,3,3);
% imagesc(x,y,n_true);
% axis image;
% xlabel('x [mm]');
% ylabel('y [mm]');
% title('True RI Distribution');
% colorbar;
%%
fprintf('\n');
fprintf('====================================================\n');
fprintf('SIMULATION SUMMARY\n');
fprintf('====================================================\n');
fprintf('Grid size                 : %d x %d\n',Nx,Ny);
fprintf('Cuvette thickness         : %.3f mm\n',d0);
fprintf('rotation angle            : %.2f rad\n',theta);
fprintf('p20/p02                   : %.16f\n',p_sample.p20/p_sample.p02);
fprintf('estimated gx/gy           : %.16f\n',ab);
fprintf('estimated nx0             : %.16f\n',mean(nx_map(:)));
fprintf('estimated ny0             : %.16f\n',mean(ny_map(:)));

%% ============================================================
% FUNCTIONS
% =============================================================


%% ------------------------------------------------------------
% Generate true SHWS wavefront
% ------------------------------------------------------------

function [W,phi,OPD] = ...
    generate_SHWS_wavefront(d,n,lambda)

    % Optical path difference
    %OPD = (n-1).*d;
    OPD = (n).*d;
    % Grid center
    [Ny,Nx] = size(OPD);

    ix0 = round(Nx/2);
    iy0 = round(Ny/2);

    % Remove piston
    W = OPD - OPD(iy0,ix0);

    % Optical phase
    phi = ...
        2*pi/lambda .* W;

end


%% ------------------------------------------------------------
% SHWS measurement simulation
% ------------------------------------------------------------

function [slope_x,slope_y,...
          spot_x,spot_y] = ...
    simulate_SHWS(W,x,y,...
                  f_ml,pixel_size,...
                  sigma_centroid_pixel)

    % Grid spacing
    dx = x(2)-x(1);
    dy = y(2)-y(1);

    % True wavefront slopes
    [dWdx,dWdy] = ...
        gradient(W,dx,dy);

    % --------------------------------------------------------
    % Ideal spot displacement
    % --------------------------------------------------------

    spot_x = f_ml .* dWdx;
    spot_y = f_ml .* dWdy;

    % Convert to CCD pixel
    spot_x_pixel = ...
        spot_x ./ pixel_size;

    spot_y_pixel = ...
        spot_y ./ pixel_size;

    % --------------------------------------------------------
    % Add centroid noise
    % --------------------------------------------------------

    spot_x_pixel = ...
        spot_x_pixel + ...
        sigma_centroid_pixel .* ...
        randn(size(spot_x_pixel));

    spot_y_pixel = ...
        spot_y_pixel + ...
        sigma_centroid_pixel .* ...
        randn(size(spot_y_pixel));

    % --------------------------------------------------------
    % Convert measured spot displacement
    % back to wavefront slopes
    % --------------------------------------------------------

    spot_x_meas = ...
        spot_x_pixel .* pixel_size;

    spot_y_meas = ...
        spot_y_pixel .* pixel_size;

    slope_x = ...
        spot_x_meas ./ f_ml;

    slope_y = ...
        spot_y_meas ./ f_ml;

end


%% ------------------------------------------------------------
% Quadratic wavefront reconstruction
% ------------------------------------------------------------

function W_fit = ...
    reconstruct_quadratic_wavefront( ...
    slope_x,slope_y,X,Y)

    % --------------------------------------------------------
    % W = a0 + a1*x + a2*y
    %       + a3*x^2 + a4*x*y + a5*y^2
    % --------------------------------------------------------

    N = numel(X);

    Ax = [ ...
        zeros(N,1), ...
        ones(N,1), ...
        zeros(N,1), ...
        2*X(:), ...
        Y(:), ...
        zeros(N,1)];

    Ay = [ ...
        zeros(N,1), ...
        zeros(N,1), ...
        ones(N,1), ...
        zeros(N,1), ...
        X(:), ...
        2*Y(:)];

    A = [Ax;Ay];

    B = [slope_x(:);slope_y(:)];

    coef = A\B;

    W_fit = ...
        coef(1) ...
        + coef(2).*X ...
        + coef(3).*Y ...
        + coef(4).*X.^2 ...
        + coef(5).*X.*Y ...
        + coef(6).*Y.^2;

end


%% ------------------------------------------------------------
% Construct polynomial wavefront
% ------------------------------------------------------------

function W = ...
    polynomial_wavefront(coef,X,Y)

    W = ...
        coef(1).*X.^2 ...
        + coef(2).*X.*Y ...
        + coef(3).*Y.^2 ...
        + coef(4).*X ...
        + coef(5).*Y ...
        + coef(6);

end
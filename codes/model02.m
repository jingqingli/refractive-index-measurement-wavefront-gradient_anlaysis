clear;clc
lambda = 632.8e-9;
k0 = 2*pi/lambda;
n_w = 1.3325;
% dx = 2.4e-6;
dx = 150e-6;
d0=0.5e-3;
% d0=0.2e-3;
% % % %% ---------- 1 读取干涉图 ----------
% II3 = double(imread('251.bmp'));
% II1 = double(imread('01.bmp'));
% II0 = double(imread('k1.bmp'));
% % % %% ---------- 1 读取干涉图 ----------
%II30 = double(readmatrix('10\4.csv'));
k1=11;
filename = sprintf('10\\%d.csv',k1);
II301 = double(readmatrix(filename));
% k2=2;
% filename = sprintf('10\\%d.csv',k2);
% II302 = double(readmatrix(filename));
II30=II301;

k1=16;
filename = sprintf('0\\%d.csv',k1);
II10=double(readmatrix(filename));
II00= double(readmatrix('空\1.csv'));

II3=II30(82:120,2:40);
II1=II10(82:120,2:40);
II0=II00(82:120,2:40);

figure;imagesc(II3);colorbar
% % 
% II3 = double(imread('B_137_x0000002_y0000005.bmp'));
% II1 = double(imread('B_1333.bmp'));
% II0 = double(imread('B_1.bmp'));
% II3 = double(imread('B_137_0.5.bmp'));
% II1 = double(imread('B_0.5_1333.bmp'));
% II0 = double(imread('B_0.5_1.bmp'));
% II3 = double(imread('B_137_0.5_x-2_y-5_10_6.bmp'));
% II1 = double(imread('B_0.5_1333.bmp'));
% II0 = double(imread('B_0.5_1.bmp'));
% 
% II3 = double(imread('20_14.bmp'));
% % II3 = double(imread('20_14.bmp'));
% II1 = double(imread('0_9.bmp'));
% II0 = double(imread('k_5.bmp'));

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

phi3=I3;
phi1=I1;
phi0=I0;
% [p23, phi3] = phase_fit(I3);
% [p21, phi1] = phase_fit(I1);
% [p20, phi0] = phase_fit(I0);

[dphix0,dphiy0]=gradient(phi0,dx);
% [dphix01,dphiy01]=gradient(phi01,dx);
[dphix1,dphiy1]=gradient(phi1,dx);
[dphix3,dphiy3]=gradient(phi3,dx);

dph21x=mean(dphix1(1))-mean(dphix0(1));
dph21y=mean(dphiy1(1))-mean(dphiy0(1));
nref=n_w;
dpx=(dphix3-dphix0)./(dphix1-dphix0);
dpy=(dphiy3-dphiy0)./(dphiy1-dphiy0);
[px0, dpxx] = phase_fit1(dpx,dx);
[py0, dpyy] = phase_fit1(dpy,dx);
% axx=(nref-1)*px0.p10;
% byy=(nref-1)*px0.p01;
% n0=(px0.p00-axx*d0/(mean(dphix1(:))-mean(dphix0(:))))*(nref-1)+1;

nxxx=(dphix3-dphix0)./(dphix1-dphix0)*(nref-1)+1;
nyyy=(dphiy3-dphiy0)./(dphiy1-dphiy0)*(nref-1)+1;
[n2x0, nxx0] = phase_fit(nxxx,dx);
[n2y0, nyy0] = phase_fit(nyyy,dx);
format long
nax=n2x0.p00
nay=n2y0.p00
naxy=(n2x0.p00+n2y0.p00)/2
% coef_new1(1)=n2x0.p20;coef_new1(2)=n2x0.p11;coef_new1(3)=n2x0.p02;
% coef_new1(4)=n2x0.p10;coef_new1(5)=n2x0.p01;coef_new1(6)=n2x0.p00;
% coef_new0(1)=n2y0.p20;coef_new0(2)=n2y0.p11;coef_new0(3)=n2y0.p02;
% coef_new0(4)=n2y0.p10;coef_new0(5)=n2y0.p01;coef_new0(6)=n2y0.p00;
% 
% % nxx=coef_new1(1)*X.^2 + coef_new1(2)*X.*Y + coef_new1(3)*Y.^2 + coef_new1(4)*X + coef_new1(5)*Y + coef_new1(6);
% % nyy=coef_new0(1)*X.^2 + coef_new0(2)*X.*Y + coef_new0(3)*Y.^2 + coef_new0(4)*X + coef_new0(5)*Y + coef_new0(6);
% nxx= coef_new1(4)*X + coef_new1(5)*Y + coef_new1(6);
% nyy= coef_new0(4)*X + coef_new0(5)*Y + coef_new0(6);
% nxx= coef_new1(5)*Y + coef_new1(6);
% nyy= coef_new0(5)*Y + coef_new0(6);

% figure;subplot(1,3,1);imagesc(nxxx);axis image;colorbar;title('nxx')
% subplot(1,3,2);imagesc(nyyy);axis image;colorbar;title('nyy')
% subplot(1,3,3);imagesc(n3);axis image;colorbar;title('GRIN')

[p23, phi3] = phase_fit(I3,dx);
% residual = phi3 - I3;
% RMS = sqrt(mean(residual(:).^2))
% [p23, phi3] = phase_fit1(I3,dx);
% residual = phi3 - I3;
% RMS = sqrt(mean(residual(:).^2))
% [p23, phi3] = phase_fit3(I3,dx);
% residual = phi3 - I3;
% RMS = sqrt(mean(residual(:).^2))

[p21, phi1] = phase_fit(I1,dx);
[p20, phi0] = phase_fit(I0,dx);

theta=(30)/180*pi;
coef_new3 = rotate_phase_coef(theta, p23);
coef_new1 = rotate_phase_coef1(theta, p21);
coef_new0 = rotate_phase_coef1(theta, p20);

Nx=length(I3);Ny=Nx;
x = (-round(Nx/2-1):round(Nx/2)-1)*dx;
y = (-round(Ny/2-1):round(Ny/2)-1)*dx;

[X,Y] = meshgrid(x,y);

phi01=coef_new1(1)*X.^2 + coef_new1(2)*X.*Y + coef_new1(3)*Y.^2 + coef_new1(4)*X + coef_new1(5)*Y + coef_new1(6);
phi00=coef_new0(1)*X.^2 + coef_new0(2)*X.*Y + coef_new0(3)*Y.^2 + coef_new0(4)*X + coef_new0(5)*Y + coef_new0(6);
phi03=coef_new3(1)*X.^2 + coef_new3(2)*X.*Y + coef_new3(3)*Y.^2 + coef_new3(4)*X + coef_new3(5)*Y + coef_new3(6);

% phi3=unwrap_phase3;
% phi1=unwrap_phase1;
% phi0=unwrap_phase0;

[dphix0,dphiy0]=gradient(phi0);
[dphix1,dphiy1]=gradient(phi1);
[dphix3,dphiy3]=gradient(phi3);

nref=n_w;
nxxx=(dphix3-dphix0)./(dphix1-dphix0)*(nref-1)+1;
nyyy=(dphiy3-dphiy0)./(dphiy1-dphiy0)*(nref-1)+1;
[n2x0, nxx0] = phase_fit1(nxxx,dx);
[n2y0, nyy0] = phase_fit1(nyyy,dx);
% figure;subplot(1,2,1);imagesc(nxxx);colorbar;
% subplot(1,2,2);imagesc(nyyy);colorbar

[dphix0,dphiy0]=gradient(phi00);
[dphix1,dphiy1]=gradient(phi01);
[dphix3,dphiy3]=gradient(phi03);


nxx=(dphix3-dphix0)./(dphix1-dphix0)*(nref-1)+1;
nyy=(dphiy3-dphiy0)./(dphiy1-dphiy0)*(nref-1)+1;
[n2x, nxx1] = phase_fit1(nxx,dx);
[n2y, nyy1] = phase_fit1(nyy,dx);
nbx=n2x.p00
nby=n2y.p00
% double(n2x.p00)-n2y.p00
% figure;subplot(1,2,1);imagesc(nxx);colorbar;
% subplot(1,2,2);imagesc(nyy);colorbar
% data4=nxx;
% [n2x, nxx1] = phase_fit(nxx);

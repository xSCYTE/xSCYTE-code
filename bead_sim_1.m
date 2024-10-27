%%% 3D split-step Fourie beam propagation method (BPM) for simulating
%%% tomographic phase microscope (TPM) by Choi et al.
%%% 1. Simulation of the TPM measurements with BPM
%%% 2. 3D reconstruction using iradon, non-negative constrain, ISTA-TV, and
%%%     FISTA-TV algorithms.
%%% U. S. Kamilov, BIG, EPFL, 2014.

%% Initialize
%%% Set up the parameters and define the computational grid.

clear; close all; clc;

lambda0 = 532e-9; % free space wavelength (m)
magnification = 100*200/180; % Mnominal*Fused/Fnominal
camera_pixel_size = 6.5e-6*2; % camera pixel size (m)

n0 = 1.56; % background refractive index (immersion oil)
k0 = 2*pi/lambda0; % free space wavenumber (1/m)
k = n0*k0; % medium wavenumber (1/m)

%%% depends on cropping
aa=1;
Nx = aa*128; % number of pixels along x
Ny = aa*128; % number of pixels along y
Nz = aa*128; % number of pixels along z

dx = 1*camera_pixel_size/magnification; % discretization step along x
dy = 1*camera_pixel_size/magnification; % discretization step along y
dz = 1*camera_pixel_size/magnification; % discretization step along z

Lx = Nx*dx; % length of computational window along x
Ly = Ny*dy; % length of computational window along y
Lz = Nz*dz; % length of computational window along z

x = dx*(-Nx/2+1:Nx/2)'; % computational grid along x (horiz)
y = dy*(-Ny/2+1:Ny/2)'; % computational grid along y (vert)
z = dz*(1:Nz)'; % computational grid along z

[XX, YY] = meshgrid(x, y); %2D meshgrid

dkx = 2*pi/Lx; % frequency discretization step along x
dky = 2*pi/Ly; % frequency discretization step along y
dkz = 2*pi/Lz; % frequency discretization step along z

kx = dkx*[0:Nx/2-1, -Nx/2:-1]'; % frequency grid along x
ky = dky*[0:Ny/2-1, -Ny/2:-1]'; % frequency grid along y
kz = dkz*[0:Nz/2-1, -Nz/2:-1]'; % frequency grid along z

[Kxx, Kyy] = meshgrid(kx, ky); % 2D frequency meshgrid

K2 = Kxx.^2+Kyy.^2; % frequency norm for all points

dphi = real(K2./(k+sqrt(k^2-K2))); % diffraction phase factor

[~, mid_index_y] = min(abs(y)); % midpoint index along x
[~, mid_index_x] = min(abs(x)); % midpoint index along y
[~, mid_index_z] = min(abs(z-Lz/2)); % midpoint index along z

%%% Create object
forwardObj = PlaneWave3D(Lx, Ly, Lz, Nx, Ny, Nz);

%% Phantom
%%% Generate a phantom consisting of a bead centered in the middle of the
%%% computation window.

[X, Y, Z] = meshgrid(x, y, z); % 3D meshgrid

dn = 0.01 + 0.003i; % refractive index variation
n = n0+dn; % refractive index of the bead
D = 2e-6; % diameter of the bead
xc = 0; % center of the bead along x
yc = 0; % center of the bead along y
zc = 0.5*Lz; % center of the bead along z
zc1 = 0.5*Lz; % center of the bead along z
zc2 = 0.5*Lz; % center of the bead along z
[~, ind_xc] = min(abs(x-xc)); % index of the bead center along x
[~, ind_yc] = min(abs(y-yc)); % index of the bead center along y
[~, ind_zc] = min(abs(z-zc1)); % index of the bead center along z
[~, ind_zc] = min(abs(z-zc2)); % index of the bead center along z

R1 = sqrt((X-xc).^2 + (Y-yc).^2 + (Z-zc1).^2);
R2 = sqrt((X-xc).^2 + (Y-yc).^2 + (Z-zc2).^2);
f = dn*double(R2 < D/2)+dn*double(R1 < D/2); % final 3D phantom object
%%% Plot the projection of the phantom along each of the spatial
%%% coordinates
figure(101);
set(101, 'Name', 'Object');
subplot(1, 3, 1);
imagesc(x, y, squeeze(abs(f(:,:,ceil(Nz/2)))));
axis square;
xlabel('x');
ylabel('y');
title('XY');
subplot(1, 3, 2);
imagesc(z, y, squeeze(abs(f(:,ceil(Nx/2),:))));
axis square;
xlabel('z');
ylabel('y');
title('YZ');
subplot(1, 3, 3);
imagesc(z, x, squeeze(abs(f(ceil(Ny/2),:,:))));
axis square;
xlabel('z');
ylabel('x');
title('XZ');
colormap gray;
clear X Y Z R;

%% Measure
%%% Scan uniquely along the axis x and store the complex field at the end
%%% of the computational window

theta = pi/3; % maximum scanning angle (radians)
Ntheta = 49; % number of scanning angles
%thetas = linspace(-Ltheta/2, Ltheta/2, Ntheta);
psais = linspace(0, 2*pi, Ntheta-1);
% thetas = linspace(-Ltheta, Ltheta, Ntheta);
% dtheta = thetas(2)-thetas(1); % angle step

% angs = thetas(end:-1:1)*180/pi-90; % for iradon in degrees

% [~, mid_index_theta] = min(abs(thetas));

gob = zeros(Ny, Nx, Ntheta); % measurements when object is removed
gem = zeros(Ny, Nx, Ntheta); % measurements when object is present
gin = zeros(Ny, Nx, Ntheta);

for l = 1:Ntheta
    if l == 1
%     angx = thetas(l);
%     angy = 0;
    
        Bx = 0.3*Lx; % beam scale (m)
%    Bx = Bx/cos(angx); % projection of the propagation plane to XY plane
    
        By = 0.3*Ly;
%    By = By/cos(angy);
    
        ain = exp(-((XX./Bx).^2 + (YY./By).^2)); % illumination beam amplitude
%    pin = k*(sin(angx)*XX + sin(angy)*YY);
        pin = 0;
    else
        Bx = 0.3*Lx/(cos(theta)*cos(psais(l-1)));
        By = 0.3*Ly/(cos(theta)*sin(psais(l-1)));
        ain = exp(-((XX./Bx).^2 + (YY./By).^2));
        pin = k*(sin(theta)*cos(psais(l-1))*XX+sin(theta)*sin(psais(l-1))*YY);
    end
        
    
%     igin = ain.*exp(1i*pin); % define a plane wave at the origin
    igin = ain.*exp(1i*pin);
%     igin = ifft2(fft2(igin));
    igin = ifft2(fft2(igin).*exp(-1i*(-Lz/2)*dphi)); % propagate to z = 0
        
    [uob, gtot_ob] = forwardObj.iforward(f, igin);
    [uem, gtot_em] = forwardObj.iforward(zeros(Ny, Nx, Nz), igin);
%     uob = ifft2(fft2(uob).*exp(-1i*(Lz/2)*dphi)); % propagate to z = 0
    figure(222);imagesc(abs(uob));colormap jet; drawnow
    
    uob = exp(1i.*(angle(uob)));
    
    gin(:,:,l) = igin;
    gem(:,:,l) = uem; % store no object measurements
    gob(:,:,l) = awgn(uob,30);% store object measurements
    
    %%% Plot fields
    figure(102);
    set(102, 'Name', sprintf('[%d/%d] Propagation Profile',...
        l, Ntheta));
    
    subplot(2, 2, 1);
    imagesc(z, x, squeeze(abs(gtot_em(mid_index_y,:,:))));
    axis square;
    xlabel('z');
    ylabel('x');
    title('Object Absent');
    
    subplot(2, 2, 2);
    imagesc(z, x, squeeze(abs(gtot_ob(mid_index_y,:,:))));
    hold on;
    plot([zc zc], [x(1) x(end)], 'r--',...
        [z(1) z(end)], [xc xc], 'r--');
    plot([z(mid_index_z) z(mid_index_z)], [x(1) x(end)], 'g:',...
        [z(1) z(end)], [x(mid_index_x) x(mid_index_x)], 'g:');
    viscircles([zc xc], D/2,...
        'DrawBackgroundCircle', false,...
        'EdgeColor', 'r',...
        'LineStyle', '--',...
        'LineWidth', 1);
    hold off;
    axis square;
    xlabel('z');
    ylabel('x');
    title('Object Present');
    colormap gray;
    
%     drawnow;    
%         figure(102);
%     set(102, 'Name', sprintf('[%d/%d] Propagation Profile',...
%         l, Ntheta));
    
    subplot(2, 2, 3);
    imagesc(x, y, squeeze(abs(gtot_em(:,:,mid_index_z))));
    axis square;
    xlabel('z');
    ylabel('x');
    title('Object Absent');
    
    subplot(2, 2, 4);
    imagesc(x, y, squeeze(abs(gtot_ob(:,:,mid_index_z))));
    hold on;
    plot([zc zc], [x(1) x(end)], 'r--',...
        [z(1) z(end)], [xc xc], 'r--');
    plot([z(mid_index_z) z(mid_index_z)], [x(1) x(end)], 'g:',...
        [z(1) z(end)], [x(mid_index_x) x(mid_index_x)], 'g:');
    viscircles([zc xc], D/2,...
        'DrawBackgroundCircle', false,...
        'EdgeColor', 'r',...
        'LineStyle', '--',...
        'LineWidth', 1);
    hold off;
    axis square;
    xlabel('x');
    ylabel('y');
    title('Object Present');
    colormap gray;
    
    drawnow;    
end

clear igin uem uob gtot_em gtot_ob;

% Perform some checks

for l = 1:Ntheta
    igin = gin(:,:,l);
    iginhat = gem(:,:,l);
    iginhat = ifft2(fft2(iginhat).*exp(-1i*(-forwardObj.Lz/2)*forwardObj.dphi)); 
    diff = norm(igin(:)-iginhat(:));
    fprintf('theta = %d: diff = %12.8e\n', l, diff);
end

%% Initialize with Radon

% roi.x = 65-32:192+32;
% roi.y = 65-32:192+32;
% 
% gem_roi = gem(roi.y, roi.x, :);
% gob_roi = gob(roi.y, roi.x, :);
% 
% fhat_roi = radonRoiEst(gem_roi, gob_roi, thetas);
% fhat_roi = fhat_roi/(k0*dz);
% 
% fhat0 = zeros(Ny, Nx, Nz);
% fhat0(roi.y, roi.x, roi.x) = fhat_roi;
% 
% figure(103);
% set(103, 'Name', 'Object');
% subplot(1, 3, 1);
% imagesc(x, y, squeeze(mean(fhat0, 3)));
% axis square;
% xlabel('x');
% ylabel('y');
% title('XY');
% subplot(1, 3, 2);
% imagesc(z, y, squeeze(mean(fhat0, 2)));
% axis square;
% xlabel('z');
% ylabel('y');
% title('YZ');
% subplot(1, 3, 3);
% imagesc(z, x, squeeze(mean(fhat0, 1)));
% axis square;
% xlabel('z');
% ylabel('x');
% title('XZ');

%% Reconstruct with BPM
fhat00= 0.*ones(128,128,128);
stepSize = 1e-3;
bounds = [0, 0.1];
numIter = 50;
lambda = 0.001;
tv_maxiter = 10;

fhat = fistaEst(-gob, -gin, forwardObj, lambda, fhat00, numIter, stepSize,...
    tv_maxiter, bounds);
% fhat = mfistaEst(gob, gin, forwardObj, lambda, 0*f, numIter, stepSize,...
%    tv_maxiter, bounds, f);

figure(446); subplot(221); imagesc(squeeze(real(fhat(:,:,64))),[0 0.03]), colorbar, colormap jet, axis equal, axis off,title('x-y');
        subplot(222); imagesc(squeeze(real(fhat(:,64,:))),[0 0.03]), colorbar, colormap jet, axis equal, axis off,title('x-z');
        subplot(223); imagesc(squeeze(real(fhat(64,:,:))),[0 0.03]), colorbar, colormap jet, axis equal, axis off,title('y-z');

figure; 
imshow3D(real(f));
colormap gray;

figure; 
imshow3D(real(fhat));
colormap gray;
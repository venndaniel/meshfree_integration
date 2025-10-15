% code for the test of Subsection 4.2.2
% compute the integral of a singular function on a paraboloid

% set plots to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')
rng('default');


sur_tol = 1e-13; % tolerance for points on surface

syms x y z

phi = @(x,y,z) z +x.^2 + y.^2 -1; % paraboloid

phixs = diff(phi, x); phiys = diff(phi, y); phizs = diff(phi, z);
phix = matlabFunction(phixs, "Vars",[x y z]);
phiy = matlabFunction(phiys, "Vars",[x y z]);
phiz = matlabFunction(phizs, "Vars",[x y z]);

nphi = sqrt(phixs.^2 + phiys.^2 + phizs.^2);
kappa = matlabFunction(diff(phixs./nphi, x) + diff(phiys./nphi, y) +...
    diff(phizs./nphi, z), "Vars",[x y z]); % compute curvature

%% boundary points and weights
rng('default');
Np = 1000;

% evenly spaced points on circle
theta = linspace(0, 2*pi, 1000+1);
theta = theta(1:1000)';
xp = cos(theta); yp = sin(theta);
zp = 0*theta;

wb = ones(Np, 1)*2*pi/Np; % boundary integral weights

%% integrate
rng('default'); % reset random number generator for replicability
Nt = 2560; % number of surface points
wmax = 11; % total number of Fourier basis functions is (2*wmax + 1)^3
q = 4; %4
T = 10; %10 
sobolev = false;

% singular function and derivatives
syms x y z
s = @(x,y,z) sqrt((x).^2 + y.^2 +(z-1).^2);
sx = diff(s, x); sy = diff(s, y); sz = diff(s, z);
sxx = diff(sx, x); syy = diff(sy, y); szz = diff(sz, z);
sxy = diff(sx, y); syz = diff(sy, z); sxz = diff(sz, x);
sx = matlabFunction(sx); sy = matlabFunction(sy); sz = matlabFunction(sz);
sxx = matlabFunction(sxx); syy = matlabFunction(syy);
szz = matlabFunction(szz);
sxy = matlabFunction(sxy); syz = matlabFunction(syz);
sxz = matlabFunction(sxz);


% place points on the surface
[xs, ys, zs] = place_b_points(phi, Nt, [-4,4], [-4, 4], [0, 2],...
    sur_tol, 25); % 1 is random

p = [0 0 1]; c = 0;

% separate points between positive and negative z
% points should all be positive for this test
pinds = [xs ys zs]*p' > c;
ninds = [xs ys zs]*p' < c;

% positive and negative points
xsp = xs(pinds); ysp = ys(pinds); zsp = zs(pinds);
xsn = xs(ninds); ysn = ys(ninds); zsn = zs(ninds);

% normal vectors to surface on boundary
[nxp, nyp, nzp] = nml(xp, yp, zp, phix, phiy, phiz);
ndot = [nxp nyp nzp]*p';
% co-normal vectors to boundary curves
tx = p(1) - ndot.*nxp; ty = p(2) - ndot.*nyp; tz = p(3) - ndot.*nzp;
tn = vecnorm([tx, ty, tz], 2, 2);
tx = tx./tn; ty = ty./tn; tz = tz./tn;

% normal vectors to surface
[nxp, nyp, nzp] = nml(xsp, ysp, zsp, phix, phiy, phiz);
[nxn, nyn, nzn] = nml(xsn, ysn, zsn, phix, phiy, phiz);

Nb = (2*wmax + 1).^3; % number of Fourier frequencies
freq = -wmax:wmax;
[fx, fy, fz] = meshgrid(freq, freq, freq);

ellx = 4; elly = 4; ellz = 4;% size of bounding box
% omega_x, omega_y, omega_z
wx = reshape(fx, [Nb, 1])*2*pi/ellx;
wy = reshape(fy, [Nb, 1])*2*pi/elly;
wz = reshape(fz, [Nb, 1])*2*pi/ellz;
wa = vecnorm([wx wy wz], 2, 2); % norm of omega vector

di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(wa))); % d_n^{-1/2}

% function values on two surface subdomains
Vpos = di'.*exp(1i*(xsp.*wx' + ysp.*wy' + zsp.*wz'));
%Vneg = di'.*exp(1i*(xsn.*wx' + ysn.*wy' + zsn.*wz'));

% singular function and derivative values
se = s(xsp, ysp, zsp);
sxe = sx(xsp, ysp, zsp);sye = sy(xsp, ysp, zsp);sze = sz(xsp, ysp, zsp);
sxxe=sxx(xsp, ysp, zsp);syye=syy(xsp, ysp, zsp);szze=szz(xsp, ysp, zsp);
sxye=sxy(xsp, ysp, zsp);syze=syz(xsp, ysp, zsp);sxze=sxz(xsp, ysp, zsp);

V = [Vpos se.*Vpos]; % augmented function values

Vlappos = [-wa.^2'.*Vpos (-wa.^2'.*se + (sxxe + syye + szze) +...
    2*(1i*wx'.*sxe + 1i*wy'.*sye + 1i*wz'.*sze)).*Vpos]; % Laplacian
% first normal derivative un-augmented
Vnpos1 = 1i*(nxp.*wx' + nyp.*wy' + nzp.*wz').*Vpos; 
% first normal derivative singular term
Vnpos2 = (1i*(nxp.*wx' + nyp.*wy' + nzp.*wz').*se +...
    (nxp.*sxe + nyp.*sye + nzp.*sze)).*Vpos; 
Vnpos = [Vnpos1 Vnpos2]; % first normal derivative

% second normal derivative unaugmented
Vnnpos1 = -(nxp.^2.*wx.^2' + nyp.^2.*wy.^2' + nzp.^2.*wz.^2' + ...
    2*nxp.*nyp.*(wx.*wy)' + ...
    2*nxp.*nzp.*(wx.*wz)' + ...
    2*nzp.*nyp.*(wz.*wy)').*Vpos;

% second normal derivative singular term
Vnnpos2 = (nxp.^2.*(-wx.^2'.*se + 2i*wx'.*sxe + sxxe)+...
    nyp.^2.*(-wy.^2'.*se + 2i*wy'.*sye + syye) +...
    nzp.^2.*(-wz.^2'.*se + 2i*wz'.*sze + szze) + ...
    2*nxp.*nyp.*(-wx'.*wy'.*se + 1i*wx'.*sye + 1i*wy'.*sxe + sxye) + ...
    2*nxp.*nzp.*(-wx'.*wz'.*se + 1i*wx'.*sze + 1i*wz'.*sxe + sxze) + ...
    2*nzp.*nyp.*(-wy'.*wz'.*se + 1i*wy'.*sze + 1i*wz'.*sye + syze)).*Vpos;
Vnnpos = [Vnnpos1 Vnnpos2]; % second normal derivative

% number of positive/negative points
Npos = length(xsp); Nneg = length(xsn); 
f = -1/4/pi./sqrt(xsp.^2 + ysp.^2 + (zsp-1).^2); % function to integrate

'solving'
% solve w/ and w/o singular terms
ap = lsqminnorm([Vlappos - Vnnpos - kappa(xsp, ysp, zsp).*Vnpos], f);
anaive = lsqminnorm([-wa.^2'.*Vpos - Vnnpos1 -...
    kappa(xsp, ysp, zsp).*Vnpos1], f);

% un-augmented boundary point values
Vb = di'.*exp(1i*(xp.*wx' + yp.*wy' + zp.*wz')); 

% singular function and derivative values on boundary
sp = s(xp, yp, zp); 
sxp = sx(xp, yp, zp); syp = sy(xp, yp, zp); szp = sz(xp, yp, zp);

% co-normal derivative values
Vtnaive = 1i*(tx.*wx' + ty.*wy' + tz.*wz').*Vb;
Vt = [Vtnaive (1i*(tx.*wx' + ty.*wy' + tz.*wz').*sp + ...
    (tx.*sxp + ty.*syp + tz.*szp)).*Vb];

% integral value (by evaluating line integral on circle)
integ = real(wb'*(Vt*(-ap)))
integnaive = real(wb'*(Vtnaive*(-anaive)))

%%
scatter3(xsp, ysp, zsp, 20, real(V*ap), 'filled');
axis('equal')
hold on;
colorbar()
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
title('Surface Poisson Solution on a Paraboloid')


function [nx, ny, nz] = nml(x, y, z, phix, phiy, phiz)
    nx = phix(x,y,z); ny = phiy(x,y,z); nz = phiz(x,y,z) + 0*x;
    nm = vecnorm([nx ny nz],2,2);
    nx = nx./nm; ny = ny./nm; nz = nz./nm;
end


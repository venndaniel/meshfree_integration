% code for the test of Subsection 4.1.3
% integrates on a genus-2 surface by splitting it into two subdomains with
% boundary

% LaTex for plots
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')
rng('default');

sur_tol = 1e-13; % tolerance for placing point on surface

% level set and derivatives
syms x y z
phi = @(x,y,z) 0.25./((x-1).^2 + y.^2) + 0.25./((x+1).^2 + y.^2)+ z.^2-...
    1 + 0.1*x.^2 + 0.25*y.^2;
phix = @(x,y,z)  -0.25*2*(x-1)./((x-1).^2 + y.^2).^2 -...
    0.25*2*(x+1)./((x+1).^2 + y.^2).^2 + 0.1*2*x;
phiy = @(x,y,z)  -0.25*2*y./((x-1).^2 + y.^2).^2 -...
    0.25*2*y./((x+1).^2 + y.^2).^2 + 0.25*2*y;
phiz = @(x,y,z)  2*z;

% need to compute curvature kappa, use symbolic differentiation
phixs = diff(phi, x); phiys = diff(phi, y); phizs = diff(phi, z);
phix = matlabFunction(phixs, "Vars",[x y z]);
phiy = matlabFunction(phiys, "Vars",[x y z]);
phiz = matlabFunction(phizs, "Vars",[x y z]);
nphi = sqrt(phixs.^2 + phiys.^2 + phizs.^2);

% curvature as MATLAB function
kappa = matlabFunction(diff(phixs./nphi, x) + diff(phiys./nphi, y) +...
    diff(phizs./nphi, z), "Vars",[x y z]);

%% Compute boundary weights
rng('default');
Np = 2000;

% boundary plane is specifed by p*[x;y;z] - c = 0
p = [0,0, 1]; % normal vector specifying boundary plane
c = 0; % constant for boundary plane

% place points on the intersection of the plane and surface
[xp, yp, zp] = place_b_points_plane(phi, Np, [-4,4], [-4, 4], [0, 2],...
    sur_tol, 25, [p, c]); % 1 is random, use 25

wb = zeros(Np, 1); % integration weights
for j = 1:Np
    inds = [1:j-1, j+1:Np]; % indices of everything except current point
    % get nearest neighbour index
    [h1, ind1] = min(vecnorm([xp(j) - xp(inds), yp(j) - yp(inds), ...
        zp(j) - zp(inds)], 2, 2));
    ind1 = inds(ind1);

    % displacement to nearest neighbour
    sep = [xp(ind1) - xp(j), yp(ind1) - yp(j), zp(ind1) - zp(j)];
    % displacement to all other points
    seps = [xp(inds) - xp(j), yp(inds) - yp(j), zp(inds) - zp(j)];

    % get next nearest point "on other side" of first nearest point
    [h2, ind2] = min(vecnorm([xp(j) - xp(inds), yp(j) - yp(inds), ...
        zp(j) - zp(inds)], 2, 2) + 100*((seps*sep')>0));
    ind2 = inds(ind2);

    % next separation vector
    sep2 = [xp(ind2) - xp(j), yp(ind2) - yp(j), zp(ind2) - zp(j)];

    % angle between two separation vectors
    alpha = acos(sep2*sep'/(norm(sep)*norm(sep2)));
    d = norm([xp(ind2) - xp(ind1), yp(ind2) - yp(ind1), ...
        zp(ind2) - zp(ind1)]);
    % linear distance
    R = d/2/sin(alpha); 
    % radius of circle connecting the three points
    theta = 2*asin(d/2/R);

    % arc angle between point and first neighbour
    alpha = 2*asin(norm(sep)/2/R); 
    % arc angle between point and second neighbour
    beta = 2*asin(norm(sep2)/2/R);

    % ensure quadratic functions of arc angle are integrated exactly
    A = [1 1 1;
         0 alpha theta;
         0 alpha^2 theta^2];
    b = [theta; theta^2/2; theta^3/3];

    w = (A\b)*R; % integration weights for 3 point arc
    % update weights by half the weights for this arc
    % (other half comes from other two trios each point is a part of)
    wb(j) = wb(j) + w(2)/2; 
    wb(ind1) = wb(ind1) + w(1)/2;
    wb(ind2) = wb(ind2) + w(3)/2;
end
'done line weights'

%% plot boundary curves
scatter3(xp, yp, zp, 30, 'k', 'filled');
title('Genus-2 Surface Subdomain Boundary')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
view(-30,50)
fontsize(18, 'pixels')
axis('equal')

%% integrate on Surface
Nt = 5120; % total number of surface points
wmax = 13; %13
q = 5; %5
T = 10; %10
sobolev = false;
lmult = 2;
rng('default'); % reset the random number generator
[xs, ys, zs] = place_b_points(phi, Nt, [-4,4], [-4, 4], [-2, 2],...
    sur_tol, 25); % 1 is random

% indices for the two subdomains: S_+ and S_-
pinds = [xs ys zs]*p' > c;
ninds = [xs ys zs]*p' < c;

% separate into two point clouds
xsp = xs(pinds); ysp = ys(pinds); zsp = zs(pinds);
xsn = xs(ninds); ysn = ys(ninds); zsn = zs(ninds);

% normal vectors to surface on boundary
[nxp, nyp, nzp] = nml(xp, yp, zp, phix, phiy, phiz);
ndot = [nxp nyp nzp]*p';
% tangent vectors to boundary curves
tx = p(1) - ndot.*nxp; ty = p(2) - ndot.*nyp; tz = p(3) - ndot.*nzp;
tn = vecnorm([tx, ty, tz], 2, 2);
tx = tx./tn; ty = ty./tn; tz = tz./tn;

% normal vectors to surface
[nxp, nyp, nzp] = nml(xsp, ysp, zsp, phix, phiy, phiz);
[nxn, nyn, nzn] = nml(xsn, ysn, zsn, phix, phiy, phiz);

Nb = (2*wmax + 1).^3; % number of Fourier frequencies
freq = -wmax:wmax;
[fx, fy, fz] = meshgrid(freq, freq, freq);

ellx = 10; elly = 6; ellz = 2; % size of bounding box
% omega_x, omega_y, omega_z
wx = reshape(fx, [Nb, 1])*2*pi/ellx;
wy = reshape(fy, [Nb, 1])*2*pi/elly;
wz = reshape(fz, [Nb, 1])*2*pi/ellz;
wa = vecnorm([wx wy wz], 2, 2); % norm of omega vector

di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(wa))); % d_n^{-1/2}

% function values on two surface subdomains
Vpos = di'.*exp(1i*(xsp.*wx' + ysp.*wy' + zsp.*wz'));
Vneg = di'.*exp(1i*(xsn.*wx' + ysn.*wy' + zsn.*wz'));

% positive subdomain
Vlappos = -wa.^2'.*Vpos; % Laplacian
% first normal derivative
Vnpos = 1i*(nxp.*wx' + nyp.*wy' + nzp.*wz').*Vpos; 
% second normal derivative
Vnnpos = -(nxp.^2.*wx.^2' + nyp.^2.*wy.^2' + nzp.^2.*wz.^2' + ...
    2*nxp.*nyp.*(wx.*wy)' + ...
    2*nxp.*nzp.*(wx.*wz)' + ...
    2*nzp.*nyp.*(wz.*wy)').*Vpos;

% negative subdomain
Vlapneg = -wa.^2'.*Vneg; 
Vnneg = 1i*(nxn.*wx' + nyn.*wy' + nzn.*wz').*Vneg; 
Vnnneg = -(nxn.^2.*wx.^2' + nyn.^2.*wy.^2' + nzn.^2.*wz.^2' + ...
    2*nxn.*nyn.*(wx.*wy)' + ...
    2*nxn.*nzn.*(wx.*wz)' + ...
    2*nzn.*nyn.*(wz.*wy)').*Vneg;

% number of points in each subdomain
Npos = length(xsp); Nneg = length(xsn); 

% solve the Poisson problem using lsqminnorm, which uses complete 
% orthogonal decomposition
% Vlappos - Vnnpos - kappa(xsp, ysp, zsp).*Vnpos is Laplace-Beltrami
'solving'
ap = lsqminnorm([Vlappos - Vnnpos - kappa(xsp, ysp, zsp).*Vnpos], ...
    [ones(Npos, 1)]);

'solving'
an = lsqminnorm([Vlapneg - Vnnneg - kappa(xsn, ysn, zsn).*Vnneg], ...
    [ones(Nneg, 1)]);

% boundary values and tangent derivatives
Vb = di'.*exp(1i*(xp.*wx' + yp.*wy' + zp.*wz'));
Vt = 1i*(tx.*wx' + ty.*wy' + tz.*wz').*Vb;

integ = real(wb'*(Vt*(an - ap))) % compute integral


function [nx, ny, nz] = nml(x, y, z, phix, phiy, phiz)
    % normal vector to surface given by a level set with partial
    % derivatives phix, phiy, phiz
    nx = phix(x,y,z); ny = phiy(x,y,z); nz = phiz(x,y,z) + 0*x;
    nm = vecnorm([nx ny nz],2,2);
    nx = nx./nm; ny = ny./nm; nz = nz./nm;
end

function [xf, yf, zf] = place_b_points_plane(lsf, N, xlim, ylim, zlim, ...
    eps, per, plane)
    % places N points on the intersection of plane and a surface given by a
    % level set specified by lsf
    % xlim, ylim, zlim specify the bounding box for the surface
    % eps is the tolerance
    % per is the number of points to test per final point on the boundary
    %   - a higher per value will give more evenly-spaced points
    % plane consists of 4 values; the first 3 give the normal vector p to 
    % the plane, and the 4th gives the constant value. The plane is given 
    % by p*[x;y;z] - c = 0
    Ntest = N*per; % total number of points to generate
    x = zeros(Ntest, 1);
    y = zeros(Ntest, 1);
    z = zeros(Ntest, 1);
    ellx = xlim(2) - xlim(1);
    elly = ylim(2) - ylim(1);
    ellz = zlim(2) - zlim(1);
    % create Ntest points on the boundary via bisection
    for found = 1:Ntest
        % find a point with a positive level set value, and a point with a
        % negative level set value
        xtp = rand()*ellx + xlim(1);
        ytp = rand()*elly + ylim(1);
        ztp = rand()*ellz + zlim(1);
        xtn = rand()*ellx + xlim(1);
        ytn = rand()*elly + ylim(1);
        ztn = rand()*ellz + zlim(1);
        % move the points to the plane
        [xtp, ytp, ztp] = proj_plane(xtp, ytp, ztp, plane);
        [xtn, ytn, ztn] = proj_plane(xtn, ytn, ztn, plane);
        % keep generating new points until a negative point is found
        while lsf(xtn, ytn, ztn) > 0
            xtn = rand()*ellx + xlim(1);
            ytn = rand()*elly + ylim(1);
            ztn = rand()*ellz + zlim(1);
            [xtn, ytn, ztn] = proj_plane(xtn, ytn, ztn, plane);
        end
        % keep generating new points until a postive point is found
        while lsf(xtp, ytp, ztp) < 0
            xtp = rand()*ellx + xlim(1);
            ytp = rand()*elly + ylim(1);
            ztp = rand()*ellz + zlim(1);
            [xtp, ytp, ztp] = proj_plane(xtp, ytp, ztp, plane);
        end
        % bisect until within tolerance
        while vecnorm([xtp - xtn, ytp - ytn, ztp - ztn], 2, 2) > eps
            xtt = (xtp + xtn)/2;
            ytt = (ytp + ytn)/2;
            ztt = (ztp + ztn)/2;
            if lsf(xtt, ytt, ztt) > 0
                xtp = xtt;
                ytp = ytt;
                ztp = ztt;
            else
                xtn = xtt;
                ytn = ytt;
                ztn = ztt;
            end
        end
        % add new point on boundary
        x(found) = (xtp + xtn)/2;
        y(found) = (ytp + ytn)/2;
        z(found) = (ztp + ztn)/2;
    end
    % add points that are the farthest away from current points on boundary
    found = 1;
    xf = zeros(N, 1);
    yf = zeros(N, 1);
    zf = zeros(N, 1);
    % start with one point
    xf(1) = x(1);
    yf(1) = y(1);
    zf(1) = z(1);
    while found < N
        % indices of batch to test
        ran = (found*per+1):(found*per + per);
        % get point farthest from current points on boundary
        tests = min(((xf(1:found) - x(ran)').^2 + (yf(1:found) - y(ran)').^2 +...
            (zf(1:found) - z(ran)').^2), [], 1);
        [~, ind] = max(tests);
        found = found + 1;
        % add point to boundary
        xf(found) = x((found-1)*per+ind);
        yf(found) = y((found-1)*per+ind);
        zf(found) = z((found-1)*per+ind);
    end
end

function [xf, yf, zf] = proj_plane(x, y, z, plane)
    % projects a point (x,y,z) to a plane
    % plane consists of 4 values; the first 3 give the normal vector p to 
    % the plane, and the 4th gives the constant value. The plane is given 
    % by p*[x;y;z] - c = 0
   p = plane(1:3);
   v = plane(4);
   c = (dot([x,y,z],p)-v)/norm(p)^2;
   xf = x - c*p(1); yf = y - c*p(2); zf = z - c*p(3);
end



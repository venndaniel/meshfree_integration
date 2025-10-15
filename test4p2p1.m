% code for the test of Subsection 4.1.1
% integrates a logarithm on the unit circle
% numbers are recorded and plotted in fig4_code.m
% The number of points (Nt) must be changed to produce values in Figure 4

% set plots to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')
rng('default');

lsf = @(x,y) x.^2 + y.^2; % 2D level set describing domain
Nt = 1280; % number of interior points
omax = 30; % total number of Fourier basis functions is (2*omax + 1)^2
Nb = (2*omax + 1)^2;
ell = 4; % side length of bounding box

% shape parameters
q = 4; 
T = 10;

% place scattered points in circle with some processing to get roughly
% evenly spaced points
[xs, ys] = place_points(lsf, Nt, [-2, 2], [-2, 2], 20, 0);

% form omega_x, omega_y arrays
freq = -omax:omax;
[wx, wy] = meshgrid(freq);
wx = 2*pi/ell*reshape(wx, [Nb, 1]);
wy = 2*pi/ell*reshape(wy, [Nb, 1]);

wa = vecnorm([wx, wy], 2, 2); % norm of omega vector

di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(wa))); % d^{-1/2}

% singular function s(x) and derivatives
syms x y z
s = @(x,y) ((x).^2 + y.^2).*log((x).^2 + y.^2);
%s = @(x,y) sin(x + y);
sx = diff(s, x); sy = diff(s, y);
sxx = diff(sx, x); syy = diff(sy, y);
sx = matlabFunction(sx); sy = matlabFunction(sy);
sxx = matlabFunction(sxx); syy = matlabFunction(syy);

V = di'.*exp(1i*(wx'.*xs + wy'.*ys)); % funcation values

V2 = [V, s(xs, ys).*V]; % function values for augmented series

Vlap = -wa.^2'.*V; % Laplacian of un-augmented series

% Laplacian for augmented series
Vlap2 = [Vlap, (sxx(xs, ys) + syy(xs, ys)).*V +... 
    2i*(wx'.*sx(xs, ys) + wy'.*sy(xs, ys)).*V + s(xs, ys).*Vlap];

% compute solution coefficients
a = lsqminnorm(Vlap2, -1/2/pi*log(sqrt(xs.^2 + ys.^2)), 1e-20);

% create boundary points
theta = linspace(0, 2*pi, 1000+1);
theta = theta(1:1000)';
xb = cos(theta); yb = sin(theta);

% matrix to evaluate boundary values
Vbn2 = [1i*(wx'.*xb + wy'.*yb).*di'.*exp(1i*(wx'.*xb + wy'.*yb)),...
    1i*(wx'.*xb + wy'.*yb).*di'.*s(xb,yb).*exp(1i*(wx'.*xb + wy'.*yb)) +...
    di'.*(sx(xb,yb).*xb + sy(xb,yb).*yb).*exp(1i*(wx'.*xb + wy'.*yb))];

abs(sum(Vbn2*a)/1000*2*pi -1/4)*4 % compute relative error



function [x, y] = place_points(lsf, N, xlim, ylim, Ntest, ew)
    % places N points inside the 2D domain given by lsf < 1
    % xlim, ylim, specify a bounding box for the domain
    % Ntrest is the number of points to test per point (higher = more
    % evenly spaced
    % ew is outer weight (higher = more points near boundary, 0 for even)

    x = zeros(N, 1);
    y = zeros(N, 1);
    found = 0; % keep track of how many points have been added
    ellx = xlim(2) - xlim(1);
    elly = ylim(2) - ylim(1);
    while found < N
        foundt = 0;
        xtt = zeros(Ntest, 1);
        ytt = zeros(Ntest, 1);
        while foundt < Ntest % create Ntest points
            xt = rand()*ellx + xlim(1);
            yt = rand()*elly + ylim(1);
            % find new point if xt, yt is outside domain
            if lsf(xt, yt) < 1 
                foundt = foundt + 1;
                xtt(foundt) = xt;
                ytt(foundt) = yt;
            end
        end
       
        if found > 0    
            tests = (ew*lsf(xtt', ytt')+1).*min((x(1:found) - xtt').^2 + ...
                (y(1:found) - ytt').^2, [], 1);
            % pick point with the highest (weighted) distance from points
            % that have already been added
            [~, ind] = max(tests); 
            found = found + 1;
            x(found) = xtt(ind);
            y(found) = ytt(ind);
        else
            % if point cloud is empty, just add the first point
            found = found + 1;
            x(found) = xtt(1);
            y(found) = ytt(1);
        end
    end
end

% code for the test of Subsection 4.1.2
% integrates on a genus-2 surface with Voronoi cell subdomains
% also includes code for a comparison to integrating via triangulation
% The number of points (Nt) must be changed to produce values in Table 2

% LaTex for plots
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')
rng('default');

sur_tol = 1e-13; % tolerance for placing point on surface

% level set and derivatives
phi = @(x,y,z) 0.25./((x-1).^2 + y.^2) + 0.25./((x+1).^2 + y.^2)+ z.^2-...
    1 + 0.1*x.^2 + 0.25*y.^2;
phix = @(x,y,z)  -0.25*2*(x-1)./((x-1).^2 + y.^2).^2 -...
    0.25*2*(x+1)./((x+1).^2 + y.^2).^2 + 0.1*2*x;
phiy = @(x,y,z)  -0.25*2*y./((x-1).^2 + y.^2).^2 -...
    0.25*2*y./((x+1).^2 + y.^2).^2 + 0.25*2*y;
phiz = @(x,y,z)  2*z;

plot_patch = false; % show surface patch while computing


Npp = 100; % number of Voronoi cells
Nt = 8000; % total number of points
wmax = 11; % total number of Fourier basis functions per patch is (2*wmax + 1)^2
q = 5; %5
T = 5; %5 (1 is osc width of domain)
ell_mult = 2;
sobolev = false;
bpm = 2; % boundary point multiplier
minbp = 5; % minimum # of boundary points

neigh = 30; % number of Voronoi nodes "near" current to consider
% a Voronoi cell with 30 sides would be unlucky
corners = 30; % max number of corners to prepare to store

neigh = min(neigh, Npp-1); % readjust number of nodes to consider if Npp<30


"Creating Point Cloud"
[xpp, ypp, zpp] = place_b_points(phi, Npp, [-5,5], [-4, 4], [-2, 2],...
    sur_tol, 25); % Voronoi nodes

[xs, ys, zs] = place_b_points(phi, Nt, [-5,5], [-4, 4], [-2, 2],...
    sur_tol, 25); % Point cloud



"Done Point Cloud"

cind = zeros(Nt, 1); % store index of closest Voronoi node to each point
in_patch = zeros(Npp, 1);
for j = 1:Nt
    % find closest Voronoi node
    [~, ci] = mink(vecnorm([xs(j)-xpp,ys(j)-ypp,zs(j)-zpp], 2, 2),1);
    cind(j) = ci;
    in_patch(ci) = in_patch(ci) + 1; % count number of points in patch
end

%% Plot Voronoi Cells
scatter3(xs, ys, zs, 10, cind, 'filled');
%scatter3(xs, ys, zs, 10, 'k', 'filled');
axis('equal')
colormap('jet')
title('Voronoi Cells on Genus-2 Surface')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
view(-30,50)
fontsize(18, 'pixels')

%% For Triangulation Comparison

"Meshing"

% Create mesh of point cloud
ptCloud = pointCloud([xs ys zs]);
[mesh,~] = pc2surfacemesh(ptCloud,"ball-pivot");
area = 0;
verts = mesh.Vertices;
faces = mesh.Faces;
xa = 0;
"Integrating"
for j = 1:length(faces(:,1))
    v1 = verts(faces(j,1),:);
    v2 = verts(faces(j,2),:);
    v3 = verts(faces(j,3),:);
    avg = 1;
    % add triangle area to total area
    area = area + avg*norm(cross(v2-v1,v3-v1))/2;
    
end

area

%% Meshfree method with Voronoi Cells

weights = zeros(Nt, 1); % store integration weights

% keep track of the cell with the most corners (used for ensuring the
% choice of "neigh" is sufficient
max_corners = 0; 
for ppind = 1:Npp
    "Voronoi Cell: " + ppind + " of " + Npp % current Voronoi cell index
    inds1 = find(cind==ppind); % get indices of points in this cell
    if plot_patch % show patch if plot_patch has been set to true
        scatter3(xs(inds1), ys(inds1), zs(inds1), 10, 'b', 'filled');
        hold on;
        axis('equal');
    end

    % find neighbouring Voronoi cells
    dists = vecnorm([xpp(ppind) - xpp, ypp(ppind) - ypp, ...
        zpp(ppind) - zpp], 2, 2);
    % get list of nearby Voronoi node indices
    [~, test] = mink(dists, min(neigh+1,Npp));
    test = test(2:(neigh+1)); % ignore self
   

    % we need to find the corners of the Voronoi cell now
    % corners occur at a point equidistant from 3 Voronoi nodes

    % keep track of which potential corners are on the boundary of the cell
    in_shape = false(neigh^2/2, 1); 
    % to store location of potential corners
    xc = zeros(neigh, 1);
    yc = xc; zc = xc;
    cpp = zeros(neigh^2/2, 2); % to store indices of corners in
    c = 0; % index for in_shape, xc, yc, zc, cpp

    for n = 1:(neigh-1)
        for m = (n+1):neigh
            % indices of potential neibouring cells
            j = test(n); k = test(m); 
            
            c = c + 1; % update corner index

            % get Vornoi node locations to test
            p1 = [xpp(ppind); ypp(ppind); zpp(ppind)];
            p2 = [xpp(j); ypp(j); zpp(j)];
            p3 = [xpp(k); ypp(k); zpp(k)];

            % compute a point equidistant to all three nodes and on the
            % surface
            cen = cent_surf2(p1, p2, p3, phi, phix, phiy, phiz, sur_tol);

            % consider the potential corner to be an actual corner of the
            % cell if two conditions are met
            % 1. its distance to the nearest Voronoi node is within 
            % tolerance of its distance to the current node (it is on the 
            % boundary of the cell in Euclidean space
            % 2. it is at least 20% closer to a surface point in the cell
            % than it is to the node under consideration (it is "connected"
            % to the surface cell, and not on another side of the surface)
            if abs(min(vecnorm(cen' - [xpp ypp zpp], 2, 2) -...
                    (vecnorm(cen' - p1', 2, 2)))) < 2*sur_tol &&...
                min(vecnorm(cen' - [xs(inds1), ys(inds1), zs(inds1)], 2, 2)) <... 
                0.8*vecnorm(p1' - cen', 2, 2)
                % add point to list of corners if these conditions are met
                xc(c) = cen(1); yc(c) = cen(2); zc(c) = cen(3);
                % record that this cornser is on the cell boundary
                in_shape(c) = true; 
                % indices of other two cells connected to corner
                cpp(c,1) = j; cpp(c, 2) = k; 
            end
        end
    end
    % shrink corner location arrays and neighbouring cell arrays to only 
    % include corners that are actually on the cell boundary
    xc = xc(in_shape); yc = yc(in_shape); zc = zc(in_shape); 
    cpp = cpp(in_shape,:);
    if plot_patch
        scatter3(xc, yc, zc, 40, 'k', 'filled');
    end
   
    % set up integration procedure
    
    % get points in patch
    x = xs(inds1); y = ys(inds1); z = zs(inds1);

    % set bounding box for Fourier series
    ellx = ell_mult*(max([x;xc]) - min([x;xc]));
    elly = ell_mult*(max([y;yc]) - min([y;yc]));
    ellz = ell_mult*(max([z;zc]) - min([z;zc]));

    Nb = (2*wmax + 1).^3; % number of basis functions
    freq = -wmax:wmax;
    [wx, wy, wz] = meshgrid(freq, freq, freq); % omega_x,y,z
    wx = reshape(wx, [Nb, 1])*2*pi/ellx;
    wy = reshape(wy, [Nb, 1])*2*pi/elly;
    wz = reshape(wz, [Nb, 1])*2*pi/ellz;
    
    wa = vecnorm([wx wy wz], 2, 2); % norm of omega
    waw = vecnorm([wx*ellx wy*elly wz*ellz], 2, 2); % rescale omega norm
    
    di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(waw))); % d_n^{-1/2}
    if sobolev
        di = 1./((2*pi/T).^q + waw.^q);
    end
    
    
    h = sqrt(4*pi/Nt); % point spacing for boundary integral
    Vban = []; % store outward conormal vector to patch for later use
    wba = []; % boundary integral weights
    max_corners = max(length(xc), max_corners);
    % cycle through edges by checking pairs of corners
    for j = 1:length(xc)-1
        for k = (j+1):length(xc)
            % check if corners share an edge by checking if they share a
            % neighbouring Voronoi cell
            if min(min(abs(cpp(j, :) - cpp(k, :)'))) == 0 

                % distance between corners
                dist = vecnorm([xc(j) - xc(k), yc(j) - yc(k),...
                    zc(j) - zc(k)], 2, 2);
                % number of points on this edge
                Nbound = max(ceil(bpm*dist/h), minbp);
                th = dist/Nbound; % point spacing
                % parametrize line between corners
                t = linspace(1/Nbound/2, 1-1/Nbound/2, Nbound)';
                % set initial boundary points on this line
                xb = (1-t)*xc(j) + t*xc(k);
                yb = (1-t)*yc(j) + t*yc(k);
                zb = (1-t)*zc(j) + t*zc(k);

                % get index of shared Voronoi point
                % candidates (copy cpp(j,:) to two rows)
                S = cpp(j, :) - [0;0]; 
                shared = S((cpp(j, :) - cpp(k, :)')==0); % shared index

                % normal vector to plane separating cells
                pnx = (xpp(shared) - xpp(ppind))*ones(Nbound, 1);
                pny = (ypp(shared) - ypp(ppind))*ones(Nbound, 1);
                pnz = (zpp(shared) - zpp(ppind))*ones(Nbound, 1);
                pnn = vecnorm([pnx, pny, pnz], 2, 2);
                pnx = pnx./pnn; pny = pny./pnn; pnz = pnz./pnn;

                % move points along the plane back to the surface
                for n = 1:Nbound
                    [xbb, ybb, zbb] = plane_inter(xb(n), yb(n), zb(n), ...
                        phi, phix, phiy, phiz, sur_tol, 1, ...
                        pnx(1), pny(1), pnz(1));
                    xb(n) = xbb; yb(n) = ybb; zb(n) = zbb;
                end
                
                
                % surface norm vector at boundary
                [bnx, bny, bnz] = nml(xb, yb, zb, phix, phiy, phiz);
                % compute tangent vector to boundary curve
                Ta = cross([pnx, pny, pnz], [bnx, bny, bnz], 2);
                Tan = vecnorm(Ta, 2, 2);
                tx = Ta(:, 1)./Tan; ty = Ta(:, 2)./Tan; tz = Ta(:, 3)./Tan;
                % compute conormal vector to boundary
                Nu = cross([bnx, bny, bnz], [tx, ty, tz], 2);
                nun = vecnorm(Nu, 2, 2);
                nux = Nu(:, 1)./nun; nuy = Nu(:, 2)./nun; nuz = Nu(:, 3)./nun;
    
                % solve first order problem to compute boundary weights
                % boundary values
                Vb = di'.*exp(1i*(wx'.*xb + wy'.*yb + wz'.*zb));
                % tangent derivative
                Vbt = 1i*(wx'.*tx + wy'.*ty + wz'.*tz).*Vb;
                % initial value
                Vi = di'.*exp(1i*(wx'.*xc(j) + wy'.*yc(j) + wz'.*zc(j)));
                % end value
                Ve = di'.*exp(1i*(wx'.*xc(k) + wy'.*yc(k) + wz'.*zc(k)));
                % matrix for solving diff. eq. with initial value 0
                Vs = [Vbt; Vi];
                % compute boundary weights
                wb = Ve*lsqminnorm(Vs,[eye(Nbound);zeros(1, Nbound)]);
                wb = wb'./sign(sum(wb)); % correct sign
                if plot_patch % plot boundaries if plot_patch
                    quiver3(xb, yb, zb, tx/20, ty/20, tz/20, 0);
                    quiver3(xb, yb, zb, nux/20, nuy/20, nuz/20, 0);
                end
                
                % conormal derivative for this edge
                Vbn = 1i*(wx'.*nux + wy'.*nuy + wz'.*nuz).*Vb;
                Vban = [Vban; Vbn]; % store conormal derivatives
                wba = [wba; wb]; % store boundary weights
            end
        end
    end
    
    % surface Poisson problem on a patch

    wba = real(wba); % force real weights
    N = length(x); % number of points in patch
    Nbound = length(wba); % number of points on boundary
    V = di'.*exp(1i*(x.*wx' + y.*wy' + z.*wz')); % function value
    [nx, ny, nz] = nml(x, y, z, phix, phiy, phiz); % surface normal
    Vlap = -wa.^2'.*V; % Laplacian
    Vn = 1i*(nx.*wx' + ny.*wy' + nz.*wz').*V; % first normal derivative
    Vnn = -(nx.^2.*wx.^2' + ny.^2.*wy.^2' + nz.^2.*wz.^2' + ...
    2*nx.*ny.*(wx.*wy)' + ...
    2*nx.*nz.*(wx.*wz)' + ...
    2*nz.*ny.*(wz.*wy)').*V; % second normal derivative
    
    Vs = [(Vlap - Vnn); Vn]; % matrix for Eq. (3.17)
    % compute weights by evaluating conormal derivative and using the
    % boundary weights to compute the line integrals
    w = wba'*(Vban*lsqminnorm(Vs, [eye(N); zeros(N)]));
    w = real(w'); % force real weights
    
    weights(inds1) = w; % set weights for this patch to w
    
    if plot_patch % reset the patch plot if plot_patch
        hold off;
        pause(0.001);
    end
end

sum(weights) % surface area!

    

function [cx, cy, cz] = cyl_inter(x, y, z, nx, ny, nz, phi, phix, phiy, ...
    phiz, tol)
    % moves a point (x,y,z) to a surface specified by the phi level set 
    % moving along a line specified by the direction [nx, ny, nz]
    % phix, phiy, phiz are partial derivatives of phi used for Newton's
    % method, tol is tolerance for Newton's method
    cx = x; cy = y; cz = z;
    iters = 4000;
    iter = 0;
    while abs(phi(cx, cy, cz)) > tol && iter < iters
        iter = iter + 1;
        px = phix(cx,cy,cz); py = phiy(cx,cy,cz); pz = phiz(cx,cy,cz);
        p = phi(cx,cy,cz);
        dn = dot([nx, ny, nz], [px, py, pz]);
        cx = cx - p*nx/dn;
        cy = cy - p*ny/dn;
        cz = cz - p*nz/dn;
        iter = iter + 1;
    end
end

function [cx, cy, cz] = plane_inter(x, y, z, phi, phix, phiy, phiz, tol, ...
    w, plx, ply, plz)
    % moves a point (x,y,z) to a surface specified by the phi level set 
    % moving along a plane specified by the normal vector [plx, ply, plz]
    % phix, phiy, phiz are partial derivatives of phi used for Newton's
    % method, tol is tolerance for Newton's method
    % plx, ply, plz must be normalized
    % w is a weight for Newton's method that will decrease automatically
    % if Newton's method fails
    err = 2*tol; 
    cx = x; cy = y; cz = z; % current x,y,z
    iters = round(4000/w^2); % max iters to attempt
    iter = 1; % current iter
    % continue moving point until on surface or max iters exceeded
    while abs(phi(cx, cy, cz)) > tol && iter < iters
        iter = iter + 1;
        % partial derivatives of level set
        px = phix(cx,cy,cz);
        py = phiy(cx,cy,cz);
        pz = phiz(cx,cy,cz);
        p = phi(cx,cy,cz); % current level set value
        ng = vecnorm([px, py, pz], 2, 2);
        dt = dot([px py pz], [plx ply plz]);
        % project level set gradient onto plane
        px = px - dt*plx; py = py - dt*ply; pz = pz - dt*plz;

        % update point
        cx = cx - w*p*px/ng^2;
        cy = cy - w*p*py/ng^2;
        cz = cz - w*p*pz/ng^2;
    end
     
    % if Newton's method failed, try again with a smaller weight
    if iter == iters
        'fail'
        [cx, cy, cz] = plane_inter(x, y, z, phi, phix, phiy, phiz, tol, ...
            0.9*w, plx, ply, plz);
    end
end


function [nx, ny, nz] = nml(x, y, z, phix, phiy, phiz)
    % normal vector to surface given by a level set with partial
    % derivatives phix, phiy, phiz
    nx = phix(x,y,z); ny = phiy(x,y,z); nz = phiz(x,y,z);
    nm = vecnorm([nx ny nz],2,2);
    nx = nx./nm; ny = ny./nm; nz = nz./nm;
end


function cen = cent(p1, p2, p3)
    % computes centre of three points p1, p2, p3 lying in the same plane
    v12 = p2 - p1;
    v13 = p3 - p1;
    v = ([v12'; v13']*[v12 v13])\([v12'*v12/2; v13'*v13/2]);
    cen = p1 + [v12 v13]*v;
end


function cen = cent_surf2(p1, p2, p3, phi, phix, phiy, phiz, tol)
    % finds a point on the surface specified by phi that is equidistant
    % to p1, p2, p3
    % phix, phiy, phiz are the partial derivatives of phi needed for
    % Newton's method

    % displacements between points
    v12 = p2 - p1;
    v13 = p3 - p1;

    % normal vector to the triangle formed by the three points
    N = cross(v12', v13'); 
    nx = N(1); ny = N(2); nz = N(3);
    nn = vecnorm([nx ny nz], 2, 2);
    nx = nx/nn; ny = ny/nn; nz = nz/nn;

    cen = cent(p1, p2, p3); % point in the centre of the three points
    % move point to the surface along the direction given by the normal to
    % the triangle
    [cenx, ceny, cenz] = cyl_inter(cen(1), cen(2), cen(3), nx, ny, nz, phi, ...
        phix, phiy, phiz, tol);
    cen = [cenx; ceny; cenz];
end


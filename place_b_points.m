function [xf, yf, zf] = place_b_points(lsf, N, xlim, ylim, zlim, eps, per)
    % creates N points on a surface specified by the level set function lsf
    % xlim, ylim, zlim specify the bounding box for the surface
    % eps is the tolerance for points on the surface
    % per is the number of points to test for each point in the final point
    % cloud

    Ntest = N*per; % total number of points to test
    x = zeros(Ntest, 1);
    y = zeros(Ntest, 1);
    z = zeros(Ntest, 1);
    ellx = xlim(2) - xlim(1);
    elly = ylim(2) - ylim(1);
    ellz = zlim(2) - zlim(1);
    % create Ntest points on the surface via bisection
    for found = 1:Ntest
        % find a point with a positive level set value, and a point with a
        % negative level set value
        xtp = rand()*ellx + xlim(1);
        ytp = rand()*elly + ylim(1);
        ztp = rand()*ellz + zlim(1);
        xtn = rand()*ellx + xlim(1);
        ytn = rand()*elly + ylim(1);
        ztn = rand()*ellz + zlim(1);
        % keep generating new points until a negative point is found
        while lsf(xtn, ytn, ztn) > 0
            xtn = rand()*ellx + xlim(1);
            ytn = rand()*elly + ylim(1);
            ztn = rand()*ellz + zlim(1);
        end
        % keep generating new points until a positive point is found
        while lsf(xtp, ytp, ztp) < 0
            xtp = rand()*ellx + xlim(1);
            ytp = rand()*elly + ylim(1);
            ztp = rand()*ellz + zlim(1);
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
        % add new point on the surface
        x(found) = (xtp + xtn)/2;
        y(found) = (ytp + ytn)/2;
        z(found) = (ztp + ztn)/2;
    end
    % add points that are the farthest away from the current point cloud
    found = 1;
    xf = zeros(N, 1);
    yf = zeros(N, 1);
    zf = zeros(N, 1);
    xf(1) = x(1);
    yf(1) = y(1);
    zf(1) = z(1);
    while found < N
        % indices of batch to test
        ran = (found*per+1):(found*per + per);
        % get point farthest from current points on surface
        tests = min(((xf(1:found) - x(ran)').^2 + (yf(1:found) - y(ran)').^2 +...
            (zf(1:found) - z(ran)').^2), [], 1);
        [~, ind] = max(tests);
        found = found + 1;
        % add point to point cloud
        xf(found) = x((found-1)*per+ind);
        yf(found) = y((found-1)*per+ind);
        zf(found) = z((found-1)*per+ind);
    end
end
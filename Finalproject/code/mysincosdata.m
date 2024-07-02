function pde = mysincosdata
%% SINCOSDATA trigonometric data for Poisson equation with mass term
%
%     f = (2*pi^2 + 2) * sin(pi*x) * sin(pi*y);
%     u = sin(pi*x) * sin(pi*y);
%     Du = (pi*cos(pi*x)*sin(pi*y), pi*sin(pi*x)*cos(pi*y));
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f', @f, 'exactu', @exactu, 'g_D', @g_D, 'Du', @Du);

    % load data (right hand side function)
    function rhs = f(p)
        x = p(:,1); y = p(:,2);
        rhs = (2*pi^2 + 2) * sin(pi*x) .* sin(pi*y);
    end
    % exact solution
    function u = exactu(p)
        x = p(:,1); y = p(:,2);
        u = sin(pi*x) .* sin(pi*y);
    end
    % Dirichlet boundary condition
    function u = g_D(p)
        u = exactu(p);
        % ?
        % u = 0;
    end
    % Derivative of the exact solution
    function uprime = Du(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = pi * cos(pi*x) .* sin(pi*y);
        uprime(:,2) = pi * sin(pi*x) .* cos(pi*y);
    end
end

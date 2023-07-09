% Test the discretization of the fractional gradient
% Created July 1, 2022

%f = @(x, y, z) x .* y .* z + y.^3 .* z - x;
f = @(x, y, z) (x + y.^2) .* y .* z - x;

x1 = 0.0; x2 = 2.0;
y1 = 1.0; y2 = 3.0;
z1 = 0.5; z2 = 5.0;

% fractional order of fractional gradient
alpha = 0.6;

% s = n-alpha; order of fractional integration
s = 1 - alpha;

subdivisions_vec = [1 2 4 8 16];
n_samples = length(subdivisions_vec);
errors_1_center = zeros(n_samples, 1);
errors_1_integrated = zeros(n_samples, 1);
errors_2_center = zeros(n_samples, 1);
x_errors_1_center = zeros(n_samples, 1);
y_errors_1_center = zeros(n_samples, 1);
y_errors_line_1_center = zeros(n_samples, 1);
y_errors_line_1_integrated = zeros(n_samples, 1);
y_errors_line_center_vs_integrated = zeros(n_samples, 1);

for iter = 1:n_samples
    nx = subdivisions_vec(iter);
    ny = subdivisions_vec(iter);
    nz = subdivisions_vec(iter);
%       nz = 0;

    [D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);
    
    ni = nx + 1;
    nj = ny + 1;
    nk = nz + 1;
    
    dx = (x2 - x1) / nx;
    dy = (y2 - y1) / ny;
    dz = (z2 - z1) / nz;
    
    n_nodes = size(D0, 2);
    values_on_nodes = zeros(n_nodes, 1);

    
    % Fill values on nodes
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                x = x1 + (i-1) * dx;
                y = y1 + (j-1) * dy;
                z = z1 + (k-1) * dz;
                
                node = node_index(i, j, k);
                values_on_nodes(node) = f(x, y, z);
            end
        end
    end
    
    
    B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);

    M1_1ps = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1+s);
    M1_1ms = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1-s);
    
    I1 = (1 ./ edge_lengths) .* (B1 * M1_1ps) .* (1 ./ edge_lengths)';

    discrete_frac_gradient_v1 =               I1 * D0 * values_on_nodes;
    discrete_frac_gradient_v2 = (B1 * M1_1ms)^-1 * D0 * values_on_nodes;


    n_edges = size(D0, 1);
    continuous_frac_gradient_center = zeros(n_edges, 1);
    continuous_frac_gradient_integrated = zeros(n_edges, 1);

    x_err_v1_center = 0.0;
    y_err_v1_center = 0.0;
    y_err_line_1_center = 0.0;
    y_err_line_1_integrated = 0.0;
    y_err_line_center_vs_integrated = 0.0;

    % Fill in continuous values at centers of edges
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                x = x1 + (i-1) * dx;
                y = y1 + (j-1) * dy;
                z = z1 + (k-1) * dz;

                % x edge
                if i ~= ni
                    edge = x_edge_index(i, j, k);
                    x = edge_coords(edge, 1);
                    y = edge_coords(edge, 2);
                    z = edge_coords(edge, 3);
                    xa = x - dx/2;
                    xb = x + dx/2;

                    continuous_frac_gradient_center(edge) = (x+(-1).*x1).^(1+(-1).*alpha).*((-1)+y.*z).*gamma(2+(-1).*alpha).^(-1);
                    if i == 1
                        continuous_frac_gradient_integrated(edge) = (((-1).*x1+xb).^(2+(-1).*alpha).*((-1)+y.*z).*gamma(3+(-1).*alpha).^(-1)) / dx;
                    else
                        continuous_frac_gradient_integrated(edge) = (((x1+(-1).*xa).*(x1+(-1).*xb)).^((-1).*alpha).*(((-1).*x1+xa) ...
      .^alpha.*(x1+(-1).*xb).^2+(-1).*(x1+(-1).*xa).^2.*((-1).*x1+xb) ...
      .^alpha).*((-1)+y.*z).*gamma(3+(-1).*alpha).^(-1)) / dx;
                    end
                    
                    x_err_v1_center = x_err_v1_center + abs(continuous_frac_gradient_center(edge) - discrete_frac_gradient_v1(edge));
                end
    
                % y edge
                if j ~= nj
                    edge = y_edge_index(i, j, k);
                    x = edge_coords(edge, 1);
                    y = edge_coords(edge, 2);
                    z = edge_coords(edge, 3);
                    ya = y - dy/2;
                    yb = y + dy/2;

                    continuous_frac_gradient_center(edge) = (y+(-1).*y1).^(1+(-1).*alpha).*(((-3)+alpha).*((-2)+alpha).*x+6.* ...
  y.*(y+y1)+3.*y1.*((-2).*alpha.*y+((-2)+alpha).*((-1)+alpha).*y1)) ...
  .*z.*gamma(4+(-1).*alpha).^(-1);
                    if j == 1
                        continuous_frac_gradient_integrated(edge) = (((-1).*y1+yb).^(2+(-1).*alpha).*(((-3)+alpha).*(((-4)+alpha).*x+ ...
  3.*((-2)+alpha).*y1.^2)+(-6).*((-2)+alpha).*y1.*yb+6.*yb.^2).*z.* ...
  gamma(5+(-1).*alpha).^(-1)) / dy;
                    else
                        continuous_frac_gradient_integrated(edge) = ((-1).*((y1+(-1).*ya).*(y1+(-1).*yb)).^((-1).*alpha).*((y1+(-1).* ...
      ya).^2.*(((-3)+alpha).*(((-4)+alpha).*x+3.*((-2)+alpha).*y1.^2)+( ...
      -6).*((-2)+alpha).*y1.*ya+6.*ya.^2).*((-1).*y1+yb).^alpha+(-1).*(( ...
      -1).*y1+ya).^alpha.*(y1+(-1).*yb).^2.*(((-3)+alpha).*(((-4)+alpha) ...
      .*x+3.*((-2)+alpha).*y1.^2)+(-6).*((-2)+alpha).*y1.*yb+6.*yb.^2)) ...
      .*z.*gamma(5+(-1).*alpha).^(-1)) / dy;
                    end
                    y_err_v1_center = y_err_v1_center + abs(continuous_frac_gradient_center(edge) - discrete_frac_gradient_v1(edge));
                    if i == ni && k == nk
                        y_err_line_1_center = y_err_line_1_center + abs(continuous_frac_gradient_center(edge) - discrete_frac_gradient_v1(edge));
                        y_err_line_1_integrated = y_err_line_1_integrated + abs(continuous_frac_gradient_integrated(edge) - discrete_frac_gradient_v1(edge));
                        y_err_line_center_vs_integrated = y_err_line_center_vs_integrated + abs(continuous_frac_gradient_center(edge) - continuous_frac_gradient_integrated(edge));
                    end
                end
    
                % z edge
                if k ~= nk
                    edge = z_edge_index(i, j, k);
                    x = edge_coords(edge, 1);
                    y = edge_coords(edge, 2);
                    z = edge_coords(edge, 3);
                    za = z - dz/2;
                    zb = z + dz/2;

                    continuous_frac_gradient_center(edge) = y.*(x+y.^2).*(z+(-1).*z1).^(1+(-1).*alpha).*gamma(2+(-1).*alpha).^(-1);
                    if k == 1
                        continuous_frac_gradient_integrated(edge) = (y.*(x+y.^2).*((-1).*z1+zb).^(2+(-1).*alpha).*gamma(3+(-1).*alpha) ...
  .^(-1)) / dz;
                    else
                        continuous_frac_gradient_integrated(edge) = (y.*(x+y.^2).*((z1+(-1).*za).*(z1+(-1).*zb)).^((-1).*alpha).*(((-1) ...
      .*z1+za).^alpha.*(z1+(-1).*zb).^2+(-1).*(z1+(-1).*za).^2.*((-1).* ...
      z1+zb).^alpha).*gamma(3+(-1).*alpha).^(-1)) / dz;
                    end
                end

            end
        end
    end
    
    n_x_edges = nx * nj * nk;
    n_y_edges = ni * ny * nk;
    x_err_v1_center = x_err_v1_center / n_x_edges;
    y_err_v1_center = x_err_v1_center / n_y_edges;
    y_err_line_1_center = y_err_line_1_center / ny;
    y_err_line_1_integrated = y_err_line_1_integrated / ny;
    y_err_line_center_vs_integrated = y_err_line_center_vs_integrated / ny;

    errors_1_center(iter) = mean(abs(discrete_frac_gradient_v1 - continuous_frac_gradient_center));
    errors_1_integrated(iter) = mean(abs(discrete_frac_gradient_v1 - continuous_frac_gradient_integrated));
    errors_2_center(iter) = mean(abs(discrete_frac_gradient_v2 - continuous_frac_gradient_center));
    x_errors_1_center(iter) = x_err_v1_center;
    y_errors_1_center(iter) = y_err_v1_center;
    y_errors_line_1_center(iter) = y_err_line_1_center;
    y_errors_line_1_integrated(iter) = y_err_line_1_integrated;
    y_errors_line_center_vs_integrated(iter) = y_err_line_center_vs_integrated;
end

% Compare methods
% loglog(subdivisions_vec, errors_1, subdivisions_vec, errors_2);
% hold on
% loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-1., '--', 'Color', '#7E2F8E');
% loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
% ylim([1e-5 1e2]);
% legend('error v1 and v2', 'error v3', '1 order convergence', '2 order convergence');

% Look at just version 2
loglog(subdivisions_vec, y_errors_line_1_integrated);
hold on
loglog(subdivisions_vec, 1.0 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
legend('error', '2nd order convergence');


xlabel('subdivisions');
ylabel('mean absolute error');
title(sprintf('Error of discrete fractional gradient for alpha=%g', alpha));

% Test the discretization of the fractional curl -- version 3
% Created July 5, 2023

% Vector field we are using:
% F[{x_, y_, z_}] := {
%    y Sin[5 x y + z] + 3 (x z - 1/2)^2 - 3 (y - 1/2)^2,
%    -y^3 + x z + z Cos[10 x y z],
%    x y (z - 1/4) + y^3 z - x + Cos[2 x y z] + 2 Sin[5 x^3 y]
%    };

err_order = 2.0;

x1 = 0.0; x2 = 1.0;
y1 = 0.0; y2 = 1.0;
z1 = 0.0; z2 = 1.0;

% This controls which plot to generate
change_alpha = true;

if change_alpha % change alpha
    s_vec = [logspace(-8, -1, 40) linspace(1e-1, 0.9999, 40)];
    alpha_vec = 1 - s_vec;
    n_samples = length(alpha_vec);
    const_subdivs = 16;
    subdivisions_vec = const_subdivs * ones(n_samples, 1);
else % Chnage number of subdivisions
    subdivisions_vec = [2 4 8 16 32];
    n_samples = length(subdivisions_vec);
    const_alpha = 0.25;
    alpha_vec = const_alpha * ones(n_samples, 1);
end

errors = zeros(n_samples, 1);

for iter = 1:n_samples
    tic;

    alpha = alpha_vec(iter);
    s = 1 - alpha;

    calculate_nonalpha_things = ~change_alpha | (iter == 1);
    
    if calculate_nonalpha_things
        nx = subdivisions_vec(iter);
        ny = subdivisions_vec(iter);
        nz = subdivisions_vec(iter);
    
        [D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);
        
        xs = node_coords(:, 1);
        ys = node_coords(:, 2);
        zs = node_coords(:, 3);
        
        %% Compute values needed for discrete fractional curl
        values_x_from_0 = (-3).*xs.*((-1)+ys).*ys+(-3/2).*xs.^2.*zs+xs.^3.*zs.^2+(1/5).*cos( ...
      zs)+(-1/5).*cos(5.*xs.*ys+zs);
        values_y_from_0 = (-1/4).*ys.^4+xs.*ys.*zs+(1/10).*xs.^(-1).*sin(10.*xs.*ys.*zs);
        values_z_from_0 = (1/4).*(zs.*((-1).*xs.*(4+ys)+2.*ys.*(xs+ys.^2).*zs+8.*sin(5.* ...
      xs.^3.*ys))+2.*xs.^(-1).*ys.^(-1).*sin(2.*xs.*ys.*zs));
    end
    
    %% Compute values needed for continuous fractional curl
    continuous_frac_curl_integrated_to_nodes_yz = (1/4).*ys.*zs.*(4.*((-1)+xs).*ys.^((-1).*alpha).*gamma(2+(-1).* ...
  alpha).^(-1)+(-1).*ys.^3.*zs.^((-1).*alpha).*gamma(2+(-1).*alpha) ...
  .^(-1)+(-1).*((-2)+alpha).^(-1).*((-1)+alpha).^(-1).*zs.^((-1).* ...
  alpha).*gamma(1+(-1).*alpha).^(-1).*(((-2)+alpha).*ys.^3+4.*xs.* ...
  zs+4.*zs.*hypergeom([1],[(3/2)+(-1/2).*alpha,2+(-1/2).*alpha],( ...
  -25).*xs.^2.*ys.^2.*zs.^2))+(-1).*((-1)+alpha).^(-1).*ys.^((-1).* ...
  alpha).*gamma(1+(-1).*alpha).^(-1).*(((-2)+alpha).^(-1).*((-12).*( ...
  12+(-7).*alpha+alpha.^2).^(-1).*ys.^3.*zs+xs.*(8+(-4).*alpha+ys+( ...
  -2).*ys.*zs))+(-40).*((-2)+alpha).^(-1).*xs.^3.*ys.*hypergeom([1], ...
  [(3/2)+(-1/2).*alpha,2+(-1/2).*alpha],(-25/4).*xs.^6.*ys.^2)+4.* ...
  hypergeom([(1/2),1],[(3/2),1+(-1/2).*alpha,(3/2)+(-1/2).*alpha],( ...
  -1).*xs.^2.*ys.^2.*zs.^2)));

    continuous_frac_curl_integrated_to_nodes_xz = (1/60).*zs.*(30.*xs.^(1+(-1).*alpha).*(2+ys.^3.*zs).*gamma(2+(-1) ...
  .*alpha).^(-1)+12.*zs.^((-1).*alpha).*((-1)+15.*xs.*((-1)+ys).*ys+ ...
  cos(5.*xs.*ys)).*gamma(2+(-1).*alpha).^(-1)+5.*xs.^(1+(-1).*alpha) ...
  .*(3.*((-2)+alpha).^(-1).*((-1)+alpha).^(-1).*gamma(1+(-1).*alpha) ...
  .^(-1).*(2.*((-2)+alpha).*ys.^3.*zs+xs.*(4+ys+(-2).*ys.*zs)+4.*(( ...
  -2)+alpha).*hypergeom([(1/2),1],[(3/2),1+(-1/2).*alpha,(3/2)+( ...
  -1/2).*alpha],(-1).*xs.^2.*ys.^2.*zs.^2))+(-5).*2.^(1+alpha).*3.^( ...
  (1/2).*alpha).*pi.^(3/2).*xs.^3.*ys.*gamma(3+(-1/2).*alpha).^(-1) ...
  .*gamma((5/6)+(-1/6).*alpha).^(-1).*gamma((7/6)+(-1/6).*alpha).^( ...
  -1).*gamma((3/2)+(-1/6).*alpha).^(-1).*hypergeom([(2/3),(5/6),1,( ...
  7/6),(4/3)],[(5/6)+(-1/6).*alpha,1+(-1/6).*alpha,(7/6)+(-1/6).* ...
  alpha,(4/3)+(-1/6).*alpha,(3/2)+(-1/6).*alpha,(5/3)+(-1/6).* ...
  alpha],(-25/4).*xs.^6.*ys.^2))+(-6).*((-1)+alpha).^(-1).*zs.^((-1) ...
  .*alpha).*gamma(1+(-1).*alpha).^(-1).*(5.*xs.*(6.*ys+(-6).*ys.^2+( ...
  (-3)+alpha).^(-1).*((-2)+alpha).^(-1).*xs.*zs.*((-9)+3.*alpha+4.* ...
  xs.*zs))+4.*hypergeom([1],[1+(-1/2).*alpha,(3/2)+(-1/2).*alpha],( ...
  -1/4).*zs.^2).*sin((5/2).*xs.*ys).^2+(-2).*((-2)+alpha).^(-1).* ...
  zs.*hypergeom([1],[(3/2)+(-1/2).*alpha,2+(-1/2).*alpha],(-1/4).* ...
  zs.^2).*sin(5.*xs.*ys)));

    continuous_frac_curl_integrated_to_nodes_xy = (1/20).*ys.*(5.*xs.^(1+(-1).*alpha).*(ys.^3+(-4).*zs).*gamma(2+( ...
  -1).*alpha).^(-1)+10.*xs.^2.*ys.^((-1).*alpha).*zs.*((-3)+2.*xs.* ...
  zs).*gamma(2+(-1).*alpha).^(-1)+5.*((-2)+alpha).^(-1).*((-1)+ ...
  alpha).^(-1).*xs.^(1+(-1).*alpha).*gamma(1+(-1).*alpha).^(-1).*((( ...
  -2)+alpha).*ys.^3+4.*xs.*zs+(-4).*((-2)+alpha).*zs.*hypergeom([( ...
  1/2),1],[(3/2),1+(-1/2).*alpha,(3/2)+(-1/2).*alpha],(-25).*xs.^2.* ...
  ys.^2.*zs.^2))+2.*((-1)+alpha).^(-1).*ys.^((-1).*alpha).*gamma(1+( ...
  -1).*alpha).^(-1).*((-30).*((-2)+alpha).^(-1).*xs.*ys+(-60).*(6+( ...
  -5).*alpha+alpha.^2).^(-1).*xs.*ys.^2+(-15).*xs.^2.*zs+10.*xs.^3.* ...
  zs.^2+2.*cos(zs)+(-2).*cos(zs).*hypergeom([1],[1+(-1/2).*alpha,( ...
  3/2)+(-1/2).*alpha],(-25/4).*xs.^2.*ys.^2)+(-10).*((-2)+alpha).^( ...
  -1).*xs.*ys.*hypergeom([1],[(3/2)+(-1/2).*alpha,2+(-1/2).*alpha],( ...
  -25/4).*xs.^2.*ys.^2).*sin(zs)));

    %% Fix the cases where x, y, or z is 0 and there could be a NaN
    ni = nx + 1;
    nj = ny + 1;
    nk = nz + 1;
    
    if x1 == 0.0
        % x == 0 (i == 1)
        for j = 1:nj     % for all y
            for k = 1:nk % for all z
                node = node_index(1, j, k);
                continuous_frac_curl_integrated_to_nodes_xz(node) = 0.0;
                continuous_frac_curl_integrated_to_nodes_xy(node) = 0.0;
                
                pos = node_coords(node, :);
                x = pos(1);
                y = pos(2);
                z = pos(3);
                values_y_from_0(node) = (-1/4).*y.^4+y.*z;
                values_z_from_0(node) = (1/2).*z.*(2+y.^3.*z);
            end
        end
    end
    
    if y1 == 0.0
        % y == 0 (j == 1)
        for i = 1:ni     % for all x
            for k = 1:nk % for all z
                node = node_index(i, 1, k);
                continuous_frac_curl_integrated_to_nodes_yz(node) = 0.0;
                continuous_frac_curl_integrated_to_nodes_xy(node) = 0.0;
                
                pos = node_coords(node, :);
                x = pos(1);
                y = pos(2);
                z = pos(3);
                values_z_from_0(node) = z+(-1).*x.*z;
            end
        end
    end

    if z1 == 0.0
        % z == 0 (k == 1)
        for i = 1:ni     % for all x
            for j = 1:nj % for all y
                node = node_index(i, j, 1);
                continuous_frac_curl_integrated_to_nodes_yz(node) = 0.0;
                continuous_frac_curl_integrated_to_nodes_xz(node) = 0.0;
            end
        end
    end

    
    %% Calculate discrete fractional curl
    if calculate_nonalpha_things
        % Create values_on_edges
        F_integrated_to_nodes_on_edges = [values_x_from_0; values_y_from_0; values_z_from_0];
        B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);
        values_on_edges = B1 * F_integrated_to_nodes_on_edges;
    
        B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
    end

    % Create I1
    M1_1ps = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1+s);
    I1 = B1 * M1_1ps .* (1 ./ edge_lengths)';
    
    % Create I2
    M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
    I2 = B2 * M2_1ps .* (1 ./ face_areas)';

    discrete_frac_curl_unintegrated = (1 ./ face_areas) .* (I2 * D1 * I1^-1 * values_on_edges);

    %% Calculate continuous fractional curl
    continuous_frac_curl_integrated_to_nodes = [continuous_frac_curl_integrated_to_nodes_yz; continuous_frac_curl_integrated_to_nodes_xz; continuous_frac_curl_integrated_to_nodes_xy];
    continuous_frac_curl_unintegrated = (1 ./ face_areas) .* (B2 * continuous_frac_curl_integrated_to_nodes);
    
    %% Calculate error
    errors(iter) = mean(abs(discrete_frac_curl_unintegrated - continuous_frac_curl_unintegrated).^err_order)^(1/err_order);

    T = toc;
    fprintf("%d,%d, %f\n", iter, nx, T);
end

if ~change_alpha
    loglog(subdivisions_vec, errors);
    hold on
    loglog(subdivisions_vec, 4 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
    xlim([subdivisions_vec(1) subdivisions_vec(end)]);
    xlabel('subdivisions');
    ylabel(sprintf('l^%g error', err_order));
    set(gca, 'XTick', subdivisions_vec);
    title(sprintf('\\alpha=%g', const_alpha));
    legend('error', '2nd order convergence');

    alpha_str = strrep(sprintf('%.2f', const_alpha), '0.', '.');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D1%s.fig', alpha_str));
    
    convergence_rate = (log(errors(end-1)) - log(errors(end))) / (log(subdivisions_vec(end)) - log(subdivisions_vec(end-1)));
    fprintf('Convergence rate: %d\n', convergence_rate);
else
    plot(alpha_vec, errors);
    xlabel('alpha');
    ylabel(sprintf('l^%g error', err_order));
    title('error vs. \alpha (linear-linear)');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D1_alpha_%d.fig', const_subdivs));
    
    figure

    loglog(1 - alpha_vec, errors);
    xlabel('1-alpha');
    ylabel(sprintf('l^%g error', err_order));
    title('error vs. 1-\alpha (log-log)');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D1_loglog_s_%d.fig', const_subdivs));
end

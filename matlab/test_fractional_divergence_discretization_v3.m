% Test the discretization of the fractional divergence -- version 3
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
change_alpha = false;

if change_alpha % change alpha
    s_vec = [logspace(-8, -1, 40) linspace(1e-1, 0.9999, 40)];
    alpha_vec = 1 - s_vec;
    n_samples = length(alpha_vec);
    const_subdivs = 16;
    subdivisions_vec = const_subdivs * ones(n_samples, 1);
else % Chnage number of subdivisions
    subdivisions_vec = [2 4 8 16];
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
    
        %% Compute values needed for discrete fractional divergence
        values_yz_from_0 = (1/50).*xs.^(-2).*((-2)+25.*xs.^2.*ys.*zs.*((-1).*ys+xs.*zs).*(( ...
      -3)+2.*ys+2.*xs.*zs)+2.*cos(5.*xs.*ys)+2.*cos(zs)+(-2).*cos(5.* ...
      xs.*ys+zs)+10.*xs.*ys.*(sin(5.*xs.*ys)+(-1).*sin(5.*xs.*ys+zs)));
        
        values_xz_from_0 = (1/100).*xs.^(-1).*ys.^(-2).*(1+25.*xs.^2.*ys.^2.*zs.*((-4).* ...
      ys.^3+xs.*zs)+(-1).*cos(10.*xs.*ys.*zs));
        
        values_xy_from_0 = (-1/5).*xs.^(-2)+(-1/2).*xs.^2.*ys+(-1/16).*xs.^2.*ys.^2+(1/4).* ...
      xs.^2.*ys.^2.*zs+(1/4).*xs.*ys.^4.*zs+(1/5).*xs.^(-2).*cos(5.* ...
      xs.^3.*ys)+(1/2).*5.^(-1/3).*ys.^(2/3).*gamma((1/3))+(1/2).*(-1) ...
      .^(2/3).*5.^(-1/3).*ys.^(2/3).*igamma((1/3),(sqrt(-1)*(-5)).* ...
      xs.^3.*ys)+(-1/2).*5.^(-1/3).*(sqrt(-1).*ys).^(2/3).*igamma((1/3), ...
      (sqrt(-1)*5).*xs.^3.*ys)+(1/2).*zs.^(-1).*sinint(2.*xs.*ys.*zs);
    end

    %% Compute values needed for continuous fractional divergence
    continuous_frac_divergence_integrated_to_nodes = (1/2).*xs.*ys.^(1+(-1).*alpha).*zs.*gamma(2+(-1).*alpha).^(-1).*( ...
  12.*((-24)+26.*alpha+(-9).*alpha.^2+alpha.^3).^(-1).*ys.^3+(-1).* ...
  zs+2.*zs.*hypergeom([(1/2),1],[(3/2),1+(-1/2).*alpha,(3/2)+(-1/2) ...
  .*alpha],(-25).*xs.^2.*ys.^2.*zs.^2)+(-1).*zs.*hypergeom([1,1],[2, ...
  1+(-1/2).*alpha,(3/2)+(-1/2).*alpha],(-25).*xs.^2.*ys.^2.*zs.^2))+ ...
  (1/4).*xs.*ys.*zs.^(1+(-1).*alpha).*gamma(3+(-1).*alpha).^(-1).*(( ...
  -8)+4.*alpha+xs.*ys.*zs+ys.^3.*zs+(-4).*((-2)+alpha).*hypergeom([( ...
  1/2),(1/2),1],[(3/2),(3/2),1+(-1/2).*alpha,(3/2)+(-1/2).*alpha],( ...
  -1).*xs.^2.*ys.^2.*zs.^2))+(-1/6).*xs.^(1+(-1).*alpha).*ys.*gamma( ...
  4+(-1).*alpha).^(-1).*(27.*xs.*zs.^2+(-9).*alpha.*xs.*zs.^2+(-12) ...
  .*xs.^2.*zs.^3+((-3)+alpha).*ys.*(3.*((-2)+alpha).*((-1)+cos(zs)) ...
  .*((-1)+hypergeom([1,1],[2,1+(-1/2).*alpha,(3/2)+(-1/2).*alpha],( ...
  -25/4).*xs.^2.*ys.^2))+10.*xs.*ys.*hypergeom([1,(3/2)],[(5/2),( ...
  3/2)+(-1/2).*alpha,2+(-1/2).*alpha],(-25/4).*xs.^2.*ys.^2).*sin( ...
  zs)));

    %% Fix the cases where x, y, or z is 0 and there could be a NaN
    if calculate_nonalpha_things
        ni = nx + 1;
        nj = ny + 1;
        nk = nz + 1;
        
        if x1 == 0.0
            % x == 0 (i == 1)
            for j = 1:nj     % for all y
                for k = 1:nk % for all z
                    node = node_index(1, j, k);
                    pos = node_coords(node, :);
                    x = pos(1);
                    y = pos(2);
                    z = pos(3);
                    
                    values_yz_from_0(node) = (-1/2).*y.^2.*((-1)+((-3)+2.*y).*z+cos(z));
                    values_xz_from_0(node) = 0.0;
                    values_xy_from_0(node) = 0.0;
                end
            end
        end
        
        if y1 == 0.0
            % y == 0 (j == 1)
            for i = 1:ni     % for all x
                for k = 1:nk % for all z
                    node = node_index(i, 1, k);
                    pos = node_coords(node, :);
                    x = pos(1);
                    y = pos(2);
                    z = pos(3);
                    
                    values_xz_from_0(node) = (1/4).*x.*(2+x).*z.^2;
                end
            end
        end
    
        if z1 == 0.0
            % z == 0 (k == 1)
            for i = 1:ni     % for all x
                for j = 1:nj % for all y
                    node = node_index(i, j, 1);
                    pos = node_coords(node, :);
                    x = pos(1);
                    y = pos(2);
                    z = pos(3);
                    
                    if x ~= 0.0 % Avoid introducing a NaN again
                        values_xy_from_0(node) = (-1/5).*x.^(-2)+x.*y+(-1/2).*x.^2.*y+(-1/16).*x.^2.*y.^2+(1/5).* ...
      x.^(-2).*cos(5.*x.^3.*y)+(1/2).*5.^(-1/3).*y.^(2/3).*gamma((1/3))+ ...
      (-1/2).*5.^(-1/3).*((sqrt(-1)*(-1)).*y).^(2/3).*igamma((1/3),( ...
      sqrt(-1)*(-5)).*x.^3.*y)+(-1/2).*5.^(-1/3).*(sqrt(-1).*y).^(2/3).* ...
      igamma((1/3),(sqrt(-1)*5).*x.^3.*y);
                    end
                end
            end
        end
    end

    %% Calculate discrete fractional divergence
    
    if calculate_nonalpha_things
        % Create values_on_faces
        F_integrated_to_nodes_on_faces = [values_yz_from_0; values_xz_from_0; values_xy_from_0];
        B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
        values_on_faces = B2 * F_integrated_to_nodes_on_faces;

        B3 = create_B3(nx, ny, nz, node_index, volume_index);
    end
    
    % Create I2
    M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
    I2 = B2 * M2_1ps .* (1 ./ face_areas)';
    
    % Create I3
    M3_1ps = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, 1+s);
    I3 = B3 * M3_1ps .* (1 ./ volumes)';

    discrete_frac_divergence_unintegrated = (1./ volumes) .* (I3 * (D2 * I2^-1) * values_on_faces);
    
    %% Calculate continuous fractional divergence
    continuous_frac_divergence_unintegrated = (1./ volumes) .* (B3 * continuous_frac_divergence_integrated_to_nodes);
    
    %% Calculate error
    diff = discrete_frac_divergence_unintegrated - continuous_frac_divergence_unintegrated;
%     errors(iter) = mean(abs(diff).^err_order)^(1/err_order);
    errors(iter) = ((volumes' * abs(diff).^err_order) / sum(volumes))^(1/err_order);

    T = toc;
    fprintf("%d,%d, %f\n", iter, nx, T);
end


if ~change_alpha
    loglog(subdivisions_vec, errors);
    hold on
    loglog(subdivisions_vec, subdivisions_vec.^-2., '--', 'Color', '#77AC30');
    xlabel('subdivisions');
    ylabel(sprintf('l^%g error', err_order));
    set(gca, 'XTick', subdivisions_vec);
    title(sprintf('\\alpha=%g', const_alpha));
    legend('error', '2nd order convergence');

    alpha_str = strrep(sprintf('%.2f', const_alpha), '0.', '.');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D2%s.fig', alpha_str));

    convergence_rate = (log(errors(end-1)) - log(errors(end))) / (log(subdivisions_vec(end)) - log(subdivisions_vec(end-1)));
    fprintf('Convergence rate: %d\n', convergence_rate);
else
    plot(alpha_vec, errors);
    xlabel('alpha');
    ylabel(sprintf('l^%g error', err_order));
    title('error vs. \alpha (linear-linear)');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D2_alpha_%d.fig', const_subdivs));
    
    figure

    loglog(1 - alpha_vec, errors);
    xlabel('1-alpha');
    ylabel(sprintf('l^%g error', err_order));
    title('error vs. 1-\alpha (log-log)');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D2_loglog_s_%d.fig', const_subdivs));
end

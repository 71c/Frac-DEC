% Test the discretization of the fractional gradient - version 3
% Created July 4, 2023

% f[x_, y_, z_] = -8 x y^2 z + 4 (x - 1/2)^2 + (y - 1/2)^2 - 
%    3 Cos[20 (x - 1/2) (y - 1) z];
f = @(x, y, z) 4.*((-1/2)+x).^2+((-1/2)+y).^2+(-8).*x.*y.^2.*z+(-3).*cos(20.*(( ...
  -1/2)+x).*((-1)+y).*z);


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
    subdivisions_vec = [2 4 8 16 32];
    n_samples = length(subdivisions_vec);
    const_alpha = 0.5;
    alpha_vec = const_alpha * ones(n_samples, 1);
end

errors = zeros(n_samples, 1);

for iter = 1:n_samples
    tic;

    nx = subdivisions_vec(iter);
    ny = subdivisions_vec(iter);
    nz = subdivisions_vec(iter);

    alpha = alpha_vec(iter);
    s = 1 - alpha;

    [D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);

    xs = node_coords(:, 1);
    ys = node_coords(:, 2);
    zs = node_coords(:, 3);
    
    %% Calculate discrete fractional gradient
    values_on_nodes = f(xs, ys, zs);
    B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);
    M1_1ps = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1+s);
    I1 = (B1 * M1_1ps) .* (1 ./ edge_lengths)';
    discrete_frac_gradient = I1 * D0 * values_on_nodes;
    discrete_frac_gradient_unintegrated = (1 ./ edge_lengths) .* discrete_frac_gradient;

    %% Calculate continuous fractional gradient
    continuous_frac_gradient_integrated_to_nodes_xs = (-4).*((-2)+alpha).^(-1).*xs.^(2+(-1).*alpha).*gamma(2+(-1).* ...
  alpha).^(-1).*((-1)+(-2).*((-3)+alpha).^(-1).*xs+(-2).*ys.^2.*zs+ ...
  15.*((-3)+alpha).^(-1).*((-1)+ys).*zs.*((-20).*xs.*((-1)+ys).*zs.* ...
  cos(10.*((-1)+ys).*zs).*hypergeom([1],[2+(-1/2).*alpha,(5/2)+( ...
  -1/2).*alpha],(-100).*xs.^2.*((-1)+ys).^2.*zs.^2)+((-3)+alpha).* ...
  hypergeom([1],[(3/2)+(-1/2).*alpha,2+(-1/2).*alpha],(-100).* ...
  xs.^2.*((-1)+ys).^2.*zs.^2).*sin(10.*zs+(-10).*ys.*zs)));

    continuous_frac_gradient_integrated_to_nodes_ys = ((-3)+alpha).^(-1).*((-2)+alpha).^(-1).*ys.^(2+(-1).*alpha).* ...
  gamma(2+(-1).*alpha).^(-1).*((-3)+alpha+2.*ys+(-16).*xs.*ys.*zs+( ...
  -60).*((1/2)+(-1).*xs).*zs.*(10.*((-1)+2.*xs).*ys.*zs.*cos(10.*(1+ ...
  (-2).*xs).*zs).*hypergeom([1],[2+(-1/2).*alpha,(5/2)+(-1/2).* ...
  alpha],(-25).*(1+(-2).*xs).^2.*ys.^2.*zs.^2)+(-1).*((-3)+alpha).* ...
  hypergeom([1],[(3/2)+(-1/2).*alpha,2+(-1/2).*alpha],(-25).*(1+(-2) ...
  .*xs).^2.*ys.^2.*zs.^2).*sin(10.*(1+(-2).*xs).*zs)));

    continuous_frac_gradient_integrated_to_nodes_zs = (-4).*((-2)+alpha).^(-1).*zs.^(2+(-1).*alpha).*gamma(2+(-1).* ...
  alpha).^(-1).*((-2).*xs.*ys.^2+(-75).*((-3)+alpha).^(-1).*(1+(-2) ...
  .*xs).^2.*((-1)+ys).^2.*zs.*hypergeom([1],[2+(-1/2).*alpha,(5/2)+( ...
  -1/2).*alpha],(-25).*(1+(-2).*xs).^2.*((-1)+ys).^2.*zs.^2));

%     % Fill in continuous_frac_gradient_integrated_to_nodes
%     n_nodes = size(D0, 2);
%     continuous_frac_gradient_integrated_to_nodes = zeros(3 * n_nodes, 1);
%     for node = 1:n_nodes
%         continuous_frac_gradient_integrated_to_nodes(node_direction_index(node, 1)) = continuous_frac_gradient_integrated_to_nodes_xs(node);
%         continuous_frac_gradient_integrated_to_nodes(node_direction_index(node, 2)) = continuous_frac_gradient_integrated_to_nodes_ys(node);
%         continuous_frac_gradient_integrated_to_nodes(node_direction_index(node, 3)) = continuous_frac_gradient_integrated_to_nodes_zs(node);
%     end
    continuous_frac_gradient_integrated_to_nodes = [continuous_frac_gradient_integrated_to_nodes_xs; continuous_frac_gradient_integrated_to_nodes_ys; continuous_frac_gradient_integrated_to_nodes_zs];
    
    continuous_frac_gradient_integrated = B1 * continuous_frac_gradient_integrated_to_nodes;
    continuous_frac_gradient_unintegrated = (1 ./ edge_lengths) .* continuous_frac_gradient_integrated;
    
    %% Calculate error
    errors(iter) = mean(abs(discrete_frac_gradient_unintegrated - continuous_frac_gradient_unintegrated).^err_order)^(1 / err_order);

    T = toc;
    fprintf("%d,%d, %f\n", iter, nx, T);
end


if ~change_alpha
    figure
    loglog(subdivisions_vec, errors);
    hold on
    loglog(subdivisions_vec, 4 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
    xlim([subdivisions_vec(1) subdivisions_vec(end)]);
    xlabel('subdivisions');
    ylabel(sprintf('l^%g error', err_order));
    set(gca, 'XTick', subdivisions_vec);
    title(sprintf('\\alpha=%g', const_alpha));
    legend('error', '2nd order convergence reference line');
    
    alpha_str = strrep(sprintf('%.2f', const_alpha), '0.', '.');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D0%s.fig', alpha_str));
    
    convergence_rate = (log(errors(end-1)) - log(errors(end))) / (log(subdivisions_vec(end)) - log(subdivisions_vec(end-1)));
    fprintf('Convergence rate: %d\n', convergence_rate);
else
    plot(alpha_vec, errors);
    xlabel('alpha');
    ylabel(sprintf('l^%g error', err_order));
    title('error vs. \alpha (linear-linear)');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D0_alpha_%d.fig', const_subdivs));
    
    figure

    loglog(1 - alpha_vec, errors);
    xlabel('1-alpha');
    ylabel(sprintf('l^%g error', err_order));
    title('error vs. 1-\alpha (log-log)');
%     savefig(sprintf('../../../papers/paper_FDEC/plots/2023-07-05/D0_loglog_s_%d.fig', const_subdivs));
end

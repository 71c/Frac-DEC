% Test the discretization of the fractional gradient - version 2
% Created September 18, 2022

f = @(x, y, z) (x + y.^2) .* y .* z - x;

err_order = 2.0;

x1 = 0.0; x2 = 1.0;
y1 = 0.0; y2 = 1.0;
z1 = 0.0; z2 = 1.0;

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
    const_alpha = 0.9;
    alpha_vec = const_alpha * ones(n_samples, 1);
end

errors = zeros(n_samples, 1);

for iter = 1:n_samples
    nx = subdivisions_vec(iter);
    ny = subdivisions_vec(iter);
    nz = subdivisions_vec(iter);

    alpha = alpha_vec(iter);
    s = 1 - alpha;

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
    
    I1 = (B1 * M1_1ps) .* (1 ./ edge_lengths)';

    discrete_frac_gradient = I1 * D0 * values_on_nodes;

    n_edges = size(D0, 1);
    continuous_frac_gradient_integrated_to_nodes = zeros(3 * n_nodes, 1);

    tic;

    % Calculate the continuous version of M_1^1 diag(1/v_1) D_1^\alpha c^0
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                x = x1 + (i-1) * dx;
                y = y1 + (j-1) * dy;
                z = z1 + (k-1) * dz;

                for l = 1:3
                    row_number = node_direction_index(node_index(i, j, k), l);

                    if l == 1 % x direction
                        continuous_frac_gradient_integrated_to_nodes(row_number) = x.^(2+(-1).*alpha).*((-1)+y.*z).*gamma(3+(-1).*alpha).^(-1);
                    elseif l == 2 % y direction
                        continuous_frac_gradient_integrated_to_nodes(row_number) = y.^(2+(-1).*alpha).*(((-4)+alpha).*((-3)+alpha).*x+6.*y.^2).*z.*gamma(5+(-1).*alpha).^(-1);
                    elseif l == 3 % z direction
                        continuous_frac_gradient_integrated_to_nodes(row_number) = y.*(x+y.^2).*z.^(2+(-1).*alpha).*gamma(3+(-1).*alpha).^(-1);
                    end
                end
            end
        end
    end

    T = toc;
    fprintf("%d,%d, %f\n", iter, nx, T);



    continuous_frac_gradient_integrated = B1 * continuous_frac_gradient_integrated_to_nodes;

    discrete_frac_gradient_unintegrated = (1 ./ edge_lengths) .* discrete_frac_gradient;
    continuous_frac_gradient_unintegrated = (1 ./ edge_lengths) .* continuous_frac_gradient_integrated;
    
    errors(iter) = mean(abs(discrete_frac_gradient_unintegrated - continuous_frac_gradient_unintegrated).^err_order)^(1 / err_order);
end


if ~change_alpha
    loglog(subdivisions_vec, errors);
    hold on
    loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
    xlim([subdivisions_vec(1) subdivisions_vec(end)]);
    xlabel('subdivisions');
    ylabel(sprintf('l^%g error', err_order));
    set(gca, 'XTick', subdivisions_vec);
    title(sprintf('Error of discrete fractional gradient for alpha=%g', const_alpha));
    legend('error', '2nd order convergence reference line');
    
    convergence_rate = (log(errors(end-1)) - log(errors(end))) / (log(subdivisions_vec(end)) - log(subdivisions_vec(end-1)));
    fprintf('Convergence rate: %d\n', convergence_rate);
else
    plot(alpha_vec, errors);
    xlabel('alpha');
    ylabel(sprintf('l^%g error', err_order));
    title(sprintf('Error of discrete fractional gradient for varying vaues of alpha at %d subdivisions', const_subdivs));
    
    figure

    loglog(1 - alpha_vec, errors);
    xlabel('1-alpha');
    ylabel(sprintf('l^%g error', err_order));
    title(sprintf('Error of discrete fractional gradient for varying vaues of 1-alpha at %d subdivisions (log scale)', const_subdivs));
end

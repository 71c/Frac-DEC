% Test the discretization of the fractional curl -- version 2
% Created September 18, 2022

% Vector field we are using:
% F(x, y, z) = [x y + (x - 1) z^2; x z - y^3; x y^2 z + z^4]

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
    
    n_edges = size(D0, 1);
    values_on_edges = zeros(n_edges, 1);
    
    n_nodes = size(D0, 2);
    continuous_frac_curl_integrated_to_nodes = zeros(3 * n_nodes, 1);
    
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                x = x1 + (i-1) * dx;
                y = y1 + (j-1) * dy;
                z = z1 + (k-1) * dz;
    
                xa = x;
                ya = y;
                za = z;
                
                xb = xa + dx;
                yb = ya + dy;
                zb = za + dz;

                % x edge
                if i ~= ni
                    edge = x_edge_index(i, j, k);
                    values_on_edges(edge) = (-1/2).*(xa+(-1).*xb).*((xa+xb).*y+((-2)+xa+xb).*z.^2);
                end
    
                % y edge
                if j ~= nj
                    edge = y_edge_index(i, j, k);
                    values_on_edges(edge) = (1/4).*ya.^4+(-1/4).*yb.^4+x.*((-1).*ya+yb).*z;
                end
    
                % z edge
                if k ~= nk
                    edge = z_edge_index(i, j, k);
                    values_on_edges(edge) = (1/10).*((-2).*za.^5+2.*zb.^5+5.*x.*y.^2.*((-1).*za.^2+zb.^2));
                end

                for l = 1:3
                    row_number = node_direction_index(node_index(i, j, k), l);

                    if l == 1 % YZ
                        if j == 1 || k == 1 % y == 0 or z == 0
                            continuous_frac_curl_integrated_to_nodes(row_number) = 0.0;
                        else
                            continuous_frac_curl_integrated_to_nodes(row_number) = x.*y.*z.^2.*(y.^(2+(-1).*alpha)+((-3)+alpha).*z.^((-1).*alpha)).*gamma(4+(-1).*alpha).^(-1);
                        end
                    elseif l == 2 % XZ
                        if i == 1 || k == 1 % x == 0 or z == 0
                            continuous_frac_curl_integrated_to_nodes(row_number) = 0.0;
                        else
                            continuous_frac_curl_integrated_to_nodes(row_number) = (-1/2).*x.^(2+(-1).*alpha).*y.^2.*z.^2.*gamma(3+(-1).*alpha).^(-1)+((-2)+x).*x.*z.^(3+(-1).*alpha).*gamma(4+(-1).*alpha).^(-1);
                        end
                    elseif l == 3 % XY
                        if i == 1 || j == 1 % x == 0 or y == 0
                            continuous_frac_curl_integrated_to_nodes(row_number) = 0.0;
                        else
                            continuous_frac_curl_integrated_to_nodes(row_number) = (-1/2).*x.^2.*y.^(1+(-1).*alpha).*(y+(-2).*(x.^(-1).*y).^alpha.*z).*gamma(3+(-1).*alpha).^(-1);
                        end
                    end
                end

            end
        end
    end
    
    B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);
    B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
    
    M1_1ps = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1+s);
    M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
    
    I1 = B1 * M1_1ps .* (1 ./ edge_lengths)';
    I2 = B2 * M2_1ps .* (1 ./ face_areas)';

    discrete_frac_curl_unintegrated = (1 ./ face_areas) .* (I2 * D1 * I1^-1 * values_on_edges);
    continuous_frac_curl_unintegrated = (1 ./ face_areas) .* (B2 * continuous_frac_curl_integrated_to_nodes);

    errors(iter) = mean(abs(discrete_frac_curl_unintegrated - continuous_frac_curl_unintegrated).^err_order)^(1/err_order);
end

if ~change_alpha
    loglog(subdivisions_vec, errors);
    hold on
    loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
    ylim([1e-5 1e-1]);
    xlim([subdivisions_vec(1) subdivisions_vec(end)]);
    xlabel('subdivisions');
    ylabel(sprintf('l^%g error', err_order));
    set(gca, 'XTick', subdivisions_vec);
    title(sprintf('Error of discrete fractional curl for alpha=%g', const_alpha));
    legend('error', '2nd order convergence');
    
    convergence_rate = (log(errors(end-1)) - log(errors(end))) / (log(subdivisions_vec(end)) - log(subdivisions_vec(end-1)));
    fprintf('Convergence rate: %d\n', convergence_rate);
else
    plot(alpha_vec, errors);
    xlabel('alpha');
    ylabel(sprintf('l^%g error', err_order));
    title(sprintf('Error of discrete fractional curl for varying vaues of alpha at %d subdivisions', const_subdivs));
    
    figure

    loglog(1 - alpha_vec, errors);
    xlabel('1-alpha');
    ylabel(sprintf('l^%g error', err_order));
    title(sprintf('Error of discrete fractional curl for varying vaues of 1-alpha at %d subdivisions (log scale)', const_subdivs));
end


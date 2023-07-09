% Tests that D_{p+1}^\alpha D_p^\alpha = 0.
% Created July 5, 2023

err_order = 2.0;
x1 = 0.0; x2 = 1.0;
y1 = 0.0; y2 = 1.0;
z1 = 0.0; z2 = 1.0;

f = @(x, y, z) 4.*((-1/2)+x).^2+((-1/2)+y).^2+(-8).*x.*y.^2.*z+(-3).*cos(20.*(( ...
  -1/2)+x).*((-1)+y).*z);


alpha_vec = [0.2 0.4 0.6 0.8];
subdivisions_vec = [2 3 4 5 6 7 8 9 10];

errors_D1D0_method1 = zeros(length(subdivisions_vec), length(alpha_vec));
errors_D1D0_method2 = zeros(length(subdivisions_vec), length(alpha_vec));

errors_D2D1_method1 = zeros(length(subdivisions_vec), length(alpha_vec));
errors_D2D1_method2 = zeros(length(subdivisions_vec), length(alpha_vec));

for subdiv_num = 1:length(subdivisions_vec)
    subdivisions = subdivisions_vec(subdiv_num);

    nx = subdivisions;
    ny = subdivisions;
    nz = subdivisions;
    
    % Create mesh
    [D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);
    
    % Create B_p matrices
    B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);
    B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
    B3 = create_B3(nx, ny, nz, node_index, volume_index);

    xs = node_coords(:, 1);
    ys = node_coords(:, 2);
    zs = node_coords(:, 3);
    
    %% Create c0 and c1

    % Create c0
    c0 = f(xs, ys, zs);
    
    values_x_from_0 = (-3).*xs.*((-1)+ys).*ys+(-3/2).*xs.^2.*zs+xs.^3.*zs.^2+(1/5).*cos(zs)+(-1/5).*cos(5.*xs.*ys+zs);
    values_y_from_0 = (-1/4).*ys.^4+xs.*ys.*zs+(1/10).*xs.^(-1).*sin(10.*xs.*ys.*zs);
    values_z_from_0 = (1/4).*(zs.*((-1).*xs.*(4+ys)+2.*ys.*(xs+ys.^2).*zs+8.*sin(5.*xs.^3.*ys))+2.*xs.^(-1).*ys.^(-1).*sin(2.*xs.*ys.*zs));
    
    % Fix the cases where x or y is 0 and there could be a NaN
    ni = nx + 1;
    nj = ny + 1;
    nk = nz + 1;
    % x == 0 (i == 1)
    for j = 1:nj     % for all y
        for k = 1:nk % for all z
            node = node_index(1, j, k);
            pos = node_coords(node, :);
            x = pos(1);
            y = pos(2);
            z = pos(3);
            values_y_from_0(node) = (-1/4).*y.^4+y.*z;
            values_z_from_0(node) = (1/2).*z.*(2+y.^3.*z);
        end
    end
    % y == 0 (j == 1)
    for i = 1:ni     % for all x
        for k = 1:nk % for all z
            node = node_index(i, 1, k);
            pos = node_coords(node, :);
            x = pos(1);
            y = pos(2);
            z = pos(3);
            values_z_from_0(node) = z+(-1).*x.*z;
        end
    end
    
    % Create c1
    F_integrated_to_nodes_on_edges = [values_x_from_0; values_y_from_0; values_z_from_0];
    c1 = B1 * F_integrated_to_nodes_on_edges;

    for alpha_num = 1:length(alpha_vec)
        alpha = alpha_vec(alpha_num);
        s = 1 - alpha;

        fprintf('%d %g\n', subdivisions, alpha);

        %% Create D_p^\alpha matrices
        % Create I1
        M1_1ps = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1+s);
        I1 = (B1 * M1_1ps) .* (1 ./ edge_lengths)';
        
        % Create I2
        M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
        I2 = B2 * M2_1ps .* (1 ./ face_areas)';
        
        % Create I3
        M3_1ps = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, 1+s);
        I3 = B3 * M3_1ps .* (1 ./ volumes)';
        
        % Create D_p^\alpha
        D0_alpha = I1 * D0;
        D1_alpha = I2 * D1 * I1^-1;
        D2_alpha = I3 * D2 * I2^-1;
        
        %% Test D D = 0
        c2_method1 = (D1_alpha * D0_alpha) * c0;
        c3_method1 = (D2_alpha * D1_alpha) * c1;

        c2_method2 = D1_alpha * (D0_alpha * c0);
        c3_method2 = D2_alpha * (D1_alpha * c1);

        errors_D1D0_method1(subdiv_num, alpha_num) = ((face_areas' * abs((1 ./ face_areas) .* c2_method1).^err_order) / sum(face_areas)) ^ (1/err_order);
        errors_D2D1_method1(subdiv_num, alpha_num) = ((volumes' * abs((1 ./ volumes) .* c3_method1).^err_order) / sum(volumes)) ^ (1/err_order);

        errors_D1D0_method2(subdiv_num, alpha_num) = ((face_areas' * abs((1 ./ face_areas) .* c2_method2).^err_order) / sum(face_areas)) ^ (1/err_order);
        errors_D2D1_method2(subdiv_num, alpha_num) = ((volumes' * abs((1 ./ volumes) .* c3_method2).^err_order) / sum(volumes)) ^ (1/err_order);
    end
end

matrices = {errors_D1D0_method1, errors_D1D0_method2, errors_D2D1_method1, errors_D2D1_method2};
titles = {'error D_1^\alpha D_0^\alpha method 1', 'error D_1^\alpha D_0^\alpha method 2', 'error D_2^\alpha D_1^\alpha method 1', 'error D_2^\alpha D_1^\alpha method 2'};
for i = 1:length(matrices)
    title_i = titles{i};
    errors = matrices{i};

    figure
    hold on
    for j = 1:size(errors, 2)
        plots(j) = plot(subdivisions_vec, errors(:, j));
        names{j} = sprintf('\\alpha=%g', alpha_vec(j));
    end
    legend(plots, names)
    xlabel('subdivisions')
    ylabel('error')
    if weighted
        title(sprintf('%s (weighted)', title_i));
    else
        title(sprintf('%s (unweighted)', title_i));
    end
end


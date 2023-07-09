% Test the discretization of the fractional divergence -- version 2
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
    const_subdivs = 8;
    subdivisions_vec = const_subdivs * ones(n_samples, 1);
else % Chnage number of subdivisions
    subdivisions_vec = [2 4 8 16];
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
    
    n_faces = size(D1, 1);
    values_on_faces = zeros(n_faces, 1);
    
    n_nodes = size(D0, 2);
    continuous_frac_divergence_integrated_to_nodes = zeros(n_nodes, 1);
    
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
    
                % yz face
                if j ~= nj && k ~= nk
                    face = yz_face_index(i, j, k);
                    values_on_faces(face) = x.*((-1/2).*ya.^2+(1/2).*yb.^2).*((-1).*za+zb)+((-1)+x).*((-1).* ...
      ya+yb).*((-1/3).*za.^3+(1/3).*zb.^3);
                end
    
                % xz face
                if i ~= ni && k ~= nk
                    face = xz_face_index(i, j, k);
                    values_on_faces(face) = (-1).*((-1).*xa+xb).*y.^3.*((-1).*za+zb)+((-1/2).*xa.^2+(1/2).* ...
      xb.^2).*((-1/2).*za.^2+(1/2).*zb.^2);
                end
    
                % xy face
                if i ~= ni && j ~= nj
                    face = xy_face_index(i, j, k);
                    values_on_faces(face) = ((-1/2).*xa.^2+(1/2).*xb.^2).*((-1/3).*ya.^3+(1/3).*yb.^3).*z+(( ...
      -1).*xa+xb).*((-1).*ya+yb).*z.^4;
                end
    
                % Fill in fractional divergence integrated to node
                node = node_index(i, j, k);
                if i == 0 || j == 0 || k == 0
                    continuous_frac_divergence_integrated_to_nodes(node) = 0.0;
                else
                    continuous_frac_divergence_integrated_to_nodes(node) = (1/6).*(x.*y.*z).^(1+(-1).*alpha).*(36.*(4.*(x.*y).^alpha.*z.^4+(( ...
  -5)+alpha).*y.^3.*(x.*z).^alpha)+(-1).*((-5)+alpha).*((-4)+alpha) ...
  .*((-3)+alpha).*x.*y.^alpha.*(x.^alpha.*y.^2.*z+z.^alpha.*(3.*y+ ...
  2.*z.^2))).*gamma(6+(-1).*alpha).^(-1);
                end
                
            end
        end
    end
    
    B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
    B3 = create_B3(nx, ny, nz, node_index, volume_index);
    
    M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
    M3_1ps = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, 1+s);

    I2 = B2 * M2_1ps .* (1 ./ face_areas)';
    I3 = B3 * M3_1ps .* (1 ./ volumes)';

    discrete_frac_divergence_unintegrated = (1./ volumes) .* (I3 * (D2 * I2^-1) * values_on_faces);
    continuous_frac_divergence_unintegrated = (1./ volumes) .* (B3 * continuous_frac_divergence_integrated_to_nodes);
    errors(iter) = mean(abs(discrete_frac_divergence_unintegrated - continuous_frac_divergence_unintegrated).^err_order)^(1/err_order);
end


if ~change_alpha
    loglog(subdivisions_vec, errors);
    hold on
    loglog(subdivisions_vec, subdivisions_vec.^-2., '--', 'Color', '#77AC30');
    xlabel('subdivisions');
    ylabel(sprintf('l^%g error', err_order));
    set(gca, 'XTick', subdivisions_vec);
    title(sprintf('Error of discrete fractional divergence for alpha=%g', const_alpha));
    legend('error', '2nd order convergence');

    convergence_rate = (log(errors(end-1)) - log(errors(end))) / (log(subdivisions_vec(end)) - log(subdivisions_vec(end-1)));
    fprintf('Convergence rate: %d\n', convergence_rate);
else
    plot(alpha_vec, errors);
    xlabel('alpha');
    ylabel(sprintf('l^%g error', err_order));
    title(sprintf('Error of discrete fractional divergence for varying vaues of alpha at %d subdivisions', const_subdivs));
    
    figure

    loglog(1 - alpha_vec, errors);
    xlabel('1-alpha');
    ylabel(sprintf('l^%g error', err_order));
    title(sprintf('Error of discrete fractional divergence for varying vaues of 1-alpha at %d subdivisions (log scale)', const_subdivs));
end

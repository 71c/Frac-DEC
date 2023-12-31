% Test the discretization of the fractional curl
% Created July 1, 2022

% Vector field we are using:
% F(x, y, z) = [x y + (x - 1) z^2; x z - y^3; x y^2 z + z^4]

x1 = 0.0; x2 = 2.0;
y1 = 1.0; y2 = 3.0;
z1 = 0.5; z2 = 5.0;

% fractional order of fractional curl
alpha = 0.9;

% s = n-alpha; order of fractional integration
s = 1 - alpha;

subdivisions_vec = [1 2 4 8 16];
n_samples = length(subdivisions_vec);
errors_1 = zeros(n_samples, 1);
errors_2 = zeros(n_samples, 1);
errors_3 = zeros(n_samples, 1);

for iter = 1:n_samples
    nx = subdivisions_vec(iter);
    ny = subdivisions_vec(iter);
    nz = subdivisions_vec(iter);

    [D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);
    
    ni = nx + 1;
    nj = ny + 1;
    nk = nz + 1;
    
    dx = (x2 - x1) / nx;
    dy = (y2 - y1) / ny;
    dz = (z2 - z1) / nz;
    
    n_edges = size(D0, 1);
    values_on_edges = zeros(n_edges, 1);
    
    n_faces = size(D1, 1);
    continuous_frac_curl = zeros(n_faces, 1);
    
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


                % yz face
                if j ~= nj && k ~= nk
                    face = yz_face_index(i, j, k);
                    x = face_coords(face, 1);
                    y = face_coords(face, 2);
                    z = face_coords(face, 3);
                    continuous_frac_curl(face) = (-1).*x.*(2.*((-2)+alpha).^(-1).*(y+(-1).*y1).^(1+(-1).*alpha).*( ...
  y+y1+(-1).*alpha.*y1).*z+(z+(-1).*z1).^(1+(-1).*alpha)).*gamma(2+( ...
  -1).*alpha).^(-1);
                end
    
                % xz face
                if i ~= ni && k ~= nk
                    face = xz_face_index(i, j, k);
                    x = face_coords(face, 1);
                    y = face_coords(face, 2);
                    z = face_coords(face, 3);
                    continuous_frac_curl(face) = (-1).*((x+(-1).*x1).^(1+(-1).*alpha).*y.^2.*z+2.*((-2)+alpha).^( ...
  -1).*((-1)+x).*(z+(-1).*z1).^(1+(-1).*alpha).*(z+z1+(-1).*alpha.* ...
  z1)).*gamma(2+(-1).*alpha).^(-1);
                end
    
                % xy face
                if i ~= ni && j ~= nj
                    face = xy_face_index(i, j, k);
                    x = face_coords(face, 1);
                    y = face_coords(face, 2);
                    z = face_coords(face, 3);
                    continuous_frac_curl(face) = ((-1).*x.*(y+(-1).*y1).^(1+(-1).*alpha)+(x+(-1).*x1).^(1+(-1).* ...
  alpha).*z).*gamma(2+(-1).*alpha).^(-1);
                end

            end
        end
    end
    
    B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);
    B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
    
    M1_1ps = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1+s);
    M1_1ms = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1-s);
    
    M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
    M2_1ms = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1-s);
    
    I1 = (1 ./ edge_lengths) .* B1 * M1_1ps .* (1 ./ edge_lengths)';
    I2 = (1 ./ face_areas) .* B2 * M2_1ps .* (1 ./ face_areas)';

    discrete_frac_curl_v1 = I2 * D1 * B1 * M1_1ms               * diag(1 ./ edge_lengths) * values_on_edges;
    discrete_frac_curl_v2 = I2 * D1 * I1^-1                     * diag(1 ./ edge_lengths) * values_on_edges;
    discrete_frac_curl_v3 = (B2 * M2_1ms)^-1 * D1 * B1 * M1_1ms * diag(1 ./ edge_lengths) * values_on_edges;
    errors_1(iter) = mean(abs(discrete_frac_curl_v1 - continuous_frac_curl));
    errors_2(iter) = mean(abs(discrete_frac_curl_v2 - continuous_frac_curl));
    errors_3(iter) = mean(abs(discrete_frac_curl_v3 - continuous_frac_curl));
end

% https://www.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html

% Compare methods
% loglog(subdivisions_vec, errors_1, subdivisions_vec, errors_2, subdivisions_vec, errors_3);
% hold on;
% loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-1., '--', 'Color', '#7E2F8E');
% loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
% ylim([1e-4 1e1]);
% legend('error v1', 'error v2', 'error v3', '1 order convergence', '2 order convergence');

% Look at just version 2
loglog(subdivisions_vec, errors_2);
hold on;
loglog(subdivisions_vec, 0.1 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
legend('error', '2nd order convergence');


xlabel('subdivisions');
ylabel('mean absolute error');
title(sprintf('Error of discrete fractional curl for alpha=%g', alpha));

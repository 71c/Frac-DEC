% Test the discretization of the fractional divergence
% Created June 30, 2022

% Vector field we are using:
% F(x, y, z) = [x y + (x - 1) z^2; x z - y^3; x y^2 z + z^4]

% x1 = 0.0; x2 = 2.0;
% y1 = 1.0; y2 = 3.0;
% z1 = 0.5; z2 = 5.0;

x1 = 0.0; x2 = 1.0;
y1 = 0.01; y2 = 1.0;
z1 = 0.01; z2 = 1.0;

% fractional order of fractional divergence
alpha = 0.8;

% s = n-alpha; order of fractional integration
s = 1 - alpha;

subdivisions_vec = [1 2 4 8 12 16];
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
    
    n_faces = size(D1, 1);
    values_on_faces = zeros(n_faces, 1);
    
    n_volumes = size(D2, 1);
    continuous_frac_divergence = zeros(n_volumes, 1);
    
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
    
                % volume
                if i ~= ni && j ~= nj && k ~= nk
                    volume = volume_index(i, j, k);
                    x = volume_coords(volume, 1);
                    y = volume_coords(volume, 2);
                    z = volume_coords(volume, 3);
                    continuous_frac_divergence(volume) = 3.*((-3)+alpha).^(-1).*((-2)+alpha).^(-1).*((-1)+alpha).^(-1).*(y+(-1).*y1).^( ...
  1+(-1).*alpha).*(2.*y.^2+(-2).*((-1)+alpha).*y.*y1+((-2)+alpha).*((-1)+alpha).* ...
  y1.^2).*gamma(1+(-1).*alpha).^(-1)+(-1).*((-1)+alpha).^(-1).*(x+(-1).*x1) ...
  .^(1+(-1).*alpha).*(y+z.^2).*gamma(1+(-1).*alpha).^(-1)+(z+(-1).*z1).^(1+( ...
  -1).*alpha).*((1+(-1).*alpha).^(-1).*x.*y.^2+4.*((-4)+alpha).^(-1).*((-3)+alpha) ...
  .^(-1).*((-2)+alpha).^(-1).*((-1)+alpha).^(-1).*(6.*z.^3+(-6).*((-1)+alpha).* ...
  z.^2.*z1+3.*((-2)+alpha).*((-1)+alpha).*z.*z1.^2+(-1).*((-3)+alpha).*((-2)+alpha) ...
  .*((-1)+alpha).*z1.^3)).*gamma(1+(-1).*alpha).^(-1);
                end
                
            end
        end
    end
    
    B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
    B3 = create_B3(nx, ny, nz, node_index, volume_index);
    
    M2_1ps = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1+s);
    M2_1ms = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1-s);
    M3_1ps = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, 1+s);
    M3_1ms = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, 1-s);
    
    I3 = (1 ./ volumes) .* B3 * M3_1ps .* (1 ./ volumes)';
    I2 = (1 ./ face_areas) .* B2 * M2_1ps .* (1 ./ face_areas)';

    discrete_frac_divergence_v1 = I3 * D2 * B2 * (M2_1ms .* (1 ./ face_areas)') * values_on_faces;
    discrete_frac_divergence_v2 = I3 * D2 * (I2^-1 .* (1 ./ face_areas)') * values_on_faces;
    discrete_frac_divergence_v3 = (B3 * M3_1ms)^-1 * D2 * B2 * (M2_1ms .* (1 ./ face_areas)') * values_on_faces;
    errors_1(iter) = mean(abs(discrete_frac_divergence_v1 - continuous_frac_divergence));
    errors_2(iter) = mean(abs(discrete_frac_divergence_v2 - continuous_frac_divergence));
    errors_3(iter) = mean(abs(discrete_frac_divergence_v3 - continuous_frac_divergence));
end

% Compare methods
loglog(subdivisions_vec, errors_1, subdivisions_vec, errors_2, subdivisions_vec, errors_3);
hold on
loglog(subdivisions_vec, 10 * subdivisions_vec.^-1., '--', 'Color', '#7E2F8E');
loglog(subdivisions_vec, 10 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
ylim([1e-2 1e3]);
legend('error v1', 'error v2', 'error v3', '1 order convergence', '2 order convergence');

% Look at just version 2
% loglog(subdivisions_vec, errors_2);
% hold on
% loglog(subdivisions_vec, 10 * subdivisions_vec.^-2., '--', 'Color', '#77AC30');
% legend('error', '2 order convergence');
% 
% xlabel('subdivisions');
% ylabel('mean absolute error');
% title(sprintf('Error of discrete fractional divergence for alpha=%g', alpha));

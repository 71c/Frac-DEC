function [M2] = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, s)
% Returns the 3D M2 matrix of order s.

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;

n_yz_faces = ni * (nj - 1) * (nk - 1);
n_xz_faces = (ni - 1) * nj * (nk - 1);
n_xy_faces = (ni - 1) * (nj - 1) * nk;
n_faces = n_yz_faces + n_xz_faces + n_xy_faces;

dx = (x2 - x1) / nx;
dy = (y2 - y1) / ny;
dz = (z2 - z1) / nz;

nzeros = ni * nj * nk * (ny * nz + nx * nz + nx * ny) / 4;

rows = zeros(nzeros, 1);
cols = zeros(nzeros, 1);
vals = zeros(nzeros, 1);
index = 1;
function set_entry(row, col, val)
    rows(index) = row;
    cols(index) = col;
    vals(index) = val;
    index = index + 1;
end

%%%%%% Outer loop: row of matrix to fill in %%%%%%%
for i = 1:ni         % x
    for j = 1:nj     % y
        for k = 1:nk % z
            for l = 1:3
                % row of matrix
                row_number = node_direction_index(node_index(i, j, k), l);

                x = x1 + (i-1) * dx;
                y = y1 + (j-1) * dy;
                z = z1 + (k-1) * dz;
                
                %%%%%%% Inner loop: column of matrix to fill in %%%%%%%
                if l == 1
                    % yz
                    for j2 = 1:j-1     % y
                        for k2 = 1:k-1 % z
                            ya = y1 + (j2-1) * dy;
                            za = z1 + (k2-1) * dz;
                            yb = ya + dy;
                            zb = za + dz;
                            
                            face = yz_face_index(i, j2, k2);
                            set_entry(row_number, face, ((y - ya)^s - (y - yb)^s) * ((z - za)^s - (z - zb)^s) / gamma(1 + s)^2);
                        end
                    end
                elseif l == 2
                    % xz
                    for i2 = 1:i-1     % x
                        for k2 = 1:k-1 % z
                            xa = x1 + (i2-1) * dx;
                            za = z1 + (k2-1) * dz;
                            xb = xa + dx;
                            zb = za + dz;

                            face = xz_face_index(i2, j, k2);
                            set_entry(row_number, face, ((x - xa)^s - (x - xb)^s) * ((z - za)^s - (z - zb)^s) / gamma(1 + s)^2);
                        end
                    end
                elseif l == 3
                    % xy
                    for i2 = 1:i-1     % x
                        for j2 = 1:j-1 % y
                            xa = x1 + (i2-1) * dx;
                            ya = y1 + (j2-1) * dy;
                            xb = xa + dx;
                            yb = ya + dy;

                            face = xy_face_index(i2, j2, k);
                            set_entry(row_number, face, ((x - xa)^s - (x - xb)^s) * ((y - ya)^s - (y - yb)^s) / gamma(1 + s)^2);
                        end
                    end

                end

            end
        end
    end
end

M2 = sparse(rows, cols, vals, 3*n_nodes, n_faces, nzeros);

end

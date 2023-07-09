function [M1] = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, s)
% Returns the 3D M1 matrix of order s.

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;

n_x_edges = (ni - 1) * nj * nk;
n_y_edges = ni * (nj - 1) * nk;
n_z_edges = ni * nj * (nk - 1);
n_edges = n_x_edges + n_y_edges + n_z_edges;

dx = (x2 - x1) / nx;
dy = (y2 - y1) / ny;
dz = (z2 - z1) / nz;

nzeros = ni * nj * nk * (nx + ny + nz) / 2;

rows = zeros(nzeros, 1);
cols = zeros(nzeros, 1);
vals = zeros(nzeros, 1, 'like', 0.0);
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
                    for i2 = 1:i-1 % x
                        xa = x1 + (i2-1) * dx;
                        xb = xa + dx;
                        
                        edge = x_edge_index(i2, j, k);
                        set_entry(row_number, edge, ((x - xa)^s - (x - xb)^s) / gamma(1 + s));
                    end
                elseif l == 2
                    for j2 = 1:j-1 % y
                        ya = y1 + (j2-1) * dy;
                        yb = ya + dy;

                        edge = y_edge_index(i, j2, k);
                        set_entry(row_number, edge, ((y - ya)^s - (y - yb)^s) / gamma(1 + s));
                    end
                elseif l == 3
                    for k2 = 1:k-1 % z
                        za = z1 + (k2-1) * dz;
                        zb = za + dz;

                        edge = z_edge_index(i, j, k2);
                        set_entry(row_number, edge, ((z - za)^s - (z - zb)^s) / gamma(1 + s));
                    end
                end
                
            end
        end
    end
end

M1 = sparse(rows, cols, vals, 3*n_nodes, n_edges, nzeros);

end

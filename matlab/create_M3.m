function [M3] = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, s)
% Returns the 3D M3 matrix of order s.

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;
n_volumes = (ni - 1) * (nj - 1) * (nk - 1);

dx = (x2 - x1) / nx;
dy = (y2 - y1) / ny;
dz = (z2 - z1) / nz;

nzeros = ni * nj * nk * nx * ny * nz / 8;

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
            x = x1 + (i-1) * dx;
            y = y1 + (j-1) * dy;
            z = z1 + (k-1) * dz;

            node = node_index(i, j, k);

            %%%%%% Inner loop: column of matrix to fill in %%%%%%%
            for i2 = 1:i-1         % x
                for j2 = 1:j-1     % y
                    for k2 = 1:k-1 % z
                        xa = x1 + (i2-1) * dx;
                        ya = y1 + (j2-1) * dy;
                        za = z1 + (k2-1) * dz;
                        
                        xb = xa + dx;
                        yb = ya + dy;
                        zb = za + dz;
                        
                        volume = volume_index(i2, j2, k2);
                        set_entry(node, volume, ((x - xa)^s - (x - xb)^s) * ((y - ya)^s - (y - yb)^s) * ((z - za)^s - (z - zb)^s) / gamma(1 + s)^3);
                    end
                end
            end


        end
    end
end

M3 = sparse(rows, cols, vals, n_nodes, n_volumes, nzeros);

end



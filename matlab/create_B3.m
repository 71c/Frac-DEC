function [B3] = create_B3(nx, ny, nz, node_index, volume_index)
% Returns the 3D B3 matrix

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;
n_volumes = (ni - 1) * (nj - 1) * (nk - 1);

nzeros = 8 * n_volumes;

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

for i = 1:ni-1         % x
    for j = 1:nj-1     % y
        for k = 1:nk-1 % z
            volume = volume_index(i,   j,   k);

            node_aaa = node_index(i,   j,   k);
            node_aab = node_index(i,   j,   k+1);
            node_aba = node_index(i,   j+1, k);
            node_abb = node_index(i,   j+1, k+1);
            node_baa = node_index(i+1, j,   k);
            node_bab = node_index(i+1, j,   k+1);
            node_bba = node_index(i+1, j+1, k);
            node_bbb = node_index(i+1, j+1, k+1);

            set_entry(volume, node_aaa, -1);
            set_entry(volume, node_aab, 1);
            set_entry(volume, node_aba, 1);
            set_entry(volume, node_abb, -1);
            set_entry(volume, node_baa, 1);
            set_entry(volume, node_bab, -1);
            set_entry(volume, node_bba, -1);
            set_entry(volume, node_bbb, 1);
        end
    end
end

B3 = sparse(rows, cols, vals, n_volumes, n_nodes, nzeros);

end
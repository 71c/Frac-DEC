function [B1] = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index)
% Returns the 3D B1 matrix.

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;

n_x_edges = (ni - 1) * nj * nk;
n_y_edges = ni * (nj - 1) * nk;
n_z_edges = ni * nj * (nk - 1);
n_edges = n_x_edges + n_y_edges + n_z_edges;

nzeros = 2 * n_edges; % number of nonzero elements

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

for i = 1:ni         % x
    for j = 1:nj     % y
        for k = 1:nk % z
            src_node = node_index(i, j, k);

            % x edge
            if i ~= ni
                edge = x_edge_index(i, j, k);
                dst_node = node_index(i+1, j, k);
                l = 1;
                set_entry(edge, node_direction_index(src_node, l), -1);
                set_entry(edge, node_direction_index(dst_node, l), 1);
            end

            % y edge
            if j ~= nj
                edge = y_edge_index(i, j, k);
                dst_node = node_index(i, j+1, k);
                l = 2;
                set_entry(edge, node_direction_index(src_node, l), -1);
                set_entry(edge, node_direction_index(dst_node, l), 1);
            end

            % z edge
            if k ~= nk
                edge = z_edge_index(i, j, k);
                dst_node = node_index(i, j, k+1);
                l = 3;
                set_entry(edge, node_direction_index(src_node, l), -1);
                set_entry(edge, node_direction_index(dst_node, l), 1);
            end
        end
    end
end

B1 = sparse(rows, cols, vals, n_edges, 3*n_nodes, nzeros);

end
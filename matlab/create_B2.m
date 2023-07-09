function [B2] = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index)
% Returns the 3D B2 matrix

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;

n_yz_faces = ni * (nj - 1) * (nk - 1);
n_xz_faces = (ni - 1) * nj * (nk - 1);
n_xy_faces = (ni - 1) * (nj - 1) * nk;
n_faces = n_yz_faces + n_xz_faces + n_xy_faces;

nzeros = 4 * n_faces;

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
            bottom_left_node = node_index(i, j, k);
    
            % yz face
            if j ~= nj && k ~= nk
                face = yz_face_index(i, j, k);
                top_left_node =     node_index(i, j, k+1);
                bottom_right_node = node_index(i, j+1, k);
                top_right_node =    node_index(i, j+1, k+1);
                l = 1;
                set_entry(face, node_direction_index(bottom_left_node, l), 1);
                set_entry(face, node_direction_index(top_right_node, l), 1);
                set_entry(face, node_direction_index(top_left_node, l), -1);
                set_entry(face, node_direction_index(bottom_right_node, l), -1);
            end
    
            % xz face
            if i ~= ni && k ~= nk
                face = xz_face_index(i, j, k);
                top_left_node =     node_index(i, j, k+1);
                bottom_right_node = node_index(i+1, j, k);
                top_right_node =    node_index(i+1, j, k+1);
                l = 2;
                set_entry(face, node_direction_index(bottom_left_node, l), 1);
                set_entry(face, node_direction_index(top_right_node, l), 1);
                set_entry(face, node_direction_index(top_left_node, l), -1);
                set_entry(face, node_direction_index(bottom_right_node, l), -1);
            end
    
            % xy face
            if i ~= ni && j ~= nj
                face = xy_face_index(i, j, k);
                top_left_node =     node_index(i, j+1, k);
                bottom_right_node = node_index(i+1, j, k);
                top_right_node =    node_index(i+1, j+1, k);
                l = 3;
                set_entry(face, node_direction_index(bottom_left_node, l), 1);
                set_entry(face, node_direction_index(top_right_node, l), 1);
                set_entry(face, node_direction_index(top_left_node, l), -1);
                set_entry(face, node_direction_index(bottom_right_node, l), -1);
            end
        end
    end
end

B2 = sparse(rows, cols, vals, n_faces, 3*n_nodes, nzeros);

end
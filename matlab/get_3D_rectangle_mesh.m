function [D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz)
% Builds a 3D rectangle mesh

ni = nx + 1;
nj = ny + 1;
nk = nz + 1;

n_nodes = ni * nj * nk;

n_x_edges = (ni - 1) * nj * nk;
n_y_edges = ni * (nj - 1) * nk;
n_z_edges = ni * nj * (nk - 1);
n_edges = n_x_edges + n_y_edges + n_z_edges;

n_yz_faces = ni * (nj - 1) * (nk - 1);
n_xz_faces = (ni - 1) * nj * (nk - 1);
n_xy_faces = (ni - 1) * (nj - 1) * nk;
n_faces = n_yz_faces + n_xz_faces + n_xy_faces;

n_volumes = (ni - 1) * (nj - 1) * (nk - 1);

node_index =   @(i, j, k) nk * (nj * (i-1) + j-1) + k;

% Ordering version 1:
% node_direction_index = @(node, l) 3 * (node - 1) + l;
% Ordering version 2:
node_direction_index = @(node, l) n_nodes * (l - 1) + node;

x_edge_index = @(i, j, k) nk * (nj * (i-1) + j-1) + k;
y_edge_index = @(i, j, k) n_x_edges + nk * ((nj-1) * (i-1) + j-1) + k;
z_edge_index = @(i, j, k) n_x_edges + n_y_edges + (nk-1) * (nj * (i-1) + j-1) + k;

yz_face_index = @(i, j, k) (nk-1) * ((nj-1) * (i-1) + j-1) + k;
xz_face_index = @(i, j, k) n_yz_faces + (nk-1) * (nj * (i-1) + j-1) + k;
xy_face_index = @(i, j, k) n_yz_faces + n_xz_faces + nk * ((nj-1) * (i-1) + j-1) + k;

volume_index = @(i, j, k) (nk-1) * ((nj-1) * (i-1) + j-1) + k;

dx = (x2 - x1) / nx;
dy = (y2 - y1) / ny;
dz = (z2 - z1) / nz;

node_coords = zeros(n_nodes, 3);
edge_coords = zeros(n_edges, 3);
face_coords = zeros(n_faces, 3);
volume_coords = zeros(n_volumes, 3);

edge_lengths = zeros(n_edges, 1);
edge_lengths(1:n_x_edges) = dx;
edge_lengths(n_x_edges+1:n_x_edges+n_y_edges) = dy;
edge_lengths(n_x_edges+n_y_edges+1:n_edges) = dz;

face_areas = zeros(n_faces, 1);
face_areas(1:n_yz_faces) = dy * dz;
face_areas(n_yz_faces+1:n_yz_faces+n_xz_faces) = dx * dz;
face_areas(n_yz_faces+n_xz_faces+1:n_faces) = dx * dy;

volumes = ones(n_volumes, 1) * (dx * dy * dz);

nzeros_D0 = 2 * n_edges;
nzeros_D1 = 4 * n_faces;
nzeros_D2 = 6 * n_volumes;

rows_D0 = zeros(nzeros_D0, 1);
cols_D0 = zeros(nzeros_D0, 1);
vals_D0 = zeros(nzeros_D0, 1);
index_D0 = 1;
function set_entry_D0(row, col, val)
    rows_D0(index_D0) = row;
    cols_D0(index_D0) = col;
    vals_D0(index_D0) = val;
    index_D0 = index_D0 + 1;
end

rows_D1 = zeros(nzeros_D1, 1);
cols_D1 = zeros(nzeros_D1, 1);
vals_D1 = zeros(nzeros_D1, 1);
index_D1 = 1;
function set_entry_D1(row, col, val)
    rows_D1(index_D1) = row;
    cols_D1(index_D1) = col;
    vals_D1(index_D1) = val;
    index_D1 = index_D1 + 1;
end

rows_D2 = zeros(nzeros_D2, 1);
cols_D2 = zeros(nzeros_D2, 1);
vals_D2 = zeros(nzeros_D2, 1);
index_D2 = 1;
function set_entry_D2(row, col, val)
    rows_D2(index_D2) = row;
    cols_D2(index_D2) = col;
    vals_D2(index_D2) = val;
    index_D2 = index_D2 + 1;
end


for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            x = x1 + (i-1) * dx;
            y = y1 + (j-1) * dy;
            z = z1 + (k-1) * dz;
            
            node = node_index(i, j, k);

            node_coords(node, :) = [x y z];
            
            % x edge
            if i ~= ni
                x_edge = x_edge_index(i, j, k);
                set_entry_D0(x_edge, node, -1);
                set_entry_D0(x_edge, node_index(i+1, j, k), 1);

                edge_coords(x_edge, :) = [x + dx/2, y, z];
            end

            % y edge
            if j ~= nj
                y_edge = y_edge_index(i, j, k);
                set_entry_D0(y_edge, node, -1);
                set_entry_D0(y_edge, node_index(i, j+1, k), 1);

                edge_coords(y_edge, :) = [x, y + dy/2, z];
            end

            % z edge
            if k ~= nk
                z_edge = z_edge_index(i, j, k);
                set_entry_D0(z_edge, node, -1);
                set_entry_D0(z_edge, node_index(i, j, k+1), 1);

                edge_coords(z_edge, :) = [x, y, z + dz/2];
            end

            % yz face
            if j ~= nj && k ~= nk
                face = yz_face_index(i, j, k);
                face_coords(face, :) = [x, y + dy/2, z + dz/2];
                close_y_edge = y_edge_index(i, j, k);
                close_z_edge = z_edge_index(i, j, k);
                far_y_edge = y_edge_index(i, j, k+1);
                far_z_edge = z_edge_index(i, j+1, k);
                set_entry_D1(face, close_y_edge, 1);
                set_entry_D1(face, close_z_edge, -1);
                set_entry_D1(face, far_y_edge, -1);
                set_entry_D1(face, far_z_edge, 1);
            end

            % xz face
            if i ~= ni && k ~= nk
                face = xz_face_index(i, j, k);
                face_coords(face, :) = [x + dx/2, y, z + dz / 2];
                close_x_edge = x_edge_index(i, j, k);
                close_z_edge = z_edge_index(i, j, k);
                far_x_edge = x_edge_index(i, j, k+1);
                far_z_edge = z_edge_index(i+1, j, k);
                set_entry_D1(face, close_x_edge, -1);
                set_entry_D1(face, close_z_edge, 1);
                set_entry_D1(face, far_x_edge, 1);
                set_entry_D1(face, far_z_edge, -1);
            end

            % xy face
            if i ~= ni && j ~= nj
                face = xy_face_index(i, j, k);
                face_coords(face, :) = [x + dx/2, y + dy/2, z];
                close_x_edge = x_edge_index(i, j, k);
                close_y_edge = y_edge_index(i, j, k);
                far_x_edge = x_edge_index(i, j+1, k);
                far_y_edge = y_edge_index(i+1, j, k);
                set_entry_D1(face, close_x_edge, 1);
                set_entry_D1(face, close_y_edge, -1);
                set_entry_D1(face, far_x_edge, -1);
                set_entry_D1(face, far_y_edge, 1);
            end

            % volume
            if i ~= ni && j ~= nj && k ~= nk
                volume = volume_index(i, j, k);
                volume_coords(volume, :) = [x + dx/2, y + dy/2, z + dz/2];
                close_yz_face = yz_face_index(i, j, k);
                close_xz_face = xz_face_index(i, j, k);
                close_xy_face = xy_face_index(i, j, k);
                far_yz_face = yz_face_index(i+1, j, k);
                far_xz_face = xz_face_index(i, j+1, k);
                far_xy_face = xy_face_index(i, j, k+1);
                set_entry_D2(volume, close_yz_face, -1);
                set_entry_D2(volume, close_xz_face, -1);
                set_entry_D2(volume, close_xy_face, -1);
                set_entry_D2(volume, far_yz_face, 1);
                set_entry_D2(volume, far_xz_face, 1);
                set_entry_D2(volume, far_xy_face, 1);
            end

        end
    end
end

D0 = sparse(rows_D0, cols_D0, vals_D0, n_edges, n_nodes, nzeros_D0);
D1 = sparse(rows_D1, cols_D1, vals_D1, n_faces, n_edges, nzeros_D1);
D2 = sparse(rows_D2, cols_D2, vals_D2, n_volumes, n_faces, nzeros_D2);

end

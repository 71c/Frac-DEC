
x1 = 0.0; x2 = 2.0;
y1 = 1.0; y2 = 3.0;
z1 = 0.5; z2 = 5.0;

nx = 3;
ny = 6;
nz = 8;
ni = nx + 1;
nj = ny + 1;
nk = nz + 1;
[D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);

M1 = create_M1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, x1, x2, y1, y2, z1, z2, 1.6);

nzeros_pred = ni * nj * nk * (nx + ny + nz) / 2;
nzeros = nnz(M1);
fprintf("M1: nzeros: %d, nzeros pred: %d\n", nzeros, nzeros_pred);


M2 = create_M2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index, x1, x2, y1, y2, z1, z2, 1.6);
nzeros_pred = ni * nj * nk * (ny * nz + nx * nz + nx * ny) / 4;
nzeros = nnz(M2);
fprintf("M2: nzeros: %d, nzeros pred: %d\n", nzeros, nzeros_pred);


M3 = create_M3(nx, ny, nz, node_index, volume_index, x1, x2, y1, y2, z1, z2, 1.6);
nzeros_pred = ni * nj * nk * nx * ny * nz / 8;
nzeros = nnz(M3);
fprintf("M3: nzeros: %d, nzeros pred: %d\n", nzeros, nzeros_pred);

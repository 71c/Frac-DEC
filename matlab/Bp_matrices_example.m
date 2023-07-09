x1 = 0.0; x2 = 1.0;
y1 = 0.0; y2 = 1.0;
z1 = 0.0; z2 = 1.0;
nx = 1;
ny = 1;
nz = 1;

[D0, D1, D2, node_coords, edge_coords, face_coords, volume_coords, edge_lengths, face_areas, volumes, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index, yz_face_index, xz_face_index, xy_face_index, volume_index] = get_3D_rectangle_mesh(x1, x2, y1, y2, z1, z2, nx, ny, nz);

B1 = create_B1(nx, ny, nz, node_index, node_direction_index, x_edge_index, y_edge_index, z_edge_index);
B2 = create_B2(nx, ny, nz, node_index, node_direction_index, yz_face_index, xz_face_index, xy_face_index);
B3 = create_B3(nx, ny, nz, node_index, volume_index);

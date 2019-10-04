function show_func(mesh,f)
col = 'w';
if length(f) == mesh.nv
    col = 'interp';
elseif length(f) == mesh.nf
    col = 'flat';
elseif length(f) == mesh.ne
    mesh = nc_mesh(mesh);
    patch('faces',mesh.om.triangles,'vertices',mesh.vertices,'FaceColor','w','EdgeColor','none'); colorbar
    col = 'interp';
end
if size(f,2) == 1
    cw = 'CData';
elseif size(f,2) == 3
    cw = 'FaceVertexCData';
end

patch('faces',mesh.triangles,'vertices',mesh.vertices,cw,f,'FaceColor',col,'EdgeColor','none'); colorbar
axis equal; cameratoolbar; cameratoolbar('SetCoordSys','none'); axis off;
function w = vortices_func( mesh, vi, U, a )

X = mesh.vertices;
v = X(vi,:);

w = zeros(mesh.nv,1);
for i = 1:size(v,1)

    vc = repmat( v(i,:), mesh.nv, 1 );
    r = mesh.normv( X - vc );
    w = w + ( U(i)/a )*exp( .5*( 1 - r.^2/a^2 ) );

end
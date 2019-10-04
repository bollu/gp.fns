function w = vortices_surf( tnl, taylor, vi, U, a, use_euclidean )

if nargin < 6
    use_euclidean = false;
end

X = tnl.mesh.vertices;
T = tnl.mesh.triangles;

v = X(vi,:);

% default parameters
if nargin < 4
    U = repmat(0.1592, size(v,1), 1);     % maximum tangential velocity
end

if nargin < 5
    a = 0.0477;     % vortex' core size
end

% initial distribution of vorticity
w = zeros(tnl.mesh.nv,1);

for i = 1:size(v,1)
    
    if ( use_euclidean || ...
         isprop( tnl.mesh, 'IVid' ) == 1 || ...
         isprop( tnl.mesh, 'map2i' ) == 1 )
        % on the plane
        vc = repmat( v(i,:), tnl.mesh.nv, 1 );
        r = tnl.mesh.normv( X - vc );
    else
        % distance of a point from vortex' center (use Geodesics in heat)
        r = geodesic_heat( tnl.mesh, vi(i) );
    end
    
    % add distribution of the current vortex
    if taylor == 1
        w = w + ( U(i)/a )*( 2 - r.^2/a^2 ).*exp( .5*( 1 - r.^2/a^2 ) );
    else
        w = w + ( U(i)/a )*exp( .5*( 1 - r.^2/a^2 ) );
    end
end


% % parameters for domain [-pi,pi]x[-pi,pi]
% U = 1;                              % maximum tangential velocity
% a = .3;                             % vortex' core size
% v = [ [-.4, 0, 0]; [.4, 0, 0] ];    % two vortices
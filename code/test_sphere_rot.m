%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

paths;

meshname = 'sphere_s3';
sphere = STREAM_TNL( meshname );
mesh = sphere.mesh;

run_sim = 1;
run_vis = 1;

%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%
global S W;

h = 0.002;
steps = 2;
NU = 0;
filename = sprintf('%s_h=%g_NU=%g_rot.mat', meshname, h, NU);

if run_sim == 1
    
S = 10*sum(sphere.mesh.LB_basis(:,[3,4]),2)+5*sphere.mesh.LB_basis(:,10);

for k = 1:1/(h*steps)
    S2 = sphere.run_sim(h, S(:,end), steps, NU);
    S = [ S S2(:,end) ];

    W = sphere.stream2vort( S );

    pause(0.001);
    save( filename, 'sphere', 'W' );
end

end

%%%%%%%%%%%%%%%%%
% Visualization %
%%%%%%%%%%%%%%%%%
if run_vis == 1
    
load( filename );
mesh = sphere.mesh;
k = size(W,2);

figure;    
for i = 1:k
    clf; show_func(mesh,W(:,i));
    pause(0.01);    
end
    

end
%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

paths;

meshname = 'sphere_s4';
sphere = STREAM_TNL( meshname );
mesh = sphere.mesh;

run_sim = 1;
run_vis = 1;

%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%
global S W;

vi = [ 37362; 37340 ];
U =  -0.1; a = 0.02;

w1 = mesh.project( vortices_func( mesh, vi(1), U, a ) );
w2 = mesh.project( vortices_func( mesh, vi(2), -U, a ) );
w = w1+w2;
% figure; show_func( mesh, w );

S = - mesh.laplaceInverseF( w );
s1 = - mesh.laplaceInverseF( w1 );
s2 = - mesh.laplaceInverseF( w2 );

h = 0.1;
steps = 200;
NU = 0.7e-4;
alpha = 50;
blendin = 0.5;

filename = sprintf( '%s_h=%g_NU=%g_jet.mat', meshname, h, NU);

if run_sim == 1
    
for k = 1:steps
    S2 = sphere.run_sim(h, S(:,end), 1, NU);
    S = [ S S2(:,end) ];
    
    % we inject a portion of the initial vorticity
    blendscale = min(h*k, blendin)/blendin;
        
    S(:,end) = S(:,end) + .5*alpha*h*blendscale*(s1+s2);
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
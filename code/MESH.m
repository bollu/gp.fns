classdef MESH < handle
    
    properties (Access='public')
        
        name
        
        vertices
        triangles
        
        N
        
        nv
        nf
        
        va
        va_inv
        
        A
        A_INV
        
        ta
        
        M % mass_matrix_barycentric
        
        W % cotan laplacian without mass matrix
        W_precon % preconditioner for W
        
        L % the full laplacian
        LL % L squared
        
        % for the smallest eigenpairs
        LB_basis
        LB_basisi
        LB_evals
        
        % triangle indices (123, 231)
        II
        JJ
    end
    
    properties (Access='protected')
                
        % stuff needed to compute vf2op:
        JJc1
        JJc2
        JJc3
        IIn
        JJn
        
        % remember for iterative solve
        last_laplace_inverse = [];
    end
    
    methods
        
        function [ mesh ] = MESH(meshname)
            
            mesh.name = meshname;
            [mesh.vertices, mesh.triangles] = MESH_READER.readOff([meshname '.off']);
                                   
            mesh.nv = length(mesh.vertices);
            mesh.nf = length(mesh.triangles);
            
            % triangle areas and vertex normals
            X = mesh.vertices; T = mesh.triangles;
            NN = cross(X(T(:,1),:)-X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:));
            mesh.ta = MESH.normv(NN)/2;
            mesh.N = NN./repmat(MESH.normv(NN),1,3);
            
            % initialize indices for spare matrices:
            mesh.II = [T(:,1);T(:,2);T(:,3)];
            mesh.JJ = [T(:,2);T(:,3);T(:,1)];
            
            mesh.IIn = [mesh.II;mesh.JJ;mesh.II;mesh.JJ];
            mesh.JJn = [mesh.JJ;mesh.II;mesh.II;mesh.JJ];
            
            % Use the Barycentric areas
            mesh.M = mesh.mass_matrix_barycentric();
            AA = sum(mesh.M,2);
            mesh.va = full(AA);
            mesh.va_inv = ones(size(mesh.va)) ./ mesh.va;
            mesh.A = spdiags(mesh.va, 0, mesh.nv, mesh.nv);
            mesh.A_INV = spdiags(mesh.va_inv, 0, mesh.nv, mesh.nv);
            
            % compute cotan Laplacian
            mesh.W = mesh.cotLaplacian();
            mesh.L = mesh.A_INV * mesh.W;
            mesh.LL = mesh.L * mesh.L;
            
            % compute some eigenvectors of L
            [evecs, evals] = eigs(mesh.W, mesh.A, min(50, mesh.nv), 1e-5);
            [~,ii] = sort(diag(evals));
            evals = evals(ii,ii); evecs = evecs(:,ii);
            mesh.LB_basis = evecs;
            mesh.LB_basisi = evecs'*mesh.A;
            mesh.LB_evals = diag(evals);
            
            % create preconditioner for W
            opts.diagcomp = 1E-5;
            opts.type = 'ict';
            opts.droptol = 1E-5;
            mesh.W_precon = ichol(mesh.W, opts);
            
            % initialize data needed for vf2op
            [ mesh.JJc1, mesh.JJc2, mesh.JJc3 ] = mesh.vf2op_initialization();
            
        end
        
        function [ l2 ] = dot( mesh, f, g )
            l2 = f'*(mesh.va.*g);
        end
        
        function [ ke ] = vort2energy( mesh, w )
            s = - mesh.laplaceInverseF( w );
            ke = - mesh.dot(s, w);
        end
        
        function [ l2 ] = dotV( mesh, U, V )
            a = U.*V;
            sa = sum(a, 2);
            l2 = dot(mesh.ta, sa);
        end
        
        function [ nrm ] = norm( mesh, f )
            n2 = mesh.dot(f,f);
            nrm = sqrt(n2);
        end
        
        function [ nrm ] = normV( mesh, V )
            n2 = mesh.dotV(V,V);
            nrm = sqrt(n2);
        end
        
        function [ dp, dp2 ] = trace_dot( mesh, U, V )
            dg = full(diag(U'*mesh.A*V));
            dp = 0.5*sum( dg );
            dp2 = 0.5*dot( mesh.va, dg );
        end

        function [ JV ] = J( mesh, V )
            JV = cross(mesh.N, V);
        end
        
        function [ gradF ] = grad( mesh, f )
            % input function is defined on vertices (nv x 1)
            T = mesh.triangles;
            Ar = repmat(mesh.ta,1,3);
            G = repmat(f(T(:,1)),1,3) .* mesh.JJc1 + ...
                repmat(f(T(:,2)),1,3) .* mesh.JJc2 + ...
                repmat(f(T(:,3)),1,3) .* mesh.JJc3;
            gradF = G./(2*Ar);
        end
        
        function [ V ] = curlF( mesh, s )
            V = - mesh.J( mesh.grad( s ));
        end
        
        function [ curlw ] = curlV( mesh, V )
            X = mesh.vertices;
            T = mesh.triangles;
            
            % check if this is the same as for vf2op...
            V1 = X(T(:,1),:);
            V2 = X(T(:,2),:);
            V3 = X(T(:,3),:);
            
            C1 = V3 - V2;
            C2 = V1 - V3;
            C3 = V2 - V1;
            
            I = [T(:,1); T(:,2); T(:,3)];
            J = ones(length(I),1);
            S = dot([C1; C2; C3], [V; V; V], 2);
            curlw = sparse(I,J,S, mesh.nv, 1)/2;
            
            curlw = curlw ./ mesh.va;
            curlw = full(curlw);
        end
        
        function [ fp ] = project( mesh, f )
            fp = f - (mesh.LB_basisi(1,:)*f)*mesh.LB_basis(:,1);
        end
        
        function [ lf ] = laplaceF( mesh, f )
            lf = mesh.L * f;
        end
        
        function [ lif ] = laplaceInverseF( mesh, f, guess, prec )
            if nargin < 3
                guess = mesh.last_laplace_inverse;
            end
            if nargin < 4
                prec = 1e-12;
            end
            if mesh.nv < 200
                lif = mesh.L \ f;
            else
                WL = mesh.W;
                PC = mesh.W_precon;
                b = mesh.project(mesh.A*f);
                [lif,flag,~,~,~] = pcg(WL,b,prec,500,PC,PC',guess);
                if flag
                    lif = mesh.L \ f;
                end
            end
            lif = mesh.project( lif );
            mesh.last_laplace_inverse = lif;
        end
        
        function M = mass_matrix_barycentric(mesh)
            Ar = mesh.ta;
            Mij = 1/12*[Ar; Ar; Ar];
            Mji = Mij;
            Mii = 1/6*[Ar; Ar; Ar];
            Mn = [Mij;Mji;Mii];
            In_MM = [mesh.II;mesh.JJ;mesh.II];
            Jn_MM = [mesh.JJ;mesh.II;mesh.II];
            M = sparse(In_MM,Jn_MM,Mn,mesh.nv,mesh.nv);
        end
        
        function [ W ] = cotLaplacian(mesh)
            X = mesh.vertices;
            T = mesh.triangles;
            
            % Find orig edge lengths and angles
            L1 = MESH.normv(X(T(:,2),:)-X(T(:,3),:));
            L2 = MESH.normv(X(T(:,1),:)-X(T(:,3),:));
            L3 = MESH.normv(X(T(:,1),:)-X(T(:,2),:));
            
            A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
            A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
            A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
            AA = [A1,A2,A3];
            AA = acos(AA);
            
            % The Cot Laplacian
            S = 0.5*cot([AA(:,3);AA(:,1);AA(:,2)]);
            Sn = [-S;-S;S;S];
            W = sparse(mesh.IIn,mesh.JJn,Sn,mesh.nv,mesh.nv);
        end
        
        function [JC1, JC2, JC3 ] = vf2op_initialization(mesh)
            X = mesh.vertices;
            T = mesh.triangles;
            Nf = mesh.N;
            v1 = X(T(:,1),:);
            v2 = X(T(:,2),:);
            v3 = X(T(:,3),:);
            C1 = v3 - v2;
            C2 = v1 - v3;
            C3 = v2 - v1;
            JC1 = cross(Nf, C1);
            JC2 = cross(Nf, C2);
            JC3 = cross(Nf, C3);
        end
        
        function op = DV(mesh, Vf)
            WW = mesh.vf2opWW(Vf);
            op = mesh.A_INV*WW;
        end
        
        function op = adDV(mesh, Vf)
            WW = mesh.vf2opWW(Vf);
            op = mesh.A_INV*WW';
        end
        
        function WW = vf2opWW(mesh, Vf)
            Sij = 1/6*[dot(mesh.JJc2,Vf,2); dot(mesh.JJc3,Vf,2); dot(mesh.JJc1,Vf,2)];
            Sji = 1/6*[dot(mesh.JJc1,Vf,2); dot(mesh.JJc2,Vf,2); dot(mesh.JJc3,Vf,2)];
            Sn = [Sij;Sji;-Sij;-Sji];
            WW = sparse(mesh.IIn,mesh.JJn,Sn,mesh.nv,mesh.nv);
        end
    end
    
    methods (Static)
        
        function [ nv ] = normv( v )
            nv = sqrt(sum(v.^2,2));
        end
        
    end
end
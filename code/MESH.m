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
            % T : FX3 -> Ix
            % X : Vx3 -> R
            X = mesh.vertices; T = mesh.triangles;
            % x21 = x[t[n, 1]] - x[t[n,2]]
            % x31 = x[t[n, 1]] - x[t[n,3]]
            NN = cross(X(T(:,1),:)-X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:));
            % mesh.ta: area of the triangle is half the area of the parallelogram.
            mesh.ta = MESH.normv(NN)/2; % Vx1
            % mesh.N: normalized normal vectors.
            mesh.N = NN./repmat(MESH.normv(NN),1,3);

            % initialize indices for spare matrices:
            mesh.II = [T(:,1);T(:,2);T(:,3)]; % 3fx1 column vector given by (T[:, 1] <> T[:, 2] <> T[:, 3])
            mesh.JJ = [T(:,2);T(:,3);T(:,1)]; % 3fx1 column vector given by (T[:, 2] <> T[:, 3] <> T[:, 1])
            % see that (mesh.II[i], mesh.JJ[i]) gives us the adjacent vertices: 1 → 2, 2 → 3, 3 → 1

            % This makes no sense, because it should be vx1.
            % see that IIn, JJn gives us all the edges: I → J, J → I, I → I, J → J
            mesh.IIn = [mesh.II;mesh.JJ;mesh.II;mesh.JJ]; %12fx1 column vector given by (II <> JJ <> II <> JJ)
            mesh.JJn = [mesh.JJ;mesh.II;mesh.II;mesh.JJ]; %12fx1 column vector given by (JJ <> II <> JJ <> II)

            % compute mass matrix.
            mesh.M = mesh.mass_matrix_barycentric();
            % what is AA? sum of each of the rows. AA[i] = Σj M[i][j].
            % AA[i] will be the total mass "adjacent" to `i`.
            AA = sum(mesh.M,2);
            % mesh.va = [1xN].
            % mesh.va is the mass matrix for each vertex 'a'. I guess 'a' is 'area' for the vertices.
            mesh.va = full(AA); % convert the sparse matrix 'AA' to a dense matrix 'va'.
            % pointwise inverse of mesh.va [1xN].
            mesh.va_inv = ones(size(mesh.va)) ./ mesh.va;
            % sparse matrix of size NxN whose diagonal is given by `mesh.va`.
            % spdiags(_, 0, _, _): place on main diagonal.
            % spdiags(_, δ, _, _): place relative to main diagonal.
            %    +ve → above main diagonal, -ve → below main diagonal
            % ChatGPT
            % =======
            %   In particular, the entries of the mass matrix A are given by the
            %   integrals of the product of pairs of basis functions over the
            %   domain of each finite element. The diagonal entries of the mass
            %   matrix correspond to the masses associated with each node in the
            %   mesh.
            %   For example, in the context of a surface mesh, the mass matrix A is
            %   typically constructed from the surface area of each triangle in the
            %   mesh. Specifically, the diagonal entries of the mass matrix are
            %   proportional to the areas of the triangles adjacent to each vertex
            %   in the mesh.
            mesh.A = spdiags(mesh.va, 0, mesh.nv, mesh.nv);
            % inverse sparse matrix of 'A', since [NxN]
            mesh.A_INV = spdiags(mesh.va_inv, 0, mesh.nv, mesh.nv);

            % compute cotan Laplacian
            mesh.W = mesh.cotLaplacian();
            mesh.L = mesh.A_INV * mesh.W;
            mesh.LL = mesh.L * mesh.L;

            % compute some eigenvectors of L.
            %  solve generalized eigenvalue problem `Wx = Ax`. gives at most 50 eigenvalues
            %  at a tolerance of relative eigenvalue error `1e-5` [ie, (Wx-Ax / λAx) < 1e-5]/
            % each column of 'evecs' in an eigenvector.
            [evecs, evals] = eigs(mesh.W, mesh.A, min(50, mesh.nv), 1e-5);
            % sort eigenvalues in ascending order, and return indexes that puts it in ascending order.
            %   TODO: I guess because the eigenvalues are negative, ascending order means that the
            %     highest magnitude eigenvalue will succeed?
            [~,ii] = sort(diag(evals));
            % top eigenvalues and eigenvectors
            evals = evals(ii,ii); evecs = evecs(:,ii);
            % the basis of the possion problem is LB? why B? 'Laplacian basis'?
            mesh.LB_basis = evecs;
            % product of transpose of eigenvectors with the mass matrix.
            % chatGPT:
            % =======
            % Overall, the code computes a matrix that maps the
            % Laplace-Beltrami eigenbasis onto the standard basis of Euclidean
            % space. This can be useful for representing functions on the
            % surface in terms of their Laplace-Beltrami eigenfunctions, or for
            % computing geometric quantities that depend on the eigenfunctions,
            % such as the curvature or heat flow on the surface.
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

        % https://graphics.stanford.edu/courses/cs468-13-spring/assets/lecture12-lu.pdf
        function M = mass_matrix_barycentric(mesh)
            Ar = mesh.ta; % take areas of all triangles: nx1
            % https://libigl.github.io/libigl-python-bindings/tut-chapter1/#mass-matrix
            Mij = 1/12*[Ar; Ar; Ar]; % Mij = 1/12 times the vector (area, area, rea)
            Mji = Mij; % the Mji is symmetric. TODO
            Mii = 1/6*[Ar; Ar; Ar]; %derivation of mass matrix.
            Mn = [Mij;Mji;Mii]; % concatenation of Mij, Mji, Mii
            % we want to define (IxJ), (JxI), (IxI) values.
            In_MM = [mesh.II;mesh.JJ;mesh.II]; %rows where values are defined in M
            Jn_MM = [mesh.JJ;mesh.II;mesh.II]; %columns where values are defined in M
            M = sparse(In_MM,Jn_MM,Mn,mesh.nv,mesh.nv);
        end

        function [ W ] = cotLaplacian(mesh)
            X = mesh.vertices;
            T = mesh.triangles;

            % Find orig edge lengths and angles
            % L1 = v[t[2]] - v[t[3]]
            L1 = MESH.normv(X(T(:,2),:)-X(T(:,3),:));
            % L2 = v[t[1]] - v[t[3]]
            L2 = MESH.normv(X(T(:,1),:)-X(T(:,3),:));
            % L3 = v[t[1]] - v[t[2]]
            L3 = MESH.normv(X(T(:,1),:)-X(T(:,2),:));

            % cosine rule of angles:
            % a^2 = b^2 + c^2 - 2bc cos(α)
            %    b^2 + c^2 - a^2 =  2bc cos(α)
            A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3); % cos(a1)
            A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3); % cos(a2)
            A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2); % cos(a3)
            AA = [A1,A2,A3];
            AA = acos(AA); % actual angles of all [a1, a2, a3]

            % http://www.geometry.caltech.edu/pubs/DGPDEC_SIG13.pdf Page 70
            % The Cot Laplacian.
            S = 0.5*cot([AA(:,3);AA(:,1);AA(:,2)]); %
            Sn = [-S;-S;S;S]; % (II → JJ), (JJ → II), (II → II), (JJ → JJ) adjacency information.

            % W[IIn[k], JJn[k]] = Sn[k], where W is a [VxV] matrix.
            % Recall that IIn, JJn, contains all possible adjacencies.
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

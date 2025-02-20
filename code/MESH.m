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

            % compute mass matrix in barycentric coords?
            % in (https://graphics.stanford.edu/courses/cs468-13-spring/assets/lecture12-lu.pdf), this is called A^triangle[i][j]
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
            % in (https://graphics.stanford.edu/courses/cs468-13-spring/assets/lecture12-lu.pdf), this is called
            % as the lumped mass matrix, with a(ii) = Area(dual cell i) ?
            % TODO: why is this legal? Why do we phrase laplacian in terms of this?
            mesh.A = spdiags(mesh.va, 0, mesh.nv, mesh.nv);
            % inverse sparse matrix of 'A', since [NxN]
            mesh.A_INV = spdiags(mesh.va_inv, 0, mesh.nv, mesh.nv);

            % compute cotan Laplacian
            mesh.W = mesh.cotLaplacian(); % why W?
            mesh.L = mesh.A_INV * mesh.W; % why L?
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
            % TODO: understand preconditioning. Why does it improve performance?
            % compute W_precon such that W_precon@W_precon = mesh.W
            mesh.W_precon = ichol(mesh.W, opts);

            % initialize data needed for vf2op
            % ChatGPT:
            % ========
            %  These matrices are used to represent the gradient and divergence
            %  operators on the mesh in the finite element method. Specifically,
            %  they are used to compute the vector field Laplacian operator, which
            %  is a generalization of the Laplace-Beltrami operator to vector
            %  fields on the surface.
            % Bollu
            % =====
            % mesh.JJc1 is a [Fx3] matrix, where mesh.JC1[k] is
            % JC1 has vectors. Let tri[k] = (v1, v2, v3). Then JC1[k] is the vector in the triangle
            % that is orthogonal to (v2 - v3) and has length (v2 - v3). It is given by: (v2 - v3) x unit-normal(tri[k])
            [ mesh.JJc1, mesh.JJc2, mesh.JJc3 ] = mesh.vf2op_initialization();

        end

        function [ l2 ] = dot( mesh, f, g )
            % inner product where 'mesh.va' is the mass matrix, since
            % the mass matrix `mesh.va` represents the inner product of basis vecotrs: A[i, j] := int_A hi . hj.
            % where hi is the hat function: 1 at the vertex 'i', 0 at all other vertices.
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
        %       f(a) = 1
        %      /|
        %     / |
        %    /  |
        %   /   H  ^
        %  /    |  | JJc1 = yhat * |c - b| = yhat * B
        % /     |
        %*---B--*
        %f(b) = f(c) = 0
        %
        % ([1, 1, 1] .* yhat * B) / (2 * 1/2BH)
        % = (yhat * B) / (B H)
        % = yhat / H
        % which is exactly what the gradient should be, if we assume that the
        %   gradient changes by 1/H over the full 'H' distance to get from '1' to
        %   '0'.
        % TODO: but this is different from the paper about using vector fields,
        %  since thee, we are supposed to have functions on faces, NOT vertices!
        %  In that theory, we were supposd to use an averaging operator in the worst case,
        %  something I do not see here.
        function [ gradF ] = grad( mesh, f )
            % input function is defined on vertices (nv x 1)
            T = mesh.triangles;
            Ar = repmat(mesh.ta,1,3); % make a vector [area[tri], area[tri], area[tri]]
            G = repmat(f(T(:,1)),1,3) .* mesh.JJc1 + ... % for all faces, compute f(vertex '1' of face), which rescales the vector that points towards vector 1
                repmat(f(T(:,2)),1,3) .* mesh.JJc2 + ...
                repmat(f(T(:,3)),1,3) .* mesh.JJc3;
            gradF = G./(2*Ar); % divide each gradient by the area of the triangle
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
          % ChatGPT:
          % =======
          %   project to space of harmonic functions on the mesh.
          %   we project to space of the laplace beltrami eigenvalues.
          %   y = x - A^{-1}xA
          %   TODO: I am still very confused about this, because it seems to be REMOVING the components that
          %   are harmonic? it's doing f - <stuff> ?
          %   unless the operator below in the <stuff> is projecting to the subspace that is *orthogonal* to harmonics?
          %   seems weird.
            fp = f - (mesh.LB_basisi(1,:)*f)*mesh.LB_basis(:,1);
        end

        function [ lf ] = laplaceF( mesh, f )
            lf = mesh.L * f;
        end

        function [ lif ] = laplaceInverseF( mesh, f, guess, prec )
            if nargin < 3
                guess = mesh.last_laplace_inverse; % laplace inverse cache.
            end
            if nargin < 4
                prec = 1e-12;
            end
            if mesh.nv < 200
                lif = mesh.L \ f; % perform inverse directly
            else
                WL = mesh.W;
                PC = mesh.W_precon;
                b = mesh.project(mesh.A*f);
                % solve WL u = b with preconditioner 'prec' that is an approximate inverse for WL.
                % TODO: see why 'prec' is an approximate inverse for 'WL'.
                [lif,flag,~,~,~] = pcg(WL,b,prec,500,PC,PC',guess);
                if flag
                    lif = mesh.L \ f; % If algorithm does not converge, directly perform inverse.
                end
            end
            lif = mesh.project( lif );
            mesh.last_laplace_inverse = lif;
        end

        % https://graphics.stanford.edu/courses/cs468-13-spring/assets/lecture12-lu.pdf
        % consider a triangle ABC
        % t=0             A ha(A) = 1
        %                 /\
        %                /  \
        %             --=----=-- ha(t) = t ; width of this =
        %              /      \
        % t=1         C--------B ha(B) = ha(C) = 0
        %
        % suppose i = j. Then:
        %  \int_A ha . ha
        % = \int_t ha(t) ha(t)
        function M = mass_matrix_barycentric(mesh) % [VxV] matrix. M[i][j] = integral_{triangle} h[i] h[j] dA
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

        % JC1 has vectors such that JC1[k] is the vector that is orthogonal to the normal of the triangle
        % as well as the edge (v2 - v3), and has length that if (v2 - v3).
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
        % chatgpt
        % ======
        % To compute the operator corresponding to the divergence of Vf, the function
        % first calls the vf2opWW method of the mesh structure to obtain a sparse
        % matrix WW that approximates the dot product of the gradient operator
        % and Vf. The vf2opWW method constructs a weighted sum of the gradients
        % of Vf at each mesh vertex, using the local weights of the Voronoi cell
        % of each vertex.
        % The function then multiplies WW by the inverse of the mass matrix A to
        % obtain an approximation of the divergence operator:
        % op = A_INV * WW
            WW = mesh.vf2opWW(Vf);
            op = mesh.A_INV*WW;
        end

        function op = adDV(mesh, Vf)
           % chatgpt
           % ======
           % div(Vf)' = - grad(Vf)
           % div(Vf)' ~ -W' * Vf

           % To compute the operator corresponding to the adjoint of the divergence
           % of Vf, the function first calls the vf2opWW method of the mesh
           % structure to obtain a sparse matrix WW that approximates the
           % dot product of the gradient operator and Vf. The vf2opWW method
           % constructs a weighted sum of the gradients of Vf at each mesh
           % vertex, using the local weights of the Voronoi cell of each
           % vertex.

           % The function then multiplies the transpose of WW by the inverse
           % of the mass matrix A to obtain an approximation of the
           % adjoint of the divergence operator:
            WW = mesh.vf2opWW(Vf);
            op = mesh.A_INV*WW';
        end

        function WW = vf2opWW(mesh, Vf)
          % chatgpt
          % ======
          % The Matlab function vf2opWW computes a matrix WW that corresponds to the
          % operator W in the finite element method, based on the vector field Vf defined
          % on the mesh.
          % First, it calculates the values of the components of the stress
          % tensor Sij and Sji at each mesh vertex using the dot product between
          % the corresponding coordinate vectors and Vf. Then, it assembles the
          % stress tensor components into a single vector Sn and uses the sparse
          % matrix function sparse to create a sparse matrix WW. The matrix WW is
          % assembled based on the indices in mesh.IIn and mesh.JJn, which are
          % the row and column indices, respectively, of the non-zero entries of
          % WW. Finally, WW is returned as the output of the function.
          % In summary, the function vf2opWW maps a vector field
          % Vf to a matrix WW that represents the corresponding
          % operator W in the finite element method.


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

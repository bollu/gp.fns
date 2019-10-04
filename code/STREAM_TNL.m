classdef STREAM_TNL < handle
    %STREAM_TNL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
    end
    
    properties (Access='private')
        ID
    end
    
    methods
        
        function tnl = STREAM_TNL( meshname )
            tnl.mesh = MESH( meshname );
            tnl.ID = speye(tnl.mesh.nv);
        end
        
        function [ Dv ] = adDV( tnl, s )
            curlS = tnl.mesh.curlF( s );
            Dv = tnl.mesh.adDV( curlS );
        end
        
        function [ Dv ] = DV( tnl, s )
            curlS = tnl.mesh.curlF( s );
            Dv = tnl.mesh.DV( curlS );
        end
        
        function [ CAY ] = cay( tnl, h, Dv, w )
            CAY = ( tnl.ID - 0.5*h*Dv ) \ (w + 0.5*h*Dv*w);
        end
        
        function [ stream ] = advect( tnl, h, stream_old, nu )
                        
            L = tnl.mesh.L;
            LL = tnl.mesh.LL;
            
            w_old = - L*stream_old;
            Dv0 = tnl.adDV( stream_old );
            LM = (tnl.ID - 0.5*h*Dv0)*w_old;
            LF = nu*h*L*LM/2;

            function [ diff, jac_sparse ] = DEL( stream )
                Dv = tnl.adDV( stream );
                w = - L*stream;
                RM = (tnl.ID + 0.5*h*Dv)*w;
                RF = nu*h*L*RM/2;
                
                diff = RM - LM + RF + LF;
                
                if nargout > 1
                    Dvm = tnl.adDV( - w );
                    DvL = Dv*L;
                    jac_sparse = -( 0.5*nu*h*LL + 0.25*nu*h^2*L*DvL + 0.25*nu*h^2*L*Dvm ...
                                    + L + (0.5*h)*DvL - (0.5*h)*Dvm );    
                end
            end
            
            w_guess = tnl.cay( - h, Dv0, w_old ) - 2*LF;
            guess = - tnl.mesh.laplaceInverseF( w_guess );
            
            tol = 1E-7;
            stream = tnl.newton_solve( @DEL, guess, tol );

        end
        
        function [ sol ] = newton_solve( tnl, FUNC, s, tol )
            
            MAX_IT = 20;
            
            cnt = 0;
            res = 1;
            
            [ Fs, Js ] = FUNC( s );
            
            fprintf('\n NEWTON SOLVE:\n');
                        
            while res > tol
                
                iterstart = tic;
                
                cnt = cnt+1;
                Fs = tnl.mesh.project( Fs );
                
                JF_INV_F = Js \ Fs;
                
                s = tnl.mesh.project(s - JF_INV_F);
                [ Fs, Js ] = FUNC( s );
                res = norm( Fs );
                
                itertime = toc(iterstart);
                
                if cnt == MAX_IT
                    break;
                end
                fprintf( '\tIT %3d, TIME=[ %4.2f s ]  res=%e\n', cnt, itertime, res );
            end
            
            if cnt == MAX_IT
                fprintf('\nNEWTON DID NOT CONVERGE!!! res=%e\n', res);
            end
            
            sol = s;
        end
         
        function S = run_sim( tnl, h, stream0, steps, nu )
            
            S = zeros( length(stream0), steps+1 );
            S(:,1) = stream0;
            
            for k = 1:steps
                fprintf( 'computing step %d of %d ', k, steps );
                dtstep = tic;
                S_next = tnl.advect( h, S(:, k), nu );
                S(:, k+1) = S_next;
                fprintf( '[%6.3f s]\n', toc(dtstep) );
            end
        end
        
        function [ W ] = stream2vort( tnl, S )
            W = - tnl.mesh.L*S;
        end
    end
    
end


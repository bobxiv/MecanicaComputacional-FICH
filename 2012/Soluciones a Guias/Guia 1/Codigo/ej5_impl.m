function [u, Nx, Ny, Nt, X, Y] = ej5_impl(Lx, Ly, Lt, dx, dy, dt)
%TP1EJ5 Diferencias finitas en el tiempo, backward Euler

    x = 0 : dx : Lx;
    y = 0 : dy : Ly;
    t = 0 : dt : Lt;
    
    Nx = length(x);
    Ny = length(y);
    Nt = length(t);
    N  = Nx * Ny * Nt;
    
    Q = inline('100*(x+y)', 'x', 'y');
    
    Tizq = 10;
    Tder = 10;
    Tsup = 10;
    qinf = 0;
    
    if Nx <= Ny
        orden = 'x';
    else
        orden = 'y';
    end
    
    K = sparse(N);
    b = zeros(N, 1);
    
    for i = 1 : Nx
        for j = 1 : Ny
            p = pos(i, j, 1, Nx, Ny, orden);
            if i == 1
                K(p, p) = 1;
                b(p) = Tizq;
            elseif i == Nx
                K(p, p) = 1;
                b(p) = Tder;
            elseif j == 1
                K(p, p) = -3/(2*dy);
                p2 = pos(i, j+1, 1, Nx, Ny, orden);
                K(p, p2) = 2/dy;
                p2 = pos(i, j+2, 1, Nx, Ny, orden);
                K(p, p2) = -1/(2*dy);
                b(p) = qinf;
            elseif j == Ny
                K(p, p) = 1;
                b(p) = Tsup;
            else
                K(p, p) = 1;
                b(p) =   0;
            end
        end
    end
    
    for k = 2 : Nt
        for i = 1 : Nx
            for j = 1 : Ny
                p = pos(i, j, k, Nx, Ny, orden);
                if i == 1
                    K(p, p) = 1;
                    b(p) = Tizq;
                elseif i == Nx
                    K(p, p) = 1;
                    b(p) = Tder;
                elseif j == 1
                    K(p, p) = -3/(2*dy);
                    p2 = pos(i, j+1, k, Nx, Ny, orden);
                    K(p, p2) = 2/dy;
                    p2 = pos(i, j+2, k, Nx, Ny, orden);
                    K(p, p2) = -1/(2*dy);
                    b(p) = qinf;
                elseif j == Ny
                    K(p, p) = 1;
                    b(p) = Tsup;
                else
                    K(p,p) = 1/dt + 2/dx^2 + 2/dy^2;
                    p2 = pos(i, j, k-1, Nx, Ny, orden);
                    K(p, p2) = -1/dt;
                    p2 = pos(i-1, j, k, Nx, Ny, orden);
                    K(p, p2) = -1/dx^2;
                    p2 = pos(i+1, j, k, Nx, Ny, orden);
                    K(p, p2) = -1/dx^2;
                    p2 = pos(i, j-1, k, Nx, Ny, orden);
                    K(p, p2) = -1/dy^2;
                    p2 = pos(i, j+1, k, Nx, Ny, orden);
                    K(p, p2) = -1/dy^2;
                    b(p) = Q(x(i), y(j));
                end
            end
        end
    end
    %u=1; return;
    
    u_vec = K \ b;
    u = [];
    for k = 1 : Nt
        u = [u; reshape(u_vec((k-1)*Nx*Ny+1:k*Nx*Ny), Ny, Nx)];
    end
    [X, Y] = meshgrid(x, y);
end
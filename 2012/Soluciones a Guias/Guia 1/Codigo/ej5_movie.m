%Crea los frames de un video simulando el desarrollo temporal de una ED
%          dT/dt - ( k*d^2T/dx^2 + k*d^2T/dy^2 ) = 100*(x+y)
%                                          EN 2D
%Requiere de un estado inicial. Los parametros son:
%
%                      - X  el vector de puntos en el espacio
%                      - T0 como estado inicial
%                      - tf el tiempo final
%                      - dx el intervalo espacial
%                      - dt el intervalo de tiempo
function [F] = ej5_movie(Tipo, Lx, Ly, Lt, dx, dy, dt)

if( strcmpi(Tipo,'Explicito') )
    [T, Nx, Ny, Nt, X, Y] = ej5_expl(Lx, Ly, Lt, dx, dy, dt);
elseif( strcmpi(Tipo,'Implicito') )
    [T, Nx, Ny, Nt, X, Y] = ej5_impl(Lx, Ly, Lt, dx, dy, dt);
elseif( strcmpi(Tipo,'CranckNicholson') )
    [T, Nx, Ny, Nt, X, Y] = ej5_cn(Lx, Ly, Lt, dx, dy, dt);
end

F(Nt) = struct('cdata', [], 'colormap', []);
figure;

for t = 1 : Nt
    surf(X, Y, T([1:Ny]+Ny*(t-1),1:Nx));

    axis([0,Lx,0,Ly,0,50]);
    grid on;
    set(gca, 'xtick', [0:dx:Lx])
    set(gca, 'ytick', [0:dy:Ly])
    set(gca, 'ztick', [0:10:50])
    xlabel('x')
    ylabel('y')
    zlabel('T')
    if( strcmpi(Tipo,'Explicito') )
        title('Explicito')
    elseif( strcmpi(Tipo,'Implicito') )
        title('Implicito')
    elseif( strcmpi(Tipo,'CranckNicholson') )
        title('Crank-Nicholson')
    end

    F(t) = getframe(gcf);
end

end
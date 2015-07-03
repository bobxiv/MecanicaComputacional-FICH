% Ensambla la contribucion nodal de (Ny,j) a la matriz global y al vector 
% independiente. Solo para nodos de CC Superior!!!
function CCSuperior(j)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global Tipoccs;%tipo de cc
    global qccs;%conductividad del borde superior
    global Tccs;%Temperatura impuesta en el borde derecho
    
    global vecX;
    global vecY;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(Ny  ,j  );
    Tim1j = WorldToMat(Ny-1,j  );
    
    dy_m = vecY(Ny-1) - vecY(Ny-2);%atras
    dy_p = vecY(Ny  ) - vecY(Ny-1);%adelante

    if( strcmp(Tipoccs, 'Neumann') )
        K(Tij,Tij  ) = K(Tij,Tij  ) + -subs(k)/(dy_m+dy_p);
        K(Tij,Tim1j) = K(Tij,Tim1j) + -subs(k)/(dy_m+dy_p);
        b(Tij      ) = b(Tij      ) + subs(qccs);
    end
    
    if( strcmp(Tipoccs, 'Dirichlet') )
        K(Tij,Tij  ) = K(Tij,Tij  ) + 1;
        b(Tij      ) = b(Tij      ) + subs(Tccs);
    end
    
end
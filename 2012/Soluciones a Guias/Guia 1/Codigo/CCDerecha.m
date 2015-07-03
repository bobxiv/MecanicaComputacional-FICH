% Ensambla la contribucion nodal de (i,Nx) a la matriz global y al vector 
% independiente. Solo para nodos de CC Derecha!!!
function CCDerecha(i)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global Tipoccd;%tipo de cc
    global qccd;%conductividad del borde derecho
    global Tccd;%Temperatura impuesta en el borde derecho
    
    global vecX;
    global vecY;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,Nx  );
    Tijm1 = WorldToMat(i  ,Nx-1);
    Tijm2 = WorldToMat(i  ,Nx-2);
    
    dx_m = vecX(2) - vecX(1);%atras
    dx_p = vecX(3) - vecX(2);%adelante

     if( EqEspaciado )%estan los nodos estan equiespaciados?
                    %si estan equiespaciados usamos una mejor aproximacion
        dx = vecX(Nx)-vecX(Nx-1);
        
        if( strcmp(Tipoccd, 'Neumann') )
            K(Tij,Tij  ) = K(Tij,Tij  ) - subs(k)*3;
            K(Tij,Tijm1) = K(Tij,Tijm1) + subs(k)*4;
            K(Tij,Tijm2) = K(Tij,Tijm2) - subs(k);
            b(Tij      ) = b(Tij      ) + subs(qccd)*2*dx;
        end
        
     else
                    %si no lo estan una peor
        dx_mm = vecX(Nx-1) - vecX(Nx-2);%entre jm1 y jm2
        
        if( strcmp(Tipoccd, 'Neumann') )
            K(Tij,Tijm1) = K(Tij,Tijm1) + subs(k);
            K(Tij,Tijm2) = K(Tij,Tijm2) - subs(k);
            b(Tij      ) = b(Tij      ) + subs(qccd)*dx_mm;
        end
         
     end
    
    if( strcmp(Tipoccd, 'Dirichlet') )
        K(Tij,Tij  ) = K(Tij,Tij  ) + 1;
        b(Tij      ) = b(Tij      ) + subs(Tccd);
    end
    
end
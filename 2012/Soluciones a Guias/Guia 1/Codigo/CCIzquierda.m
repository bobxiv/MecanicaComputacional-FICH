% Ensambla la contribucion nodal de (i,1) a la matriz global y al vector 
% independiente. Solo para nodos de CC Izquierda!!!
function CCIzquierda(i)

    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global Tipocci;%tipo de cc
    global qcci;%conductividad del borde izquierdo
    global Tcci;%Temperatura impuesta en el borde derecho
    
    global vecX;
    global vecY;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,1  );
    Tijp1 = WorldToMat(i  ,1+1  );
    Tijp2 = WorldToMat(i  ,1+2  );
    
    

    if( EqEspaciado )%estan los nodos estan equiespaciados?
                    %si estan equiespaciados usamos una mejor aproximacion
        dx = vecX(2)-vecX(1);
        
        if( strcmp(Tipocci, 'Neumann') )
            K(Tij,Tij  ) = K(Tij,Tij  ) + subs(k)*3;
            K(Tij,Tijp1) = K(Tij,Tijp1) - subs(k)*4;
            K(Tij,Tijp2) = K(Tij,Tijp2) + subs(k);
            b(Tij      ) = b(Tij      ) + subs(qcci)*2*dx;
        end
        
    else
                    %si no lo estan una peor
        dx_pp = vecX(3) - vecX(2);%entre jp1 y jp2
                    
        if( strcmp(Tipocci, 'Neumann') )
            K(Tij,Tijp1) = K(Tij,Tijp1) - subs(k);
            K(Tij,Tijp2) = K(Tij,Tijp2) + subs(k);
            b(Tij      ) = b(Tij      ) + subs(qcci)*dx_pp;
        end
        
    end
    
    if( strcmp(Tipocci, 'Dirichlet') )
        K(Tij,Tij  ) = K(Tij,Tij  ) + 1;
        b(Tij      ) = b(Tij      ) + subs(Tcci);
    end
    
end
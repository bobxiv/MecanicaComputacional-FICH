% Ensambla la contribucion nodal de (1,j) a la matriz global y al vector 
% independiente. Solo para nodos de CC Inferior!!!
function CCInferior(j)
    
    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global Tipoccin;%tipo de cc
    global qccin;%conductividad del borde inferior
    global Tccin;%Temperatura impuesta en el borde derecho
    
    global vecX;
    global vecY;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(1  ,j  );
    Tip1j = WorldToMat(1+1,j  );
    
    dy_m = vecY(2  ) - vecY(1);%atras
    dy_p = vecY(3  ) - vecY(2);%adelante

    if( strcmp(Tipoccin, 'Neumann') )
        K(Tij,Tij  ) = K(Tij,Tij  ) + -subs(k)/(dy_m+dy_p);
        K(Tij,Tip1j) = K(Tij,Tip1j) + -subs(k)/(dy_m+dy_p);
        b(Tij      ) = b(Tij      ) + subs(qccin);
    end
    
    if( strcmp(Tipoccin, 'Dirichlet') )
        K(Tij,Tij  ) = K(Tij,Tij  ) + 1;
        b(Tij      ) = b(Tij      ) + subs(Tccin);
    end
    
end
% Ensambla la contribucion nodal de (1,j) a la matriz global y al vector 
% independiente. Solo para nodos de CC Inferior!!!
function CCInferior(j)
    
    global K;
    global Nx;
    global Ny;
    global b;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor

    global Tipoccin;%tipo de cc
    global qccin;%conductividad del borde inferior
    global Tccin;%Temperatura impuesta en el borde derecho
    
    global vecR;
    global vecTheta;
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(1  ,j  );
    Tip1j = WorldToMat(1+1,j  );
    Tip2j = WorldToMat(1+2,j  );

    if( strcmp(Tipoccin, 'Neumann') )
        dtheta = vecTheta(2  ) - vecTheta(1);%adelante
    
        K(Tij,Tij  ) = - 3/(2*dtheta);
        K(Tij,Tip1j) = 4/(2*dtheta);
        K(Tij,Tip2j) = - 1/(2*dtheta);
        b(Tij      ) = -subs(qccin);
    end
    
    if( strcmp(Tipoccin, 'Dirichlet') )
        K(Tij,Tij  ) = 1;
        b(Tij      ) = subs(Tccin);
    end
    
end
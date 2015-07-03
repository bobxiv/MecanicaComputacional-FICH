% Ensambla la contribucion nodal de (i,j) a la matriz global
% solo para elementos (i,j) internos!!! Para CC usar las funciones de CC
function ContribucionNodal(i,j)

    global K;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    
    global Dimension;%Dimension del problema
    
    global vecX;
    global vecY;
    
    global EqEspaciado;%si los nodos estan equiespaciados
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    Tij   = WorldToMat(i  ,j  );
    Tijp1 = WorldToMat(i  ,j+1);
    Tijm1 = WorldToMat(i  ,j-1);

    if( EqEspaciado )%estan los nodos estan equiespaciados?
                    %si estan equiespaciados usamos una mejor aproximacion
        dx = vecX(j)-vecX(j-1);
        
        if( strcmp(Dimension, '2D') )
            Tip1j = WorldToMat(i+1,j  );
            Tim1j = WorldToMat(i-1,j  );
            
            dy = vecY(i)-vecY(i-1);

            K(Tij,Tij  ) = K(Tij,Tij  ) - subs(k)*2*dy^2 - subs(k)*2*dx^2 + subs(c)*dx^2*dy^2;
            K(Tij,Tim1j) = K(Tij,Tim1j) + subs(k)*dx^2;
            K(Tij,Tip1j) = K(Tij,Tip1j) + subs(k)*dx^2;
            K(Tij,Tijm1) = K(Tij,Tijm1) + subs(k)*dy^2;
            K(Tij,Tijp1) = K(Tij,Tijp1) + subs(k)*dy^2;
        else
            K(Tij,Tij  ) = K(Tij,Tij  ) - (subs(k)*2)/dx^2 + subs(c);
            K(Tij,Tijm1) = K(Tij,Tijm1) + subs(k)/dx^2;
            K(Tij,Tijp1) = K(Tij,Tijp1) + subs(k)/dx^2;
        end
        
    else
                    %si no lo estan una peor
        dx_m = vecX(j  ) - vecX(j-1);%atras
        dx_p = vecX(j+1) - vecX(j  );%adelante
        
        if( strcmp(Dimension, '2D') )
            Tip1j = WorldToMat(i+1,j  );
            Tim1j = WorldToMat(i-1,j  );

            dy_m = vecY(i  ) - vecY(i-1);%atras
            dy_p = vecY(i+1) - vecY(i  );%adelante

            K(Tij,Tij  ) = K(Tij,Tij  ) - subs(k)*2*(dy_m*dy_p) - subs(k)*2*(dx_m*dx_p) + subs(c)*(dx_m*dx_p)*(dy_m*dy_p);
            K(Tij,Tim1j) = K(Tij,Tim1j) + subs(k)/dy_m^2;
            K(Tij,Tip1j) = K(Tij,Tip1j) + subs(k)/dy_p^2;
            K(Tij,Tijm1) = K(Tij,Tijm1) + subs(k)/dx_m^2;
            K(Tij,Tijp1) = K(Tij,Tijp1) + subs(k)/dx_p^2;
        else
            K(Tij,Tij  ) = K(Tij,Tij  ) - (subs(k)*2)/(dx_m*dx_p) + subs(c);
            K(Tij,Tijm1) = K(Tij,Tijm1) + (subs(k)*2)/(dx_m*(dx_m+dx_p));
            K(Tij,Tijp1) = K(Tij,Tijp1) + (subs(k)*2)/(dx_p*(dx_m+dx_p));
        end
        
    end
    
end
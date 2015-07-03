% Ensambla la contribucion nodal de (l,m) a la matriz global
% solo para elementos (l,m) internos!!! Para CC usar las funciones de CC
function ContribucionNodal(l,m)

    global K;
    
    %Parametros de Ecuacion de Calor
    global Q;%Fuente de Calor
    global k;%Difusividad
    global c;%Perdida de calor al ambiente
    global Tamb;%Temperatura ambiente
    global h;
    
    global Tipocci; %tipo de cc izquierda
    global Tipoccd; %tipo de cc derecha
    global Tipoccin;%tipo de cc inferior
    global Tipoccs; %tipo de cc superior
    
    global Dimension;%Dimension del problema
    
    global minX;
    global maxX;
    global minY;
    global maxY;
    
    global N;%Funciones de forma
    global W;%Funciones de peso
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij

    %Tij   = WorldToMat(l  ,m  );
    
    syms symm syml x y;
    Wl = sym(subs(W,syml,l,0));
    Nm = subs(N,symm,m,0);
    
    dMixta = -Wl*h*Nm;
    
    if( strcmp(Dimension, '2D') )
        
        dInterno = -k*diff(Wl,x)*diff(Nm,x) - k*diff(Wl,y)*diff(Nm,y) + Wl*c*Nm;

        FH = GetFHandle2D(dInterno);
        K(l,m  ) = K(l,m  ) + dblquad(FH,minX,maxX,minY,maxY);%quad2d es otra opcion
%          if( isempty( symvar(dInterno) ) )
%              dInternoH = @(x,y) repmat(subs(dInterno), size(x));
%              K(l,m  ) = K(l,m  ) + dblquad(dInternoH,minX,maxX,minY,maxY);
%          else
%              dInternoH = matlabFunction(dInterno);
%              K(l,m  ) = K(l,m  ) + dblquad(dInternoH,minX,maxX,minY,maxY);%quad2d es otra opcion
%          end

        %Mixta
        if( strcmpi( Tipocci, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,x,minX,0));
            K(l,m  ) = K(l,m  ) + quad(FH,minY,maxY);
            %dMixtaH = matlabFunction(subs(dMixta,x,0));
            %K(l,m  ) = K(l,m  ) + quad(dMixtaH,minY,maxY);
        end
        if( strcmpi( Tipoccd, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,x,maxX,0));
            K(l,m  ) = K(l,m  ) + quad(FH,minY,maxY);
            %dMixtaH = matlabFunction(subs(dMixta,x,1));
            %K(l,m  ) = K(l,m  ) + quad(dMixtaH,minY,maxY);
        end
        if( strcmpi( Tipoccin, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,y,minY,0));
            K(l,m  ) = K(l,m  ) + quad(FH,minX,maxX);
            %dMixtaH = matlabFunction(subs(dMixta,y,0));
            %K(l,m  ) = K(l,m  ) + quad(dMixtaH,minX,maxX);
        end
        if( strcmpi( Tipoccs, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,y,maxY,0));
            K(l,m  ) = K(l,m  ) + quad(FH,minX,maxX);
            %dMixtaH = matlabFunction(subs(dMixta,y,1));
            %K(l,m  ) = K(l,m  ) + quad(dMixtaH,minX,maxX);
        end
        
        %Dirichlet
         if( strcmpi( Tipocci, 'Dirichlet' ) )
            dDirichlet = -Wl*Nm + k*Wl*diff(Nm,x)*-1;%subs(diff(Nm,x),{x},{-x},0)
            FH = GetFHandle1D(subs(dDirichlet,x,minX,0));
            %dDirichletH = matlabFunction(subs(dDirichlet,x,minX));
            K(l,m  ) = K(l,m  ) + quad(FH,minY,maxY);
        end
        if( strcmpi( Tipoccd, 'Dirichlet' ) )
            dDirichlet = -Wl*Nm + k*Wl*diff(Nm,x);
            FH = GetFHandle1D(subs(dDirichlet,x,maxX,0));
            %dDirichletH = matlabFunction(subs(dDirichlet,x,maxX));
            K(l,m  ) = K(l,m  ) + quad(FH,minY,maxY);
        end
        if( strcmpi( Tipoccin, 'Dirichlet' ) )
            dDirichlet = -Wl*Nm + k*Wl*diff(Nm,y)*-1;%subs(diff(Nm,y),{y},{-y},0)
            FH = GetFHandle1D(subs(dDirichlet,y,minY,0));
            %dDirichletH = matlabFunction(subs(dDirichlet,y,minY));
            K(l,m  ) = K(l,m  ) + quad(FH,minX,maxX);
        end
        if( strcmpi( Tipoccs, 'Dirichlet' ) )
            dDirichlet = -Wl*Nm + k*Wl*diff(Nm,y);
            FH = GetFHandle1D(subs(dDirichlet,y,maxY,0));
            %dDirichletH = matlabFunction(subs(dDirichlet,y,maxY));
            K(l,m  ) = K(l,m  ) + quad(FH,minX,maxX);
        end
        
    else
        
        dInterno = -k*diff(Wl,x)*diff(Nm,x) + Wl*c*Nm;
        
        FH = GetFHandle1D(dInterno);
        K(l,m  ) = K(l,m  ) + quad(FH,minX,maxX);
%         if( isempty( symvar(dInterno) ) )
%             dInternoH = @(x) repmat(subs(dInterno), size(x));
%             K(l,m  ) = K(l,m  ) + quad(dInternoH,minX,maxX);
%         else
%             dInternoH = matlabFunction(dInterno);
%             K(l,m  ) = K(l,m  ) + quad(dInternoH,minX,maxX);
%         end
        
        %Mixta
        if( strcmpi( Tipocci, 'Mixta' ) )
            K(l,m  ) = K(l,m  ) + subs(dMixta,x,minX,0);
        end
        if( strcmpi( Tipoccd, 'Mixta' ) )
            K(l,m  ) = K(l,m  ) + subs(dMixta,x,maxX,0);
        end
        
        %Dirichlet
         if( strcmpi( Tipocci, 'Dirichlet' ) )
            dDirichlet = -Wl*Nm + k*Wl*diff(Nm,x);
            K(l,m  ) = K(l,m  ) + subs(dDirichlet,x,minX,0);
        end
        if( strcmpi( Tipoccd, 'Dirichlet' ) )
            dDirichlet = -Wl*Nm + k*Wl*diff(Nm,x);
            K(l,m  ) = K(l,m  ) + subs(dDirichlet,x,maxX,0);
        end
        
    end
    
end
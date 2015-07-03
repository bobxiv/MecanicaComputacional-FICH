% Ensambla la contribucion nodal de (i,j) al vector independiente
% solo para elementos (i,j) internos!!! Para CC usar las funciones de CC
function ContribucionIndependiente(l)

    global b;
    
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
    
    global qcci;%conductividad del borde izquierdo
    global qccd;%conductividad del borde derecho
    global qccin;%conductividad del borde inferior
    global qccs;%conductividad del borde superior
    
    global Tcci;%dirichlet del borde izquierdo
    global Tccd;%dirichlet del borde derecho
    global Tccin;%dirichlet del borde inferior
    global Tccs;%dirichlet del borde superior
    
    global Dimension;%Dimension del problema
    
    global minX;
    global maxX;
    global minY;
    global maxY;
    
    global N;%Funciones de forma
    global W;%Funciones de peso
    
    % Matriz K se accede ( indice de ecuacion a modificar, indice de termino en la ecuacion )
    % Un nodo (i,j) en el mundo se corresponde con un (solo) indice Tij
    
    %Tij   = WorldToMat(l  ,m);
    
    syms symm syml x y;
    Wl = subs(W,syml,l);
    dInterno = -Wl*(-c*Tamb+Q);
    
    %dNeumann = Wl*q; pero cada CC tiene su q propia!
    
    dMixta = -Wl*h*Tamb;
    
    if( strcmp(Dimension, '2D') )
        
        %estos if isempty son para el caso donde el simbolico es una cte
        %y quad se rompe si lo intenta evaluar
        
        FH = GetFHandle2D(dInterno);
        b(l) = b(l) + dblquad(FH,minX,maxX,minY,maxY);%quad2d es otra opcion
        %dInternoH = @(valX,valY) subs(dInterno,{x,y},{valX,valY});
        %b(l) = b(l) + dblquad(dInternoH,minX,maxX,minY,maxY);%quad2d es otra opcion
        
        %Termino en CC Mixta
        if( strcmpi( Tipocci, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,x,minX,0));
            b(l) = b(l) + quad(FH,minY,maxY);
            %dMixtaH = matlabFunction(subs(dMixta,x,0));
            %b(l) = b(l) + quad(dMixtaH,minY,maxY);
        end
        if( strcmpi( Tipoccd, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,x,maxX,0));
            b(l) = b(l) + quad(FH,minY,maxY);
            %dMixtaH = matlabFunction(subs(dMixta,x,1));
            %b(l) = b(l) + quad(dMixtaH,minY,maxY);
        end
        if( strcmpi( Tipoccin, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,y,minY,0));
            b(l) = b(l) + quad(FH,minX,maxX);
            %dMixtaH = matlabFunction(subs(dMixta,y,0));
            %b(l) = b(l) + quad(dMixtaH,minX,maxX);
        end
        if( strcmpi( Tipoccs, 'Mixta' ) )
            FH = GetFHandle1D(subs(dMixta,y,maxY,0));
            b(l) = b(l) + quad(FH,minX,maxX);
            %dMixtaH = matlabFunction(subs(dMixta,y,1));
            %b(l) = b(l) + quad(dMixtaH,minX,maxX);
        end
        
        %Termino en CC Neumann
        if( strcmpi( Tipocci, 'Neumann' ) )
            dNeumann = Wl*qcci;
            FH = GetFHandle1D(subs(dNeumann,x,minX,0));
            b(l) = b(l) + quad(FH,minY,maxY);
%             if( isempty( symvar(dNeumann) ) )
%                 dNeumannH = @(x) repmat(subs(qcci), size(x));
%                 b(l) = b(l) + quad(dNeumannH,minY,maxY);
%             else
%                 dNeumannH = matlabFunction(subs(dNeumann,x,0));%aca suponemos que posee solo x
%                 %dNeumannH = matlabFunction(dNeumann);
%                 b(l) = b(l) + quad(dNeumannH,minY,maxY);
%             end
        end
        if( strcmpi( Tipoccd, 'Neumann' ) )
            dNeumann = Wl*qccd;
            FH = GetFHandle1D(subs(dNeumann,x,maxX,0));
            b(l) = b(l) + quad(FH,minY,maxY);    
%             if( isempty( symvar(dNeumann) ) )
%                 dNeumannH = @(x) repmat(subs(qccd), size(x));
%                 b(l) = b(l) + quad(dNeumannH,minY,maxY);
%             else
%                 dNeumannH = matlabFunction(subs(dNeumann,x,1));%aca suponemos que posee solo x
%                 %dNeumannH = matlabFunction(dNeumann);
%                 b(l) = b(l) + quad(dNeumannH,minY,maxY);
%             end
        end
        if( strcmpi( Tipoccin, 'Neumann' ) )
            dNeumann = Wl*qccin;
            FH = GetFHandle1D(subs(dNeumann,y,minY,0));
            b(l) = b(l) + quad(FH,minX,maxX);    
%             if( isempty( symvar(dNeumann) ) )
%                 dNeumannH = @(x) repmat(subs(qccin), size(x));
%                 b(l) = b(l) + quad(dNeumannH,minX,maxX);
%             else
%                 dNeumannH = matlabFunction(subs(dNeumann,y,0));%aca suponemos que posee solo y
%                 %dNeumannH = matlabFunction(dNeumann);
%                 b(l) = b(l) + quad(dNeumannH,minX,maxX);
%             end
        end
        if( strcmpi( Tipoccs, 'Neumann' ) )
            dNeumann = Wl*qccs;
            FH = GetFHandle1D(subs(dNeumann,y,maxY,0));
            b(l) = b(l) + quad(FH,minX,maxX);   
%             if( isempty( symvar(dNeumann) ) )
%                 dNeumannH = @(x) repmat(subs(qccs), size(x));
%                 b(l) = b(l) + quad(dNeumannH,minX,maxX);
%             else
%                 dNeumannH = matlabFunction(subs(dNeumann,y,1));%aca suponemos que posee solo y
%                 %dNeumannH = matlabFunction(dNeumann);
%                 b(l) = b(l) + quad(dNeumannH,minX,maxX);
%             end
        end
        
        %Termino en CC Dirichlet
        if( strcmpi( Tipocci, 'Dirichlet' ) )
            dDirichlet = -Wl*Tcci;
            FH = GetFHandle1D(subs(dDirichlet,{x},{minX},0));
            b(l) = b(l) + quad(FH,minY,maxY);            
%             if( isempty( symvar(dDirichlet) ) )
%                 dDirichletH = @(x) repmat(subs(dDirichlet), size(x));
%                 b(l) = b(l) + quad(dDirichletH,minY,maxY);
%             else
%                 dDirichletH = matlabFunction(subs(dDirichlet,x,minX));
%                 b(l) = b(l) + quad(dDirichletH,minY,maxY);
%             end
        end
        if( strcmpi( Tipoccd, 'Dirichlet' ) )
            dDirichlet = -Wl*Tccd;
            FH = GetFHandle1D(subs(dDirichlet,x,maxX,0));
            b(l) = b(l) + quad(FH,minY,maxY);
%             if( isempty( symvar(dDirichlet) ) )
%                 dDirichletH = @(x) repmat(subs(dDirichlet), size(x));
%                 b(l) = b(l) + quad(dDirichletH,minY,maxY);
%             else
%                 dDirichletH = matlabFunction(subs(dDirichlet,x,maxX));
%                 b(l) = b(l) + quad(dDirichletH,minY,maxY);
%             end
        end
        if( strcmpi( Tipoccin, 'Dirichlet' ) )
            dDirichlet = -Wl*Tccin;
            FH = GetFHandle1D(subs(dDirichlet,y,minY,0));
            b(l) = b(l) + quad(FH,minX,maxX);
%             if( isempty( symvar(dDirichlet) ) )
%                 dDirichletH = @(x) repmat(subs(dDirichlet), size(x));
%                 b(l) = b(l) + quad(dDirichletH,minX,maxX);
%             else
%                 dDirichletH = matlabFunction(subs(dDirichlet,y,minY));
%                 b(l) = b(l) + quad(dDirichletH,minX,maxX);
%             end
        end
        if( strcmpi( Tipoccs, 'Dirichlet' ) )
            dDirichlet = -Wl*Tccs;
            FH = GetFHandle1D(subs(dDirichlet,y,maxY,0));
            b(l) = b(l) + quad(FH,minX,maxX);
%             if( isempty( symvar(dDirichlet) ) )
%                 dDirichletH = @(x) repmat(subs(dDirichlet), size(x));
%                 b(l) = b(l) + quad(dDirichletH,minX,maxX);
%             else
%                 dDirichletH = matlabFunction(subs(dDirichlet,y,maxY));
%                 b(l) = b(l) + quad(dDirichletH,minX,maxX);
%             end
        end
        
    else
        
        %estos if isempty son para el caso donde el simbolico es una cte
        %y quad se rompe si lo intenta evaluar
        
        FH = GetFHandle1D(dInterno);
        b(l) = b(l) + quad(FH,minX,maxX);
%         if( isempty( symvar(dInterno) ) )%si c*Tamb es 0 o Q es 0 se da
%             dInternoH = @(x) repmat(subs(dInterno), size(x));
%             b(l) = b(l) + quad(dInternoH,minX,maxX);
%         else
%             dInternoH = matlabFunction(dInterno);
%             b(l) = b(l) + quad(dInternoH,minX,maxX);
%         end
        
        %Termino en CC Neumann
        if( strcmpi( Tipocci, 'Neumann' ) )
            dNeumann = Wl*qcci;
            if( isempty( symvar(dNeumann) ) )
                b(l) = b(l) + subs(dNeumann);
            else
                b(l) = b(l) + subs(dNeumann,minX);
            end
        end
        if( strcmpi( Tipoccd, 'Neumann' ) )
            dNeumann = Wl*qccd;
            if( isempty( symvar(dNeumann) ) )
                b(l) = b(l) + subs(dNeumann);
            else
                b(l) = b(l) + subs(dNeumann,maxX);
            end
        end
        
        %Termino en CC Dirichlet
        if( strcmpi( Tipocci, 'Dirichlet' ) )
            dDirichlet = -Wl*Tcci;
            if( isempty( symvar(dDirichlet) ) )
                b(l) = b(l) + subs(dDirichlet);
            else
                b(l) = b(l) + subs(dDirichlet,minY);
            end
        end
        if( strcmpi( Tipoccd, 'Dirichlet' ) )
            dDirichlet = -Wl*Tccd;
            if( isempty( symvar(dDirichlet) ) )
                b(l) = b(l) + subs(dDirichlet);
            else
                b(l) = b(l) + subs(dDirichlet,maxX);
            end
        end
        
    end
    
end


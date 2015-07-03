%Genera una malla estructurada en el entorno x en [0,Lx],
%con Nx nodos en el espacio Lx
%
%Parametros:
%           Nx numero de nodos en el eje x
%           Lx largo del dominio en x sobre el que se tiraran puntos
function [malla] = GeneradorMallaUnidimensional(Nx, Lx)
	
	malla.Nn=Nx;              %Numero de nodos en total
	malla.Ne=(Nx-1);          %Numero de elementos en total
	malla.n=zeros(malla.Nn,1);%Vector de Nodos
	malla.e=zeros(malla.Ne,2);%Matriz de elementos
                              %Cada fila tiene los indices de 1 elemento
                              %RECORDAR QUE USAMOS ELEMENTOS LINEALES
                              %CON 2 NODOS
	malla.Q=zeros(malla.Nn,1);%Vector de la fuente de temperatura 
                              %en cada nodo
	malla.phi=0*ones(malla.Nn,1);%coeficientes de aproximacion
    
    %Posicionar cada nodo en el espacio
    malla.n(:, 1) = 0:Lx/(malla.Nn-1):Lx;
    
	for i=1:(malla.Nn-1)
        malla.e(i,1) = i;
        malla.e(i,2) = i+1;
    end
    
    malla.h = Lx/(malla.Nn-1);
    
end

	


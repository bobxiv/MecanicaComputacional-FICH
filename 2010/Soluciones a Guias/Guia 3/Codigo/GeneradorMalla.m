%Genera una malla estructurada en el entorno x en [0,Lx] y y en [0 Ly],
%con Nx nodos en el espacio Lx y Ny nodos en el espacio Ly
%
%Parametros:
%           Nx numero de nodos en el eje x
%           Ny numero de nodos en el eje y
%           Lx largo del dominio en x sobre el que se tiraran puntos
%           Ly largo del dominio en y sobre el que se tiraran puntos
function [malla] = GeneradorMalla(Nx, Ny, Lx, Ly)
	dx=Lx/(Nx-1);
	dy=Ly/(Ny-1);
	
	malla.Nn=Nx*Ny;           %Numero de nodos en total
	malla.Nx=Nx;              %Numero de nodos en x
	malla.Ny=Ny;              %Numero de nodos en y
	malla.Ne=(Nx-1)*(Ny-1)*2; %Numero de elementos en total
	malla.n=zeros(malla.Nn,2);%Vector de Nodos
                              %Cada fila tiene primero la posicion en x
                              %y segundo la posicion en y de un nodo
	malla.e=zeros(malla.Ne,3);%Matriz de elementos
                              %Cada fila tiene los indices de 1 elemento
                              %RECORDAR QUE USAMOS ELEMENTOS TRIANGULARES
	malla.Q=zeros(malla.Nn,1);%Vector de la fuente de temperatura 
                              %en cada nodo
	malla.phi=0*ones(malla.Nn,1);%coeficientes de aproximacion

    %Posicionar cada nodo en el espacio
	Cn=1;
    for j=0:(Ny-1)
        for i=0:(Nx-1)
            malla.n(Cn, :) = [ i*dx, j*dy ];
            Cn = Cn + 1;
        end
    end
    
	Ce=1;
	for i=0:(Ny-2)
		n0=i*Nx+1;
		n1=n0+1;
		for j=1:(Nx-1)
            malla.e(Ce, :) = [n0,n1,n1+Nx];
            Ce = Ce + 1;

            malla.e(Ce, :) = [n0,n1+Nx,n0+Nx];
            Ce = Ce + 1;
            
            n0 = n0 + 1;
            n1 = n1 + 1;
        end
    end
    
    %Establecer la temperatura en todo nodo igual a 0
	for i=1:malla.Nn
		malla.Q(i)=0;
    end
    
end
	
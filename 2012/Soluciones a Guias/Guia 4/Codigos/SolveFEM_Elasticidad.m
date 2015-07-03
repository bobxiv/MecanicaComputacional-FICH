% Resuelve la ecuacion de elasticidad estacionaria 
%   | dsig_x/dx + dtau/dy + F_x = 0
%   | dtau_x/dx + dsig/dy + F_y = 0
%                                          EN 2D
% Utiliza:
%       * Elementos triangulares lineales genericos
%       * Galerkin
%       * Tension Plana o Deformacion Plana
%
% En particular se considera una Fuerza F=[0,0]
% Se considera que solo hay CC Dirichlet y Newmann Natural en el resto
%
% Los parametros son:
%
%                 - ProbElasticidadViga   especifica algunas propiedades
%                 - NombreMallaCargar     nombre del archivo de malla
%                 - Plot        estructura con informacion de ploteo
%
%               Para ayuda ver el script de ejemplo MainScript.m
function [ a ] = SolveFEM_Elasticidad( ProbElasticidadViga, NombreMallaCargar, Plot )
    
    eval(NombreMallaCargar);

    %Parametros del Problema y Malla
    %-------------------------------
    global coordinates;%Buffer de nodos
    global elements;   %Buffer de indices
    global fixnodes;   %Indice nodo - Desplazamiento
    global sideload;   %Indices de nodos - Fuerza (no hay en este caso)
    
    global F;          %Fuente de Peso
    global nu;         %
    global E;          %Tension Plana
    global ro;         %Densidad
    global FlexDerecha;%Desplazamiento del lado derecho de la viga
    
    F           = ProbElasticidadViga.F;
    nu          = poiss;
    E           = young;
    ro          = denss;
    FlexDerecha = ProbElasticidadViga.FlexDerecha;
    
    %Creamos la matriz global y el vector global
    %-------------------------------------------
    
    % 2 DoF por nodo - Intercalados los desplazamientos en [ux, uy]
    global K;
    K = sparse(2*length(coordinates(:,1)),2*length(coordinates(:,1)));
    global b;
    b = zeros(2*length(coordinates(:,1)),1);
    
    %Matriz de Hook
    global D;
    %Tension Plana
    D = (E/(1-nu^2)) *  [ 1 , nu, 0;
                         nu, 1 , 0;
                         0 , 0 , (1-nu)/2 ];
  
     %Deformacion Plana
%     D = (E/((1+nu)*(1-2*nu))) * [ 1+nu, 1-nu, 0;
%                                   1-nu, 1+nu, 0;
%                                   0   , 0   ,(1-2*nu)/2 ];
    
    %Ensamblado de la matriz global para nodos internos sin CC
    %---------------------------------------------------------
    
    %Recorre por cada elemento y acumula su contribucion
    for l= 1:length(elements(:,1))
        if( length(elements(1,:)) == 3 )
            ContribucionElementalTriangulo(l);%Agrega la contribucion del elemento l a K
        elseif( length(elements(1,:)) == 4 )
            ContribucionElementalCuadrilatero(l);%Agrega la contribucion del elemento l a K
        end
    end
    
    %Ensamblado de las CC en matriz global y vector global
    %-----------------------------------------------------
    
    %Calculo de Condiciones Dirichlet
    for l= 1:length(fixnodes(:,1))
        ContribucionDirichlet(l);%Agrega la contribucion de CC Dirichlet l
    end
    
    %La CC Newmann no importan en este caso en particular porque son CC
    %naturales por lo que se anulan al debilitar
    %Calculo de Condiciones Newmann por cara
    %     for l= 1:length(sideload(:,1))
    %         ContribucionNewmannCara(l);%Agrega la contribucion de CC Newmann l
    %     end
    
    %Resolucion de la matriz global
    %------------------------------
    
    %Soluciona el sistema 
    %   K*a=b -> a = inv(K)*b
    a = K\b;
    
    %Plotear si se lo pide
    %---------------------
    
    if( Plot.hacer )
        %dibujar:
        quiver(coordinates(:,1),coordinates(:,2),a(1:2:length(a)),a(2:2:length(a)));
        
        hold on;
        color = zeros(length(coordinates(:,1)),1);
        for k=1:length(color)
            color(k) = norm([a(2*k-1), a(2*k)]);
        end
        trisurf(elements,coordinates(:,1),coordinates(:,2), color);
        
        %hold on;
        trimesh(elements,coordinates(:,1)+a(1:2:length(a)),coordinates(:,2)+a(2:2:length(a)));
        
        xlabel('x');
        ylabel('y');
        zlabel('Desplazamiento');
            
        title (Plot.titulo);
    end

end
% Ensambla la contribucion elemental sin CC del elemento l en la matriz global
% para un elemento Triangular Lineal Generico
function ContribucionElementalTriangulo(l)

    global K;
    global D;
    
    %Parametros de Ecuacion
    global ro;         %Densidad
    
    %Malla
    global coordinates;%Buffer de nodos
    global elements;   %Buffer de indices
    
    nodos = coordinates(elements(l,:),:);%vector de los nodos del elemento
    
    % Componentes de los nodos (esta asignacion es por claridad nomas)
    xi = nodos(1,1);
    yi = nodos(1,2);
    xj = nodos(2,1);
    yj = nodos(2,2);
    xk = nodos(3,1);
    yk = nodos(3,2);
    
    % Calculo del area del elemento
    area = (1/2)* det( [xi, yi, 1; 
                        xj, yj, 1; 
                        xk, yk, 1] );
    
    % Matriz Bi
    %ai = (xj*yk-xk*yj)/(2*area);
    bi = (yj-yk)/(2*area);
    ci = (xk-xj)/(2*area);
    
    Bi =   [bi, 0; 
            0, ci; 
            ci, bi];
   % Matriz Bj
    %aj = (xk*yi-xi*yk)/(2*area);
    bj = (yk-yi)/(2*area);
    cj = (xi-xk)/(2*area);
    
    Bj =   [bj, 0; 
            0, cj; 
            cj, bj];
    % Matriz Bk
    %ak = (xi*yj-xj*yi)/(2*area);
    bk = (yi-yj)/(2*area);
    ck = (xj-xi)/(2*area);
    
    Bk =   [bk, 0; 
            0, ck; 
            ck, bk];
    
    Kelemental = area * ro *  [ Bi'*D*Bi , Bi'*D*Bj , Bi'*D*Bk;
                                Bj'*D*Bi , Bj'*D*Bj , Bj'*D*Bk;
                                Bk'*D*Bi , Bk'*D*Bj , Bk'*D*Bk ];
          
    iME = ones(6,1);
    iME([1,3,5,2,4,6]) = [2*elements(l,:)-1,2*elements(l,:)];
    
    K(iME,iME) = K(iME,iME) + Kelemental(:,:);
    
end
% Ensambla la contribucion elemental sin CC del elemento l en la matriz global
% para un elemento Cuadrangular (alineado con los ejes) Bilineal
function ContribucionElementalCuadrilatero(l)

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
    xq = nodos(4,1);
    yq = nodos(4,2);
    
    % Calculo los deltas
    %(en realidad no son necesarias tantas comprobaciones pero con 
    % esto me cubro de cualquier cosa rara del mallador)
    dx = max( [abs(xj-xi), abs(xk-xj), abs(xq-xk), abs(xk-xi), abs(xq-xi), abs(xj-xq)] );
    dy = max( [abs(yj-yi), abs(yk-yj), abs(yq-yk), abs(yk-yi), abs(yq-yi), abs(yj-yq)] );
    
    syms x;
    syms y;
    
    % Matriz Bi
    %ai = (xi+dx)*(yi+dy)/(dx*dy);
    bi = -(yi+dy)/(dx*dy);
    ci = -(xi+dx)/(dx*dy);
    di = 1/(dx*dy);

    Bi =   [bi+di*y, 0; 
            0      , ci+di*x; 
            ci+di*x, bi+di*y];
   % Matriz Bj
    %aj = -xi*(yi+dy)/(dx*dy);
    bj = (yi+dy)/(dx*dy);
    cj = xi/(dx*dy);
    dj = -1/(dx*dy);
    
    Bj =   [bj+dj*y, 0; 
            0      , cj+dj*x; 
            cj+dj*x, bj+dj*y];
    % Matriz Bk
    %ak = xi*yi/(dx*dy);
    bk = -yi/(dx*dy);
    ck = -xi/(dx*dy);
    dk = 1/(dx*dy);
    
    Bk =   [bk+dk*y, 0; 
            0      , ck+dk*x; 
            ck+dk*x, bk+dk*y];
        
    % Matriz Bq
    %aq = -(xi+dx)*yi/(dx*dy);
    bq = yi/(dx*dy);
    cq = (xi+dx)/(dx*dy);
    dq = -1/(dx*dy);
    
    Bq =   [bq+dq*y, 0; 
            0      , cq+dq*x; 
            cq+dq*x, bq+dq*y];
    
    Kelemental = double(int( int( ro *  [ Bi'*D*Bi , Bi'*D*Bj , Bi'*D*Bk , Bi'*D*Bq;
                                          Bj'*D*Bi , Bj'*D*Bj , Bj'*D*Bk , Bj'*D*Bq;
                                          Bk'*D*Bi , Bk'*D*Bj , Bk'*D*Bk , Bk'*D*Bq;
                                          Bq'*D*Bi , Bq'*D*Bj , Bq'*D*Bk , Bq'*D*Bq], x,xi,xi+dx), y,yi,yi+dy));

    iME = ones(8,1);
    iME([1,3,5,7,2,4,6,8]) = [2*elements(l,:)-1,2*elements(l,:)];
    
    K(iME,iME) = K(iME,iME) + Kelemental(:,:);
    
end

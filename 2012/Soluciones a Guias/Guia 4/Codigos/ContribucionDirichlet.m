% Ensambla la condicion Dirichlet l (nodal)
function ContribucionDirichlet(l)

    global K;
    global b;
    
    %CC Dirichlet
    global fixnodes;
    
    %Malla
    global coordinates;%Buffer de nodos
    global elements;   %Buffer de indices
 
    %vector con los indices de los nodos: los ux y luego los uy
    %iME = [2*fixnodes(l,1)-1, 2*fixnodes(l,1)];%indice Matriz Elemento
    
    %K(iME,:) = zeros(length(iME),2*length(coordinates(:,1)));
    %K(iME,iME) = ones(length(iME),length(iME));
    
    %b(iME) = fixnodes(l,3);
   
    index = 2*fixnodes(l,1)-1+fixnodes(l,2)-1;
    
    K(index,:) = zeros(1,2*length(coordinates(:,1)));
    K(index,index) = 1;
    
    b(index) = fixnodes(l,3);
end


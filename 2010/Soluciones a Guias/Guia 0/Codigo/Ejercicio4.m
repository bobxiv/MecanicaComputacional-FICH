%Nota si se corre testarea(modificado para Matlab) este llamara a este
%metodo para el calculo del area

%Calcula el area total de los triangulos que se pasan
%en los arreglos vertexs y adjacency
function [areaTot] = Ejercicio4(vertexs, adjacency )

    %Datos de prueba:
    %vertexs = [
    %            1, 1, 0;
    %            -1, -1, 0;
    %            1, -1, 0;
    %            1, 1, 0;
    %            0.5, 1, 0;
    %            -1, 1, 0
    %          ];

    %adjacency =[
    %            1, 2, 3;
    %            1, 3, 4;
    %            1, 4, 5;
    %            1, 5, 6;
    %            1, 6, 2
    %          ];
      
   areaTot = 0.0;
   for i = 1:length(adjacency)
    [P1,P2,P3]    = GetPrmitive(i,vertexs,adjacency);
    rectangleArea = norm( cross(P1-P2,P2-P3) );
    areaTot       = areaTot + rectangleArea / 2.0;
   end

end


%Obtiene la primitiva(triangulo) i de los datos formados 
%por vertexs y adjacency
function [P1,P2,P3] = GetPrmitive(i, vertexs, adjacency)
    P1 = vertexs(adjacency(i,1),:);
    P2 = vertexs(adjacency(i,2),:);
    P3 = vertexs(adjacency(i,3),:);
    return;
end

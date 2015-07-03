%Calcula el area total de los triangulos que se pasan
%en los arreglos vertexs y adjacency
function [area, vol] = Ejercicio5( x )

    %Datos de prueba:
    %x = [
    %      1,  1, 0;
    %      0,  1, 0;
    %      0,  0, 0;
    %      0,  0, 1
    %               ];
      
   %Calculo del area:
   area = [0.0 0.0 0.0;
           0.0 0.0 0.0;
           0.0 0.0 0.0;
           0.0 0.0 0.0];
   for i = 1:4
    [P1,P2,P3]    = GetPrmitive(i,x);
    area(i,:)     = cross(P1-P2,P2-P3);
    rectangleArea = norm( area(i,:) );
    area(i,:)     = area(i,:) / rectangleArea;%no es nesesario normalizar y 
                                              %luego ajustar pero asi se lee mas facil
    area(i,:)     = area(i,:) * (rectangleArea / 2.0);
   end
   
   %Calculo del volumen:
   vol = 0.0;
   
   [P1,P2,P3]    = GetPrmitive(1,x);
   P4            = x(4,:);
   
   PlaneNormal = cross(P1-P2,P2-P3);
   PlaneNormalUnit = PlaneNormal / norm(PlaneNormal);
   
   %Distancia al punto mas cercano en el plano [P1,P2,P3] desde P4
   dist = abs( PlaneNormalUnit * P4' - PlaneNormal * P1' );
   
   %Porque en un cubo hay 4 tetrahedros y en medio cubo hay 2
   %ademas la normal da el area de un rectangulo
   vol = ( (norm(PlaneNormal) / 2.0)* dist ) / 2.0;

end


%Obtiene la primitiva(triangulo) i de los datos formados 
%por vertexs y adjacency
function [P1,P2,P3] = GetPrmitive(i, x)
    P1 = x(i,:);
    if( (i+1) > 4 )
        P2 = x( mod(i+1,4)+1,:);
    else
        P2 = x( i+1,:);
    end
    if( (i+2) > 4 )
        P3 = x( mod(i+2,4)+1,:);
    else
        P3 = x( i+2,:);
    end
        
    return;
end
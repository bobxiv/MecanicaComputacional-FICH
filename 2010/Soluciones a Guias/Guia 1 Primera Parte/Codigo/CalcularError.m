%Calcula el error como la norma de 2 de la diferencia entre la solucion
%analitica y el resultado de los calculos
function [e] = ErrorNorma2(Tanalitic, T)
    
    Errors = Tanalitic - T;
    
    e = norm(Errors);
    
end
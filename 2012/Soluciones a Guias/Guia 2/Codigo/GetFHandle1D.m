function [ FH ] = GetFHandle1D( simbolico )
    syms x y;
    if( isempty( symvar(simbolico) ) )
        FH = @(x) repmat(subs(simbolico), size(x)); 
    else
        %se supone que solo puede tener ahora 1 variable
        %veamos cual es
        aux = symvar(simbolico);
        FH = @(valX) subs(simbolico,{aux(1)},{valX});
    end

end

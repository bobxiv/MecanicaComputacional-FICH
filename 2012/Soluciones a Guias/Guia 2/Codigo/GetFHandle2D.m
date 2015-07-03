function [ FH ] = GetFHandle2D( simbolico )
    syms x y;
    if( isempty( symvar(simbolico) ) )
        FH = @(x,y) repmat(subs(simbolico), size(x)); 
    else
        FH = @(valX,valY) subs(simbolico,{x,y},{valX,valY});
    end

end


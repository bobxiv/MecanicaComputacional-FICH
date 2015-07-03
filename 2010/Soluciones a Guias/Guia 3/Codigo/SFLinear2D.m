% a + b*x + c*y
% Coef = (a,b,c)
% Pos  = (x,y)
function [ res ] = SFLinear2D(a, b, c, x, y)

    %res = Coef(1) + Coef(2)*Pos(1) + Coef(3)*Pos(2);
    res = a + b*x + c*y;

end


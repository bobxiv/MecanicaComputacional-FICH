function [ e ] = ErrorCuadraticoMedio(Tanalitic, T)

    e = 0;
    for i = 1:length(T)
        e = e + (Tanalitic(i)-T(i))^2;
    end
    e = e / length(T);

end


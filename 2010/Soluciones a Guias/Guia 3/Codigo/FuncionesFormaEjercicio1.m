function [ ] = FuncionesFormaEjercicio1( h )

M = [ (-h/2)^2 (-h/2) 1 ;
      0        0      1 ;
      (h/2)^2  (h/2)  1];
  
f_i = [1;
       0;
       0];
f_j = [0;
       1;
       0];
f_k = [0;
       0;
       1];

c_i = M\f_i;
c_j = M\f_j;
c_k = M\f_k;

x=[ -h/2: h/500 :h/2];

N_i = c_i(1) .* x.^2 + c_i(2) .* x + c_i(3);
plot(x,N_i,'y');
hold on;

N_j = c_j(1) .* x.^2 + c_j(2) .* x + c_j(3);
plot(x,N_j,'r');
hold on;

N_k = c_k(1) .* x.^2 + c_k(2) .* x + c_k(3);
plot(x,N_k,'g');
hold on;

plot(x,zeros(length(x),1),'b');

xlabel('x');
legend('N_i', 'N_j', 'N_k');

end


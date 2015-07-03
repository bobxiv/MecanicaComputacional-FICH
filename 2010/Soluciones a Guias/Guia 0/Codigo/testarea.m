%$Id: testarea.m,v 1.2 2003/08/22 12:22:45 mstorti Exp $
%% Permite verificar si ls funcion `area' esta bien
%% Genera una malla de N*N quads en un plano en el espacio 
%% z=1, 0<=x,y,<=1 y lo divide en 2*N*N triangulos. Despues lo proyecta
%% sobre la esfera, fomando un sector de 1/24 de la superficie de la
%% esfera. 
N=100;

%% genera los quads sobre la superficie
xx = (0:N)'/N;
xx = xx(:,ones(N+1,1));
y = xx';
z = ones(size(xx));

xx = xx(:);
y  = y(:);
z  = z(:);

%% x es un vector de (N+1)^2 x 3
x = [xx, y, z];

%% genera las conectividades (grid)
%% Este es el primer panel cuadrangular. Los otros se
%% generan por repeticion
q1 = [1 2 N+3 N+2];

%% Va agregando paneles por repeticion
gl1 = [];
for k=1:N
  gl1 = [gl1;
	 q1+k-1];
end

%% A esta altura gl1 contiene una tira de quads.
%% Genera el resto por repeticion
grid = [];
for k=1:N
  grid = [grid;
	  gl1+(k-1)*(N+1)];
end

%% Divide cada quad en dos tri
grid=[grid(:,[1 2 3]);
      grid(:,[3 4 1])];

%% Verifica el calculo del area de un cuadrado en el espacio
%A=Area(x,grid);
A=Ejercicio4(x,grid);
display('Area del cuadrado, deberia dar 1: ');
display(A);

%% Proyecta sobre la esfera
r = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
for k=1:3
  x(:,k) = x(:,k)./r;
end

%% Verifica cuanto da para la esfera
%A=Area(x,grid);
A=Ejercicio4(x,grid);
display('Area 1/24 de la superficie de la esfera, A = pi/6 = 0.5236 : ');
display(A);
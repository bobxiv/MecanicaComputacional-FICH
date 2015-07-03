function [] = Ejercicio6( )

%malla = GeneradorMalla(60, 60, 1, 1);
malla = GeneradorMalla(8, 8, 1, 1);

%Seteo de Ctes
malla.k = 1.0;

for m=1:malla.Nn
    malla.Q(m) = 0.0;
end

%Estado inicial en el tiempo
a = zeros(malla.Nn,1);
for m=1:malla.Nn
    a(m) = 0.0;
end

a = SetaBordersToDirichlet(malla, a);
    
malla.alpha = -1.0;

dx=1/(8-1);
%nu = 1e-8;
%Fo1 = nu*dt/dx^2
%Fo2 = nu*dt/dy^2

%dt = dx^2/(4*malla.k);
dt=0.00001;          %Paso en el tiempo
iterations = 2*90000;%En el tiempo

%Seteo de Condiciones Iniciales
%malla.AbajoCI = 'Newmann';
malla.AbajoCI = 'Dirichlet';
malla.Abajoq  = -15;
malla.AbajoPhi  = 100;
%malla.ArribaCI = 'Newmann';
malla.ArribaCI = 'Dirichlet';
malla.Arribaq = 0;
malla.ArribaPhi  = 100;

%malla.IzquierdaCI = 'Newmann';
malla.IzquierdaCI = 'Dirichlet';
malla.Izquierdaq  = -15;
malla.IzquierdaPhi  = 100;
%malla.DerechaCI = 'Newmann';
malla.DerechaCI = 'Dirichlet';
malla.Derechaq = 0;
malla.DerechaPhi = 100;

%Comienzo de resolucion

[elemCount,nodesIndexs] = size(malla.e);

%Calculo de las areas de cada elemento triangular
malla.area = zeros(elemCount,1);
for m=1:elemCount
    
    xi = malla.n( malla.e(m,1) , 1 );
    yi = malla.n( malla.e(m,1) , 2 );
    xj = malla.n( malla.e(m,2) , 1 );
    yj = malla.n( malla.e(m,2) , 2 );
    xk = malla.n( malla.e(m,3) , 1 );
    yk = malla.n( malla.e(m,3) , 2 );
    Pi = [xi,yi,0];
    Pj = [xj,yj,0];
    Pk = [xk,yk,0];
     
    area = cross(Pj-Pi,Pk-Pj);
    malla.area(m) = norm(area)/ 2.0;
end

%Calculo de coeficientes de cada Funcion de Forma
% a + b*x + c*y
malla.nCoef = zeros(elemCount,3,3);
for m=1:elemCount
    %1 es el i
    %2 es el j
    %3 es el k
     xi = malla.n( malla.e(m,1) , 1 );
     xj = malla.n( malla.e(m,2) , 1 );
     xk = malla.n( malla.e(m,3) , 1 );
     yi = malla.n( malla.e(m,1) , 2 );
     yj = malla.n( malla.e(m,2) , 2 );
     yk = malla.n( malla.e(m,3) , 2 );
     
     M = [1 xi yi;
          1 xj yj;
          1 xk yk];
      
     b = [1;
          0;
          0];
     
     coef = M\b;
     malla.nCoef(m,1,:) = coef;
     
     b = [0;
          1;
          0];
      
     coef = M\b;
     malla.nCoef(m,2,:) = coef;
       
      b = [0;
           0;
           1];
        
      coef = M\b;
      malla.nCoef(m,3,:) = coef;
end

%Matrices globales
K = zeros(malla.Nn,malla.Nn);
K = sparse(K);
C = zeros(malla.Nn,malla.Nn);
C = sparse(C);
f = zeros(malla.Nn,1);
a_next = zeros(malla.Nn,1);

%Calculo de la matriz elemental
%Y Ensamblado de la matriz global
for m=1:elemCount

    Coef = malla.nCoef(m,1,:);
    ai = Coef(1);
    bi = Coef(2);
    ci = Coef(3);
    Coef = malla.nCoef(m,2,:);
    aj = Coef(1);
    bj = Coef(2);
    cj = Coef(3);
    Coef = malla.nCoef(m,3,:);
    ak = Coef(1);
    bk = Coef(2);
    ck = Coef(3);
    
    K_e = [ (bi^2+ci^2)   ,(bi*bj+ci*cj), (bi*bk+ci*ck); %i o 1
            (bj*bi+cj*ci) ,(bj^2+cj^2)  , (bj*bk+cj*ck); %j o 2
            (bk*bi+ck*ci) ,(bk*bj+ck*cj), (bk^2+ck^2)  ];%k o 3
        
    K_e = K_e * malla.k * malla.area(m);
        
    %Se agrega la contribucion elemental a la matriz global
    K( malla.e(m,1) , malla.e(m,1) ) = K( malla.e(m,1) , malla.e(m,1) ) + K_e(1,1);
    K( malla.e(m,1) , malla.e(m,2) ) = K( malla.e(m,1) , malla.e(m,2) ) + K_e(1,2);
    K( malla.e(m,1) , malla.e(m,3) ) = K( malla.e(m,1) , malla.e(m,3) ) + K_e(1,3);
    
    K( malla.e(m,2) , malla.e(m,1)) = K( malla.e(m,2) , malla.e(m,1) ) + K_e(2,1);
    K( malla.e(m,2) , malla.e(m,2)) = K( malla.e(m,2) , malla.e(m,2) ) + K_e(2,2);
    K( malla.e(m,2) , malla.e(m,3)) = K( malla.e(m,2) , malla.e(m,3) ) + K_e(2,3);
    
    K( malla.e(m,3) , malla.e(m,1)) = K( malla.e(m,3) , malla.e(m,1)) + K_e(3,1);
    K( malla.e(m,3) , malla.e(m,2)) = K( malla.e(m,3) , malla.e(m,2)) + K_e(3,2);
    K( malla.e(m,3) , malla.e(m,3)) = K( malla.e(m,3) , malla.e(m,3)) + K_e(3,3);
    
    
    [X,Y,Wx,Wy]=triquad(8,[malla.n(malla.e(m,1),1) malla.n(malla.e(m,1),2); malla.n(malla.e(m,2),1) malla.n(malla.e(m,2),2); malla.n(malla.e(m,3),1) malla.n(malla.e(m,3),2)]);
    Ni=@(x,y) ai+bi*x+ci*y;
    Nj=@(x,y) aj+bj*x+cj*y;
    Nk=@(x,y) ak+bk*x+ck*y;
    
    C_e = [ (Wx'*feval(Ni,X,Y)*feval(Ni,X,Y)*Wy) (Wx'*feval(Ni,X,Y)*feval(Nj,X,Y)*Wy) (Wx'*feval(Ni,X,Y)*feval(Nk,X,Y)*Wy); %i
            (Wx'*feval(Nj,X,Y)*feval(Ni,X,Y)*Wy) (Wx'*feval(Nj,X,Y)*feval(Nj,X,Y)*Wy) (Wx'*feval(Nj,X,Y)*feval(Nk,X,Y)*Wy); %j
            (Wx'*feval(Nk,X,Y)*feval(Ni,X,Y)*Wy) (Wx'*feval(Nk,X,Y)*feval(Nj,X,Y)*Wy) (Wx'*feval(Nk,X,Y)*feval(Nk,X,Y)*Wy)];%k
    
    C_e = C_e * (-malla.alpha);
    
    %Se agrega la contribucion elemental a la matriz global
    C( malla.e(m,1) , malla.e(m,1) ) = C( malla.e(m,1) , malla.e(m,1) ) + C_e(1,1);
    C( malla.e(m,1) , malla.e(m,2) ) = C( malla.e(m,1) , malla.e(m,2) ) + C_e(1,2);
    C( malla.e(m,1) , malla.e(m,3) ) = C( malla.e(m,1) , malla.e(m,3) ) + C_e(1,3);
    
    C( malla.e(m,2) , malla.e(m,1)) = C( malla.e(m,2) , malla.e(m,1) ) + C_e(2,1);
    C( malla.e(m,2) , malla.e(m,2)) = C( malla.e(m,2) , malla.e(m,2) ) + C_e(2,2);
    C( malla.e(m,2) , malla.e(m,3)) = C( malla.e(m,2) , malla.e(m,3) ) + C_e(2,3);
    
    C( malla.e(m,3) , malla.e(m,1)) = C( malla.e(m,3) , malla.e(m,1)) + C_e(3,1);
    C( malla.e(m,3) , malla.e(m,2)) = C( malla.e(m,3) , malla.e(m,2)) + C_e(3,2);
    C( malla.e(m,3) , malla.e(m,3)) = C( malla.e(m,3) , malla.e(m,3)) + C_e(3,3);
    
    
    %Contribucion de fe con termino fuente
    Qprom = (malla.Q(malla.e(m,1)) + malla.Q(malla.e(m,2)) + malla.Q(malla.e(m,3))) / 3.0;
    f(malla.e(m,1)) = f(malla.e(m,1)) + 1/3.0 * malla.area(m) * Qprom;
    f(malla.e(m,2)) = f(malla.e(m,2)) + 1/3.0 * malla.area(m) * Qprom;
    f(malla.e(m,3)) = f(malla.e(m,3)) + 1/3.0 * malla.area(m) * Qprom;

end


%Contribucion de Condiciones de contorno a f

%Abajo CI
if strcmp(malla.AbajoCI , 'Dirichlet' ) == 1
    %Condicion Dirichlet abajo
    for m=1:malla.Nx
        aIndex = m;
        K( aIndex , : ) = zeros(1,malla.Nn);
        K( aIndex , aIndex ) = 1;
        
        C( aIndex , : ) = zeros(1,malla.Nn);
        C( aIndex , aIndex ) = 1;
        f(aIndex) = malla.AbajoPhi;
    end
else
    %Condicion Newman abajo
    for m=1:malla.Nx
        dist = ( (malla.n(m,1)-malla.n(m+1,1))^2 + (malla.n(m,2)-malla.n(m+1,2))^2 )^(1/2);
        f(m) = f(m) - 1/2 * malla.Abajoq * dist;
    end
end

%Arriba CI
if strcmp(malla.ArribaCI , 'Dirichlet' ) == 1
    %Condicion Dirichlet arriba
    for m=1:malla.Nx
        aIndex = malla.Nn-(m-1);
        K( aIndex , : ) = zeros(1,malla.Nn);
        K( aIndex , aIndex ) = 1;
        
        C( aIndex , : ) = zeros(1,malla.Nn);
        C( aIndex , aIndex ) = 1;
        f(aIndex) = malla.ArribaPhi;
    end
else
    %Condicion Newman arriba
    for m=1:malla.Nx
        aIndex = malla.Nn-(m-1);
        dist = ( (malla.n(aIndex-1,1)-malla.n(aIndex,1))^2 + (malla.n(aIndex-1,2)-malla.n(aIndex,2))^2 )^(1/2);
        f(aIndex) = f(aIndex) - 1/2 * malla.Arribaq * dist;
    end
end

%Izquierda CI
if strcmp(malla.IzquierdaCI , 'Dirichlet' ) == 1
    %Condicion Dirichlet izquierda
    for m=1:malla.Ny
        aIndex = (m-1)*malla.Nx+1;
        K( aIndex , : ) = zeros(1,malla.Nn);
        K( aIndex , aIndex ) = 1;
        
        C( aIndex , : ) = zeros(1,malla.Nn);
        C( aIndex , aIndex ) = 1;
        f(aIndex) = malla.IzquierdaPhi;
    end
else
    %Condicion Newman izquierda
    for m=1:malla.Ny-1
        aIndex = (m-1)*malla.Nx+1;
        dist = ( (malla.n(aIndex,1)-malla.n(aIndex+malla.Nx,1))^2 + (malla.n(aIndex,2)-malla.n(aIndex+malla.Nx,2))^2 )^(1/2);
        f(aIndex) = f(aIndex) - 1/2 * malla.Izquierdaq * dist;
    end
end

%Derecha CI
if strcmp(malla.DerechaCI , 'Dirichlet' ) == 1
    %Condicion Dirichlet derecha
    for m=1:malla.Ny
        aIndex = (m-1)*malla.Nx+malla.Nx;
        K( aIndex , : ) = zeros(1,malla.Nn);
        K( aIndex , aIndex ) = 1;
        
        C( aIndex , : ) = zeros(1,malla.Nn);
        C( aIndex , aIndex ) = 1;
        f(aIndex) = malla.DerechaPhi;
    end
else
   %Condicion Newman derecha
    for m=1:malla.Ny-1
        aIndex = (m-1)*malla.Nx+malla.Nx;
        dist = ( (malla.n(aIndex,1)-malla.n(aIndex+malla.Nx,1))^2 + (malla.n(aIndex,2)-malla.n(aIndex+malla.Nx,2))^2 )^(1/2);
        f(aIndex) = f(aIndex) - 1/2 * malla.Derechaq * dist;
    end
end

salir = false;
time = 0;
for i=1:iterations
    %Solucion del problema
    a_next = dt* inv(C) * ( (C/dt - K) * a + f);
    a = a_next;
    
    for i=1:malla.Nn
        if a(i) > 50
            salir = true;
        else
            salir = false;
            break
        end
    end
    if salir == true
        time = i*dt;
        break;
    end
    
end


[X,Y] = meshgrid(0:1.0/(malla.Nx-1):1.0,0:1.0/(malla.Ny-1):1.0);

for j=0:(malla.Ny-1)
    Z(j+1,:) = a(j*malla.Nx+1 : (j+1)*malla.Nx);
end

surf(X,Y,Z);
xlabel('x');
ylabel('y');
zlabel('T');
title('k*(d^2T/dx^2+d^2T/dy^2)+Q-alpha*dT/dt=0')

end


function [newa] = SetaBordersToDirichlet(malla,a)

newa = a;

for m=1:malla.Nx
    aIndex = m;
    newa(aIndex) = 100.0;
end
for m=1:malla.Nx
    aIndex = malla.Nn-(m-1);
    newa(aIndex) = 100.0;
end
for m=1:malla.Ny
    aIndex = (m-1)*malla.Nx+1;
    newa(aIndex) = 100.0;
end
for m=1:malla.Ny
    aIndex = (m-1)*malla.Nx+malla.Nx;
    newa(aIndex) = 100.0;
end

end

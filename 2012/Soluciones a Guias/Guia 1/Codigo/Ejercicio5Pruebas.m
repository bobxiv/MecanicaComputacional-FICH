% Ejercicio 5 Guia MDF 2012

% Parametros del problema a simular
%----------------------------------
Lx = 1;
Ly = 1;
Lt = 0.6;
dx = 0.1;
dy = 0.1;

% Importante para windows
set(gcf,'Renderer','zbuffer');

% Video de metodo explicito de MDF en el tiempo
%----------------------------------------------
dt = 0.002;
Fexplicito = ej5_movie('explicito', Lx, Ly, Lt, dx, dy, dt);

Vexplicito = VideoWriter('Videos/explicito.avi');

open(Vexplicito);
for k=1:length(Fexplicito)
    writeVideo(Vexplicito,Fexplicito(k));
end
close(Vexplicito);

disp('Se a guardado el video del metodo explicito');

% Video de metodo implicito de MDF en el tiempo
%----------------------------------------------
dt = 0.002;
Fimplicito = ej5_movie('implicito', Lx, Ly, Lt, dx, dy, dt);

Vimplicito = VideoWriter('Videos/implicito.avi');

open(Vimplicito);
for k=1:length(Fimplicito)
    writeVideo(Vimplicito,Fimplicito(k));
end
close(Vimplicito);

disp('Se a guardado el video del metodo implicito');

% Video de metodo Cranck-Nicholson de MDF en el tiempo
%-----------------------------------------------------
dt = 0.002;
Fcn = ej5_movie('CranckNicholson', Lx, Ly, Lt, dx, dy, dt);

Vcn        = VideoWriter('Videos/CranckNicholson.avi');

open(Vcn);
for k=1:length(Fcn)
    writeVideo(Vcn,Fcn(k));
end
close(Vcn);

disp('Se a guardado el video del metodo Cranck-Nicholson');
function [Lx, dx, dy, x, y, kx, ky, k2, Sx, Sy, gx, gy, tsup, Eo, M, lambda_r, lambda_m, La, teta, fx, fy, rr] = parametros_globais(Ly, nx, ny, d) %  Constantes for programa
    Lx = 4*Ly;      %  Comprimento Fisico na direcao x
    dx = Lx/nx;     %  Discretizacao de Lx
    dy = Ly/ny;     %  Discretizacao de Ly
    
    %  Pontos de colocacao da malha euleriana
    x = zeros(1, nx);
    for i = 1:nx, x(i) = (i-1)*dx; end
    y = zeros(1, ny);
    for j = 1:ny, y(j) = (j - fix(0.35*ny))*dy; end
    
    % Wave vector
    kx = zeros(1, nx);
    for i = 1:nx/2+1, kx(i) = 2.0*pi*(i - 1)/(nx*dx); end
    for i = nx/2+2:nx, kx(i) = 2.0*pi*(i-1-nx)/(nx*dx); end

    ky = zeros(1, ny);
    for j=1:ny/2+1, ky(j)=2.0*pi*(j-1)/(ny*dy); end
	for j=ny/2+2:ny, ky(j)=2*pi*(j-1-ny)/(ny*dy); end 

    k2 = zeros(nx,ny);  %  k2=kx^2 + ky^2
	for i = 1:nx
        for j = 1:ny, k2(i,j) = kx(i)^2+ky(j)^2 + eps/2; end 
    end 

    % ---------------------------------------------------------------------
    %  Dados de Entrada
	Sx = 1;         %  frequencia for GT
	Sy = 1;         %  frequencia for GT

	gx = -980.0;
	gy = 0;

	tsup = 1;
	Eo = 1;              %  Numero de Eotvos
	M = 0.01;            %  Numero de Morton
	lambda_r = 0.5;      %  razao entre as fases
	lambda_m = 0.5;
	La = 120.0d0;

	rr = d/2;
	rrteta = 20.0;
	teta = rr/rrteta;

    fx = zeros(nx, ny);
	fy = zeros(nx, ny);
end
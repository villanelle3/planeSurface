function [imin, jmax, iimin, iimax, jjmax, zp] = gera_geometria(nx, ny, Ly)

% Posicao das condicoes de contorno

    imin =  2;                   %  comeco da zona de buffer
    imax = nx;                   %  nao usado
    jmax = fix(0.35*ny);         % comeco do dominio util
    jmin = 1;                    % nao usado
    
    iimin = fix(0.3333*nx);      % fim do dominio complementar
    iimax = nx;                  %  fim do dominio util
    jjmin = fix(0.1*ny);         %  nao usado
    jjmax = jmax + fix(0.114/Ly * ny);       %  limite superior do dominio util
    zp  = fix(0.019/Ly * ny);                % fix(0.05*Lx/dx)% zona de forçagem 

end

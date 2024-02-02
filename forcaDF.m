function [fx, fy, u, v, ufis, vfis] = forcaDF(u, v, ufis, vfis, nx, ny, phi, Qtx, Qty, kx, ky, k2, nit, jmax, jjmax, dt, zp, rho, nstep, nsavec)
    % Complex arrays
    F0x = zeros(nx, ny);
    F0y = zeros(nx, ny);

    % Real arrays
    sfx_local = zeros(1, nx);   % Real arrays
    Cf = zeros(1, nx);
    Cf2 = zeros(1, nx);
    Cf3 = zeros(1, nx);
    Cf4 = zeros(1, nx);

    % Scalars
    maxfx0 = 1;
    maxfy0 = 1;
    maxfx = 1;
    maxfy = 1;
    it = 0;
    epslon = 1e-4;

    % ---------------------------------------------------------------------
    fx = zeros(nx, ny);
    fy = zeros(nx, ny);

    % for ifi = 1:nfi
    %     front(ifi).SSFx = zeros(nx, ny);  % Initialize SSFx field
    %     front(ifi).SSFy = zeros(nx, ny);  % Initialize SSFy field
    % end
    % ---------------------------------------------------------------------
    % Zona de amortecimento (Buffer Zone)
    ZBx = phi .* (ufis - Qtx); % + phib * (ufis - Qtxb)
    ZBy = phi .* (vfis - Qty); % + phib * (vfis - Qtyb)
    
    ZBx = fft2(ZBx);
    ZBy = fft2(ZBy);

    [ZBx, ZBy] = projecao(ZBx, ZBy, nx, ny, kx, ky, k2); % Assuming projecao returns two matrices ux and uy
    u = u - ZBx;
    v = v - ZBy;

    while (maxfx/maxfx0 >= epslon || maxfy/maxfy0 >= epslon) && it < nit
        fxd = zeros(nx, ny);
        fyd = zeros(nx, ny);
        for i = 1:nx
            % Loop for the top surface of the plate
            for j = jmax:jmax 
                fxd(i,j) = real((0.0 - ufis(i,j)) / dt);
                fyd(i,j) = real((0.0 - vfis(i,j)) / dt);
                sfx_local(i) = sfx_local(i) + fxd(i,j);
            end
            
            % Loop for the bottom surface of the plate
            for j = jmax-zp:jmax-zp
                fxd(i,j) = real((0.0 - ufis(i,j)) / dt);
                fyd(i,j) = real((0.0 - vfis(i,j)) / dt);
            end
        end
        
        for i = 1:nx
            % Casca das condicoes de topo
            for j = jjmax:jjmax
                fyd(i,j) = real((vfis(i,jjmax-1)-vfis(i,j))/dt);
                fxd(i,j) = real((ufis(i,jjmax-1)-ufis(i,j))/dt);
            end
            for j = jjmax+zp:jjmax+zp
                fyd(i,j) = real((vfis(i,jjmax-1)-vfis(i,j))/dt);
                fxd(i,j) = real((ufis(i,jjmax-1)-ufis(i,j))/dt);
            end
        end

        
        maxfx = max(max(real(fxd) - real(F0x)));
        maxfy = max(max(real(fyd) - real(F0y)));
        
        if it == 1 && nstep >= 1
            maxfx0 = maxfx;
            maxfy0 = maxfy;
        end
        
        F0x = real(fxd);
        F0y = real(fyd);
        it = it + 1;

        % Impose Eulerian force on the estimated field
        fxd = fft2(fxd);
        fyd = fft2(fyd);
        fx = fx + fxd;
        fy = fy + fyd;
        
        [fxd, fyd] = projecao(fxd, fyd, nx, ny, kx, ky, k2); 
        
        u = u + dt / rho * fxd;
        v = v + dt / rho * fyd;
        
        ufis = u;
        vfis = v;
        
        ufis = real(ifft2(ufis));
        vfis = real(ifft2(vfis));
    end
    % Add ZBx and ZBy to fx and fy
    fx = fx + ZBx;
    fy = fy + ZBy;
    
    % Calculating spectral derivative
    dudye = ufis;
    
    dudye = fft2(dudye);  % 2-D FFT
    dudye = ifft2(1i * ky .* dudye);  % Multiply by 1i * ky for spectral derivative

    if mod(nstep, nsavec) == 0 && nstep > 1
    % Calculo dos coeficientes de atrito
        % Pela forca lagrangeana 
        for i = iimin+1:iimax-1
            Cf(i) = (-2.0 * sfx_local(i) * dy) / (rho * umax^2);
        end
        
        % Por aproximacao de diferencas finitas
        for i = iimin+1:iimax-1
            Cf2(i) = (2.0 * ni * rho * real(ufis(i,jmax+1) - ufis(i,jmax)) / dy) / (rho * umax^2);
        end
        
        % Por derivada espectral
        for i = iimin+1:iimax-1
            Cf3(i) = (2.0 * ni * rho * dudye(i,jmax)) / (rho * umax^2);
        end
        
        % Por derivada espectral
        for i = iimin+1:iimax-1
            Cf4(i) = (2.0 * ni * rho * dudye(i,jmax+1)) / (rho * umax^2);
        end
    end

end

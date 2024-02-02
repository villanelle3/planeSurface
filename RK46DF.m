function [ufis, vfis, u, v] = RK46DF(dt, ni, nx, ny, ufis, vfis, kx, ky, k2, u, v)
    alfa = [0.0, -0.691750960670, -1.727127405211, -0.694890150986, ...
        -1.039942756197, -1.531977447611];
    beta = [0.122, 0.477263056358, 0.381941220320, 0.447757195744, ...
        0.498614246822, 0.186648570846];
    c = [0.0, 0.122, 0.269115878630, 0.447717183551, 0.749979795490, ...
        0.898555413085];

    k1x = zeros(nx, ny);    k1y = zeros(nx, ny);

    for ii = 1:6
        viscx = ni .* k2 .* u; % Calculation of the viscous term in the x-direction
        viscy = ni .* k2 .* v; % Calculation of the viscous term in the y-direction
        ut = u;
        vt = v;
        utfis = ufis;
        vtfis = vfis;
        [tnlxe, tnlye] = tnlinear(ut, vt, utfis, vtfis, nx, ny, kx, ky);

        RHSx = -tnlxe - viscx; % + fontex; % Calculation of the right-hand side of the equation in x
        RHSy = -tnlye - viscy; % + fontey; % Calculation of the right-hand side of the equation in y
        
        [RHSx, RHSy] = projecao(RHSx,RHSy, nx, ny, kx, ky, k2);
        
        % Coefficients of RK4
        k1x = alfa(ii) * k1x + dt * RHSx;
        k1y = alfa(ii) * k1y + dt * RHSy;

        % Update velocity prediction
        u = u + k1x * beta(ii);
        v = v + k1y * beta(ii);
        
        ufis = u;
        vfis = v;
        ufis = real(ifft2(ufis));
        vfis = real(ifft2(vfis));
    end
end

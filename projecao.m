function [ax, ay] = projecao(ax, ay, nx, ny, kx, ky, k2)
    tmp = zeros(nx, ny);
    for j = 1:ny
        for i = 1:nx
            tmp(i,j) = (ax(i,j) * kx(i) + ay(i,j) * ky(j)) / k2(i,j);
            ax(i,j) = ax(i,j) - tmp(i,j) * kx(i);
            ay(i,j) = ay(i,j) - tmp(i,j) * ky(j);
        end
    end
end

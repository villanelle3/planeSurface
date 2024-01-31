function [u, v, p] = ci(t_ci, nx, ny, Sx, Sy, ni, t, umax, ucf, rr, teta, jmax, jjmax)
    u = zeros(nx, ny);
    v = zeros(nx, ny);
    p = zeros(nx, ny);

switch t_ci
    % ---------------------------------------------------------------------
    case 0
        return
    % ---------------------------------------------------------------------
    case 1
        u = umax.*ones(nx, ny);
    % ---------------------------------------------------------------------
    case 2
        for i = 1:nx
		    for j = 1:ny
			    u(i,j) = sin(Sx*x(i))*cos(Sy*y(j))*exp(-ni*t*(Sx^2+Sy^2));
			    v(i,j) = -Sx/Sy*sin(Sy*y(j))*cos(Sx*x(i))*...
                    exp(-ni*t*(Sx^2+Sy^2));
			    p(i,j) = ((((Sx^2)/(Sy^2))*cos(Sy*y(j)))^2-...
                    (sin(Sx*x(i)))^2)*0.5*exp(-2.0*ni*t*(Sx^2+Sy^2));
		    end 
        end
    % ---------------------------------------------------------------------
    case 3
        for i = 1:nx
		    for j = 1:ny
			    u(i,j)=sin(Sx*x(i))*cos(Sy*y(j))*cos(pi*ni*t*(Sx^2+Sy^2));
			    v(i,j)=-Sx/Sy*sin(Sy*y(j))*cos(Sx*x(i))*...
                    cos(ni*t*pi*(Sx^2+Sy^2));
			    p(i,j)=0.5*cos(2.0*pi*ni*t)^2*(cos(x(i))^2+cos(y(j))^2);
            end
        end
    % ---------------------------------------------------------------------
    case 4  % jato espacial
        for i = 1:nx
		    for j = 1:ny
			    u(i, j) = (umax + ucf)/2.0 - (umax-ucf)/2.0*...
                        tanh(.25*rr/teta*((abs(y(j))/rr)-(rr/abs(y(j)))));
		    end 
        end 
        for j = 2:fix(ny/2.0d0)
	        if (real(u(1,j) - u(1,j-1)) > (ucf)), jmin = j-1; end 
        end 
        
        for j=int(ny/2.0d0):ny
	        if (real(u(1,j)-u(1,j-1))<-(ucf)), jmax=j-1; end 
        end 
    % ---------------------------------------------------------------------
    case 5  %   Placa plana com condicao de entrada imposta
        u = umax.*ones(nx, ny);
        u_bla = zeros(1, ny); 
        v_bla = zeros(1, ny);
        for i=1:nx 
		    for j = jmax:jjmax
		        u(i, j) = u_bla(j);
			    v(i ,j) = v_bla(j);
            end
        end
    % ---------------------------------------------------------------------
    otherwise
        disp('Unknown method.')
end

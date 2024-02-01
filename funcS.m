function [S] = funcS(xe)
    if (xe <= 0)
	    S = 0;
    elseif (xe > 0 && xe < 1)
	    S =  1/ (1 + exp(1/(xe-1) + 1/xe));
    else
	    S = 1;
    end
end
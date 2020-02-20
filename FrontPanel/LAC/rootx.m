function [rootx, status] = rootx(x)
    
    status = 'ok';
    if x>0
        rootx = sqrt(x);
    else
        rootx = [-1];
        status = 'negative values have no square root';
    end
end
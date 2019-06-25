function [tf] = isoidx_check_bad(i, badcyc)

    tf = false;
    if intersect(i, abs(badcyc))
        tf = true;
    end 

end

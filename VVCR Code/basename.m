% basename.m;
function new_name = basename (old_name, len)

% returns a smaller segment of the provided pathname.

nslash = 0;
oldlen = length(old_name);

new_name = old_name;

for i = oldlen : -1 : 1
    if strcmp (old_name(i), '\')
        nslash = nslash + 1;
    end
    if nslash == len
        new_name = old_name(i:oldlen);
        return;
    end
end

return;

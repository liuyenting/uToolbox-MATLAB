function list = listproduct
%LISTPRODUCT List installed product on current computer.

info = ver;
% create a cell array of the field 'Name' in 'info' struct.
list = {info.Name};

end


function colindex = findcol(header, colname)
%FINDCOL Find specified column in the header.
%
%   TBA

if isrow(colname) 
    colname = colname.';
end

nquery = size(colname, 1);

colindex = zeros([nquery, 1]);
for i = 1:nquery
    filter = startsWith(header, colname(i));
    colindex(i) = find(filter, 1, 'first');
end

end


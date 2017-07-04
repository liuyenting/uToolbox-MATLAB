function b = islicensed(tbname)
%ISLICENSED Check whether the toolbox of interest is licensed.
%
%   B = ISLICENSED(tbname) verifies the toolbox TBNAME. TBNAME is one of 
%   the name enlisted by function LISTPRODUCT.
%
%   See also: LISTPRODUCT

tblist = util.license.listproduct;
for n = tblist
    if strcmp(n, tbname)
        flex = name2flex(tbname);
        b = license('test', flex);
        return;
    end
end

warning('util:license:islicensed', 'Toolbox is not installed.');
b = false;

end

function flex = name2flex(name)
%NAME2FLEX Convert general product name to internal flex name.

pid = com.mathworks.product.util.ProductIdentifier.get(name);
flex = char(pid.getFlexName());

end

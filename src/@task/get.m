function val = get(t, propName)
% GET Get properties from the specified object
% and return the value
switch propName
    case 'Rii'
        val = t.Rii;
    case 'Ri_i'
        val = t.Ri_i;
    case 'R_i_i'
        val = t.R_i_i;
    case 'ri'
        val = t.ri;
    case 'r_i'
        val = t.r_i;
    case 'z'
        val = t.z;
    case 'usize'
        val = t.usize;
    case 'xsize'
        val = t.xsize;
    case 'nt'
        val = t.nt;
    case 'type'
        val = t.type;
    otherwise
        error([propName,' is not a property of this class'])
end

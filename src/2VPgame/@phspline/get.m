function val = get(e,propName)
% GET Get properties from the specified object
% and return the value

switch propName
    case 't'
        val = e.t;
    case 'M'
        val = e.M;
    case 'S'
        val = e.S;
    otherwise
        error([propName,' is not a property of this class'])
end




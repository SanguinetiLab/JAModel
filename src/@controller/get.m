function val = get(t, propName)
% GET Get properties from the specified object
% and return the value
switch propName
    case 'task'
        val = t.task;
    case 'lambda'
        val = t.lambda0;
    case 'a'
        val = t.a;
    case 'prior'
        val = t.prior;
    case 'umin'
        val = t.umin;
    case 'umax'
        val = t.umax;
    otherwise
        error([propName,' is not a property of this class'])
end

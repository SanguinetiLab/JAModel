function t=set(t, propName,val)
% GET Get properties from the specified object
% and return the value
switch propName
    case 'task'
        t.task=val;
    case 'lambda0'
        t.lambda0=val;
    case 'a'
        t.a = val;
    case 'prior'
        t.prior=val;
    case 'umin'
        t.umin = val;
    case 'umax'
        t.umax = val;
    otherwise
        error([propName,' is not a property of this class'])
end

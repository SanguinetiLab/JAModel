function val = get(g, propName)
% GET Get properties from the specified object
% and return the value
switch propName
    case 'task1'
        val = g.task1;
    case 'task2'
        val = g.task2;
    otherwise
        error([propName,' is not a property of this class'])
end

function val = get(e,propName)
% GET Get properties from the specified object
% and return the value

switch propName
    case 'A'
        val = e.A;     
    case 'B'
        val = e.B;
    case 'C'
        val = e.C;
    case 'D'
        val = e.D;
    case 'SigmaX'
        val = e.SigmaX;
    case 'SigmaY';
        val = e.SigmaY;

    case 'x0'
        val = e.x0;
    case 'P0'
        val = e.P0;
    case 'xsize'
        val = e.xsize;
    case 'ysize'
        val = e.ysize;
    otherwise
        error([propName,' is not a property of this class'])
end




function ps = parsetfield(ps,gh,fn)

f = getfield(ps,fn);
c = class(f);
switch c
    case 'double'
        s = get(gh,'string');
        v = str2num(s);
        if ~isempty(v)
            ps=setfield(ps,fn,v);
        else
            errordlg(['You must enter a numeric value'],'Bad Input','modal');
            txl = num2str(f);
            set(gh,'string',txl);
        end
    case 'char'
        s = get(gh,'string');
        ps=setfield(ps,fn,s);
    case 'logical'
        v  = logical(get(gh,'value'));
        ps = setfield(ps,fn,v);
    otherwise
        errordlg(['Invalid field class in ps: only char & double allowed'],'Bad Input','modal');
        txl = get(gh,'string');
        set(gh,'string',get(gh,'string'));
end


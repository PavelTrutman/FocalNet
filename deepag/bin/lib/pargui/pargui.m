% psg = pargui(ps,parname,gsz,xname,fname) - GUI for experiments
% 
% ps        ... structure with ps.ParName1, ps.ParName2, ... 
% parname   ... the name of the structure as a string
% gsz       ... the size of the grid for buttons [rnum cnum] 
% xname     ... callback function name
% fname     ... callback button label
% psg       ... pargui info structure
%
% The function creates a window with buttons for each field in ps named by
% the respective field name. The type of the button is determined by the
% type of the variable:
%                       logical ... radio button
%                       cell    ... {1} is printed   
%                       other   ... edit field
% Editing the fields changes status of ps. Pressing the single execution
% button calls function xname.
%
%
% See also PARGUIDEMO, CLEARALL, 

% T. Pajdla, pajdla@neovision.cz
% 24 Oct 2005
function ps = pargui(ps,gsz,xname,fname)

    if nargin<4
        xname = [];
    end
    if nargin<5
        fname = xname;
    end   
    
    psg.os = computer;
    psg.os = psg.os(1:3);
    
    switch psg.os
        case 'MAC'
            psg.par.th = 20; % set the grid of fields
            psg.par.adj = [0 0 0 0];
        otherwise
            psg.par.th = 20; % set the grid of fields
            psg.par.adj = [0 0 0 2*23];
    end     
    
    psg.par.bcolor  = [0.9 0.9 0.9];        % parameters    
    psg.par.fcolor  = [1 1 1];             
    psg.par.halign  = 'left';
    
    if ~isempty(ps)
        psg.par.fn = fieldnames(ps);
    else
        psg.par.fn = [];
    end
    psg.par.gn = gsz;                                  % the grid  
    
    for i=1:length(psg.par.fn)
        psg.par.p(i).t.text = psg.par.fn{i};           % fields
        psg.par.p(i).p.val  = getfield(ps,psg.par.fn{i});
    end
    
    psg.par.h = gcf;                        % get the current window
    set(psg.par.h,'menubar','none');        % no toolbar    
    set(psg.par.h,'numbertitle','off');
    set(psg.par.h,'color',psg.par.bcolor);  % set background color
    drawnow;
    psg.par.sz = get(psg.par.h,'position');  % get its size   
    %psg.par.sz = psg.par.sz + psg.par.adj;
    %set(psg.par.h,'position',psg.par.sz);
    if 0
        stck = dbstack;                         % get the caller's name         
        if length(stck)>1
            xname = stck(2).name;                            
        else
            xname = [];
        end
    end
    set(gcf,'name',fname);

    psg.par.sz = psg.par.sz([3 4]);
    if psg.par.gn(1)<2
        psg.par.dx = (psg.par.sz-[0 psg.par.th])./(psg.par.gn([2 1]));    % grid position increments
        [psg.par.gx psg.par.gy] =  meshgrid(0:psg.par.dx(1):psg.par.sz(1)-psg.par.dx(1),psg.par.sz(2)); 
    else
        % grid position increments
        psg.par.dx = psg.par.sz./psg.par.gn([2 1]);
        [psg.par.gx psg.par.gy] =  meshgrid(0:psg.par.dx(1):psg.par.sz(1)-psg.par.dx(1),psg.par.sz(2):-psg.par.dx(2):0); 
    end
    psg.par.gx = psg.par.gx+2;
    psg.par.gy = psg.par.gy;
    psg.par.xx = psg.par.gx(1);
    psg.par.xy = psg.par.gy(1);
    psg.par.gx = psg.par.gx(2:end);
    psg.par.gy = psg.par.gy(2:end);
    % execute button
    if ~isempty(xname)
        psg.par.xb.pos = [psg.par.xx psg.par.xy-psg.par.th+4 psg.par.dx(1) psg.par.th-2];
        psg.par.xb.h   = uicontrol('Style','pushbutton','string',xname,'position',psg.par.xb.pos,'BackgroundColor','red','HorizontalAlignment','center');
        set(psg.par.xb.h,'fontname','Arial');    
        %set(psg.par.xb.h,'Callback',sprintf('GuiExec=1;%s;',xname));
        set(psg.par.xb.h,'Callback',[xname '(ps, pdAll);']);
    end
    % parameter fields
    for i=1:min([length(psg.par.fn) length(psg.par.gx)]);
        psg.par.p(i).t.pos = [psg.par.gx(i) psg.par.gy(i)-psg.par.th psg.par.dx(1)/2 psg.par.th];  % a field        
        psg.par.p(i).p.pos = [psg.par.p(i).t.pos(1)+psg.par.p(i).t.pos(3),...
                              psg.par.p(i).t.pos(2),...
                              psg.par.p(i).t.pos(3:4)];  
        cls = class(psg.par.p(i).p.val);
        if ~strcmpi(cls,'cell')
            psg.par.p(i).t.h   = uicontrol('Style','text','String',psg.par.p(i).t.text,'position',psg.par.p(i).t.pos,'BackgroundColor',psg.par.bcolor,'HorizontalAlignment',psg.par.halign);
        end
        switch cls
            case {'double', 'char'}
                set(psg.par.p(i).t.h,'fontname','Arial');
                switch class(psg.par.p(i).p.val)
                    case 'double'
                        txl = num2str(psg.par.p(i).p.val);
                    case 'char'
                        txl = psg.par.p(i).p.val;
                end
                psg.par.p(i).p.h   = uicontrol('Style','edit','string',txl,'position',psg.par.p(i).p.pos,'BackgroundColor',psg.par.bcolor,'HorizontalAlignment',psg.par.halign);        
                set(psg.par.p(i).p.h,'Callback',@(gh,callbackdata)parSetField(psg.par.p(i).p.h, psg.par.p(i).t.text));
            case 'logical'
                psg.par.p(i).p.h   = uicontrol('Style','radiobutton','position',psg.par.p(i).p.pos,'BackgroundColor',psg.par.bcolor,'HorizontalAlignment',psg.par.halign);        
                set(psg.par.p(i).p.h,'Callback',@(gh,callbackdata)parSetField(psg.par.p(i).p.h, psg.par.p(i).t.text));
            case 'cell'
                psg.par.p(i).p.h   = uicontrol('Style','text','position',psg.par.p(i).p.pos,'BackgroundColor',psg.par.bcolor,'HorizontalAlignment',psg.par.halign);                        
                set(psg.par.p(i).p.h,'string',psg.par.p(i).p.val);
            otherwise
                psg.par.p(i).p.h   = uicontrol('Style','text','position',psg.par.p(i).p.pos,'BackgroundColor',psg.par.bcolor,'HorizontalAlignment',psg.par.halign);                        
                set(psg.par.p(i).p.h,'string',cls);
        end
    end
    
    function parSetField(gh,fn)
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
        assignin('base', 'ps', ps);
    end
end

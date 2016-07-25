% thrgui(t) - Thresholding GUI
%  
% Assumes a figure with an image. Modifies colormap set to 256 levels
% between the minimal and maximal intensities. 
%
% Uparrow  ... threshold = threshold + 1
% Downrrow ... threshold = threshold + 1
% Ctrl+    ... 10 * increment

% T. Pajdla, pajdla@neovision.cz
% 18 May 2005
function thrgui(t)
    if nargin<1
        t = 0.5;
    end
    ch=get(gca,'children');
    ch=get(ch(end));
%     if max(size(ch))>1
%         error('thrgui: to many children on axes');
%         error('        probably not a figure with single image, quitting ...');
%     end
%     cm = get(gcf,'colormap');
    im = ch(1).CData; % get image    
    if isinteger(im) % make double and normalize to [0,1]
        im = double(im);
    end
    mn = min(im(:)); % minimal imag evalue
    mx = max(im(:)); % maximal image value    
    ct = [0:1/255:1]'*[1 1 1]; % 255 values in between
    if isempty(get(gcf,'userdata')); % is userdata empty, set it there
        ud.thi = round(t*255);
        ud.cm  = ct;
        ud.ct  = ct*(mx-mn)+mn;
        ud.ctm = ct;
        set(gcf,'userdata',ud);
    else
        error('thrgui: userdata not empty, quitting ...');
    end
    set(gcf,'colormap',ud.ctm);
    set(gcf,'KeyPressFcn',@thrclbk);
    title(sprintf('\\^v\\pm1,Ctrl=\\pm10: threshold = %.3f',ud.ct(ud.thi)));
return

% callback function
function thrclbk(src,ed)
    ud=get(gcf,'userdata');
    s = 0;
    switch ed.Key
        case 'uparrow'
            s = +1;
        case 'downarrow'
            s = -1;
    end
    if s
        if ~isempty(ed.Modifier)
            if strcmpi('control',ed.Modifier{1})
                s = s * 10;
            end
        end
        ud.thi = ud.thi+s;
        if ud.thi < 1
            ud.thi = 1;
        end
        if ud.thi > size(ud.ct,1)
            ud.thi = size(ud.ct,1);
        end        
        set(gcf,'userdata',ud);
        ct = ud.ctm;
        ct(ud.thi:end,:)=ones(size(ct,1)-ud.thi+1,1)*[1 0 0];
        set(gcf,'colormap',ct);
        title(sprintf('\\^v\\pm1,Ctrl=\\pm10: threshold = %.3f',ud.ct(ud.thi)));
    end
return




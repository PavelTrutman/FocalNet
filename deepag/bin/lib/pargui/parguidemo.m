% PARGUIDEMO - GUI for experiments demonstration
%
% See also PARGUI

% T. Pajdla, pajdla@neovision.cz
% 25 Oct 2005

% GUI structures
% Clear, close all if for the first time from the command line.
if ~exist('GuiExec')
    closeall            % hardcore close all
    clear               % clear variables
    %% Initialize parameters structure ps
    ps.ShowFields       = logical(0);  
    ps.StringField      = 'string'; 
    ps.DoubleField      = 1;
    ps.Exit             = logical(0);
    %% Costruct GUI for the ps parameters
    gf(1) = subfig(3,6,1);             % GUI window
    set(gcf,'closerequestfcn',[]);     % make GUI unclosable
    pargui(ps,'ps',[size(fieldnames(ps),1)+1 1],'parguidemo','pargui xmple');   % create the GUI 
    parcutscreen(gf(1));               % shrink the screen 
    subfig(3,4,1);
    title('subfig(3,4,1)');
    GuiExec = 1;                       % GUI is running     
else %close all except the unclosable GUI and keep existing parameters
    close all
end
% Close and clear all and exit
if ps.Exit
    closeall
    clear all
    return
end
if ps.ShowFields
    disp(sprintf('ps.SelectionButton1 = %d',ps.ShowFields));
    disp(sprintf('ps.StringField = %s',ps.StringField));    
    disp(sprintf('ps.StringField = %f',ps.DoubleField));    
end




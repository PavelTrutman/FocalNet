function imagehandle = imaqplot(obj, event, himage)
%IMAQPLOT Create and update an image acquisition plot display.
% 
%    IMAQPLOT(OBJ, EVENT) creates and updates a live image 
%    acquisition plot display using video input object OBJ. EVENT is used 
%    to allow IMAQPLOT to be used as a function handle callback for OBJ's 
%    callback properties.
%
%    Example:
%       % Construct a video input object and a live image plot display.
%       obj = videoinput('winvideo', 1);
%       triggerconfig(obj, 'manual');
%       set(obj, 'TimerFcn', @imaqplot, 'TimerPeriod', 0.1)
%       start(obj)
%
%    IMAQPLOT(OBJ) creates an image acquisition plot display 
%    using video input object OBJ. The plot display will not update 
%    automatically. Calling IMAQPLOT repeatedly will update the plot
%    display with the latest image frame provided by GETSNAPSHOT.
%
%    Example:
%       % Construct a video input object and plot the current image.
%       obj = videoinput('matrox', 1);
%       imaqplot(obj)
%
%    IMAQPLOT(OBJ, HIMAGE)
%    IMAQPLOT(OBJ, EVENT, HIMAGE) creates a plot display using HIMAGE, an 
%    image object handle. This is useful if you wish to integrate an image
%    plot display into an existing GUI.
%
%    Example:
%       % Use an existing image object handle for a live plot display.
%       obj = videoinput('matrox', 1);
%       himage = imagesc( getsnapshot(obj) );
%       triggerconfig(obj, 'manual');
%       set(obj, 'TimerFcn', {@imaqplot, himage}, 'TimerPeriod', 0.1)
%       start(obj)
%
%       % Use an existing image object handle for a plot display.
%       obj = videoinput('matrox', 1);
%       himage = imagesc( getsnapshot(obj) );
%       imaqplot(obj, himage)
%
%    HIMAGE = IMAQPLOT(...) returns HIMAGE, the handle to the image object 
%    used for the image plot display.
%
%    Note, if an image object handle is not provided, one is created, in 
%    which case the video input object's Name property is used to configure 
%    the image object handle's Tag property. This allows the image object 
%    handle to be uniquely identified for the video input object.
%
%    See also IMAQDEVICE/GETSNAPSHOT, IMAQDEVICE/SET, IMAQDEVICE/START, 
%    VIDEOINPUT, IMAGESC.

%    CP 4-13-04
%    Copyright 2004 The MathWorks, Inc. 

%
% imaqplot(obj, himage)
% - set himage limits correctly.
%
% imaqplot(obj, event himage)
% - set himage limits correctly.
%
% vid = videoinput('winvideo');
% vid.TimerFcn = @imaqplot;
% 

try
    % Get a handle to the image object.
    if nargin==1
        % IMAQPLOT(OBJ)
        himage = localGetImageHandle(obj);
    elseif nargin==2
        if isstruct(event)
            % IMAQPLOT(OBJ, EVENT)
            himage = localGetImageHandle(obj);
        else
            % IMAQPLOT(OBJ, HIMAGE)
            himage = event;
        end
    end
    
    % Get the latest data to plot.
    if isrunning(obj)
        set(himage, 'CData', peekdata(obj, 1) );
    else
        set(himage, 'CData', getsnapshot(obj) );
    end
    
    % Adjust the axis limits in case the image resolution changed.
    res = get(obj, 'ROIPosition');
    ax = get(himage, 'Parent');
    set(ax, 'XLim', [0.5 res(3)+0.5], 'YLim', [0.5 res(4)+0.5], ...
        'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
    
    % Return the image object handle if requested.
    if nargout==1,
        imagehandle = himage;
    end
catch
    % Error gracefully.
    error('MATLAB:imaqplot:error', ...
        sprintf('IMAQPLOT is unable to plot correctly.\n%s', lasterr))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function himage = localGetImageHandle(obj)
% Locate the image handle with this object's name. 
% If one does not exist, create a new one.

% Find the image object.
name = get(obj, 'Name');
himage = findobj(findall(0), 'Type', 'image', 'Tag', name);
if ~isempty(himage)
    % FINDOBJ sometimes finds duplicate handles.
    himage = himage(1);
else
    % Create a new image object with the right dimensions.
    nbands = get(obj, 'NumberOfBands');
    himage = imagesc(rand(100, 100, nbands));
    set(himage, 'Tag', name);
    
    % Hide the figure so no one plots into it.
    set(gcf, 'HandleVisibility', 'off');
end


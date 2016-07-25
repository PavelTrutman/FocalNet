% h = imaqcam(varargin) - Live camera GUI
%
%  h = imaqcam(cam[,mode,cmd]) ... creates the live camera gui and opens
%                                  the camera in a given 'mode' 
%                                  [implicit if missing or ''] 
%                              cmd ... {'preview','close'}  
%
%      cam ... {'Sony1394', 'LgtxUSB'}
%      h   ... live camera gui structure
%
%  imaqcam(h,'close')  ... closes and cleares the existing gui
%
% See also IMAQHELP

% (c) T.Pajdla, www.neovision.cz, Sep 14, 2005
function cam = imaqcam(varargin)

if nargin<2
    varargin{2} = '';
end
if nargin<3
    varargin{3} = '';
end
cam  = varargin{1};

%% create the camera
if ischar(cam) 
    mode = varargin{2};
    cmd  = varargin{3};
    switch cam
        case 'Sony1394'
            try 
                h.vid = videoinput('dcam',1,mode);
                set(h.vid,'SelectedSourceName','input1')
                h.src = getselectedsource(h.vid);
                % setup camera
                set(h.src,'GainMode','manual');
                set(h.src,'ShutterMode','manual');
                set(h.src,'Gain',800);    % in [384 1023]
                set(h.src,'Shutter',400); % in [3 1150]
                h.imSize = get(h.vid,'VideoResolution');
                h.imSize = h.imSize([2 1]);
                cam = h;
            catch
                hw = imaqhwinfo('dcam','DeviceInfo');
                disp(['Supported video formats:' sprintf('\n') sprintf('%s\n',hw.SupportedFormats{:})]);
                error('imaqcam: invalid video format');
            end
        case 'LgtxUSB'
            try 
                h.vid = videoinput('winvideo',1,mode);
                set(h.vid,'SelectedSourceName', 'input1')
                h.src = getselectedsource(h.vid);
                % setup camera
                %set(h.src,'ColorEnable','off');
                set(h.src,'BacklightCompensation','off');
                set(h.src,'ExposureMode','manual');
                set(h.src,'PanMode','manual');
                set(h.src,'TiltMode','manual');
                set(h.src,'ZoomMode','manual');
                set(h.src,'Zoom',50);       % [50 150]
                set(h.src,'Brightness',0);  % [-10000 10000]
                set(h.src,'Gamma',500);     % [1 500]
                set(h.src,'Exposure',1);    % [1 10]
                h.imSize = get(h.vid,'VideoResolution');
                h.imSize = h.imSize([2 1]);
                cam = h;
            catch
                hw = imaqhwinfo('winvideo','DeviceInfo');
                disp(['Supported video formats:' sprintf('\n') sprintf('%s\n',hw.SupportedFormats{:})]);
                error('imaqcam: invalid video format');
            end
        otherwise
            error('imaqcam: unknown camera');
    end
else
    cmd  = varargin{2};
end

%% Exectute commands
switch cmd
    case 'close' % close the gui
        closepreview(cam.vid);
        delete(cam.vid);
        cam.vid = [];
        cam.src = [];
    case 'preview'
        preview(cam.vid);
    case ''
    otherwise
        warning(sprintf('imacqcam: unknown command: ''%s''',cmd));
end
return

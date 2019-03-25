function [op]=openexperimentwindow(cfg)
%% cfg.mode = 'openwindow';
%   cfg.skipsynctests = 1;
%   cfg.debugrect = 0;
%   cfg.backgroundcolor = [127 127 127 255];
%   cfg.blendfunction = 'no';
%   cfg.sourcefactor
%   cfg.destinationfactor
%       % e.g. % Screen('BlendFunction', window, 'GL_ONE', 'GL_ONE');
%              % Screen('BlendFunction', window, 'GL_DST_ALPHA', 'GL_ONE_MINUS_DST_ALPHA');
%              % Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
%   cfg.viewingdistance
%     op.windowptr = window;
%     op.windowrect = windowRect;
%     op.widthinpix = wid_inpix;
%     op.heightinpix = height_inpix;
%     op.widthincm = wid_incm;
%     op.heightincm = height_incm;
%     op.xcenter = xCenter;
%     op.ycenter = yCenter;
%     op.framerate = framerate;
%     op.pixperdeg_h
%     op.pixperdeg_v
%%
switch cfg.mode
    case 'openwindow'
        PsychDefaultSetup(1);
        commandwindow;
        Screen('Preference', 'SkipSyncTests',cfg.skipsynctests);
        
        
        screens = Screen('Screens');
        screenNumber = max(screens);
        
        if cfg.debugrect == 1
            rect =  [800 400 1440 850];
        else
            rect = [];
        end
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        
        
        [window, windowRect] = PsychImaging('OpenWindow',...
            screenNumber, cfg.backgroundcolor,rect);
        [wid_inpix, height_inpix] = Screen('WindowSize', window); % sz in px
        [xCenter, yCenter] = RectCenter(windowRect);
        [wid_incm, height_incm] = Screen('DisplaySize', screenNumber); % phys sz in mm
        
        framerate = FrameRate(window);
        
        switch cfg.blendfunction
            case 'yes'
                try
                    Screen('BlendFunction', window, cfg.sourcefactor, cfg.destinationfactor);
                catch
                    error('Define source and destination factors for the blending function');
                end
            otherwise
                disp('No blending function is being used.');
        end
        
        [ppd_h,ppd_v] = deg2pix(1,cfg.viewingdistance,window);
        
        op.windowptr = window;
        op.windowrect = windowRect;
        op.widthinpix = wid_inpix;
        op.heightinpix = height_inpix;
        op.widthincm = wid_incm;
        op.heightincm = height_incm;
        op.xcenter = xCenter;
        op.ycenter = yCenter;
        op.framerate = framerate;
        op.pixperdeg_h = ppd_h;
        op.pixperdeg_v = ppd_v;
    otherwise
        disp('Not developed yet!');
end
end
%% Function that converts degrees to pixels
function  [qPix_h,qPix_v] = deg2pix(deg,dist,windowpointer,roundthenumber)
if nargin < 3
    windowpointer = 0;
    warning('IF THE WINDOW POINTER IS OTHER THAN "0" PLEASE SPECIFY!');
end
if nargin < 4
    roundthenumber = 1;
    warning('The number of pixels are rounded to an integer. If you want otherwise the 4th input should be "0"');
end

screens = Screen('Screens');
screen_no = max(screens);

[wid_inpix,height_inpix] = Screen('WindowSize', windowpointer);
[wid_incm, height_incm] = Screen('DisplaySize', screen_no);

qPix_h = dist*tan(deg*pi/180)/...
    (wid_incm/(10*wid_inpix));

qPix_v = dist*tan(deg*pi/180)/...
    (height_incm/(10*height_inpix));

if roundthenumber
    qPix_h = round(qPix_h);
    qPix_v = round(qPix_v);
end
end
function op = createdoubledriftstimulus(cfg)
%%
%% cfg.mode = 'linearrect'
%   cfg.pixelsperdegree
%   cfg.framerate
%   cfg.distancefromthefixation %(degrees)
%   cfg.xcenter
%   cfg.ycenter
%   cfg.angleoftheinitialpointvector %the angle of the vector drawn from
%                                    %the origin to the initial position of
%                                    %the stimulus.
%   cfg.quadrantlocation %(i.e. in which quadrant of the coordinate system
%                        %the stimulus will appear)
%   cfg.trialduration %(seconds)
%   cfg.reversal %(if 1, the rect ends at the initial location)
%   cfg.horizontaldirection %('left','right','none')
%   cfg.verticaldirection %('up','down','none')
%   cfg.envelopesize %(degrees)
%   cfg.externalspeed
%   cfg.externalmotionorientation;
%
%   "op" is the (number of frames x 4) destination rect for each frame;
%% cfg.mode = 'pinknoise' or '1/f'
%   cfg.pixelsperdegree
%   cfg.framerate
%   cfg.internalspeed
%   cfg.trialduration %(seconds)
%   cfg.envelopesize %(degrees)
%   cfg.backgroundcolour
%   cfg.tar_frequency
%   cfg.numberofoctaves
%   cfg.randomnoise % = 1 only for the control stimulus without illusion
%   cfg.reversal %(if 1, the stimulus ends at the initial location)
%
%   "op.mat" is the stimulus matrix (pixels x pixels x frames)
%   "op.rect' is the stimulus rect

%% cfg.mode = 'gabor'
%  cfg.function = 'create'
%   cfg.pixelsperdegree
%   cfg.spatialfrequency
% 	cfg.contrast % = .5;
% 	cfg.aspectratio % = 1;
% 	cfg.bgoffset % = [0 0 0 1];
% 	cfg.gaborwidth
% 	cfg.gaborheight
% 	cfg.disablenorm % = 1;
% 	cfg.precontrasmultiplier % = 1;
% 	cfg.windowptr
% 	cfg.phase
% 	cfg.sigma
%
%   "op.ptr" gabor texture
%   "op.rect"
%   "op.propmat" properties matrix
%
% cfg.mode = 'gabor';
% cfg.function = 'updatepropmat'
%   cfg.propmat = op.propmat
%   cfg.phasefactor
% 	cfg.framerate
% 	cfg.internalcyclespersec
%
%   "op" propmat
%
% cfg.mode = 'gabor';
% cfg.function = 'creatematrix'
%     cfg.gaborsize;
%     cfg.pixelsperdegree;
%     cfg.trialduration;
%     cfg.framerate;
%     cfg.angle;
%     cfg.sigma
%     cfg.cyclesperdeg;

%% cfg.mode = 'gratingmask';
%  cfg.function = 'create';
%     cfg.pixelsperdegree 
%     cfg.gratingwidth
%     cfg.gratingheight
%     cfg.spatialfrequency
%     cfg.cyclesize
%     cfg.bgcolor
%     cfg.contrast

%  cfg.mode = 'gratingmask';
%  cfg.function = 'updaterect';
%     cfg.pixelsperdegree 
%     cfg.framenumber
%     cfg.cyclesize
%     cfg.cyclespersecond
%     cfg.framerate
%     cfg.gratingwidth
%     cfg.gratingheight

%% cfg.mode = 'perlinnoise'
%     cfg.size;
%     cfg.pixelsperdegree

%% cfg.mode = 'gaussianmask';
%     cfg.pixelsperdegree
%     cfg.stimulusmatrix
%     cfg.sigma
%%
switch cfg.mode
    case 'linearrect'
        
        ppd = cfg.pixelsperdegree;
        frmrate = cfg.framerate;
        dist_decent = cfg.distancefromthefixation;
        th_decent = cfg.angleoftheinitialpointvector;
        whQuad = cfg.quadrantlocation;
        qFrm = round(cfg.trialduration*cfg.framerate);
        reversal = cfg.reversal;
        dir_h = cfg.horizontaldirection;
        dir_v = cfg.verticaldirection;
        env_sz = round(cfg.envelopesize*cfg.pixelsperdegree);
        env_sz = env_sz+(mod(env_sz,2)==0);
        rect_env = [0 0 env_sz env_sz];
        ext_spd = cfg.externalspeed;
        ext_or = cfg.externalmotionorientation;
        cent_x = cfg.xcenter;
        cent_y = cfg.ycenter;
        
        dplc_perfrm = displaceperfrm(ext_spd,frmrate);
        [dplc_perfrm_h,dplc_perfrm_v] = length_inxy(dplc_perfrm,ext_or,ppd,ppd);
        
        
        
        % rows: 1st is left; 2nd is right; 1st is up; 2nd is down
        if reversal
            extdrift_h(1,:) = -1*[0:dplc_perfrm_h:(qFrm/2)*dplc_perfrm_h (qFrm/2)*dplc_perfrm_h-1:-1*dplc_perfrm_h:1];
            extdrift_v(1,:) = -1*[0:dplc_perfrm_v:(qFrm/2)*dplc_perfrm_v (qFrm/2)*dplc_perfrm_v-1:-1*dplc_perfrm_v:1];
        else
            extdrift_h(1,:) = -1*[0:dplc_perfrm_h:(qFrm)*dplc_perfrm_h];
            extdrift_v(1,:) = -1*[0:dplc_perfrm_v:(qFrm)*dplc_perfrm_v];
        end
        
        extdrift_h(2,:) = -1*extdrift_h(1,:);
        extdrift_v(2,:) = -1*extdrift_v(1,:);
        
        if isempty(extdrift_v)
            extdrift_v = zeros(size(extdrift_h));
        end
        
        if isempty(extdrift_h)
            extdrift_h = zeros(size(extdrift_v));
        end

        [x,y]= length_inxy(dist_decent,th_decent,ppd,ppd);
        dist_decent_x = x + max(abs(extdrift_h(1,:)));
        dist_decent_y = y + max(abs(extdrift_v(1,:)));
        
        signx = 2*(whQuad==2 || whQuad==4) - 1;
        signy = 2*(whQuad==3 || whQuad==4) - 1;
        
        whH = 1+strcmp(dir_h,'right');
        whV = 1+strcmp(dir_v,'down');
        
        for whFrm = 1:qFrm
            pos_x = cent_x +...
                signx * dist_decent_x +...
                extdrift_h(whH,whFrm);
            pos_y = cent_y +...
                signy * dist_decent_y +...
                extdrift_v(whV,whFrm);
            rect_dst(whFrm,:) = CenterRectOnPointd(...
                rect_env, ...
                pos_x,...
                pos_y);
        end
        
        op = rect_dst;
        
    case {'pinknoise','1/f'}
        ppd = cfg.pixelsperdegree;
        frmrate = cfg.framerate;
        int_spd = cfg.internalspeed;
        qFrm = round(cfg.trialduration*cfg.framerate);
        env_sz = round(cfg.envelopesize*cfg.pixelsperdegree);
        env_sz = env_sz+(mod(env_sz,2)==0);
        gray = cfg.backgroundcolour;
        tar_fq = cfg.tar_frequency*cfg.pixelsperdegree;
        qOct = cfg.numberofoctaves;
        randomnoise = cfg.randomnoise;
        isReversal = cfg.reversal;
        
        stim_mat = pinknoisetexture(ppd,int_spd,1/frmrate,qFrm,env_sz,gray(1),tar_fq,qOct,isReversal,randomnoise);
        [sz,~,~] = size(stim_mat);
        op.mat = stim_mat;
        op.rect = [0 0 sz sz];
    case 'gabor'
        switch cfg.function
            case 'create'
                ppd = cfg.pixelsperdegree;
                spat_fq = cfg.spatialfrequency/cfg.pixelsperdegree;
                cont = cfg.contrast;
                asp_rat = cfg.aspectratio;
                bg_offset = cfg.bgoffset;
                wid = cfg.gaborwidth;
                wid = round(wid*ppd);
                wid = wid + (mod(wid,2)==0);
                height = cfg.gaborheight;
                height = round(height*ppd);
                height = height + (mod(height,2)==0);
                dblnorm = cfg.disablenorm;
                precont_multplr = cfg.precontrasmultiplier;
                w = cfg.windowptr;
                phase = cfg.phase;
                sig = cfg.sigma*cfg.pixelsperdegree;

                [gabor.ptr, gabor.rect] = CreateProceduralGabor(w, ...
                    wid,...
                    height,...
                    [],...
                    bg_offset,...
                    dblnorm,...
                    precont_multplr);

                gabor.propmat = [phase,...
                    spat_fq,...
                    sig,...
                    cont,...
                    asp_rat,...
                    0, 0, 0];

                op = gabor;
            case 'creatematrix'
                % Marvin's code
                env_sz = cfg.gaborsize;
                ppd = cfg.pixelsperdegree;
                trl_dur = cfg.trialduration;
                temporal_freq = 1/trl_dur;
                RefreshRate = cfg.framerate;
                sigma = cfg.sigma*ppd;
                cycperdeg = cfg.cyclesperdeg;
                cycperdeg2 = cycperdeg+4;

                gaborSizePix = env_sz *ppd;
                if ~mod(gaborSizePix,2)
                    gaborSizePix = gaborSizePix + 1;
                end
                
                halfgaborsize=floor(gaborSizePix/2);
                spatfreq = cycperdeg/ppd; % cycles per deg in cycles per pix
                pixelsPerCycle = ceil(1/spatfreq);
                intDriftCycPerFrame = temporal_freq/RefreshRate; % this fraction of a cycle per frame 
                intDriftPixPerFrame = ceil(pixelsPerCycle * intDriftCycPerFrame);
                framesPerCyc = ceil(RefreshRate/temporal_freq);
                % Compute actual cosine grating:
                gratingNEW = repmat(sin(linspace(0,2*pi,pixelsPerCycle)),gaborSizePix,ceil(gaborSizePix/pixelsPerCycle));
                gratingNEW = [gratingNEW,gratingNEW(:,2:end),gratingNEW(:,2:end);
                              gratingNEW,gratingNEW(:,2:end),gratingNEW(:,2:end)]; % makes grating from -1 to 1
%                 gratingNEW = gratingNEW + 1; % now grating from 0 to 2
%                 gratingNEW = gratingNEW/2; % grating from 0 to 1
%                 gratingNEW = gratingNEW * (255); % grating from 0 to 255

                spatfreq2 = cycperdeg2/ppd; % cycles per deg in cycles per pix
                pixelsPerCycle2 = ceil(1/spatfreq2);
                intDriftCycPerFrame = temporal_freq/RefreshRate; % this fraction of a cycle per frame 
                intDriftPixPerFrame2 = ceil(pixelsPerCycle2 * intDriftCycPerFrame);
                framesPerCyc = ceil(RefreshRate/temporal_freq);
                % Compute actual cosine grating:
                gratingNEW2 = repmat(sin(linspace(0,2*pi,pixelsPerCycle2)),gaborSizePix,ceil(gaborSizePix/pixelsPerCycle2));
                gratingNEW2 = [gratingNEW2,gratingNEW2(:,2:end),gratingNEW2(:,2:end);
                              gratingNEW2,gratingNEW2(:,2:end),gratingNEW2(:,2:end)]; % makes grating from -1 to 1
%                 gratingNEW2 = gratingNEW2 + 1; % now grating from 0 to 2
%                 gratingNEW2 = gratingNEW2/2; % grating from 0 to 1
%                 gratingNEW = gratingNEW * (255); % grating from 0 to 255
                
                % make gabors 
                % gabor = nan(2*sourcerect_size+1, 2*sourcerect_size+1,4);
                [x,y] = meshgrid(-halfgaborsize:halfgaborsize,-halfgaborsize:halfgaborsize);
                gabor = nan(gaborSizePix, gaborSizePix,3);
%                 gaussian =round(255 * (1 - exp(-(x.^2)/(2*sigma.^2)-(y.^2)/(2*sigma.^2))));
                gaussian =exp(-(x.^2)/(2*sigma.^2)-(y.^2)/(2*sigma.^2));               
                ngabors = ceil(framesPerCyc);
                % gabortex = nan(1,ngabors);
                for texi = 1:ngabors
                    startpoint = mod(texi*intDriftPixPerFrame,pixelsPerCycle)+ 1;
%                     startpoint = (texi-1)*intDriftPixPerFrame + 1;
                    grating = gratingNEW(1: size(gaussian,1),...
                        startpoint:startpoint -1 + size(gaussian,1));
                    
                    startpoint2 = mod(texi*intDriftPixPerFrame2,pixelsPerCycle2)+ 1;
%                     startpoint = (texi-1)*intDriftPixPerFrame + 1;
                    grating2 = gratingNEW2(1: size(gaussian,1),...
                        startpoint2:startpoint2 -1 + size(gaussian,1));
                    
                    for layeri = 1:3
                        gabor(:,:,layeri) = rescale(grating + transpose(grating2));
                    end
                    
                    for layeri = 1:3
                        gabor(:,:,layeri) = gabor(:,:,layeri)*255;%.*  (gaussian * 255);%((-1*gaussian) + 255);
                    end
%                     opaqueThresh = 0;
%                     gaussian2 = ((-1*gaussian) + 255);
%                     gaussian2 = 0.*(gaussian2 < opaqueThresh) + gaussian2.*(gaussian2 >=opaqueThresh);
                    gabor(:,:,4) = gaussian*255;%gaussian2;

                %     gabors(:,:,:,texi) = gabor;
                %     gaborins(:,:,texi) = grating .*  ((-1*gaussian) + 255);
                    gaborstim(texi,:,:,:) = gabor;%Screen('MakeTexture', w, gabor);
                end
                
                op = gaborstim;
                
            case 'updatepropmat'
                propmat = cfg.propmat;
                phasefactor = cfg.phasefactor;
                cyc = cfg.internalcyclespersec;
                incperfrm = cyc*360/cfg.framerate;
                
                propmat(1) = propmat(1)...
                    + phasefactor ...
                    * incperfrm;
                
                op = propmat;
        end
    case 'gratingmask'
        switch cfg.function
            case 'create'
                ppd = cfg.pixelsperdegree;
                wid = round(cfg.gratingwidth*cfg.pixelsperdegree);
                wid = wid + (mod(wid,2)==0);
                height = round(cfg.gratingheight*cfg.pixelsperdegree);
                height = height + (mod(height,2)==0);
                spat_fq = cfg.spatialfrequency;
                pixpercyc = round(cfg.cyclesize*cfg.pixelsperdegree);
                fq = (1/pixpercyc)*2*pi;
                bgcol = cfg.bgcolor;
                cont = cfg.contrast;

                texsize = wid/2;
                x = meshgrid(-texsize:texsize + pixpercyc, 1);

                grating = bgcol * cos(fq*x) + bgcol;
                mask = ones(1, numel(x), 2) * bgcol;

                mask(:, :, 2)= grating .* cont;

                op = mask;
            case 'updaterect'
                ppd = cfg.pixelsperdegree;
                frameCounter = cfg.framenumber-1;
                pixpercyc = round(cfg.cyclesize*cfg.pixelsperdegree);
                shiftPerFrame = (cfg.cyclespersecond*pixpercyc)/cfg.framerate;
                wid = round(cfg.gratingwidth*cfg.pixelsperdegree);
                wid = wid + (mod(wid,2)==0);
                height = round(cfg.gratingheight*cfg.pixelsperdegree);
                height = height + (mod(height,2)==0);
                
                
                xoffset = mod(frameCounter * shiftPerFrame, pixpercyc);
                srcRect = [xoffset 0 xoffset+wid height];
                
                op = srcRect;
        end
    case 'perlinnoise'
        %from Liwei's ddspiral_nulling2
        
        m = round(cfg.size*cfg.pixelsperdegree);
        m = m + (mod(m,2)==0);
        
        
        s = zeros([m,m]);     % Prepare output image (size: m x m)
        w = m;
        i = 0;
        while w > 8
            i = i + 1;
            d = interp2(randn([m,m]), i-1, 'spline');
            s = s + i * d(1:m, 1:m);
            w = w - ceil(w/2 - 1);
        end
        s = (s - min(min(s(:,:)))) ./ ((max(max(s(:,:))) - min(min(s(:,:))))*4)+.375;%range .375 to .625
        
        % want a circle
        [x,y]=meshgrid(1:m,1:m);
        s((x-m/2).^2 + (y-m/2).^2 > (m/2)^2) = .5;
        
        op = s*255;
    case 'gaussianmask'
        ppd = cfg.pixelsperdegree;
        stim = cfg.stimulusmatrix;
        sig = cfg.sigma*ppd;
        [env_sz,~] = size(stim);
        
        imWidth = floor(env_sz/2);
        [gx,gy]=meshgrid(-imWidth:imWidth, -imWidth:imWidth);
        env = exp( -((gx.^2)+(gy.^2)) /(2*(sig)^2));
        
        op = stim.*env+127;
end
end
%% Draw Fixation
function drawfixation(w,x,y,sz,fix_col,flip,wait)
if nargin<6
    flip = 'dont';
    wait = 0;
end
if nargin<7
    wait = 0;
end
Screen('DrawDots', w,[x y], sz, fix_col, [], 2);
if strcmp(flip,'flip')
    Screen('Flip',w);
end
if wait
    WaitSecs(wait);
end
end
%%
function mov = recdisplay(rec,mov,command,var1,var2)
% mov.name
% mov.framerate
% mov.dir
% 1. command='create'
%       rec==true = 1
%       
% 2. command='record'
%       rec==true = 1
%       var1: window pointer
%       var2: duration of the frame to stay in the display in sec
%
% 3. command='finalize'
%       rec==true = 1

if rec
    if ~isfield(mov,'name')
        mov.name = ['genericmoviename',num2str(randi(10000)),'.avi'];
    else
        mov.name = [mov.name,'.avi'];
    end
    if ~isfield(mov,'framerate'),   mov.framerate = 60;    end
    
    switch command
        case 'create'
            if isfield(mov,'dir')
             cd(mov.dir)
            end
            mov.obj = VideoWriter(mov.name);
            mov.obj.FrameRate = mov.framerate;
        case 'record'            
            currFrame = Screen('GetImage', var1);
            if nargin<5
                var2 = 0;
                repeatframe = 1;
            else
                repeatframe = round(var2*mov.framerate);
            end
            for whFrm = 1:repeatframe
                open(mov.obj);
                writeVideo(mov.obj,currFrame);
            end
        case 'finalize'
            close(mov.obj)        
    end
end
end
%%
%%
function [stim] = pinknoisetexture(ppd,int_spd,flp_intv,qFrm,env_sz,bg_col,tar_fq,qOct,isReversal,isControl)

if nargin<9
    isControl = 0;
end

if isControl
    stepf = round(ppd*(int_spd*flp_intv)/2);
    img = zeros(env_sz, env_sz, qFrm+2);
else
    stepf = round(ppd*(int_spd*flp_intv));
    img = zeros(env_sz*2+stepf*qFrm, env_sz);
end

pink = (255 * makepinknoise(img, tar_fq, qOct,isControl,stepf))  - bg_col;
[stim] = framesIllusion(stepf, env_sz, bg_col, tar_fq, qFrm, pink,isReversal,isControl);
end
%%
function img = makepinknoise(img, w, qOct,isControl,step,persistence, lacunarity)
% creates and sum successive noise functions, each with higher
% frequency and lower amplitude. nomalize output within [0 1]
%
% input:
% - im: initial matrix
% - wl: grid size in pixels of the first noise octave (lowest spatial frequency)
% - octaves: number of noisy octaves with increasing spatial frequency to
%            be added
%
% optional:
% - lacunarity: frequency multiplier for each octave (usually set to 2 so
%               spatial frequency doubles each octave)
% - persistence: amplitude gain (usually set to 1/lacunarity)
%
% When lacunarity=2 and persistence=0.5 you get ~ 1/f noise
%
if nargin <= 5
    lacunarity = 2;
    persistence = 0.5;
end
if isControl
    [n, m, v] = size(img);
    a = 1;
    for oct = 1:qOct
        rndim = -1 +2*rand(ceil(n/w),ceil(m/w),step*ceil(v/w));   % uniform
        [Xq,Yq,Zq] = ndgrid(linspace(1,size(rndim,2),m),linspace(1,size(rndim,1),n),linspace(1,size(rndim,3),v));
        d = interp3(rndim,Xq,Yq,Zq, 'cubic');
        img = img + a*d(1:n, 1:m, 1:v);
        a = a*persistence;
        w = w/lacunarity;
    end
    img = (img - min(min(min(img(:,:,:))))) ./ (max(max(max(img(:,:,:)))) - min(min(min(img(:,:,:)))));
else
    
    [n, m] = size(img);
    a = 1;
    for oct = 1:qOct
        rndim = -1 + 2*rand(ceil(n/w),ceil(m/w)); % uniform
        [Xq,Yq] = meshgrid(linspace(1,size(rndim,2),m),linspace(1,size(rndim,1),n));
        d = interp2(rndim,Xq,Yq, 'cubic');
        img = img + a*d(1:n, 1:m);
        a = a*persistence;
        w = w/lacunarity;
    end
    img = (img - min(min(img(:,:)))) ./ (max(max(img(:,:))) - min(min(img(:,:))));
end
end
%%
function[m] = framesIllusion(stepf,env_sz,bg_col, tar_fq, nFrames, noiseimg,isReversal, control)

% here we prepare one set texture with drifting internal motion (2D noise)
% - this version works with the perceptual task left/right tilt -

if isReversal
    reversal = nFrames/2;
else
    reversal = nFrames;
end

% gaussian envelope
imWidth = floor(env_sz/2);
[gx,gy]=meshgrid(-imWidth:imWidth, -imWidth:imWidth);
env = exp( -((gx.^2)+(gy.^2)) /(2*(tar_fq)^2));
m = zeros(env_sz, env_sz, nFrames);

if control
    c = 3;
    for fi=1:nFrames
        if fi>1
            if fi<=reversal
                c = c+1; %step;
            else
                c = c-1; %step;
            end
        end
        noisePatt = noiseimg(:,:,c);
        m(:,:,fi) = uint8(bg_col + noisePatt.*env);
    end
else
    segBeg = 1;
    segEnd = env_sz;
    segBeg2 = 1 + stepf*reversal;
    segEnd2 = env_sz + stepf*reversal;
    
    % compute textures for individual frames
    cf = 0; cb = 0; fi = 0;
    for i=1:nFrames
        if i<=reversal
            aBeg = segBeg + (cf*stepf);
            aEnd = segEnd + (cf*stepf);
            cf = cf+1;
        else
            aBeg = segBeg2 - (cb*stepf);
            aEnd = segEnd2 - (cb*stepf);
            cb = cb+1;
        end
        fi = fi + 1;
        noisePatt = noiseimg(aBeg:aEnd,:);
        m(:,:,fi) = uint8(bg_col + noisePatt.*env);
    end
end
end
function draw_DDpatches(scr,coordX,whFrm,pink,int_thX,ext_thX,varargin)

isDrawDot = false;
info = size(coordX);

if ~isempty(varargin)
    for idx = 1:length(varargin)
        argN = varargin{idx};
        if isscalar(argN) || ischar(argN)
            switch argN
                case 'SuperimposeDots'
                    isDrawDot = true;
                    dot_gunX(1,:) = [1,0,0,1]*255;
                    dot_gunX(2,:) = [0,1,0,1]*255;
                    dot_sz = 3;
                    dot_clr = repmat([1;2],info(1)/2,1);
                case 'DotColor'
                    dot_gunX(1,:) = varargin{idx+1};
                    dot_gunX(2,:) = varargin{idx+2};
                case 'DotSize'
                    dot_sz = scr.pixperdeg_h*varargin{idx+1};
                case 'AssignColor'
                    dot_clr = varargin{idx+1};
                    
            end
        end
    end
end

if length(whFrm)==1
    whFrm = ones(1,info(1))*whFrm;
end

if length(info) == 4
    
    if length(int_thX) == 1
        int_thX = ones(1,info(1))*int_thX;
    end
    
    if length(ext_thX) == 1
        ext_thX = ones(1,info(1))*ext_thX;
    end
    
    for whPatch = 1:info(1)
        if whFrm(whPatch)
            int_th = int_thX(whPatch);
            ext_th = ext_thX(whPatch);

            coord = squeeze(coordX(whPatch,ext_th,whFrm(whPatch),1:2));
            rect = CenterRectOnPointd(pink.rect,coord(1),coord(2));

            Screen('DrawTexture', scr.windowptr, pink.tex(whPatch, whFrm(whPatch)), [], rect,int_th);
            if isDrawDot
                draw_dots(scr.windowptr,coord(1),coord(2),'Size',dot_sz,'Color',dot_gunX(dot_clr(whPatch),:));
            end
        end
    end
else
    int_th = int_thX;
    for whPatch = 1:info(1)
        coord = squeeze(coordX(whPatch,whFrm(whPatch),1:2));
        rect = CenterRectOnPointd(pink.rect,coord(1),coord(2));
        
        Screen('DrawTexture', scr.windowptr, pink.tex(whFrm(whPatch)), [], rect,int_th);
        if isDrawDot
            draw_dots(scr.windowptr,coord(1),coord(2),'Size',dot_sz,'Color',dot_gunX(dot_clr(whPatch),:));
        end
    end
end
end
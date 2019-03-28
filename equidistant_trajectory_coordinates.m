function coord = equidistant_trajectory_coordinates(cfg)
% cfg.numberoftrajectories
% cfg.trialduration
% cfg.framerate
% cfg.degrees_persecond
% cfg.horizontalpixels_perdegree
% cfg.verticalpixels_perdegree
% cfg.equidistant_point
% cfg.distancefromcenter
% cfg.xcenter
% cfg.ycenter
% cfg.angleoftrajectory

trj_th = deg2rad(cfg.angleoftrajectory);

trjN_cfg = cfg;
qTrj = cfg.numberoftrajectories;
deg1 = 180-(360/qTrj)/2;
thX = divide360(qTrj,deg1);

for whPos = 1:qTrj
    for whTrj = 1:length(trj_th)
        trjN_cfg.angleoftrajectory = trj_th(whTrj);
        trjN_cfg.angleofinitialpointvector = thX(whPos);
        [hor_coord, ver_coord] = create_linear_trajectory(trjN_cfg);
        coord(whPos,whTrj,:,:) = [hor_coord', ver_coord'];   
    end
end
coord = squeeze(coord);
end
%%
function th = divide360(q,deg1)
cw = nsidedpoly(q);
cw = cw.Vertices;
[th,~] = cart2pol(cw(:,1),cw(:,2));
th = mod(rad2deg(th),360);
th = mod(th+deg1-th(1),360);
th = deg2rad(th);
end
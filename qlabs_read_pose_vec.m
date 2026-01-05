function y = qlabs_read_pose_vec(reset)
% y = [x; y; th; valid]  (4x1)

hQBot = evalin('base','hQBot');

persistent x0 y0 init
if isempty(init), init = false; end

[ok, loc, rot, ~] = hQBot.get_world_transform();

if ~ok
    y = zeros(4,1);
    return;
end

xw = double(loc(1));
yw = double(loc(2));

% yaw guess: rot(3)
thw = atan2(sin(double(rot(3))), cos(double(rot(3))));

if ~init || reset ~= 0
    x0 = xw; y0 = yw; init = true;
end

x = xw - x0;
yy = yw - y0;

y = [x; yy; thw; 1];
end

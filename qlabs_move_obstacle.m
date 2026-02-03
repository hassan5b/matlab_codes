function y = qlabs_move_obstacle(u)
% u = [x_world; y_world; z_world; yaw_deg]  (yaw optional)
% Returns dummy scalar y for Simulink.

global hObs OB_SCALE OB_Z

y = 0;

if isempty(hObs) || ~isvalid(hObs)
    return;
end

xw = u(1);
yw = u(2);

if numel(u) >= 3
    zw = u(3);
else
    zw = OB_Z;
end

if numel(u) >= 4
    yaw = u(4);
else
    yaw = 0;
end

loc = [xw, yw, zw];
rot = [0, 0, yaw];          % degrees
scl = OB_SCALE;

% Non-blocking update
hObs.set_transform_degrees(loc, rot, scl, false);
end

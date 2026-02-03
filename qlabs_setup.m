caseNum = 1; % 1, 2, 3, or 4

system('quanser_host_peripheral_client.exe -q');
pause(2)
system('quanser_host_peripheral_client.exe  -uri tcpip://localhost:18444 &');

% MATLAB Path

newPathEntry = fullfile(getenv('QAL_DIR'), '0_libraries', 'matlab', 'qvl');
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(newPathEntry, pathCell));
else
  onPath = any(strcmp(newPathEntry, pathCell));
end

if onPath == 0
    path(path, newPathEntry)
    savepath
end

% Stop RT models
try
    qc_stop_model('tcpip://localhost:17000', 'qbot_platform_driver_virtual')
    pause(1)
    qc_stop_model('tcpip://localhost:17000', 'QBotPlatform_Workspace')
catch error
end
pause(1)
goal = [7; -1];      % example
obs_c = [4; 0.5]; % example
obs_D = 0.5;       % example
% QLab connection
qlabs = QuanserInteractiveLabs();
connection_established = qlabs.open('localhost');

if connection_established == false
    disp("Failed to open connection.")
    return
end
disp('Connected')
verbose = true;
num_destroyed = qlabs.destroy_all_spawned_actors();

map_offset_x = 4;

% ---- Flooring ----
hFloor0 = QLabsQBotPlatformFlooring(qlabs);

% For an 8x8 arena, we scale the previous 4x4 layout by 2
% (positions doubled and scale doubled).
tileScale = [2, 2, 1];

% center (shifted)
hFloor0.spawn_id(0, [-1.2+map_offset_x,  1.2, 0], [0,0,-pi/2], tileScale, 5, false);

% corners (shifted)
hFloor0.spawn_id(1, [ 1.2+map_offset_x,  3.6, 0], [0,0,-pi/2], tileScale, 0, false);
hFloor0.spawn_id(2, [ 3.6+map_offset_x, -1.2, 0], [0,0, pi  ], tileScale, 0, false);
hFloor0.spawn_id(3, [-1.2+map_offset_x, -3.6, 0], [0,0, pi/2], tileScale, 0, false);
hFloor0.spawn_id(4, [-3.6+map_offset_x,  1.2, 0], [0,0,    0], tileScale, 0, false);

% sides (shifted)
hFloor0.spawn_id(5, [-1.2+map_offset_x,  1.2, 0], [0,0,    0], tileScale, 5, false);
hFloor0.spawn_id(6, [ 1.2+map_offset_x,  1.2, 0], [0,0,-pi/2], tileScale, 5, false);
hFloor0.spawn_id(7, [ 1.2+map_offset_x, -1.2, 0], [0,0, pi  ], tileScale, 5, false);
hFloor0.spawn_id(8, [-1.2+map_offset_x, -1.2, 0], [0,0, pi/2], tileScale, 5, false);

% ---- Walls ----
hWall = QLabsWalls(qlabs, verbose);

% 8x8 boundaries in the unshifted frame are x = ±4, y = ±4
xL = -4 + map_offset_x;
xR =  4 + map_offset_x;
yT =  4;
yB = -4;

% Use the same segment spacing idea as before, extended to cover 8 m
ys = [ 3.6, 2.4, 1.2, 0, -1.2, -2.4, -3.6 ];
xs = [ 3.6, 2.4, 1.2, 0, -1.2, -2.4, -3.6 ] + map_offset_x;

% Left Walls (x = -4)
hWall.spawn_degrees([xL,  3.6, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xL,  2.4, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xL,  1.2, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xL,  0.0, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xL, -1.2, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xL, -2.4, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xL, -3.6, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);

% Right Walls (x = +4)
hWall.spawn_degrees([xR,  3.6, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xR,  2.4, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xR,  1.2, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xR,  0.0, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xR, -1.2, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xR, -2.4, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xR, -3.6, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);

% Top Walls (y = +4), rotate 90 deg
hWall.spawn_degrees([xs(1), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(2), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(3), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(4), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(5), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(6), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(7), yT, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);

% Bottom Walls (y = -4), rotate 90 deg
hWall.spawn_degrees([xs(1), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(2), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(3), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(4), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(5), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(6), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);
hWall.spawn_degrees([xs(7), yB, 0.1], [0, 0, 90]); hWall.set_enable_dynamics(true);


%xhWall.spawn_degrees([4, 0, 0.1], [0, 0, 0]); hWall.set_enable_dynamics(true);

% ----- obstacle----
global hObs OB_ACTOR OB_SCALE OB_Z;

hObs = QLabsBasicShape(qlabs);

OB_ACTOR = 200;                 % choose an unused actor number
OB_Z     = 0.25;                % height above ground
OB_SCALE = [0.75 0.75 0.25];     % visual diameter in x/y, height in z

% initial obstacle position in WORLD frame
obs0_world = [4, 1, OB_Z];   % set these based on your map offset handling
rot0_deg   = [0 0 0];

% Pick the shape you want to visualize (sphere is usually best for a circle)
cfg = QLabsBasicShape.SHAPE_SPHERE;        % or SHAPE_CYLINDER :contentReference[oaicite:1]{index=1}

hObs.spawn_id_degrees(OB_ACTOR, obs0_world, rot0_deg, OB_SCALE, cfg, true);


% ---- QBot ----
hQBot = QLabsQBotPlatform(qlabs, verbose);

% Spawn slightly inside from the left wall (x = -4 + map_offset_x)
start_x = (-4 + map_offset_x) + 0.5;
location = [start_x, 0, 0; -1.35, 0.3, 0; -1.5, 0, 0; -1.5, 0, 0];

rotation = [0, 0, 0;    0,   0, 0;   0, 0, 90;  0, 0, -90];


% Select case 1 for the (0,0) spot
hQBot.spawn_id_degrees(0, location(caseNum, :), rotation(caseNum, :), [1, 1, 1], 1);
hQBot.possess(hQBot.VIEWPOINT_TRAILING);

    file_workspace = fullfile(getenv('RTMODELS_DIR'), 'QBotPlatform', 'QBotPlatform_Workspace.rt-win64');
    file_driver    = fullfile(getenv('RTMODELS_DIR'), 'QBotPlatform', 'qbot_platform_driver_virtual.rt-win64');


% Start RT models
pause(2)
system(['quarc_run -D -r -t tcpip://localhost:17000 ', file_workspace]);
pause(1)
system(['quarc_run -D -r -t tcpip://localhost:17000 ', file_driver, ' -uri tcpip://localhost:17098']);
pause(3)


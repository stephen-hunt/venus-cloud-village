% Planet Constants (Venus)
rho_o = 65.52;  %kg/m^3
H     = 15800;  %m
R     = 6052E3; %m
g_v   = 8.87;   %m/s^2
c     = .00016;

% Instantiate Venus
venus = Planet(R, g_v, H, rho_o, c);

% Vehicle Constants
m              = 50000; %kg
Cd             = 1.651;
D_HIAD         = 16;    %m
epsilon        = .8;
Cx             = 1.1;
Cd_tube        = .9;
D_tube         = 53.9;  %m
D_tube_initial = 27;    %m

% Instantiate the entry vehicle
entryVehicle = Vehicle(m, D_HIAD, Cd, epsilon, Cx, D_tube, D_tube_initial, Cd_tube);

% Trajectory Constants
LtD             = .8;
FPA             = 3.5;   %degrees
h_e             = 280E3; %m
V_e             = 11200; %m/s
initFloatDepAlt = 100E3; %m
landingAlt      = 60E3;  %m
trim            = 15;    %degrees
finFloatDepAlt  = 70E3;  %m

% Instantiate the simulator
mission = EDLSim(h_e, V_e, LtD, FPA, trim, venus, entryVehicle, initFloatDepAlt, finFloatDepAlt, landingAlt);

% Simulate various pieces of the EDL sequence
mission.plotVAvsH();
mission.plotDRDvsH();
mission.plotqvsH();
mission.printPeakWallTemp();
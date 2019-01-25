clear;
%%  VARIABLE

% common variable

G=6.674e-11;                %

%   Earth and Moon
EARTH_MASS=59760e20;        %kg
MOON_MASS=737e20;           %kg
EARTH_RADIUS=6378e3;        %meter
MOON_RADIUS=1738e3;         %meter
ERATH_DELTAT=600;           %s
MOON_DELTAT=600;            %s
EARTH_LOVE_NUMBER=0.299;
MOON_LOVE_NUMBER=0.03;
EARTH_DISSIPATION_FUNCTION=12;
MOON_DISSIPATION_FUNCTION=27;


%  Pluto and caron
PLUTO_MASS=870.3/G*1e9;         %kg
CHARON_MASS=101.4/G*1e9;        %kg
PLUTO_RADIUS=1153e3;        %meter
CHARON_RADIUS=606e3;        %meter
PLUTO_DELTAT=600;           %s
CHARON_DELTAT=600;          %s
PLUTO_LOVE_NUMBER=0.058;
CHARON_LOVE_NUMBER=0.006;
PLUTO_DISSIPATION_FUNCTION=100;
CHARON_DISSIPATION_FUNCTION=100;

%% savedata


save variable.mat
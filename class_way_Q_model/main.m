clear;clf;

load('variable.mat');

%% initiate EarthMoon Class
EarthMoon=TenTai;
EarthMoon.phip=0;
EarthMoon.phis=0;
EarthMoon.wp_div_n=27;
EarthMoon.ws_div_n=1;
EarthMoon.X=60;
EarthMoon.e=0.05;
EarthMoon.orbital_velocity=1022;
EarthMoon.n=2*pi/(27.3127*24*60*60);
EarthMoon.Mp=EARTH_MASS;
EarthMoon.Ms=MOON_MASS;
EarthMoon.Rp=EARTH_RADIUS;
EarthMoon.Rs=MOON_RADIUS;
EarthMoon.Dtp=ERATH_DELTAT;
EarthMoon.Dts=MOON_DELTAT;
EarthMoon.k2p=EARTH_LOVE_NUMBER;
EarthMoon.k2s=MOON_LOVE_NUMBER;
EarthMoon.Qp=EARTH_DISSIPATION_FUNCTION;
%EarthMoon.Qs=MOON_DISSIPATION_FUNCTION;

%% initiate PlutoClaron Class

PlutoClaron=TenTai;

PlutoClaron.phip=0;
PlutoClaron.phis=0;
PlutoClaron.wp_div_n=5.5;
PlutoClaron.ws_div_n=2;
PlutoClaron.X=4;
PlutoClaron.e=0;
PlutoClaron.orbital_velocity=0.003705;
PlutoClaron.n=2*pi/(6.3867*24*60*60);
PlutoClaron.Mp=PLUTO_MASS;
PlutoClaron.Ms=CHARON_MASS;
PlutoClaron.Rp=PLUTO_RADIUS;
PlutoClaron.Rs=CHARON_RADIUS;
PlutoClaron.Dtp=PLUTO_DELTAT;
PlutoClaron.Dts=CHARON_DELTAT;
PlutoClaron.k2p=PLUTO_LOVE_NUMBER;
PlutoClaron.k2s=CHARON_LOVE_NUMBER;
PlutoClaron.Qp=PLUTO_DISSIPATION_FUNCTION;
%PlutoClaron.Qs=CHARON_DISSIPATION_FUNCTION;


%%
processing=["orbital_evolution","oribtial_evolution_consider_i","orbital_evolution_Qmodel"];

switch processing(2)

    case {"orbital_evolution"}
    %% orbital_evolution Detltat model  /Solve and Plot  
        for i=0:0.1:0.1
            
            %[t x] = ode15s(@PlutoClaron.orbital_evolution_change_w,[1e-2 1e30],PlutoClaron.ini_of_evolution);
            %[t1 x1] = ode15s(@PlutoClaron.orbital_evolution,[1e-3 1e9],[5.5 2 4 0]);
            [t1 x1] = ode15s(@EarthMoon.orbital_evolution,[1e-3 1e12],[27 1 60 0]);
            a=(x1(length(x1),:));
            %[t2 x2] = ode15s(@PlutoClaron.orbital_evolution,[3.201779e+08 1e9],a);
            %[t2 x2] = ode15s(@PlutoClaron.orbital_evolution,[3.201779e+08 1e9],[1.0 1.000 37.6504 0]);
            %t=[t1;t2];
            %x=[x1;x2];
            t=[t1];
            x=[x1];
            tx=[t x];
            %[t x] = ode15s(@EarthMoon.orbital_evolution_change_w,[1e-4 1e11],[27 1 60 i]);
            subplot(2,2,3);
            title('wp/n')
            semilogx(t,x(:,1))
            hold on
            subplot(2,2,4);
            title('ws/n')
            semilogx(t,x(:,2))
            hold on
            subplot(2,2,1);
            title('a/Rp')
            semilogx(t,x(:,3))
            hold on
            subplot(2,2,2);
            title('e')
            semilogx(t,x(:,4))
            hold on
        
        end
        
%        [t2 x2] = ode15s(@PlutoClaron.orbital_evolution,[1e8 1e9],a);
        

    case {"oribtial_evolution_consider_i"}
    %% orbital_evolution Detltat model  /Solve and Plot  
        for i=0:0.1:0.1

            %[t x] = ode15s(@PlutoClaron.orbital_evolution_change_w,[1e-2 1e30],PlutoClaron.ini_of_evolution);
            [t x] = ode15s(@PlutoClaron.oribtial_evolution_consider_i,[1e-3 1.3e9],[5.5 2 4 0 1]);
            %[t x] = ode15s(@EarthMoon.orbital_evolution_change_w,[1e-4 1e11],[27 1 60 i]);
            subplot(3,2,3);
            title('wp/n')
            semilogx(t,x(:,1))
            hold on
            subplot(3,2,4);
            title('ws/n')
            semilogx(t,x(:,2))
            hold on
            subplot(3,2,1);
            title('a/Rp')
            semilogx(t,x(:,3))
            hold on
            subplot(3,2,2);
            title('e')
            semilogx(t,x(:,4))
            hold on
            subplot(3,2,5);
            title('i')
            semilogx(t,x(:,5))
            hold on
        end
                
    case {"orbital_evolution_Qmodel"}
    %% orbital_evolution_Qmodel and just consider one body. It is like "solar system dynamics".
        %[t x] = ode45(@EarthMoon.orbital_evolution_Qmodel,[1e-2 1e15],EarthMoon.ini_of_evolution_Qmodel);
        [t x] = ode45(@EarthMoon.orbital_evolution_Qmodel,[1e-2 1e15],[27 60]);
        subplot(2,1,1);
        title('wp/n')
        semilogx(t,x(:,1))
        subplot(2,1,2);
        title("a/Rp")
        semilogx(t,x(:,2))
    case []
end








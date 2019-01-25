classdef TenTai
    %% State properties 
    properties      
        phip;   phis;
        wp_div_n;     ws_div_n;
        orbital_distance_div_planet_radius;
        eccentricity;
        orbital_velocity;
        %n;
        ini_t;  fin_t;
    end
    properties (Dependent)
        ini_of_evolution;
        ini_of_motion;
        ini_of_evolution_Qmodel;
        ini_of_motion_Dt;
        tspan;
    end
    %% General properties
    properties      
        Mp;        Ms;
        Rp;        Rs;
        Dtp;       Dts;
        k2p;       k2s;
        Qp;        Qs;
    end
    
    properties (Constant)
        G=6.674e-11;
        ADt=10;
    end
    
    properties (Dependent)
        mu;
        Cp;        Cs;
        n;
    end
    
    methods      
        %% Dependent function include State properties and General properties
        %general
        function temp=get.mu(obj)
            temp=(obj.Mp*obj.Ms)/(obj.Mp+obj.Ms);
        end

        function temp = get.Cp(obj)
            temp =0.4*obj.Mp*obj.Rp^2;
        end
        function temp = get.Cs(obj)
            temp =0.4*obj.Ms*obj.Rs^2;
        end
        function temp = get.n(obj)
            temp =(obj.G*(obj.Mp+obj.Ms)*(4*obj.Rp)^3)^0.5;;
        end
        
        %state
        function temp = get.ini_of_evolution(obj)
            temp = [obj.wp_div_n obj.ws_div_n obj.orbital_distance_div_planet_radius obj.eccentricity];
        end
        function temp = get.ini_of_motion(obj)
            temp = [obj.orbital_distance_div_planet_radius,0,0,obj.orbital_velocity/obj.Rp 2*pi/(24*60*60)];
        end
        function temp = get.ini_of_evolution_Qmodel(obj)
            temp = [obj.wp_div_n obj.orbital_distance_div_planet_radius];
        end
        function temp = get.ini_of_motion_Dt(obj)
            temp = [obj.orbital_distance_div_planet_radius 0 0 obj.orbital_velocity/obj.Rp obj.phip obj.wp_div_n obj.phis obj.ws_div_n];
        end
        function temp = get.tspan(obj)
            temp = [obj.ini_t obj.fin_t];
        end
        
        
        
        %% orbital_evolution Detltat model
        %solve and plot
        function evolution(obj)
            
            [t,x] = ode15s(@obj.orbital_evolution,obj.tspan,obj.ini_of_evolution);
            subplot(2,2,3);
            title('wp/n')
            semilogx(t,x(:,1))          
            subplot(2,2,4);
            title('ws/n')
            semilogx(t,x(:,2))          
            subplot(2,2,1);
            title('a/Rp')
            semilogx(t,x(:,3))          
            subplot(2,2,2);
            title('e')
            semilogx(t,x(:,4))
          
        end

        %odefun
        
        function dydt=oribtial_evolution_consider_i(obj,~,x)
            n=2*pi/(6.3867*24*60*60);
            dydt(1)=(-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*cos(x(5))*x(1)-(1+(27/2)*x(4)^2))...
                    +3/2*x(1)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(2)=(-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
                    +3/2*x(2)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(3)=(x(3)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(4)=(x(4)*(27*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
            %dydt(5)=-3/2*((n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.mu/obj.Ms)^(1/2)*(x(1)*cos(2)+obj.ADt*x(2))*sin(x(5))*(1+7/2*x(4)^2);
            dydt(5)=-sin(x(5));
            
            dydt = dydt';
        end
        
        function dydt=orbital_evolution(obj,~,x)

            dydt(1)=-3*obj.G/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*x(1)-(1+(27/2)*x(4)^2));
            dydt(2)=-3*obj.G/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2));
            dydt(3)=x(3)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt));
            dydt(4)=x(4)*(27*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)+obj.ADt*x(2))-(1+obj.ADt)));

            dydt = dydt';
        end
        function dydt=orbital_evolution_change_w(obj,~,x)
            
            %{
            dydt(1)=(-3*(obj.G)/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*x(1)-(1+(27/2)*x(4)^2))...
                    +1/1.2*3/2*x(1)*(6*(obj.G)/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(2)=(-3*(obj.G)/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
                    +1/1.2*3/2*x(2)*(6*(obj.G)/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(3)=1/1.2*(x(3)*(6*(obj.G)/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(4)=(x(4)*(27*(obj.G)/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
            %}
            
            
            
            n=2*pi/(6.3867*24*60*60);
            dydt(1)=(-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*x(1)-(1+(27/2)*x(4)^2))...
                    +3/2*x(1)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(2)=(-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
                    +3/2*x(2)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(3)=(x(3)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(4)=(x(4)*(27*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
            
            
            
            %{
            f1=((1+3*x(4)^2+3/8*x(4)^4)/(1-x(4)^2)^(9/2));
            f2=((1+15/2*x(4)^2+45/8*x(4)^4+5/16*x(4)^6)/(1-x(4)^2)^6);
            f3=((1+31/2*x(4)^2+255/8*x(4)^4+185/16*x(4)^6+25/64*x(4)^8)/(1-x(4)^2)^(15/2));
            f4=((1+3/2*x(4)^2+1/8*x(4)^4)/(1-x(4)^2)^5);
            f5=((1+15/4*x(4)^2+15/8*x(4)^4+5/64*x(4)^6)/(1-x(4)^2)^(13/2));
            %}
            %{
            dydt(1)=-3*obj.G/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*(f1*x(1)-f2)...
                    +3/2*x(1)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(f2*(x(1)+obj.ADt*x(2))-f3*(1+obj.ADt));
            dydt(2)=-3*obj.G/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*(f1*x(1)-f2)...
                    +3/2*x(2)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(f2*(x(1)+obj.ADt*x(2))-f3*(1+obj.ADt));
            dydt(3)=x(3)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(f2*(x(1)+obj.ADt*x(2))-f3*(1+obj.ADt));
            dydt(4)=x(4)*(27*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((f4*11/18*(x(1)+obj.ADt*x(2))-f5*(1+obj.ADt)));
            %}
            %{
            n=2*pi/(6.3867*24*60*60);
            dydt(1)=-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*(f1*x(1)-f2)...
                    +3/2*x(1)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(f2*(x(1)+obj.ADt*x(2))-f3*(1+obj.ADt));
            dydt(2)=-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*(f1*x(1)-f2)...
                    +3/2*x(2)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(f2*(x(1)+obj.ADt*x(2))-f3*(1+obj.ADt));
            dydt(3)=x(3)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(f2*(x(1)+obj.ADt*x(2))-f3*(1+obj.ADt));
            dydt(4)=x(4)*(27*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((f4*11/18*(x(1)+obj.ADt*x(2))-f5*(1+obj.ADt)));
            %}
            %{
            n=2*pi/(6.3867*24*60*60);
            dydt(1)=-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*(((1+3*x(4)^2+3/8*x(4)^4)/(1-x(4)^2)^(9/2))*x(1)-((1+15/2*x(4)^2+45/8*x(4)^4+5/16*x(4)^6)/(1-x(4)^2)^6))...
                    +3/2*x(1)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(((1+15/2*x(4)^2+45/8*x(4)^4+5/16*x(4)^6)/(1-x(4)^2)^6)*(x(1)+obj.ADt*x(2))-((1+31/2*x(4)^2+255/8*x(4)^4+185/16*x(4)^6+25/64*x(4)^8)/(1-x(4)^2)^(15/2))*(1+obj.ADt));
            dydt(2)=-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*(((1+3*x(4)^2+3/8*x(4)^4)/(1-x(4)^2)^(9/2))*x(1)-((1+15/2*x(4)^2+45/8*x(4)^4+5/16*x(4)^6)/(1-x(4)^2)^6))...
                    +3/2*x(2)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(((1+15/2*x(4)^2+45/8*x(4)^4+5/16*x(4)^6)/(1-x(4)^2)^6)*(x(1)+obj.ADt*x(2))-((1+31/2*x(4)^2+255/8*x(4)^4+185/16*x(4)^6+25/64*x(4)^8)/(1-x(4)^2)^(15/2))*(1+obj.ADt));
            dydt(3)=x(3)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*(((1+15/2*x(4)^2+45/8*x(4)^4+5/16*x(4)^6)/(1-x(4)^2)^6)*(x(1)+obj.ADt*x(2))-((1+31/2*x(4)^2+255/8*x(4)^4+185/16*x(4)^6+25/64*x(4)^8)/(1-x(4)^2)^(15/2))*(1+obj.ADt));
            dydt(4)=x(4)*(27*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((((1+3/2*x(4)^2+1/8*x(4)^4)/(1-x(4)^2)^5)*11/18*(x(1)+obj.ADt*x(2))-((1+15/4*x(4)^2+15/8*x(4)^4+5/64*x(4)^6)/(1-x(4)^2)^(13/2))*(1+obj.ADt)));
            %}
            
            dydt = dydt';
        end
        %% orbital_motion Detltat model and just consider planet
        %solve and plot
        function motion_evolution(obj)
        
            [t,x] = ode15s(@obj.orbital_motion,obj.tspan,obj.ini_of_motion);
            subplot(3,2,1);
            title('x/Rp')
            semilogx(t,x(:,1))
            subplot(3,2,2);
            title("x'/Rp")
            semilogx(t,x(:,2))
            subplot(3,2,3);
            title('y/Rp')
            semilogx(t,x(:,3))
            subplot(3,2,4);
            title("y'/R")
            semilogx(t,x(:,4))
            subplot(3,2,5);
            title('wp/n')
            semilogx(t,x(:,5))
        end
        %odefun
        function dydt=orbital_motion(obj,~,x)

            r2=x(1)^2+x(3)^2;

            dydt(1) = x(2);
            dydt(2) = -3*obj.k2p*obj.G*obj.Ms*obj.Rp^-3/r2^4*(x(1)+obj.Dtp*(2*(x(1)*(1)*x(2)+x(1)*x(3)*x(4))/r2+x(5)*x(3)+x(2)));
            dydt(3) = x(4);
            dydt(4) = -3*obj.k2p*obj.G*obj.Ms*obj.Rp^-3/r2^4*(x(3)+obj.Dtp*(2*(x(1)*(2)*x(3)+x(3)*x(3)*x(4))/r2-x(5)*x(1)+x(4))); 
            dydt(5) = (-3*obj.k2p*obj.G*obj.Ms^2*obj.Rp^-1/r2^3*obj.Dtp*(x(5)+(-x(4)*x(1)+x(3)*x(2))/r2))/obj.Cp;

            dydt = dydt';
        end
        
        %% orbital_evolution_Qmodel and just consider one body. It is like "solar system dynamics".
        function evolution_Qmodel(obj)
            [t,x] = ode15s(@obj.orbital_evolution_Qmodel,obj.tspan,obj.ini_of_evolution_Qmodel);
            semilogx(t,x)
        end
        %odefun
        function dydt=orbital_evolution_Qmodel(obj,~,x)

            n=2*pi/(27*24*60*60);
            
            dydt(1)=-3*(n*(x(2)*obj.Rp)^3/(obj.Mp+obj.Ms)*obj.Mp^2/(2*obj.Cp))*(obj.k2p/obj.Qp)*(obj.Rp^-1/x(2)^6)*sign(x(1)-1)+3/2*x(1)*3*n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(2)^-5*sign(x(1)-1);
            dydt(2)=3*n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(2)^-4*sign(x(1)-1);

            dydt = dydt';
        end
        
        %% oribtital motion Dt_model without considering eccentricity.
        
        
        %odefun
        function dydt=motion_Dt(obj,~,x)

            r2=x(1)^2+x(3)^2;

            dydt(1) = x(2);
            dydt(2) = (-3*obj.k2p*obj.G*obj.Ms^2*obj.Rp^-3/r2^4*(x(1)+obj.Dtp*(2*(x(1)*(1)*x(2)+x(1)*x(3)*x(4))/r2+x(5)*x(3)+x(2)))...
                      -3*obj.k2s*obj.G*obj.Mp^2*(obj.Rs^5/obj.Rp^-8)/r2^4*(x(1)+obj.Dts*(2*(x(1)*(1)*x(2)+x(1)*x(3)*x(4))/r2+x(6)*x(3)+x(2))))/obj.mu;
            dydt(3) = x(4);
            dydt(4) = (-3*obj.k2p*obj.G*obj.Ms^2*obj.Rp^-3/r2^4*(x(3)+obj.Dtp*(2*(x(1)*(2)*x(3)+x(3)*x(3)*x(4))/r2-x(5)*x(1)+x(4)))...
                      -3*obj.k2s*obj.G*obj.Mp^2*(obj.Rs^5/obj.Rp^-8)/r2^4*(x(3)+obj.Dts*(2*(x(1)*(2)*x(3)+x(3)*x(3)*x(4))/r2-x(6)*x(1)+x(4))))/obj.mu;
            dydt(5) = (-3*obj.k2p*obj.G*obj.Ms^2*obj.Rp^-1/r2^3*obj.Dtp*(x(5)+(-x(4)*x(1)+x(3)*x(2))/r2))/obj.Cp;
            dydt(6) = (-3*obj.k2s*obj.G*obj.Mp^2*(obj.Rs^5/obj.Rp^6)/r2^3*obj.Dtp*(x(6)+(-x(4)*x(1)+x(3)*x(2))/r2))/obj.Cs;

            dydt = dydt';
        end
        
        %% oribtital motion Dt_model without considering eccentricity.
        
        
        %odefun
        function dydt=motion_Dt_myway(obj,~,x)

            r2=x(1)^2+x(3)^2;

            dydt(1) = x(2);
            dydt(2) = -3*obj.k2p*obj.G*(obj.Ms^2/obj.Mp)*obj.Rp^-3/r2^4*(obj.Dtp*(2*(x(1)*(1)*x(2)+x(1)*x(3)*x(4))/r2+x(5)*x(3)+x(2)))...
                      -3*obj.k2s*obj.G*(obj.Mp^2/obj.Ms)*(obj.Rs^5/obj.Rp^-8)/r2^4*(obj.Dts*(2*(x(1)*(1)*x(2)+x(1)*x(3)*x(4))/r2+x(6)*x(3)+x(2)));
            dydt(3) = x(4);
            dydt(4) = -3*obj.k2p*obj.G*(obj.Ms^2/obj.Mp)*obj.Rp^-3/r2^4*(obj.Dtp*(2*(x(1)*(2)*x(3)+x(3)*x(3)*x(4))/r2-x(5)*x(1)+x(4)))...
                      -3*obj.k2s*obj.G*(obj.Mp^2/obj.Ms)*(obj.Rs^5/obj.Rp^-8)/r2^4*(obj.Dts*(2*(x(1)*(2)*x(3)+x(3)*x(3)*x(4))/r2-x(6)*x(1)+x(4)));
            dydt(5) = (-3*obj.k2p*obj.G*obj.Ms^2*obj.Rp^-1/r2^3*obj.Dtp*(x(5)+(-x(4)*x(1)+x(3)*x(2))/r2))/obj.Cp;
            dydt(6) = (-3*obj.k2s*obj.G*obj.Mp^2*(obj.Rs^5/obj.Rp^6)/r2^3*obj.Dtp*(x(6)+(-x(4)*x(1)+x(3)*x(2))/r2))/obj.Cs;

            dydt = dydt';
        end
        
        
        
    end
end
classdef TenTai
    %% State properties 
    properties      
        phip;   phis;
        wp_div_n;     ws_div_n;
        X;              %orbital_distance_div_planet_radius
        e;              %eccentricity
        orbital_velocity;
        n;
    end
    %% General properties
    properties      
        Mp;        Ms;
        Rp;        Rs;
        Dtp;       %Dts;
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
        Dts;
    end
    
    methods      
        %% Dependent function include State properties and General properties
        %general
        function temp=get.mu(obj)
            temp=(obj.Mp*obj.Ms)/(obj.Mp+obj.Ms);
        end
        function temp = get.Cp(obj)
            temp =0.383*obj.Mp*obj.Rp^2;
        end
        function temp = get.Cs(obj)
            temp =0.4*obj.Ms*obj.Rs^2;
        end
        function temp = get.Dts(obj)
            temp =obj.Dtp*obj.ADt*(obj.k2p/obj.k2s)*(obj.Ms/obj.Mp)^2*(obj.Rp/obj.Rs)^5;
        end
        
        %% odefun
        % orbital_evolution Detltat model
        
         function dydt=orbital_evolution(obj,~,x)            
       
%             dydt(1)=(-3*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*x(1)-(1+(27/2)*x(4)^2))...
%                     +3/2*x(1)*(6*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
%             dydt(2)=(-3*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
%                     +3/2*x(2)*(6*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
%             dydt(3)=(x(3)*(6*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
%             dydt(4)=(x(4)*(27*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
%             
            dydt(1)=(-3*obj.G/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*x(1)-(1+(27/2)*x(4)^2))...
                    +3/2*x(1)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(2)=(-3*obj.G/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
                    +3/2*x(2)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(3)=(x(3)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(4)=(x(4)*(27*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
            
            
            %n=2*pi/(6.3867*24*60*60);
            %dydt(1)=(-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*x(1)-(1+(27/2)*x(4)^2))...
            %        +3/2*x(1)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            %dydt(2)=(-3*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
            %        +3/2*x(2)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            %dydt(3)=(x(3)*(6*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            %dydt(4)=(x(4)*(27*(n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
            
            dydt = dydt';
            
        end
        
        function dydt=oribtial_evolution_consider_i(obj,~,x)

            dydt(1)=(-3*obj.G/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*cos(x(5))*x(1)-(1+(27/2)*x(4)^2))...
                    +3/2*x(1)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(2)=(-3*obj.G/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
                    +3/2*x(2)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(3)=(x(3)*(6*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
            dydt(4)=(x(4)*(27*obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
            dydt(5)=-3/2*(obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.mu/obj.Ms)^(1/2)*(x(1)*cos(2)+obj.ADt*x(2))*sin(x(5))*(1+7/2*x(4)^2);
            
            dydt = dydt';
        end
        
        % orbital_evolution_Qmodel and just consider one body. It is like "solar system dynamics".
        function dydt=orbital_evolution_Qmodel(obj,~,x)
                   
            dydt(1)=-3*(obj.n*(x(2)*obj.Rp)^3/(obj.Mp+obj.Ms)*obj.Mp^2/(2*obj.Cp))*(obj.k2p/obj.Qp)*(obj.Rp^-1/x(2)^6)*sign(x(1)-1)+3/2*x(1)*3*n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(2)^-5*sign(x(1)-1);
            dydt(2)=3*obj.n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(2)^-4*sign(x(1)-1);

            dydt = dydt';
        end
        %orbital evolution Qmdel and consider i. it's an planet produce tidal force.
        function dydt=oribtal_evolution_earth(obj,~,x)
              
           dydt(1)=3*obj.n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(1)^-4;
           dydt(2)=19/8*x(2)*3*obj.n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(1)^-5;
           dydt(3)=(-1/4)*sin(x(3))*3*obj.n*(obj.k2p/obj.Qp)*(obj.Ms/obj.Mp)*x(1)^-4;
           
        end
    end
end
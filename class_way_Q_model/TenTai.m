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
        Dtp;       Dts;
        k2p;       k2s;
        Qp;        %Qs;
    end
    
    properties (Constant)
        G=6.674e-11;
        ADt=10;
        AQ=1.15;
    end
    
    properties (Dependent)
        mu;
        Cp;        Cs;
        Qs;
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
        function temp = get.Qs(obj)
            temp =obj.Qp/obj.AQ*(obj.k2s/obj.k2p)*(obj.Mp/obj.Ms)^2*(obj.Rs/obj.Rp)^5;
        end
        
        %% odefun
        % orbital_evolution Detltat model
        
         function dydt=orbital_evolution(obj,~,x)
            error=10^-4;
            if x(1)>3/2
                Dp=15/2;Ep=51/4;Fp=57/8;
                Sp=1;
            elseif x(1)==3/2
                Dp=-19/4;Ep=-45/8;Fp=-33/16;
                Sp=1;
            elseif (x(1)<3/2) && (x(1)>1+error)
                Dp=-17;Ep=-24;Fp=-45/4;
                Sp=1;
            elseif abs(x(1)-1)<=error
                Dp=-12;Ep=-19;Fp=-21/2;
                Sp=0;
            elseif (x(1)<(1-error)) && (x(1)>1/2)
                Dp=-7;Ep=-14;Fp=-39/4;
                Sp=-1;
            else
                Dp=-7.5;Ep=-14.25;Fp=-75/8;
                Sp=-1;
            end

            if x(2)>3/2
                Ds=15/2;Es=51/4;Fs=57/8;
                Ss=1;
                v=1;
            elseif x(2)==3/2
                Ds=-19/4;Es=-45/8;Fs=-33/16;
                Ss=1;
                v=1;
            elseif (x(2)<3/2) && (x(2)>1+error)
                Ds=-17;Es=-24;Fs=-45/4;
                Ss=1;
                v=1;
            elseif abs(x(2)-1)<=error
                Ds=-12;Es=-19;Fs=-21/2;
                Ss=0;
                v=0;
            elseif (x(2)<(1-error)) && (x(2)>1/2)
                Ds=-7;Es=-14;Fs=-39/4;
                Ss=-1;
                v=0;
            else
                Ds=-7.5;Es=-14.25;Fs=-75/8;
                Ss=-1;
                v=0;
                
            end
            
%              Dp=15/2;Ep=51/4;Fp=57/8;
%              Ds=15/2;Es=51/4;Fs=57/8;
            

%             dydt(1)=(-3*(obj.G)/((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))/(2*obj.Cp*x(3)^6)*obj.k2p/obj.Qp*obj.Ms^2*(obj.Rp^-1)*(sign(x(1)-1)+x(4)^2*Dp)...
%                     +3/2*x(1)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(sign(x(1)-1)+obj.AQ*sign(x(2)-1)+x(4)^2*(Ep+obj.AQ*Es))))*31536000;
%             dydt(2)=(-3*(obj.G)/((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))/(2*obj.Cs*x(3)^6)*obj.k2s/obj.Qs*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*((x(2)-1)+x(4)^2*Ds)...
%                     +3/2*x(2)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(sign(x(1)-1)+obj.AQ*sign(x(2)-1)+x(4)^2*(Ep+obj.AQ*Es))))*31536000;
%             dydt(3)=x(3)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(sign(x(1)-1)+obj.AQ*sign(x(2)-1)+x(4)^2*(Ep+obj.AQ*Es)))*31536000;
%             dydt(4)=x(4)*(((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Fp+obj.AQ*Fs))*31536000;
%             

            
            dydt(1)=(-3*(obj.G)/((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))/(2*obj.Cp*x(3)^6)*obj.k2p/obj.Qp*obj.Ms^2*(obj.Rp^-1)*(Sp+x(4)^2*Dp)...
                    +3/2*x(1)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Sp+obj.AQ*Ss+x(4)^2*(Ep+obj.AQ*Es))))*31536000;
            dydt(2)=(-3*(obj.G)/((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))/(2*obj.Cs*x(3)^6)*obj.k2s/obj.Qs*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*(Ss+x(4)^2*Ds)...
                    +3/2*x(2)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Sp+obj.AQ*Ss+x(4)^2*(Ep+obj.AQ*Es))))*31536000*v;
            dydt(3)=x(3)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Sp+obj.AQ*Ss+x(4)^2*(Ep+obj.AQ*Es)))*31536000;
            dydt(4)=x(4)*(((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Fp+obj.AQ*Fs))*31536000;

% 
% 
% 
%             if dydt(1)==0
% 
%             
%             
%                 end
            
         
         
%             dydt(1)=(-3*(obj.n*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(2*obj.Cp*x(3)^6)*obj.k2p/obj.Qp*obj.Ms^2*(obj.Rp^-1)*(sign(x(1)-1))...
%                     +3/2*x(1)*(3*obj.n*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(sign(x(1)-1)+obj.AQ*sign(x(2)-1))))*31536000;
%             dydt(2)=(-3*(obj.n*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(2*obj.Cs*x(3)^6)*obj.k2s/obj.Qs*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((x(2)-1))...
%                     +3/2*x(2)*(3*obj.n*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(sign(x(1)-1)+obj.AQ*sign(x(2)-1))))*31536000;
%             dydt(3)=x(3)*(3*obj.n*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(sign(x(1)-1)+obj.AQ*sign(x(2)-1)))*31536000;
%             dydt(4)=0;


            dydt = dydt';

        end
        
        function dydt=oribtial_evolution_consider_i(obj,~,x)

%             dydt(1)=(-3*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cp*x(3)^6)*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-1)*((1+(15/2)*x(4)^2)*cos(x(5))*x(1)-(1+(27/2)*x(4)^2))...
%                     +3/2*x(1)*(6*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
%             dydt(2)=(-3*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.Cs*x(3)^6)*obj.k2s*obj.Dts*obj.Mp^2*(obj.Rp^5/obj.Rs^6)*((1+(15/2)*x(4)^2)*x(2)-(1+(27/2)*x(4)^2))...
%                     +3/2*x(2)*(6*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
%             dydt(3)=(x(3)*(6*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((1+27/2*x(4)^2)*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+23*x(4)^2)*(1+obj.ADt)))*31536000;
%             dydt(4)=(x(4)*(27*(obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.Rp^-3)*((11/18*(x(1)*cos(x(5))+obj.ADt*x(2))-(1+obj.ADt))))*31536000;
%             dydt(5)=-3/2*((obj.n^2*(x(3)*obj.Rp)^3/(obj.Mp+obj.Ms))/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.mu/obj.Ms)^(1/2)*(x(1)*cos(2)+obj.ADt*x(2))*sin(x(5))*(1+7/2*x(4)^2);

            error=10^-4;
            if x(1)>3/2
                Dp=15/2;Ep=51/4;Fp=57/8;
                Sp=1;
            elseif x(1)==3/2
                Dp=-19/4;Ep=-45/8;Fp=-33/16;
                Sp=1;
            elseif (x(1)<3/2) && (x(1)>1+error)
                Dp=-17;Ep=-24;Fp=-45/4;
                Sp=1;
            elseif abs(x(1)-1)<=error
                Dp=-12;Ep=-19;Fp=-21/2;
                Sp=0;
            elseif (x(1)<(1-error)) && (x(1)>1/2)
                Dp=-7;Ep=-14;Fp=-39/4;
                Sp=-1;
            else
                Dp=-7.5;Ep=-14.25;Fp=-75/8;
                Sp=-1;
            end

            if x(2)>3/2
                Ds=15/2;Es=51/4;Fs=57/8;
                Ss=1;
                v=1;
            elseif x(2)==3/2
                Ds=-19/4;Es=-45/8;Fs=-33/16;
                Ss=1;
                v=1;
            elseif (x(2)<3/2) && (x(2)>1+error)
                Ds=-17;Es=-24;Fs=-45/4;
                Ss=1;
                v=1;
            elseif abs(x(2)-1)<=error
                Ds=-12;Es=-19;Fs=-21/2;
                Ss=0;
                v=1;
            elseif (x(2)<(1-error)) && (x(2)>1/2)
                Ds=-7;Es=-14;Fs=-39/4;
                Ss=-1;
                v=1;
            else
                Ds=-7.5;Es=-14.25;Fs=-75/8;
                Ss=-1;
                v=1;
                
            end
            
     
            dydt(1)=(-3*(obj.G)/((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))/(2*obj.Cp*x(3)^6)*obj.k2p/obj.Qp*obj.Ms^2*(obj.Rp^-1)*(Sp+x(4)^2*Dp)...
                    +3/2*x(1)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Sp+obj.AQ*Ss+x(4)^2*(Ep+obj.AQ*Es))))*31536000;
            dydt(2)=(-3*(obj.G)/((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))/(2*obj.Cs*x(3)^6)*obj.k2s/obj.Qs*obj.Mp^2*(obj.Rs^5/obj.Rp^6)*(Ss+x(4)^2*Ds)...
                    +3/2*x(2)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Sp+obj.AQ*Ss+x(4)^2*(Ep+obj.AQ*Es))))*31536000*v;
            dydt(3)=x(3)*(3*((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Sp+obj.AQ*Ss+x(4)^2*(Ep+obj.AQ*Es)))*31536000;
            dydt(4)=x(4)*(((obj.G*(obj.Mp+obj.Ms)/(x(3)*obj.Rp)^3)^(1/2))*obj.k2p/obj.Qp*(obj.Ms/obj.Mp)/x(3)^5*(Fp+obj.AQ*Fs))*31536000;
            dydt(5)=-3/2*(obj.G/(obj.mu*x(3)^8))*obj.k2p*obj.Dtp*obj.Ms^2*(obj.mu/obj.Ms)^(1/2)*(1*cos(2)+obj.ADt*1)*tan(x(5))*(1+7/2*x(4)^2);

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
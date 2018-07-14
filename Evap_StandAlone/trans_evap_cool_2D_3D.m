function trans_evap_cool_2D_3D
%% If you found this code helpful, please cite the following journal article.

% A. Pathak, M. Raessi, Steady-state and transient solutions to drop evaporation 
% in a finite domain: Alternative benchmarks to the d2 law, Int. J. Heat Mass Transfer (2018),
% https://doi.org/10.1016/j.ijheatmasstransfer.2018.06.071

% The details of the governing equations and numerical descritization are
% also given in the article.

% The test cases presented in the Journal article are also setup in the
% current code. Any new fluid or test case could be added by modifying
% the code. See the  following comment.
  
% This code has been tested on MATLAB (2017a) version. MATLAB is a
% recommended software to run this code. The free software Octave could
% also be used, however we found it to be slower. 

% Make sure you add AdaptiveMesh, figstate, and SteamTable to the path of
% your Matlab. In case you are unable to do it, copy all the files in those
% directories to the current directory.

%%
    % this is the main function for modeling 1D radial transient evaporation
	% change the following to switch between problems that are already setup:
	% (1) Ri       [Initial radius]
	% (2) Pamb     [Ambient pressure]
	% (3) fluid    [fluid type: 1 is for water, 2 is for n-heptane, 
    %              3 is for isothermal test case presented in the journal article ]
	% (4) cno:     [which case to run for each fluid]
    % (5) cyl      [=1 for circle and =2 for sphere]
    % (6) uniform  [=1 for a uniform grid; =0 for a non-uniform adaptive grid]
	
	% To add a new fluid / case , add the following in addition to what is 
	% listed aboveakc
    
	% (7) fine tune twfin, final simulation time
	% (8) fine tune xlim, ylim; its better to remove them altogether initially, see the evolution
	%     and then add them.
	% (9) Add a saturation pressure equation, in Yint (Ts)
	% 	
%%
    close all;
	%pkg load odepkg;
	%debug_on_warning(1);
	%global parameters
	global lN gN lc lb zeta zeta_l uniform X0 
	global rho_l rho_g D k_l k_g Cp_l Cp_g h_lg Ti_g Pamb Mf Ma
	global Yr Lr Tr beta_d t_d
	global R Tl Tg gY
	global fluid cyl isothermal zopt
    zopt = optimset('TolX',1.e-6);
	 
	% Ri = 150.0e-6; %[m] initial radius of drop; % also length scale
	% Pamb = 1.01325; % [bar]
	Ri = 2.5e-6; %[m]
	Pamb = 28.6; % [bar]
    %Ri = 200.e-6; %[m]
	%Pamb = 1.01325; % [bar]
    cyl = 2; % 1 for circular drop; 2 for spherical drop
	uniform = 0; % 0 and 1 for non-uniform and uniform meshes respectively
    isothermal = 0;  
	fluid = 2; % fluid = 1 for water, fluid = 2 for n-heptane 
	cno = zeros(2,1); cno(2)=1;
	
	% mesh
	LA = 0.0; LB = 1.0; lN = 11; lh = (LB - LA)/(lN-1); 
	GA = 1.0; GB = 2.0; gN = 31; gh = (GB - GA)/(gN-1); X0 = 4.0;
	% outer domain radius X0 is 4 times the initial drop radius
	
	if(uniform == 1)
		zeta = linspace(GA,GB,gN); % zeta is the transformed mesh
		zeta_l = linspace(LA,LB,lN);
	else
		zeta = linspace(GA,GB,gN); % zeta is the transformed mesh
		zeta_l = linspace(LA,LB,lN);
%       you could also use a non-uniform grid from the beginning. But you
%       would need a file to read from. This file would contain the zeta
%       positions of each grid point, starting from the drop center to the
%       domain boundary

% 		if( fluid==1 )
% 		 zeta_read = dlmread('mesh_read_11_31');
%          zeta_l = zeta_read(1:lN);
%          zeta   = zeta_read(lN+1:end);
% 		elseif ( fluid==2 )
% 		 zeta_read = dlmread('mesh_read_11_31_n_heptane');
%          zeta_l = zeta_read(1:lN);
% 		 zeta   = zeta_read(lN+1:end);
%         end
	end

	order=2; % currently order is always 2; one can program higher (or lower) order 
			     % discretizations and use this flag in the code	
	
	% fluid properties %taken from Tanguy 2007
	if( fluid == 1 ) 
		rho_l = 1000.0; % kg / m3
		rho_g = 1.226;  % kg / m3
		D     = 2.e-5;  % m2 / s
		Ti_l  = 353.0;  % K
		Ti_g  = 573.0*cno(1)+373.0*cno(2);  % K
		k_l   = 0.6;    % W / m-K
		k_g   = 0.046;  % W / m-K
		Cp_l  = 4180.0; % J / kg-K
		Cp_g  = 1000.0; % J / kg-K
		h_lg  = 2.3e6;  % J / kg
		Mf    = 18.0;   % kg / kmol
		Ma    = 29.0;   % kg / kmol 
	elseif( fluid==2 )
		rho_l = 626.7;  % kg / m3
		rho_g = 23.38*cno(1) + 17.51*cno(2);  % kg / m3
		D     = 4.26e-7*cno(1) + 6.77e-7*cno(2); % m2 / s
		Ti_l  = 363.0;  % K
		Ti_g  = 423.0*cno(1) + 563.0*cno(2);  % K
		k_l   = 0.1121; % W / m-K
		k_g   = 0.0356*cno(1) + 0.04428*cno(2); % W / m-K
		Cp_l  = 2505.0; % J / kg-K
		Cp_g  = 1036.0*cno(1) + 1053.0*cno(2); % J / kg-K
		h_lg  = 3.23e5; % J / kg
		Mf    = 100.2;  % kg / kmol
		Ma    = 29.0;   % kg / kmol	
    elseif( fluid == 3 ) 
		rho_l = 10.0; % kg / m3
		rho_g = 1.0;  % kg / m3
		D     = 2.e-3;  % m2 / s
		Ti_l  = 469.715557507987;  % K
		Ti_g  = 669.715557507987;  % K
		k_l   = 0.6;    % W / m-K
		k_g   = 0.295;  % W / m-K
		Cp_l  = 500.0; % J / kg-K
		Cp_g  = 287.2; % J / kg-K
		h_lg  = 7649.90219031717;  % J / kg
		Mf    = 18.0;   % kg / kmol
		Ma    = 29.0;   % kg / kmol         
	end
	% initial condition
	gY = zeros(gN       ,1); % vapor mass fraction on the proper mesh (zeta)
	Tl = zeros(lN       ,1);
	Tg = zeros(gN       ,1);
	Xi = zeros(2*gN+lN-4,1); % "variable" vector: vapor frac (gN-2), radius(1) , liq. temp (lN-1), gas temp (gN-2)
	
	lb = lN + gN - 2; % last index of "variable" Tl
	lc = lb + gN - 2; % last index of "variable" Tg
	
	% reference variables for non-dimensionalizing
	Yr  = Yint( Ti_l ) ;
	Lr  = Ri ;
	Tr  = Ti_l ;
	beta_d = max ( [ D , k_l/(rho_l*Cp_l) , k_g/(rho_g*Cp_g) ] ) ;
	t_d = Lr^2 / beta_d ;
	
	% local variables corresponding to vapor fraction and temperatures in the two domains
	gY(1        ) = Yint(Ti_l) / Yr ;
	Tl(1   :lN  ) = Ti_l / Tr ;
	Tg(1        ) = Ti_l / Tr ;
	Tg(2   :gN  ) = Ti_g / Tr ;
	
% Smooth Initial solution
% Temperature

% The following few lines show an example, if you wish to 
% setup a smooth initial condition[these lines are specific to fluid=3].
    
%     Tl(1   :lN  ) = Ti_l / Tr ;
%     Ts = Ti_l;
%     A = log( 1+Cp_g*(Ti_g-Ts)/h_lg )/log(X0);
%     r = zeros(gN,1);
%     for k=1:gN
%         r(k) = (zeta(k)-1.0)*(X0*Lr-Ri) + Ri;
%     end
%     
%     Tg = (Ts/Tr) + ( h_lg / (Cp_g*Tr) ) * ( (r/Ri).^A-1.0 );
%     
% % Vapor mass fraction
%     Ys = Yint(Ti_l); 
%     A = log( 1.0/(1-Ys) ) / log(X0);
%     gY = 1.0/Yr - ( (1-Ys)/Yr )*( r/Ri ).^A;
% % Smooth initialization ends here

	% The input variable populated using the local variables.
	Xi(1   :gN-2) = gY(2:gN-1) ;
	Xi(gN-1     ) = Ri / Lr ;
	Xi(gN  :lb  ) = Tl(1:lN-1) ;
	Xi(lb+1:lc  ) = Tg(2:gN-1) ;
		
	% ADAPTIVE MESH REFINEMENT NEAR THE INTERFACE after every few cycles
	% AMR kicks only after Init_Step; Initially the mesh is uniform; it
	% could be made non-uniform the beginning; see the comment above
    
	jj = 1;
	
	twfin = ( 2.0 * (rho_l / rho_g ) * Ri^2 / D ) * cno(1) + 5.5 * cno(2);
	if( fluid == 2 ) 
		twfin = (1.5 * (rho_l / rho_g ) * Ri^2 / D )*cno(1) + 1.5e-4 * cno(2); 
    elseif( fluid==3 )
        twfin = (rho_l / rho_g ) * Ri^2 / D ; twfin = 1.814e-4; twfin = 1.4e-4;
    end
	
	fig_handle = figure(1);
	maxfig(fig_handle,1);
    Init_Step = 1.e-4*twfin; % Initial Step
    N_amr = 11;              % number of times AMR is performed during simulation;
    a_amr = 2*Init_Step; b_amr = (1/(N_amr-1))*log(twfin/a_amr);
    tf1 = Init_Step; tf2 = a_amr*exp(b_amr*(0:N_amr-1));
    %tf = [0.0:Init_Step/t_d:Init_Step/t_d 0.1*twfin/t_d:0.1*twfin/t_d:twfin/t_d];
    tf = [tf1 tf2];
    tf = tf/t_d;
    
	%tf = 0.0:0.1*twfin/t_d:twfin/t_d;
	for j = 1:length(tf)-1
		if( j==1 )
			vopt = odeset ('InitialStep',0.0005*twfin/t_d,'MaxStep',0.05*twfin/t_d,'AbsTol',1.0e-6,'RelTol',1.e-3);
		else
			vopt = odeset ('InitialStep',0.05*twfin/t_d,'MaxStep',0.05*twfin/t_d,'AbsTol',1.0e-6,'RelTol',1.e-3);
		end
		% [t X] = ode23s_local(@func_ode,[tf(j) tf(j+1)],Xi,vopt);
		 [t X] = ode23s(@func_ode,[tf(j) tf(j+1)],Xi,vopt);
		
		D2 = zeros(length(t),1);
		intT = zeros(length(t),1);
		
		for i=1:length(t)
			gY(2:gN-1) = X(i,1   :gN-2);
			R          = X(i,gN-1     );
			Tl(1:lN-1) = X(i,gN  :lb  );
			Tg(2:gN-1) = X(i,lb+1:lc  );
			
			Tint = fzero(@func_Tint,Tl(lN-1));
			intT(i) = Tint;
			Tl(lN) = Tint;
			Tg(1 ) = Tint;
			gY(1 ) = Yint(Tint*Tr)/Yr;
			
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compute the transient term %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dr_dt = dR_dt(gY,R);
            for ii=2:gN-1 % i is mesh index
                I=ii-1; % the index for unknowns (vector x)
                ri = (zeta(ii)-1)*(X0-R) + R;
                if(ii==2)
                    x1 = [zeta(ii-1), zeta(ii), zeta(ii+1), zeta(ii+2)];
                    y1 = [gY(ii-1), gY(ii), gY(ii+1), gY(ii+2)];
                    grad2_i = dy2dx2_for(x1,y1);
                else
                    x1 = [zeta(ii-2), zeta(ii-1), zeta(ii), zeta(ii+1)];
                    y1 = [gY(ii-2), gY(ii-1), gY(ii), gY(ii+1)];
                    grad2_i = dy2dx2_back(x1,y1);
                end
                x1 = [zeta(ii-1), zeta(ii), zeta(ii+1)];
                y1 = [gY(ii-1), gY(ii), gY(ii+1)];
                grad_i = dydx(x1,y1);
                diff = (grad2_i/(X0-R)^2) + cyl*grad_i/(ri*(X0-R));
                conv = (rho_l/rho_g)*(R/ri)^cyl*dr_dt*grad_i/(X0-R);
                %trans_Y(I) = conv +  ( D / beta_d ) * diff;
                trans_r(I) = ri;
                frat(I) = abs( conv*beta_d/(D*diff+1.e-25) );
                trans_Y(I) = conv*beta_d/(D*diff+1.e-25) + 1;
            end
            ctrans_Y = cumtrapz( trans_r,abs(trans_Y) );
            ctrans_r = cumtrapz( trans_r,ones( 1,length(trans_r) ) );
            o_trans(jj) = ctrans_Y(end) / ctrans_r(end);
            o_transm(jj) = max(abs(trans_Y));
            
            cfrat_Y = cumtrapz( trans_r,frat );
            o_frat(jj) = cfrat_Y(end)/ctrans_r(end);
            o_fratm(jj) = max(frat);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%plot
			%figure(1);
			D2(i) = (R*Lr/Ri)^2;
			
			subplot(2,3,1)
			hold on;
			
			h1=plot(t(i)*t_d,D2(i),'-db','linewidth',2,'markersize',3);
			xlabel('time'); ylabel('(D/D_0)^2');
			if(fluid==1)
				ylim([0 1.0]); xlim([0 1.1*twfin]);
			elseif(fluid==2)
				ylim([0 1.0]); xlim([0 1.1*twfin]);
			end
			title(sprintf('time %.2e',t(i)*t_d));
			set(gca,'FontSize',12);
			
			r = zeros(gN,1);
			for k=1:gN
				r(k) = (zeta(k)-1)*(X0-R) + R;
			end
			
			%figure(2);
			subplot(2,3,2)
			plot(r*Lr*1.e6,gY*Yr,'-db','linewidth',2,'markersize',4);
			xlabel('r (\mu m)'); ylabel('Y (vapor mass fraction)');
			title( sprintf('time %.2e',t(i)*t_d));
			set(gca,'FontSize',12);
			if(fluid==1)
				xlim([0 600]); ylim([0 0.3]);
			elseif(fluid==2)
				xlim([0 10]); ylim([0 0.25]*cno(1)+[0 0.65]*cno(2));
            end	
            
			T2 = zeros(lN+gN-1,1);
			R2 = zeros(lN+gN-1,1);
			
			T2(1:lN) = Tl(1:lN);
			for k=1:lN
				R2(k) = zeta_l(k) * R; 
			end
			R2(lN+1:lN+gN-1) = r(2:gN);
			T2(lN+1:lN+gN-1) = Tg(2:gN);
			
			%figure(3);
            
			subplot(2,3,3)
			plot(R2*Lr*1.e6,T2*Tr,'-dm','linewidth',2,'markersize',4);
			xlabel('r (\mu m)'); ylabel('T (temperature)');
			title(sprintf('time %.2e',t(i)*t_d));
			set(gca,'FontSize',12);
			if(fluid==1)
				xlim([0 600]); ylim([300 600]*cno(1)+[310 380]*cno(2));
			elseif(fluid==2)
				xlim([0 10]); ylim([360 425]*cno(1)+[360 580]*cno(2));
			end	
			
			%figure(6);
			subplot(2,3,6)
			hold on;
			plot(t(i)*t_d,intT(i)*Tr,'-dk','linewidth',2,'markersize',3);
			xlabel('time (s)'); ylabel('Interface Temperature (K)');
			title(sprintf('time %.2e',t(i)*t_d));
			set(gca,'FontSize',12);
			if(fluid==1)
				xlim([0 1.1*twfin]); ylim([336 350]*cno(1)+[310 350]*cno(2));
			elseif(fluid==2)
				xlim([0 1.1*twfin]); ylim([365 410]*cno(1)+[375 480]*cno(2));
			end				
			
			% save the variable to be displayed later
			o_time(jj) = t(i)*t_d;
			o_D2(jj) = D2(i);
			o_intT(jj) = intT(i)*Tr;
			
			o_r(jj,:) = r(:)*Lr;
			o_gY(jj,:) = gY(:)*Yr;
			
			o_R2(jj,:) = R2(:)*Lr;
			o_T2(jj,:) = T2(:)*Tr;
			
			jj = jj + 1;
			
			pause(0.0005);
        end
	    
		%the final answer before remeshing
		L = length(t);
		gY(2:gN-1) = X(L,1   :gN-2); 
		R          = X(L,gN-1     );
		Tl(1:lN-1) = X(L,gN  :lb  );
		Tg(2:gN-1) = X(L,lb+1:lc  );
		
		Tint = fzero(@func_Tint,Tl(lN-1));
		Tl(lN) = Tint;
		Tg(1 ) = Tint;
		gY(1 ) = Yint(Tint*Tr)/Yr;
		
		if(j < length(tf) - 1)
		
			%remeshing
			if(uniform == 0)
				rho_mesh = mesh_density_fnct(2,gN,R2(lN:lN+gN-1),[gY Tg]',3); %[gY Tg]',3);
                if( fluid==3 )
                    rho_mesh = mesh_density_fnct(1,gN,R2(lN:lN+gN-1),gY',3);
                end
				r1 = meshgen_deboor(gN,R2(lN:lN+gN-1),rho_mesh); %dlmwrite('mesh_read_81_long',zeta1,'precision','%.14e');
				gY1 = pchip(R2(lN:lN+gN-1),gY,r1);
				Tg1 = pchip(R2(lN:lN+gN-1),Tg,r1);
				
				rho_mesh = mesh_density_fnct(1,lN,R2(1:lN),Tl',3);
				r2 = meshgen_deboor(lN,R2(1:lN),rho_mesh);
				Tl2 = pchip(R2(1:lN),Tl,r2);
			else
				r1 = R2(lN:lN+gN-1);
				gY1 = gY;
				Tg1 = Tg;
				
				r2 = R2(1:lN);
				Tl2 = Tl;
			end
			
			if(uniform == 0)
				%figure(4);
				subplot(2,3,4)
				plot(r1*Lr*1.e6,gY1*Yr,'-db','linewidth',2,'markersize',4); title(sprintf('time %.2e vap frac just after remeshing',t(end)*t_d));
				xlabel('r1 (\mu m)'); ylabel('Y (vapor fraction)');
				set(gca,'FontSize',12);
				
				%figure(5);
				subplot(2,3,5)
				plot(r2*Lr*1.e6,Tl2*Tr,'-dg','linewidth',2,'markersize',4); title(sprintf('time %.2e liq temp just after remeshing',t(end)*t_d));
				xlabel('r2 (\mu m)'); ylabel('Tl (K)');
				set(gca,'FontSize',12);
				pause(0.1);
			end
		
			for k=1:gN
				zeta(k) = 1.0 + (r1(k) - R) / (X0 - R);
			end
			gY = gY1;
			Tg = Tg1; 
			
			for k=1:lN
				zeta_l(k) = r2(k) / R;
			end
			Tl = Tl2;
			%dlmwrite('mesh_read_11_31_n_heptane',[zeta_l';zeta'],'precision','%.14e');
			
			Xi(1   :gN-2) = gY(2:gN-1);
			Xi(gN-1     ) = R;
			Xi(gN  :lb  ) = Tl(1:lN-1);
			Xi(lb+1:lc  ) = Tg(2:gN-1);
		end
	end
	
	%figure(2);
	subplot(2,3,2)
	hold on;
	theo_Y = zeros(gN,1);
	for k=1:gN
		r(k) = (zeta(k)-1)*(X0-R) + R;
	end
	
	for k=1:gN
		theo_Y(k) = theor_vap_frac(r(k),R,Tint);
	end
	
	plot(r*Lr*1.e6,theo_Y,'-r','linewidth',2);
	legend('Numerical 1D transient', 'Analytical (finite domain)');
	set(gca,'FontSize',12);
	legend boxoff
    
    if( cyl==1 )
        figure(2);
        plot(r*Lr*1.e6,gY*Yr,'-r','linewidth',2);
        xlabel('r (\mu m)'); ylabel('Y (vapor mass fraction)');
        hold on;
        plot(r*Lr*1.e6,theo_Y,'ok','linewidth',2);
        legend('Transient VOF', 'Steady State Analytical');
        set(gca,'FontSize',20,'fontweight','bold');
        legend boxoff
    end 
    
% finite domain steady-state diameter
    if( cyl==2 )
        lambda  = 1 / X0;
        Ys      = Yint( Tint*Tr );
        coeff_K = ( 8.0 * rho_g * D / ( rho_l * 4.0 * Ri^2 ) ) * log( 1.0/(1.0-Ys) );
        dia2 = linspace(0.1,1,21);
        t2 = dia2 - 1 - (2.0/3.0)*lambda*(dia2.^1.5-1);
        t2 = t2 * (-1/coeff_K);
        
        t3 = dia2 - 1;
        t3 = t3 * (-1/coeff_K);
        
        subplot(2,3,1); hold on;
        h2=plot(t2,dia2,'-r','linewidth',2);
        h3=plot(t3,dia2,'-g','linewidth',2);
        %ylim([0 1.0]); xlim([0 2.0*twfin]);
        legend([h1,h2,h3],{'Numerical 1D transient', 'Analytical (finite domain)', 'Analytical (infinite domain)'});
        set(gca,'FontSize',12);
        legend boxoff
        
    end

	%save('result_variable','o_time','o_D2','o_intT','o_r','o_gY','o_R2','o_T2');
	%save heptane2.mat o_time o_D2 o_intT o_r o_gY o_R2 o_T2
%     save fluid3a.mat o_time o_D2 o_r o_gY o_intT
%     save transient2D.mat o_time o_trans o_transm Yr t_d o_frat o_fratm
end

function y=u_g(r,R,dr_dt)
	global rho_l rho_g cyl
	y = -( rho_l/rho_g-1.0 ) * ( R/r )^cyl * ( dr_dt ) ; %taken from Lee and Son 2015
end

function y=u_l(r,R,dr_dt)
	y = 0.0 ;
end

function y = Yint(Ts)
	global h_lg
	global fluid Pamb Mf Ma isothermal
	X_A = 0;
	if( fluid == 1 ) 
		X_A = XSteam('psat_T',Ts-273.15); %from steam table
	elseif (fluid ==2 )
		if( Ts > 295.0 )
			A = 4.02832; 
			B = 1268.636;
			C = -56.199;
		else
			A = 4.81803; 
			B = 1635.409;
			C = -27.338;
		end
		X_A = 10.0 ^ ( A - B/(Ts+C) );
	end
	X_A = X_A / Pamb; %divided by atmospheric pressure in bar
	%X_A = exp(-h_lg * (18.0 / 8314.0) * (1.0/Ts - 1.0/373.0));
	M_mix = Mf*(X_A) + Ma*(1.0-X_A);
	y = X_A * Mf / M_mix; 
    if( fluid==3 ) 
        if(isothermal==1) 
            y = 0.667;
        else
           d = 500/log(0.99/0.42);
           y = 0.42 * exp((Ts-200)/d);
        end
    end 
end

function y = dR_dt(Y,R)
	global rho_l rho_g D 
	global zeta X0
	global beta_d Yr
	
	x1 = [zeta(1) zeta(2) zeta(3)];
	y1 = [Y(1) Y(2) Y(3)];
	grad_Y = ( dydx_for(x1,y1) ) / (X0 - R);
	%grad_Y = (y1(2) - y1(1)) / ( (x1(2) - x1(1)) * (X0 - R)); %1st order
	y = ( rho_g * D / ( rho_l * beta_d ) ) * grad_Y / ( 1.0/Yr - Y(1) ); 
end

function xdot = func_ode(t,x)
	% x is of length 2*gN+lN-4
	% first gN-2 entries are the vapor fraction in the gas phase
	% (gN-1)th is the radius of the drop
	% (lN-1) entries of liquid region temperature
	% (gN-2) entries of gas temperature
	
	global lN gN lc lb zeta zeta_l X0           %mesh related
	global D k_l k_g rho_l rho_g Cp_l Cp_g h_lg %fluid properties
	global R Tl Tg gY                           %local variables
	global Tr Yr beta_d
    global cyl zopt
    
	% allocation 
	xdot = zeros(2*gN+lN-4,1); 
	
	% transfer values from x to local variables
	gY(2:gN-1) = x(1   :gN-2);
	R          = x(gN-1     );  
	Tl(1:lN-1) = x(gN  :lb  );
	Tg(2:gN-1) = x(lb+1:lc  );
	 
	% let us compute the interface temperature
	Tint = fzero( @func_Tint,Tl(lN-1) ); % find by doing fixed point iteration
	Tl(lN) = Tint; 
	Tg(1 ) = Tint;
	gY(1 ) = Yint( Tint*Tr ) / Yr ;
	dr_dt = dR_dt(gY,R);
		
	for i=2:gN-1 % i is mesh index
		I=i-1; % the index for unknowns (vector x)
		ri = (zeta(i)-1)*(X0-R) + R;
		ug = u_g(ri,R,dr_dt);
		
		if(i==2)
			x1 = [zeta(i-1), zeta(i), zeta(i+1), zeta(i+2)];
			y1 = [gY(i-1), gY(i), gY(i+1), gY(i+2)];
			grad2_i = dy2dx2_for(x1,y1);
		else
			x1 = [zeta(i-2), zeta(i-1), zeta(i), zeta(i+1)];
			y1 = [gY(i-2), gY(i-1), gY(i), gY(i+1)];
			grad2_i = dy2dx2_back(x1,y1);
		end
		x1 = [zeta(i-1), zeta(i), zeta(i+1)];
		y1 = [gY(i-1), gY(i), gY(i+1)];
		grad_i = dydx(x1,y1);
		diff = (grad2_i/(X0-R)^2) + cyl*grad_i/(ri*(X0-R));
		conv = -( ug - ( 2.0-zeta(i) ) * dr_dt ) * grad_i / (X0-R);
		rhs = conv +  ( D / beta_d ) * diff;
		
		xdot(I) = rhs;
	end
	xdot(gN-1) = dr_dt;
	h = zeta_l(2) - zeta_l(1);
	
	for i=1:lN-1
		I = gN + i - 1;
		if(i==1)
			grad_i = 0.0;
			x1 = [-h, zeta_l(1), zeta_l(2), zeta_l(3)];
			y1 = [Tl(2), Tl(1), Tl(2), Tl(3)];
			grad2_i = dy2dx2_for(x1,y1);
		elseif(i==2)
			x1 = [zeta_l(i-1), zeta_l(i), zeta_l(i+1)];
			y1 = [Tl(i-1), Tl(i), Tl(i+1)];
			grad_i = dydx(x1,y1);
			x1 = [-h, zeta_l(i-1), zeta_l(i), zeta_l(i+1)];
			y1 = [Tl(i), Tl(i-1), Tl(i), Tl(i+1)];
			grad2_i = dy2dx2_back(x1,y1);
		else
			x1 = [zeta_l(i-1), zeta_l(i), zeta_l(i+1)];
			y1 = [Tl(i-1), Tl(i), Tl(i+1)];
			grad_i = dydx(x1,y1);
			x1 = [zeta_l(i-2), zeta_l(i-1), zeta_l(i), zeta_l(i+1)];
			y1 = [Tl(i-2), Tl(i-1), Tl(i), Tl(i+1)];
			grad2_i = dy2dx2_back(x1,y1);
		end
		rhs = ( k_l / (rho_l*Cp_l*beta_d) ) * ( (1.0/R^2)*grad2_i + cyl*grad_i/( (zeta_l(i)+h*1.e-20) * R^2 ) ) ...
			  - ( -zeta_l(i) * dr_dt * grad_i/(R) ) ;
		xdot(I) = rhs;
	end
	
	for i=2:gN-1
		I = lb + i -1;
		if(i==2)
			x1 = [zeta(i-1), zeta(i), zeta(i+1), zeta(i+2)];
			y1 = [Tg(i-1), Tg(i), Tg(i+1), Tg(i+2)];
			grad2_i = dy2dx2_for(x1,y1);
		else
			x1 = [zeta(i-2), zeta(i-1), zeta(i), zeta(i+1)];
			y1 = [Tg(i-2), Tg(i-1), Tg(i), Tg(i+1)];
			grad2_i = dy2dx2_back(x1,y1);
		end
		x1 = [zeta(i-1), zeta(i), zeta(i+1)];
		y1 = [Tg(i-1), Tg(i), Tg(i+1)];
		grad_i = dydx(x1,y1);
		
		ri = (zeta(i)-1)*(X0-R) + R;
		ug = u_g(ri,R,dr_dt);
		
		diff = ( k_g / ( rho_g*Cp_g*beta_d) ) * ( grad2_i / (X0 - R)^2 + cyl*grad_i / ( ri*(X0-R) ) );
		conv = - ( ug - ( 2.0-zeta(i) )*dr_dt ) * grad_i / ( X0-R );
		rhs = diff + conv;
		xdot(I) = rhs;
	end
end

%%%%%%%% second order accurate derivatives on non-uniform grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derived from Taylor's series expansion 
function der = dy2dx2_back(x,y)
	h = x(4) - x(3);
	alpha = (x(3) - x(2)) / h;
	beta  = (x(2) - x(1)) / h;
	
	M = [-alpha,   -(alpha+beta); ...
		 -alpha^3, -(alpha+beta)^3];
	r = [-1;-1];
	X = M \ r;
	A = X(1); B = X(2);
	der = (2/h^2)*(1/(1+A*alpha^2+B*(alpha+beta)^2))*...
	      (y(4)+A*y(2)+B*y(1) - (1+A+B)*y(3));
end

function der = dy2dx2_for(x,y)
	h = x(2) - x(1);
	alpha = (x(3) - x(2)) / h;
	beta  = (x(4) - x(3)) / h;
	
	M = [alpha,   (alpha+beta); ...
		 alpha^3, (alpha+beta)^3];
	r = [1;1];
	X = M \ r;
	A = X(1); B = X(2);
	der = (2/h^2)*(1/(1+A*alpha^2+B*(alpha+beta)^2))*...
	      (y(1)+A*y(3)+B*y(4) - (1+A+B)*y(2));
end

function der = dydx(x,y)
	h = x(3) - x(2);
	alpha = (x(2) - x(1)) / h;
	A = -(1/alpha^2);
	
	der = (A*y(1) + y(3) - (1+A)*y(2)) / (h*(1-A*alpha));
end

function der = dydx_back(x,y)
	h = x(3) - x(2);
	alpha = (x(2) - x(1)) / h;
	A = -1/(1+alpha)^2;
	
	der = -(y(2)+A*y(1)-(1+A)*y(3)) / (h*(1.0+A*(1.0+alpha)));
end

function der = dydx_for(x,y)
	h = x(2) - x(1);
	alpha = (x(3) - x(2))/h;
	A = -1/(1+alpha)^2;
	
	der = (y(2)+A*y(3)-(1+A)*y(1)) / (h*(1.0+A*(1.0+alpha)));
end
%%%%%%%%%%%%%%%% derivatives' definitions end here %%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = theor_vap_frac(r,R,Tint)
	global rho_g D X0
	global Tr Lr 
    global cyl
	
	Ys = Yint(Tint*Tr);
	Y_inf = 0.0;
	R_eff = Lr * R * X0 / (X0 - R); %R_eff = R;
	m_dot = (4.0*pi*R_eff*rho_g*D) * log( (1-Y_inf) / (1-Ys) );
	y = 1 - (1-Ys) * exp(-m_dot/ (4*pi*rho_g*D*r*Lr) ) / ( exp(-m_dot/ ( 4*pi*rho_g*D*R*Lr ) ) );
    if( cyl==1 )
        m_dot = 2.0*pi*rho_g*D*log( (1-Y_inf)/(1-Ys) ) / log(X0/R);
        y = 1 - (1-Ys) * (r/R)^( m_dot / (2.0*pi*rho_g*D) );
    end 
end

function y = func_Tint(Tint)
	global lN zeta zeta_l X0   % mesh related
	global k_l k_g rho_l h_lg  % fluid properties
	global R Tl Tg gY          % local variables
	global Tr Yr beta_d fluid isothermal
	
	xl = [zeta_l(lN-2) ,zeta_l(lN-1), zeta_l(lN)];
	yl = [Tl(lN-2) , Tl(lN-1), Tint];
	xg = [zeta(1) , zeta(2), zeta(3)];
	yg = [Tint  , Tg(2), Tg(3)];
	
	gY(1)  = Yint( Tint*Tr ) / Yr ;
	% grad_l = (yl(3) - yl(2)) / ( (xl(3) - xl(2)) * R );
	% grad_g = (yg(2) - yg(1)) / ( (xg(2) - xg(1)) *(X0-R));
	grad_l = ( dydx_back(xl,yl) ) / R;
	grad_g = ( dydx_for (xg,yg) ) / (X0 - R);
	
	% per unit area m_dot
	m_dot = -( rho_l*beta_d ) * dR_dt(gY,R); % mass per unit area per second
	
	y = k_l * Tr * grad_l - k_g * Tr * grad_g + m_dot * h_lg;
	
    if( fluid==3 && isothermal==1 ) 
        y = 0.0;
    end 
	% fprintf("Tint y grad_l grad_g m_dot Yint= (% e % e % e % e % e % e)\n",Tint,y,k_l*grad_l, -k_g*grad_g, m_dot*h_lg, Yint(Tint)); 
	
end

function xdot = func_ode2(x)
	% wrapper for Jacobian computation
	t = 0;
	xdot = func_ode(t,x);
end

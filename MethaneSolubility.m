function NetCH4Solubility = MethaneSolubility(P,T,Salinity,SalinityUnits,Datatype)
% This script controls a set of functions for calculating the solubility of methane in water as a function of pressure, temperature and salinity.
% Pressure, P, should be a column vector in MPa.
% Temperature, T, should be a column vector in Kelvin.
% Salinity can be either a column or row vector, in either molality (moles per kilogram pure water) or parts per thousand (ppt) units.
% SalinityUnits should be written as either 'mol/kg' or 'ppt' (include the single quotes).
% The output, NetCH4Solubility, is a matrix of solubilities (in moles of methane per kilogram of water) of size [(Length(P), Length(T), Length(Salinity)]
% To plot a surface plot at a set salinity, input surf(T,P,NetCH4Solubility(:,:,x)), where x is the index of the salinity of choice.
% Datatype should be written as either 'range' or 'points,' and distinguishes between two cases:
    %'range' describes a case in which a solubilty result is required for every possible P, T, S combination over the entire range of data given.
    %'points' describes a case in which a solubility result is required for a dataset of P,T,S points, meaning solubility is calculated only for P(1), T(1), S(1), then P(2), T(2), S(2), then P(i), T(i), S(i)...

% In the presence of hydrate, the solubility is calculated according to the equations in: 
% Tishchenko, P., C. Hensen, K. Wallmann, and C. S. Wong (2005), Calculation of the stability and solubility of methane hydrate in seawater, Chemical Geology, 219, 37-52
% At pressures and temepratures for which hydrate does not form, the methane solubilities are calculated based on two works by Duan et al:
% Duan, Z. H., N. Moller, J. Greenberg, and J. H. Weare (1992), The Prediction of Methane Solubility in Natural-Waters to High Ionic-Strength from O°C to 250°C and from 0 to 1600 Bar, Geochimica et Cosmochimica Acta, 56(4), 1451-1460
% Duan, Z. H., and S. D. Mao (2006), A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar, Geochimica et Cosmochimica Acta, 70(13), 3369-3386

% This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior.
% For more information, see the official USGS copyright policy at http://www.usgs.gov/visual-id/credit_usgs.html#copyright
% The functions XSteam and Progress bar used in this script are open source code available from Mathworks (links given at the function calls below).
% Though this code has been tested, no guarantee of its accuracy can be provided, and the results should be tested against known values at the conditions of the user's interest.

% Demonstration:
% To calculate solubilities across a range of pressures, temperatures and salinities:
    %P = [5:.5:15]';
    %T = [273:.5:293]';
    %S = [0:10:30];
    %SalinityUnits = 'ppt';
    %Datatype = 'range';
    %NetCH4Solubility = MethaneSolubility(P,T,S,SalinityUnits,Datatype);
    %surf(T,P,NetCH4Solubility(:,:,1))
    %xlabel('Temperature (K)')
    %ylabel('Pressure (MPa)')
    %zlabel('Methane Solubility (moles CH4 per kg water)')

% To calculate solubilities for a collection of linked pressure, temperature, salinity datapoints:
    %P = [5:.5:15]';
    %T = [273:1:293]';
    %S = [0:1:20]';
    %SalinityUnits = 'ppt';
    %Datatype = 'points';
    %NetCH4Solubility = MethaneSolubility(P,T,S,SalinityUnits,Datatype);
    %plot(T,NetCH4Solubility,'bo')
    %xlabel('Temperature (K)')
    %ylabel('Methane Solubility (moles CH4 per kg water)')
% End demonstration

%Input Parameters


[Smolality, Sppt] = SalinityConverter(Salinity,SalinityUnits);      % Provides the program with molality and ppt units for salinity, regardless of which unit the user has chosen

    
    
    if strcmp(Datatype,'range')
        progressbar('Salinity', 'Pressure')
        %Preallocate array sizes:
            CH4Solubility_hydrate = zeros(length(P),length(T),length(Salinity));
            CH4Solubility_nohydrate = zeros(length(P),length(T),length(Salinity));
        
        for j = 1:length(Salinity)
            

            CH4Solubility_hydrate(:,:,j) = Tishchenko_2005(P,T,Sppt(j),Datatype);
            % This function calculates the solubility of methane in water in the presence of hydrate using the formulation from Tishchenko_2005 (in moles of methane per kilogram of water)
                % P = pressure in MPa
                % T = temperature in K
                % Sppt = salinitity in parts per thousand
        
            for k = 1: length(P)
                CH4Solubility_nohydrate(k,:,j) = Duan_1992_2006(P(k)*10,T,Smolality(j));
                progressbar([],k/length(P))

                % The Duan approach requires pressure in bar, where 10 bar = 1MPa
                % T is in Kelvin
                % Smolality = salinity in moles of NaCl per kilogram of pure water
            end
            progressbar(j/length(Salinity), [])
        end
    else
        progressbar('Datapoints')
        %Preallocate array sizes:
            CH4Solubility_hydrate = zeros(length(P));
            CH4Solubility_nohydrate_temp = zeros(length(P));
        
        CH4Solubility_hydrate = Tishchenko_2005(P,T,Sppt,Datatype);%CH4Solubility_hydrate_temp(k) = Tishchenko_2005(P(k),T(k),Sppt(k),Datatype);
            % This function calculates the solubility of methane in water in the presence of hydrate using the formulation from Tishchenko_2005 (in moles of methane per kilogram of water)
                % P = pressure in MPa
                % T = temperature in K
                % Sppt = salinitity in parts per thousand
        for k = 1: length(P)
 
            CH4Solubility_nohydrate_temp(k) = Duan_1992_2006(P(k)*10,T(k),Smolality(k));
            progressbar(k/length(P))

            % The Duan approach requires pressure in bar, where 10 bar = 1MPa
            % T is in Kelvin
            % Smolality = salinity in moles of NaCl per kilogram of pure water
        end        
        CH4Solubility_nohydrate = CH4Solubility_nohydrate_temp(:,1);
    end

NetCH4Solubility = min(CH4Solubility_hydrate, CH4Solubility_nohydrate);        % At any given P, T and Salinity point, the methane solubility is given by the smaller of the two intersecting solubility curves.  The peak solubility occurs at the equilibrium conditions for methane hydrate.
    

% ----------------------------------------
% This section contains the support functions:
        % function [Smolality, Sppt] = SalinityConverter(S,units)
            % Converts between molality (Smolality, in moles of salt per kilogram water) and salinity (Sppt, in grams of salt per kg of solution (water + salt)) assuming only NaCl is present in the water.
        % function [CH4SolAtPress] = Tishchenko_2005(P,T,S)
            % This function uses the the equations in Tishchenko_2005 to calculate the solubility of methane in the presence of methane hydrate.
        % function CH4solubility_nohydrate = Duan_1992_2006(P,T,mNa)
            % This function uses the the equations in Duan_1992 and 2006 to calculate the solubility of methane in water in the absence of gas hydrate.
        % function PhiCH4=fugacitycoefficientCH4(P,T)
            % This uses the equation of state described in the appendix of Duan_1992 to obtain the fugacity coefficient needed for the solubility calculation.
        % function [muCH4RT,lambdaCH4Na,chiCH4NaCl] = interactionparameters_2006(P,T)
            % This function uses the parameters and equations in Duan_2006 to calculate interaction parameters needed for the solubility calculation.
        % function p4_T2 = XSteam_p4_T(T)
            % This is a steam table function needed for obtaining the partial pressure of water.  XSteam is a function by Magnus Holmgren, based on the standard IAPWS IF-97.  Only a portion of the full code is needed here.  XSteam can be downloaded in its entirity from Mathworks: http://www.mathworks.com/matlabcentral/fileexchange/9817
   

% ------

function [Smolality, Sppt] = SalinityConverter(S,units)
% Converts between molality (Smolality, in moles of salt per kilogram water) and salinity (Sppt, in grams of salt per kg of solution (water + salt)) assuming only NaCl is present in the water.
    Msalt = 58.4428;           % Grams of salt (NaCl) per mole
    if strcmp(units,'mol/kg')            %if units == 'mol/kg'
        Smolality = S;
        Sppt = (S.*Msalt)./(1 + (S.*Msalt/1000));   % Molality to Salinity (grams of salt per kilogram of salty water)
    else
        Smolality = S./(Msalt.*(1-S./1000));             % Salinity to Molality (moles of salt per kg water)
        Sppt = S;
    end
end

% ------

function [CH4SolAtPress] = Tishchenko_2005(P,T,S,Datatype)
    % This function uses the the equations in:
    % Tishchenko, P., C. Hensen, K. Wallmann, and C. S. Wong (2005), Calculation of the stability and solubility of methane hydrate in seawater, Chemical Geology, 219(1-4), 37-52.
    % The purpose of this function is to calculate the following three
    % parameters:
    % Pdiss = Methane hydrate dissociation pressure (MPa)
    % Chyd = Methane solubility in water in the presence of hydrate (moles of methane per kilogram of water)
    % Cnohyd = Methane solubility in water with no hydrate (moles of methane per kilogram of water)
    % The inputs are:
    % P = Pressure (MPa)
    % T = Temperature (K)
    % S = Salinity (Where the salinity of seawater is ~35 ppt)

    % Preallocate array size:
        %CH4SolAtPress = zeros(length(P),length(T));

    % Calculate the methane hydrate dissociation pressure, Pdiss (MPa), using equation 24

    lnPdiss = -1.6444866e3 - 0.1374178.*T + (5.4979866e4)./T + (2.64118188e2).*log(T) + ...
        S.*(1.1178266e4 + 7.67420344.*T - (4.515213e-3).*(T.^2) - (2.04872879e5)./T - ...
        (2.17246046e3).*(log(T))) + (S.^2).*(1.70484431e2 + 0.118594073.*T - ...
        (7.0581304e-5).*(T.^2) - (3.09796169e3)./T - 33.2031996.*log(T));
    Pdiss = exp(lnPdiss);

    % Calculate the methane solubility in the presence of methane hydrate, Chyd (mol/kg), using equation 25
    lnChyd = -2.5640213e5 - (1.6448053e2).*T + (9.1089042e-2).*(T.^2) + (4.90352929e6)./T + ...
        (4.93009113e4).*log(T) + S.*(-5.16285134e2 - 0.33622376.*T + (1.88199047e-4).*(T.^2) + ...
        (9.76525718e3)./T + 9.9523354e1.*log(T));
    
    % Calculate the methane solubility in the presence of methane hydrate at arbitrary pressure, Chydpress (mol/kg), using equation 28
    if strcmp(Datatype,'range')
        % Preallocate array size:
            CH4SolAtPress = zeros(length(P),length(T));
        for i = 1: length(P)
        lnChydpress = lnChyd + ((5.04597e-2) + (7.64415e-4).*S - ((3.90236e-4) + (5.48947e-6).*S).*T + ...
            ((7.06154e-7) + (9.87742e-9).*S).*(T.^2)).*(P(i) - Pdiss) + ((7.57285e-5) - (1.90867e-8).*S - ...
            (1.4483e-10).*(S.^2) - ((1.96207e-7) - (6.67456e-11).*S).*T).*((P(i)-Pdiss).^2);
        CH4SolAtPress(i,:) = exp(lnChydpress);
        end
    else
        lnChydpress = lnChyd + ((5.04597e-2) + (7.64415e-4).*S - ((3.90236e-4) + (5.48947e-6).*S).*T + ...
            ((7.06154e-7) + (9.87742e-9).*S).*(T.^2)).*(P - Pdiss) + ((7.57285e-5) - (1.90867e-8).*S - ...
            (1.4483e-10).*(S.^2) - ((1.96207e-7) - (6.67456e-11).*S).*T).*((P - Pdiss).^2);
        CH4SolAtPress = exp(lnChydpress);
    end
        
end

% ------

function CH4solubility_nohydrate = Duan_1992_2006(P,T,mNa)
    % This function calculates the solubility of methane in water in the absence of gas hydrate.
    % The equations are given in:
    % Duan, Z. H., N. Moller, J. Greenberg, and J. H. Weare (1992), The Prediction of Methane Solubility in Natural-Waters to High Ionic-Strength from O°C to 250°C and from 0 to 1600 Bar, Geochimica et Cosmochimica Acta, 56(4), 1451-1460
    % Duan, Z. H., and S. D. Mao (2006), A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar, Geochimica et Cosmochimica Acta, 70(13), 3369-3386
    % P = Pressure (bar)
    % T = Temperature (K)
    % mNa = molality of salt (moles of salt per kilogram of water, mol/kg)
    
    % The methane solubility is given by equation 10:
    % ln(mCH4) - ln(xCH4·phiCH4·P) - muCH4/RT - 2·lamdaCH4,Na·(mNa + mK + 2·mCa + 2·mMg) - 0.06·mSO4 + 0.00624·mNa·mCl
    % Each "m" is the concentration in molality (moles of solute/kg water)
        % We only consider pure or salt water here, so mK, mCa, mMg and mSO4 are all zero.
    % xCH4 is the mole fraction of methane in the gas phase, which will be calculated below according to equation 4 of Duan_2006.
    % phiCH4 is the fugacity coefficient for methane, which is calculated according to equations 10 and 11 in Duan_1992.
    % muCH4 is the chemical potential of methane (at a reference state). The ratio or muCH4/RT is calculated below according to equation 7 in Duan_1992.
    % lambdaCH4,Na is an "interaction parameter," also to be calculated using equation 7 of Duan_1992.
    
    % Calculate xCH4, the mole fraction of methane in the gas phase (Equation 4 in Duan_1992).
    PartialpressureH2O = zeros(length(T),1);      % Preallocate array size:
    for i = 1:length(T)
        if T(i) > 273.15
            PartialpressureH2O(i) = XSteam_p4_T(T(i));

        else
            PartialpressureH2O(i) = .0061;      % 0.0061 is the result for XSteam at .01°C.  For such small partial pressures, particularly given the methane pressures generally in use when considering methane hydrate, this value could simply be replaced with zero and ignored.
        end
    end
    xCH4 = (P - PartialpressureH2O)./P;                 % This is equation 4 in Duan_1992.  XSteam provides pressure in bar, which is what Duan also uses, so no conversion is necessary
    
    % Calculate the fugacity coefficient phiCH4
    phiCH4 = fugacitycoefficientCH4(P,T);                   % PhiCH4 is the fugacity of methane. P in bar, T in Kelvin
    
    % Calculate the interaction parameters needed for the solubility equation
    [muCH4RT,lambdaCH4Na,chiCH4NaCl] = interactionparameters_2006(P,T);         % Equation 7, with parameters from Table 2, in Duan, Z. H., and S. D. Mao (2006), A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar, Geochimica et Cosmochimica Acta, 70(13), 3369-3386.
    % Calculate the methane solubility with equation 10 in Duan
        %These extra molalities are used when handling more complex mixtures.  For water containing only NaCl, mNa = mCL and the other molalities are not needed.
            mK = 0;
            mCa = 0;
            mMg = 0;
            mSO4 = 0;
            mCl = mNa;
            lambdaCH4SO4 = 0.0332;      % This is from Duan_2006
    lnmCH4 = log(xCH4.*phiCH4.*P) - muCH4RT - 2*lambdaCH4Na.*(mNa + mK + 2*mCa + 2*mMg)- chiCH4NaCl.*(mNa + mK + 2*mCa + 2*mMg).*(mCl + 2.*mSO4) - 4.*lambdaCH4SO4*mSO4;  %This is from Duan_1992
    CH4solubility_nohydrate = exp(lnmCH4);
end

% ------

function PhiCH4=fugacitycoefficientCH4(P,T)
    % This uses the equation of state described in the appendix of:
    % Duan, Z. H., N. Moller, J. Greenberg, and J. H. Weare (1992), The Prediction of Methane Solubility in Natural-Waters to High Ionic-Strength from O°C to 250°C and from 0 to 1600 Bar, Geochimica et Cosmochimica Acta, 56(4), 1451-1460
    % T must be in K, P must be in bars.


    % Fit parameters from Table A1 of Duan_1992

    a1 = 8.72553928e-2;
    a2 = -7.52599476e-1;
    a3 = 3.75419887e-1;
    a4 = 1.07291342e-2;
    a5 = 5.49626360e-3;
    a6 = -1.84772802e-2;
    a7 = 3.18993183e-4;
    a8 = 2.11079375e-4;
    a9 = 2.01682801e-5;
    a10 = -1.65606189e-5;
    a11 = 1.19614546e-4;
    a12 = -1.08087289e-4;
    alpha = 4.48262295e-2;
    beta = 7.5397e-1;
    gamma = 7.7167e-2;

    % Critical Temperature (Tc) and Pressure (Pc) from Appendix of Duan_1992

    Tc = 190.6;         % In K
    Pc = 46.41;         % In bar

    % Convenient parameter clusters for the equation of state (Used in both Duan 1992a and b)

    Pr = (P)/Pc;         % Pressure needs to be in bar for this calculation
    Tr = (T)/Tc;         % T needs to be in Kelvin for the equation of state calculation
    VrSeed = 1;         % This is a total guess, and is needed to obtain a first guess of Vr for the fit function.  Because the two terms we are attempting to match are essentialy monotonic, they have only one crossover point and the intial guess can be just about anything positive.

    B = a1 + a2./(Tr.^2) + a3./(Tr.^3);
    C = a4 + a5./(Tr.^2) + a6./(Tr.^3);
    D = a7 + a8./(Tr.^2) + a9./(Tr.^3);
    E = a10 + a11./(Tr.^2) + a12./(Tr.^3);
    F = alpha./(Tr.^3);

    % Obtain the appropriate Vr and compressibility factor, Z by equating the left and right-hand sides of equation A1 in Duan_1992
    Vr = zeros(length(T),1);
    for i = 1: length(T)
        Vrtemp = VrSeed;
        Vrtemp = fminsearch(@Vrfinder,Vrtemp);
        Vr(i) = Vrtemp;
    end
    Z = 1 + B./Vr + C./(Vr.^2) + D./(Vr.^4) + E./(Vr.^5) + (F./(Vr.^2)).*(beta + gamma./(Vr.^2)).*exp(-gamma./(Vr.^2));

    % Obtain the methane fugacity, PhiCH4 from equation A2 in Duan_1992
    lnPhiCH4 = Z - 1 - log(Z) + B./Vr + C./(2*(Vr.^2)) + D./(4*(Vr.^4)) + E./(5*(Vr.^5)) + (F./(2*gamma)).*(beta + 1 - (beta + 1 + gamma./(Vr.^2)).*exp(-gamma./(Vr.^2)));
    PhiCH4 = exp(lnPhiCH4);

    % Numerically calculate Vr
    function minerrorVr = Vrfinder(Vrtemp)    
        Z = 1 + B(i)./Vrtemp + C(i)./(Vrtemp.^2) + D(i)./(Vrtemp.^4) + E(i)./(Vrtemp.^5) + (F(i)./(Vrtemp.^2)).*(beta + gamma./(Vrtemp.^2)).*exp(-gamma./(Vrtemp.^2));
        minerrorVr = (((Pr.*Vrtemp./Tr(i)) - Z)*1000).^2;
    end
end

% ------
    
function [muCH4RT,lambdaCH4Na,chiCH4NaCl] = interactionparameters_2006(P,T)
    % This function uses the fit parameters in Table 3 of Duan_2006 to calculate interaction parameters using Equation 9 in Duan_2006
    c1 = [0.83143711e1, -0.81222036, -0.29903571e-02];
    c2 = [-0.72772168e-3, 0.10635172e-2, 0];
    c3 = [0.21489858e4, 0.18894036e3, 0];
    c4 = [-0.14019672e-4, 0, 0];
    c5 = [-0.66743449e6, 0, 0];
    c6 = [0.76985890e-2, 0.44105635e-4, 0];
    c7 = [-0.50253331e-5, 0, 0];
    c8 = [-0.30092013e1, 0, 0];
    c9 = [0.48468502e3, 0, 0];
    c10 = [0, -0.46797718e-10, 0];
        
    muCH4RT = c1(1) + c2(1)*T + c3(1)./T + c4(1)*(T.^2) + c5(1)./(T.^2) + c6(1)*P + c7(1)*P.*T + c8(1)*P./T + c9(1)*P./(T.^2) + c10(1).*(P.^2).*T;
    lambdaCH4Na = c1(2) + c2(2)*T + c3(2)./T + c4(2)*(T.^2) + c5(2)./(T.^2) + c6(2)*P + c7(2)*P.*T + c8(2)*P./T + c9(2)*P./(T.^2) + c10(2).*(P.^2).*T;
    chiCH4NaCl = c1(3);
end

% ------
function p4_T2 = XSteam_p4_T(T)
        % XSteam is a function by Magnus Holmgren, based on the standard IAPWS IF-97.  XSteam can be downloaded from Mathworks: http://www.mathworks.com/matlabcentral/fileexchange/9817       
        % This application requires only one function (p4_T) within the comprehensive XSteam approach
    %Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    %Section 8.1 The Saturation-Pressure Equation
    %Eq 30, Page 33
    teta = T - 0.23855557567849 / (T - 650.17534844798);
    a = teta ^ 2 + 1167.0521452767 * teta - 724213.16703206;
    B = -17.073846940092 * teta ^ 2 + 12020.82470247 * teta - 3232555.0322333;
    C = 14.91510861353 * teta ^ 2 - 4823.2657361591 * teta + 405113.40542057;
    p4_T2 = ((2 * C / (-B + (B ^ 2 - 4 * a * C) ^ 0.5)) ^ 4)*10;     % XSteam uses this function to take in a Kelvin temperature and provide a pressure in MPa.  The XSteam result is multiplied by 10 here to obtain the bar units required for Duan.
end

% ------
function progressbar(varargin)
    % Progressbar is a function by Steve Hoelzer, and can be downloaded from Mathworks: http://www.mathworks.com/matlabcentral/fileexchange/6922-progressbar
% Description:
%   progressbar() provides an indication of the progress of some task using
% graphics and text. Calling progressbar repeatedly will update the figure and
% automatically estimate the amount of time remaining.
%   This implementation of progressbar is intended to be extremely simple to use
% while providing a high quality user experience.
%
% Features:
%   - Can add progressbar to existing m-files with a single line of code.
%   - Supports multiple bars in one figure to show progress of nested loops.
%   - Optional labels on bars.
%   - Figure closes automatically when task is complete.
%   - Only one figure can exist so old figures don't clutter the desktop.
%   - Remaining time estimate is accurate even if the figure gets closed.
%   - Minimal execution time. Won't slow down code.
%   - Randomized color. When a programmer gets bored...
%
% Example Function Calls For Single Bar Usage:
%   progressbar               % Initialize/reset
%   progressbar(0)            % Initialize/reset
%   progressbar('Label')      % Initialize/reset and label the bar
%   progressbar(0.5)          % Update
%   progressbar(1)            % Close
%
% Example Function Calls For Multi Bar Usage:
%   progressbar(0, 0)         % Initialize/reset two bars
%   progressbar('A', '')      % Initialize/reset two bars with one label
%   progressbar('', 'B')      % Initialize/reset two bars with one label
%   progressbar('A', 'B')     % Initialize/reset two bars with two labels
%   progressbar(0.3)          % Update 1st bar
%   progressbar(0.3, [])      % Update 1st bar
%   progressbar([], 0.3)      % Update 2nd bar
%   progressbar(0.7, 0.9)     % Update both bars
%   progressbar(1)            % Close
%   progressbar(1, [])        % Close
%   progressbar(1, 0.4)       % Close
%
% Notes:
%   For best results, call progressbar with all zero (or all string) inputs
% before any processing. This sets the proper starting time reference to
% calculate time remaining.
%   Bar color is choosen randomly when the figure is created or reset. Clicking
% the bar will cause a random color change.
%
% Demos:
%     % Single bar
%     m = 500;
%     progressbar % Init single bar
%     for i = 1:m
%       pause(0.01) % Do something important
%       progressbar(i/m) % Update progress bar
%     end
% 
%     % Simple multi bar (update one bar at a time)
%     m = 4;
%     n = 3;
%     p = 100;
%     progressbar(0,0,0) % Init 3 bars
%     for i = 1:m
%         progressbar([],0) % Reset 2nd bar
%         for j = 1:n
%             progressbar([],[],0) % Reset 3rd bar
%             for k = 1:p
%                 pause(0.01) % Do something important
%                 progressbar([],[],k/p) % Update 3rd bar
%             end
%             progressbar([],j/n) % Update 2nd bar
%         end
%         progressbar(i/m) % Update 1st bar
%     end
% 
%     % Fancy multi bar (use labels and update all bars at once)
%     m = 4;
%     n = 3;
%     p = 100;
%     progressbar('Monte Carlo Trials','Simulation','Component') % Init 3 bars
%     for i = 1:m
%         for j = 1:n
%             for k = 1:p
%                 pause(0.01) % Do something important
%                 % Update all bars
%                 frac3 = k/p;
%                 frac2 = ((j-1) + frac3) / n;
%                 frac1 = ((i-1) + frac2) / m;
%                 progressbar(frac1, frac2, frac3)
%             end
%         end
%     end
%
% Author:
%   Steve Hoelzer
%
% Revisions:
% 2002-Feb-27   Created function
% 2002-Mar-19   Updated title text order
% 2002-Apr-11   Use floor instead of round for percentdone
% 2002-Jun-06   Updated for speed using patch (Thanks to waitbar.m)
% 2002-Jun-19   Choose random patch color when a new figure is created
% 2002-Jun-24   Click on bar or axes to choose new random color
% 2002-Jun-27   Calc time left, reset progress bar when fractiondone == 0
% 2002-Jun-28   Remove extraText var, add position var
% 2002-Jul-18   fractiondone input is optional
% 2002-Jul-19   Allow position to specify screen coordinates
% 2002-Jul-22   Clear vars used in color change callback routine
% 2002-Jul-29   Position input is always specified in pixels
% 2002-Sep-09   Change order of title bar text
% 2003-Jun-13   Change 'min' to 'm' because of built in function 'min'
% 2003-Sep-08   Use callback for changing color instead of string
% 2003-Sep-10   Use persistent vars for speed, modify titlebarstr
% 2003-Sep-25   Correct titlebarstr for 0% case
% 2003-Nov-25   Clear all persistent vars when percentdone = 100
% 2004-Jan-22   Cleaner reset process, don't create figure if percentdone = 100
% 2004-Jan-27   Handle incorrect position input
% 2004-Feb-16   Minimum time interval between updates
% 2004-Apr-01   Cleaner process of enforcing minimum time interval
% 2004-Oct-08   Seperate function for timeleftstr, expand to include days
% 2004-Oct-20   Efficient if-else structure for sec2timestr
% 2006-Sep-11   Width is a multiple of height (don't stretch on widescreens)
% 2010-Sep-21   Major overhaul to support multiple bars and add labels
%

persistent progfig progdata lastupdate

% Get inputs
if nargin > 0
    input = varargin;
    ninput = nargin;
else
    % If no inputs, init with a single bar
    input = {0};
    ninput = 1;
end

% If task completed, close figure and clear vars, then exit
if input{1} == 1
    if ishandle(progfig)
        delete(progfig) % Close progress bar
    end
    clear progfig progdata lastupdate % Clear persistent vars
    drawnow
    return
end

% Init reset flag 
resetflag = false;

% Set reset flag if first input is a string
if ischar(input{1})
    resetflag = true;
end

% Set reset flag if all inputs are zero
if input{1} == 0
    % If the quick check above passes, need to check all inputs
    if all([input{:}] == 0) && (length([input{:}]) == ninput)
        resetflag = true;
    end
end

% Set reset flag if more inputs than bars
if ninput > length(progdata)
    resetflag = true;
end

% If reset needed, close figure and forget old data
if resetflag
    if ishandle(progfig)
        delete(progfig) % Close progress bar
    end
    progfig = [];
    progdata = []; % Forget obsolete data
end

% Create new progress bar if needed
if ishandle(progfig)
else % This strange if-else works when progfig is empty (~ishandle() does not)
    
    % Define figure size and axes padding for the single bar case
    height = 0.03;
    width = height * 8;
    hpad = 0.02;
    vpad = 0.25;
    
    % Figure out how many bars to draw
    nbars = max(ninput, length(progdata));
    
    % Adjust figure size and axes padding for number of bars
    heightfactor = (1 - vpad) * nbars + vpad;
    height = height * heightfactor;
    vpad = vpad / heightfactor;
    
    % Initialize progress bar figure
    left = (1 - width) / 2;
    bottom = (1 - height) / 2;
    progfig = figure(...
        'Units', 'normalized',...
        'Position', [left bottom width height],...
        'NumberTitle', 'off',...
        'Resize', 'off',...
        'MenuBar', 'none' );
    
    % Initialize axes, patch, and text for each bar
    left = hpad;
    width = 1 - 2*hpad;
    vpadtotal = vpad * (nbars + 1);
    height = (1 - vpadtotal) / nbars;
    for ndx = 1:nbars
        % Create axes, patch, and text
        bottom = vpad + (vpad + height) * (nbars - ndx);
        progdata(ndx).progaxes = axes( ...
            'Position', [left bottom width height], ...
            'XLim', [0 1], ...
            'YLim', [0 1], ...
            'Box', 'on', ...
            'ytick', [], ...
            'xtick', [] );
        progdata(ndx).progpatch = patch( ...
            'XData', [0 0 0 0], ...
            'YData', [0 0 1 1] );
        progdata(ndx).progtext = text(0.99, 0.5, '', ...
            'HorizontalAlignment', 'Right', ...
            'FontUnits', 'Normalized', ...
            'FontSize', 0.7 );
        progdata(ndx).proglabel = text(0.01, 0.5, '', ...
            'HorizontalAlignment', 'Left', ...
            'FontUnits', 'Normalized', ...
            'FontSize', 0.7 );
        if ischar(input{ndx})
            set(progdata(ndx).proglabel, 'String', input{ndx})
            input{ndx} = 0;
        end
        
        % Set callbacks to change color on mouse click
        set(progdata(ndx).progaxes, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        set(progdata(ndx).progpatch, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        set(progdata(ndx).progtext, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        set(progdata(ndx).proglabel, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        
        % Pick a random color for this patch
        changecolor([], [], progdata(ndx).progpatch)
        
        % Set starting time reference
        if ~isfield(progdata(ndx), 'starttime') || isempty(progdata(ndx).starttime)
            progdata(ndx).starttime = clock;
        end
    end
    
    % Set time of last update to ensure a redraw
    lastupdate = clock - 1;
    
end

% Process inputs and update state of progdata
for ndx = 1:ninput
    if ~isempty(input{ndx})
        progdata(ndx).fractiondone = input{ndx};
        progdata(ndx).clock = clock;
    end
end

% Enforce a minimum time interval between graphics updates
myclock = clock;
if abs(myclock(6) - lastupdate(6)) < 0.01 % Could use etime() but this is faster
    return
end

% Update progress patch
for ndx = 1:length(progdata)
    set(progdata(ndx).progpatch, 'XData', ...
        [0, progdata(ndx).fractiondone, progdata(ndx).fractiondone, 0])
end

% Update progress text if there is more than one bar
if length(progdata) > 1
    for ndx = 1:length(progdata)
        set(progdata(ndx).progtext, 'String', ...
            sprintf('%1d%%', floor(100*progdata(ndx).fractiondone)))
    end
end

% Update progress figure title bar
if progdata(1).fractiondone > 0
    runtime = etime(progdata(1).clock, progdata(1).starttime);
    timeleft = runtime / progdata(1).fractiondone - runtime;
    timeleftstr = sec2timestr(timeleft);
    titlebarstr = sprintf('%2d%%    %s remaining', ...
        floor(100*progdata(1).fractiondone), timeleftstr);
else
    titlebarstr = ' 0%';
end
set(progfig, 'Name', titlebarstr)

% Force redraw to show changes
drawnow

% Record time of this update
lastupdate = clock;
end

% ------------------------------------------------------------------------------
function changecolor(h, e, progpatch) %#ok<INUSL>
% Change the color of the progress bar patch

% Prevent color from being too dark or too light
colormin = 1.5;
colormax = 2.8;

thiscolor = rand(1, 3);
while (sum(thiscolor) < colormin) || (sum(thiscolor) > colormax)
    thiscolor = rand(1, 3);
end

set(progpatch, 'FaceColor', thiscolor)
end

% ------------------------------------------------------------------------------
function timestr = sec2timestr(sec)
% Convert a time measurement from seconds into a human readable string.

% Convert seconds to other units
w = floor(sec/604800); % Weeks
sec = sec - w*604800;
d = floor(sec/86400); % Days
sec = sec - d*86400;
h = floor(sec/3600); % Hours
sec = sec - h*3600;
m = floor(sec/60); % Minutes
sec = sec - m*60;
s = floor(sec); % Seconds

% Create time string
if w > 0
    if w > 9
        timestr = sprintf('%d week', w);
    else
        timestr = sprintf('%d week, %d day', w, d);
    end
elseif d > 0
    if d > 9
        timestr = sprintf('%d day', d);
    else
        timestr = sprintf('%d day, %d hr', d, h);
    end
elseif h > 0
    if h > 9
        timestr = sprintf('%d hr', h);
    else
        timestr = sprintf('%d hr, %d min', h, m);
    end
elseif m > 0
    if m > 9
        timestr = sprintf('%d min', m);
    else
        timestr = sprintf('%d min, %d sec', m, s);
    end
else
    timestr = sprintf('%d sec', s);
end
end
end

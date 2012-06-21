function handles = Shepard(Sand,Silt,Clay)
% This script makes a ternary diagram for sand, silt and clay.  The data
% can be input as masses or as "percent of specimen."  This function
% assumes the entire specimen is represented, and renormalizes the data so
% Sand + Silt + Clay = Total Sample.

% This plotting style was introduced by Shepard in his 1954 article:
% Shepard, F.P., 1954, Nomenclature based on sand-silt-clay ratios: Journal Sedimentary Petrology, v. 24, p. 151-158.
% This script is a modification of the TERNPLOT script by Carl Sandrock (20020827)
% http://www.mathworks.com/matlabcentral/fileexchange/2299-ternplot
% This modification simply adds the labels required for the Shepard diagram.

% NOTE: The plotted result seems to be screen dependent.  If the result is
% a tiny triangle with words crammed together, simply use the plot window
% to zoom in, and everything will fall into place.


%        clay
%         / \
%        /   \
%    Sand --- Silt 

% Author: Carl Sandrock 20020827


% This software is in the public domain because it contains materials that originally came from the United States Geological Survey,
% an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at
% http://www.usgs.gov/visual-id/credit_usgs.html#copyright
% The function TERNPLOT used in this script is an open source code available from Mathworks.


% Demonstration
    % To illustrate the plotting capacity of this script for a small
    % collection of points, enter:
    % Shepard([10 20 15 45 70],[10 65 35 05 20],[80 15 50 50 10])
    % In general, column vectors containing a complete set would be defined
    % prior to calling the Shepard function.

    
% Set up the intersection points to form the interior boundaries of the
% Shepard plot
% Boundary points for the pure phase regions
[c_sac(1), c_sac(2)] = frac2xy(.25,.75);
[c_sc(1), c_sc(2)] = frac2xy(0,.75);
[sa_csa(1), sa_csa(2)] = frac2xy(.75,.25);
[sa_ssa(1), sa_ssa(2)] = frac2xy(.75,0);
[s_cs(1), s_cs(2)] = frac2xy(0,.25);
[s_sas(1), s_sas(2)] = frac2xy(.25,0);
% Boundary points for the intersections with the pure phase regions
[c(1), c(2)] = frac2xy(.125,.75);
[sa(1), sa(2)] = frac2xy(.75, .125);
[s(1), s(2)] = frac2xy(.125,.125);
% Boundary Points for homogeneous central region
[sasc_c(1), sasc_c(2)] = frac2xy(.2,.6);
[sasc_csa(1), sasc_csa(2)] = frac2xy(.4,.4);
[sasc_sa(1), sasc_sa(2)] = frac2xy(.6,.2);
[sasc_sas(1),sasc_sas(2)] = frac2xy(.4,.2);
[sasc_s(1), sasc_s(2)] = frac2xy(.2,.2);
[sasc_cs(1), sasc_cs(2)] = frac2xy(.2,.4);
% Boundary points on the midpoints of the triangle sides
[sac_csa(1),sac_csa(2)] = frac2xy(.5,.5);
[ssa_sas(1),ssa_sas(2)] = frac2xy(.5,0);
[cs_sc(1),cs_sc(2)] = frac2xy(0,.5);
% Done picking boundary points




A = Sand;
B = Silt;
C = Clay;

xoffset = 0.04;
yoffset = 0.02;

majors = 4;

Total = (A+B+C);
fA = A./Total;
fB = B./Total;
fC = 1-(fA+fB);

[x, y] = frac2xy(fA, fC);

% Sort data points in x order
[x, i] = sort(x);
y = y(i);

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');

set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state
	%plot axis lines
	hold on;
	
    plot ([0 1 0.5 0],[0 0 sin(1/3*pi) 0], 'color', tc, 'linewidth',1,...
                   'handlevisibility','off');
	set(gca, 'visible', 'off');
    

    % plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata', [0 1 0.5 0], 'ydata', [0 0 sin(1/3*pi) 0], ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end
    
	% Generate labels
	majorticks = linspace(0, 1, majors + 1);
	majorticks = majorticks(2:end-1);
	labels = num2str(majorticks'*100);
	
	zerocomp = zeros(1, length(majorticks));
	
	% Plot right labels
	[lxc, lyc] = frac2xy(zerocomp, majorticks);
	text(lxc, lyc, [repmat('  ', length(labels), 1) labels]);
	
	% Plot bottom labels
	[lxb, lyb] = frac2xy(1-majorticks, zerocomp); % fA = 1-fB
	text(lxb, lyb, labels, 'VerticalAlignment', 'Top');
	
	% Plot left labels
	[lxa, lya] = frac2xy(majorticks, 1-majorticks);
	text(lxa-xoffset, lya, labels);
	
end;

% Reset defaults
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

%Draw the separator lines on the Shephard Diagram

plot([c_sac(1) c_sc(1)],[c_sac(2) c_sc(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sa_csa(1) sa_ssa(1)],[sa_csa(2) sa_ssa(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([s_sas(1) s_cs(1)],[s_sas(2) s_cs(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sasc_c(1) c(1)],[sasc_c(2) c(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sac_csa(1) sasc_csa(1)],[sac_csa(2) sasc_csa(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sasc_sa(1) sa(1)],[sasc_sa(2) sa(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([ssa_sas(1) sasc_sas(1)],[ssa_sas(2) sasc_sas(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sasc_s(1) s(1)],[sasc_s(2) s(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([cs_sc(1) sasc_cs(1)],[cs_sc(2) sasc_cs(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sasc_c(1) sasc_sa(1)],[sasc_c(2) sasc_sa(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sasc_sa(1) sasc_s(1)],[sasc_sa(2) sasc_s(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');
plot([sasc_s(1) sasc_c(1)],[sasc_s(2) sasc_c(2)], 'color', 'k', 'linewidth',1,...
                   'handlevisibility','off');

% Write in the Sand, Silt and Clay labels

%Sand
sandlocation = [0.168581 0.481926 12.052];
sandlabel = text(sandlocation(1), sandlocation(2), 'Sand Content (%)');
set(sandlabel,'FontName','Times','FontSize',18,'Rotation',60,'HorizontalAlignment','center');
%Silt
siltlocation = [0.506579 -0.0776316 12.052];
siltlabel = text(siltlocation(1), siltlocation(2), 'Silt Content (%)');
set(siltlabel,'FontName','Times','FontSize',18,'HorizontalAlignment','center');
%Clay
claylocation = [0.832 0.481926 12.052];
claylabel = text(claylocation(1), claylocation(2), 'Clay Content (%)');
set(claylabel,'FontName','Times','FontSize',18,'Rotation',-60,'HorizontalAlignment','center');

% Write in the individual unit names
clayonlylocation = [.5 .73];
clayonly = text(clayonlylocation(1), clayonlylocation(2),'Clay');
set(clayonly,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

sandyclaylocation = [.4 .5];
sactext(1) = {'Sandy'};
sactext(2) = {'Clay'};
sandyclay = text(sandyclaylocation(1),sandyclaylocation(2),sactext);
set(sandyclay,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

siltyclaylocation = [.6 sandyclaylocation(2)];
sctext(1) = {'Silty'};
sctext(2) = {'Clay'};
siltyclay = text(siltyclaylocation(1), siltyclaylocation(2),sctext);
set(siltyclay,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

clayeysandlocation = [.25 .25];
csatext(1) = {'Clayey'};
csatext(2) = {'Sand'};
clayeysand = text(clayeysandlocation(1), clayeysandlocation(2), csatext);
set(clayeysand,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

clayeysiltlocation = [.75 clayeysandlocation(2)];
cstext(1) = {'Clayey'};
cstext(2) = {'Silt'};
clayeysilt = text(clayeysiltlocation(1), clayeysiltlocation(2), cstext);
set(clayeysilt,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

homogenouslocation = [clayonlylocation(1) .3];
homtext(1) = {'Sand'};
homtext(2) = {'+'};
homtext(3) = {'Silt'};
homtext(4) = {'+'};
homtext(5) = {'Clay'};
homogeneous = text(homogenouslocation(1), homogenouslocation(2), homtext);
set(homogeneous,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

sandonlylocation = [.12 .08];
sandonly = text(sandonlylocation(1),sandonlylocation(2),'Sand');
set(sandonly,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

siltysandlocation = [sandyclaylocation(1), sandonlylocation(2)];
ssatext(1) = {'Silty'};
ssatext(2) = {'Sand'};
siltysand = text(siltysandlocation(1), siltysandlocation(2),ssatext);
set(siltysand,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

sandysiltlocation = [siltyclaylocation(1), sandonlylocation(2)];
sastext(1) = {'Sandy'};
sastext(2) = {'Silt'};
sandysilt = text(sandysiltlocation(1), sandysiltlocation(2),sastext);
set(sandysilt,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

siltonlylocation = [.88 sandonlylocation(2)];
siltonly = text(siltonlylocation(1),siltonlylocation(2),'Silt');
set(siltonly,'FontName','Times','FontSize',10,'HorizontalAlignment','center');

% Plot data

q = plot(x, y,'bo');

if nargout > 0
    hpol = q;
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end

end
%-----------------------------------------%

% Support Functions frac2xy and radians

function [x, y] = frac2xy(fA, fC);
y = fC*sin(radians(60));
x = 1 - fA - y*cot(radians(60));
end

function radians = radians(degrees)
% RADIANS (DEGREES)
% Conversion function from Radians to Degrees.
% Richard Medlock 12-03-2002
radians = ((2*pi)/360)*degrees;
end
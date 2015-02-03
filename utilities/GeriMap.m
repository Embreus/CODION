function map = GeriMap(varargin)
%Defines a color map that is suitable for both color and b&w printing
%(courtesy of Gergely Papp). 
% Usage: 
%       colormap(MyColorMap())
%       colors = MyColorMap(nColors);
%       colors = MyColorMap(nColors,includeLight);
%       colors = MyColorMap(nColors,includeLight,includeDark);
% 
% includeLight = 0 improves the visibility on white background by excluding
% the lightest colors of the map. includeDark = 0 similarly excludes the 
% darkest colors, for instance to avoid confusion with truly black lines.
% If nColors = 2, black and red are returned rather than black and white.


  map = [0 0 0;...
        .15 .15 .5;...
        .3 .15 .75;...
        .6 .2 .50;...
        1 .25 .15;...
        .9 .5 0;...
        .9 .75 .1;...
        .9 .9 .5;...
        1 1 1];
    
if nargin >= 1 && isnumeric(varargin{1})
    if nargin >= 2
        if isscalar(varargin{2}) && isnumeric(varargin{2})
            includeLight = varargin{2};
        else
            includeLight = 1;
            warning('Invalid value passed for the includeLight parameter.');
        end        
    else
        includeLight = 1;
    end
    if nargin == 3
        if isscalar(varargin{3}) && isnumeric(varargin{3})
            includeDark = varargin{3};
        else
            includeDark = 1;
            warning('Invalid value passed for the includeDark parameter.');
        end
    else
        includeDark = 1;
    end
    nColors = varargin{1};
    
    %Determine which part of the map to use
    if nColors == 2 %If only two colors, use black and red
        cMin = 1;
        cMax = 5;
    else
        if includeLight
            cMax = 9;     %There are 9 entries in the map  
        else %Ignore the lightest part of the map for improved visibility on white background
            cMax = 7;              
        end
        if includeDark
            cMin = 1;
        else
            cMin = 3;
        end
    end
    
    colors = linspace(cMin,cMax,nColors);
    map = interp1(map,colors); %Interpolate between them
end
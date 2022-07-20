%% Settings
%wire_dia = 0.8128; %mm 20 awg
wire_dia = 0.6426; %mm 22 awg
%wire_dia = 0.5106; %mm 24 awg
insulation_thickness = 0.022+0.80;
coil_dia = 5; %in
turns =1200; %turns
secondary_length = 30; %in

%TEMP
resonance_freq = 500e3;


%% Main Function

L = getInductance(coil_dia, secondary_length, turns)
[Rac,Rdc] = getRac(wire_dia, insulation_thickness, coil_dia, turns,resonance_freq)
Q = getQ(resonance_freq, L, Rac)



function freq = getResonantFreq(inductance, capacitance)
    freq = 1/(2 * pi * sqrt(inductance * capacitance));
end

function Z = getImpedance(inductance, capacitance)
    Z = sqrt(inductance / capacitance);
end

function Q = getQ(freq, inductance, Rac)
    Q = 2*pi*freq * inductance / Rac;
end

function [Rac,Rdc]= getRac(wire_dia, insulation_thickness, coil_dia, turns, resonance_freq)
    %% Constants
    copper_conductivity = 58.7e6; %Siemens/m
    in2m = 0.0254; %m/in
    u0 = 4e-7 * pi;
    %% Resistance Math
    wire_length = coil_dia*pi*turns*in2m; %m
    wire_area = (wire_dia/2*1e-3)^2*pi; %m^2

    Rdc = wire_length / (wire_area*copper_conductivity); %Ω
    skin_depth = getSkinDepth(resonance_freq,copper_conductivity);
    %% Prox effect Math https://en.wikipedia.org/wiki/Proximity_effect_(electromagnetism) 
    center_center = (wire_dia+insulation_thickness);
    w = resonance_freq * 2 * pi;
    N = wire_dia / (center_center + wire_dia * 0.1494) ; %This is some clusterfuck aprox based on mixing dowell, Dixon theory max packing, and bartoli 
    alpha = sqrt(i * w * u0 * copper_conductivity * N);
    x = wire_dia * 1e-3 * alpha;
    Rac_scalar = real(x * coth(x)); %single layer Dowell method
    Rac = Rdc * Rac_scalar; %Ω
end

function L = getInductance(coil_dia, secondary_length, turns)
    in2m = 0.0254; %m/in
    A = pi * coil_dia^2 / 4 * in2m^2; %m^2
    u0 = 4e-7 * pi;
    L = u0 * turns^2 * A / (secondary_length * in2m); 
end

%% Support function for support functions
function skin_depth = getSkinDepth(freq,conductivity)
    u0 = 4e-7 * pi;
    skin_depth = 1/sqrt(pi * u0 * conductivity * freq);
    skin_depth = skin_depth * 1e3; %convert to mm
end



%{
function [Rac,Rdc]= getRac(wire_dia, insulation_thickness, coil_dia, turns,resonance_freq)
    %% Constants
    copper_conductivity = 58.7e6; %Siemens/m
    in2m = 0.0254; %m/in
    %% Resistance Math
    wire_length = coil_dia*pi*turns*in2m; %m
    wire_area = (wire_dia/2*1e-3)^2*pi; %m^2

    Rdc = wire_length / (wire_area*copper_conductivity); %Ω
    skin_depth = getSkinDepth(resonance_freq,copper_conductivity);
    %% Prox effect Math https://en.wikipedia.org/wiki/Proximity_effect_(electromagnetism) 
    center_center = (wire_dia+insulation_thickness);
    d = (sqrt(pi) / 2) * wire_dia / center_center;
    Rac_scalar = (d/skin_depth/2) * coth(d/skin_depth); %single layer Dowell method
    Rac = Rdc * Rac_scalar; %Ω
end
%}
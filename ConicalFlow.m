%% AE 4610 (Spring 2022) : PROJECT
% Written by Wonhee Lee
% First Written on 27th April 2022
% Need to Extend to Flow Properties

close all
clear all
clc

%% Plot Settings
% =============================================================
% Plot Font Sizes
sizeLabel = 38;
sizeLegend = 32;
sizeGCA = 32;
  
sz1 = 2;
sz2 = 45;
sz3 = 25;
clr0 = [0.0000 0.0000 0.0000];
clr1 = [0.0000 0.4470 0.7410];
clr2 = [0.8500 0.3250 0.0980];
clr3 = [0.9290 0.6940 0.1250];
clr4 = [0.4940 0.1840 0.5560];
clr5 = [0.4660 0.6740 0.1880];
clr6 = [0.3010 0.7450 0.9330];
clr7 = [0.6350 0.0780 0.1840];

%% Flow Condition
gma = 1.4;          % [-]   Specific Heat Ratio

% Simplification
pu = (gma - 1) / 2;

M1 = 2;             % [-]   Inlet Mach Number
phiShock = 33.9304904;      % [deg] Shock Wave Angle % 39.7841

% Unit Conversion
phiShock = phiShock * (pi / 180);


%% Initial Conditions
% delta: Flow Angle right After the Shock (Oblique Shock Relations)
delta = atan2((2 * cot(phiShock) * (M1^2 * (sin(phiShock))^2 - 1)),... 
        (M1^2 * (gma + cos(2 * phiShock))) + 2);

% M2   : Flow Speed right After the Shock (Oblique Shock Relations)
Mn1 = M1 * sin(phiShock);
Mn2 = sqrt((1 + pu * Mn1^2) / (gma * Mn1^2 - pu));
M2  = Mn2 / sin(phiShock - delta);

% Normalized Speed
V_2    = MtoV_(M2, pu);
V_r2   =   V_2 * cos(phiShock - delta);
V_phi2 = - V_2 * sin(phiShock - delta);


%% ODE Info.
% ODE: <Taylor-Maccoll>
% x is phi
% y is V_r
% z is dV_r / dphi 
syms x y z v
v = [y z]';
f = @(x, v) [v(2)
             ((1 / pu) * v(1) * v(2)^2 - (1 - v(1)^2 - v(2)^2) * (2 * v(1) + v(2) * cot(x))) /...
             (1 - v(1)^2 - v(2)^2 * (1 + (1 / pu)))];

% Bounds
a = phiShock;      v_a = [V_r2 V_phi2]';


%% Numerical Method Settings
% Step Size: phi
h = - 10 ^ (- 1) * (pi / 180);


%% Main
% -------------------------------------------------------------------------
% <R-K> 4th Order
[phiConeRK4, xListRK4, vListRK4, numIterRK4] =...
rk4(h, v, f, a, v_a);

phiConeRK4

% V_ & M
V_ListRK4   = sqrt(vListRK4(:, 1).^2 + vListRK4(:, 2).^2);
MListRK4    = V_toM(V_ListRK4, pu);


% <R-K> 4th Order, Adaptive
epsTol = 10 ^ (- 10); % 12.6
[phiConeRK4A, xListRK4A, vListRK4A, hListRK4A, numIterRK4A] =...
rk4Adapt(h, v, f, a, v_a, epsTol);

phiConeRK4A

% V_ & M
V_ListRK4A  = sqrt(vListRK4A(:, 1).^2 + vListRK4A(:, 2).^2);
MListRK4A   = V_toM(V_ListRK4A, pu);


% % <R-K> 4th Order, Adaptive [NEW]
% h_start = h * 10 ^ 2;
% [phiConeRK4AN, xListRK4AN, vListRK4AN, hListRK4AN, numIterRK4AN] =...
% rk4AdaptNew(h_start, v, f, a, v_a);
% 
% phiConeRK4AN
% 
% % V_
% V_ListRK4AN = sqrt(vListRK4AN(:, 1).^2 + vListRK4AN(:, 2).^2);
% MListRK4AN  = V_toM(V_ListRK4AN, pu);


%% Plot: y
hold on
plot(xListRK4  , vListRK4  (:, 1), 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr1, 'MarkerFaceColor', clr1)
plot(xListRK4A , vListRK4A (:, 1), 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr2, 'MarkerFaceColor', clr2)
% plot(xListRK4AN, vListRK4AN(:, 1), 'LineStyle', '-', 'Marker', 'o',...
%     'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr3, 'MarkerFaceColor', clr3)


%% Plot Settings
%     xlabel(['\phi [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    xlabel(['\theta [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    ylabel('V_{normalized, r} [-]', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    
    set(gca,'fontsize', sizeGCA)
    % title('Channel 0')

% ----------------------------------------------------------------------------------------------------------- %
    % Legend Settings
    [~, objh] = legend({['<R-K>: 4th'], ['<R-K>: 4th, Adaptive']}, 'Location', 'southwest', 'FontSize', sizeLegend, 'Orientation', 'vertical');
    % title(lgd, 'Kd')
    
    % Change Legend Marker Size
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12);          %// set marker size as desired
    set(objhl, 'linewidth', 12); 
    
% ----------------------------------------------------------------------------------------------------------- %
    % Minor Ticks
    set(gca,'XMinorTick','on')
    set(gca,'yMinorTick','on')
    % xtickformat('%1.0f')
    % ytickformat('%.4f')
    
    set(gca,'linewidth', 1);
    
    % Axes Thickness
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.XAxis.Color = 'k';
    ax.YAxis.LineWidth = 2;
    ax.YAxis.Color = 'k';
    
    % Set Axes Limits
%     axis([xListRK4(1) xListRK4(end) -inf +inf])
    
    % Log Scale
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    
    grid on 
    box  on
    hold off

%% Plot: z
figure
hold on
plot(xListRK4  , vListRK4  (:, 2), 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr1, 'MarkerFaceColor', clr1)
plot(xListRK4A , vListRK4A (:, 2), 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr2, 'MarkerFaceColor', clr2)
% plot(xListRK4AN, vListRK4AN(:, 2), 'LineStyle', '-', 'Marker', 'o',...
%     'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr3, 'MarkerFaceColor', clr3)


%% Plot Settings
%     xlabel(['\phi [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    xlabel(['\theta [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    ylabel('V_{normalized, \phi} [-]', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    
    set(gca,'fontsize', sizeGCA)
    % title('Channel 0')

% ----------------------------------------------------------------------------------------------------------- %
    % Legend Settings
    [~, objh] = legend({['<R-K>: 4th'], ['<R-K>: 4th, Adaptive']}, 'Location', 'southwest', 'FontSize', sizeLegend, 'Orientation', 'vertical');
    % title(lgd, 'Kd')
    
    % Change Legend Marker Size
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12);          %// set marker size as desired
    set(objhl, 'linewidth', 12); 
    
% ----------------------------------------------------------------------------------------------------------- %
    % Minor Ticks
    set(gca,'XMinorTick','on')
    set(gca,'yMinorTick','on')
    % xtickformat('%1.0f')
    % ytickformat('%.4f')
    
    set(gca,'linewidth', 1);
    
    % Axes Thickness
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.XAxis.Color = 'k';
    ax.YAxis.LineWidth = 2;
    ax.YAxis.Color = 'k';
    
    % Set Axes Limits
%     axis([xListRK4(1) xListRK4(end) -inf +inf])
    
    % Log Scale
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    
    grid on 
    box  on
    hold off

%% Plot: V_
figure
hold on
plot(xListRK4  , V_ListRK4  , 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr1, 'MarkerFaceColor', clr1)
plot(xListRK4A , V_ListRK4A , 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr2, 'MarkerFaceColor', clr2)
% plot(xListRK4AN, V_ListRK4AN, 'LineStyle', '-', 'Marker', 'o',...
%     'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr3, 'MarkerFaceColor', clr3)


%% Plot Settings
%     xlabel(['\phi [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    xlabel(['\theta [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    ylabel('V_{normalized} [-]', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    
    set(gca,'fontsize', sizeGCA)
    % title('Channel 0')

% ----------------------------------------------------------------------------------------------------------- %
    % Legend Settings
    [~, objh] = legend({['<R-K>: 4th'], ['<R-K>: 4th, Adaptive']}, 'Location', 'northwest', 'FontSize', sizeLegend, 'Orientation', 'vertical');
    % title(lgd, 'Kd')
    
    % Change Legend Marker Size
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12);          %// set marker size as desired
    set(objhl, 'linewidth', 12); 
    
% ----------------------------------------------------------------------------------------------------------- %
    % Minor Ticks
    set(gca,'XMinorTick','on')
    set(gca,'yMinorTick','on')
    % xtickformat('%1.0f')
    % ytickformat('%.4f')
    
    set(gca,'linewidth', 1);
    
    % Axes Thickness
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.XAxis.Color = 'k';
    ax.YAxis.LineWidth = 2;
    ax.YAxis.Color = 'k';
    
    % Set Axes Limits
%     axis([xListRK4(1) xListRK4(end) -inf +inf])
    
    % Log Scale
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    
    grid on 
    box  on
    hold off

%% Plot: M
figure
hold on
plot(xListRK4  , MListRK4  , 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr1, 'MarkerFaceColor', clr1)
plot(xListRK4A , MListRK4A , 'LineStyle', '-', 'Marker', 'o',...
    'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr2, 'MarkerFaceColor', clr2)
% plot(xListRK4AN, MListRK4AN, 'LineStyle', '-', 'Marker', 'o',...
%     'LineWidth', 1, 'MarkerSize', sz1, 'Color', clr3, 'MarkerFaceColor', clr3)


%% Plot Settings
%     xlabel(['\phi [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    xlabel(['\theta [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    ylabel('Mach Number _{downstream} [-]', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    
    set(gca,'fontsize', sizeGCA)
    % title('Channel 0')

% ----------------------------------------------------------------------------------------------------------- %
    % Legend Settings
    [~, objh] = legend({['<R-K>: 4th'], ['<R-K>: 4th, Adaptive']}, 'Location', 'northwest', 'FontSize', sizeLegend, 'Orientation', 'vertical');
    % title(lgd, 'Kd')
    
    % Change Legend Marker Size
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12);          %// set marker size as desired
    set(objhl, 'linewidth', 12); 
    
% ----------------------------------------------------------------------------------------------------------- %
    % Minor Ticks
    set(gca,'XMinorTick','on')
    set(gca,'yMinorTick','on')
    % xtickformat('%1.0f')
    % ytickformat('%.4f')
    
    set(gca,'linewidth', 1);
    
    % Axes Thickness
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.XAxis.Color = 'k';
    ax.YAxis.LineWidth = 2;
    ax.YAxis.Color = 'k';
    
    % Set Axes Limits
%     axis([xListRK4(1) xListRK4(end) -inf +inf])
    
    % Log Scale
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    
    grid on 
    box  on
    hold off

%% Plot: h
figure
hold on
plot(1:numIterRK4  , abs(h) * (180 / pi).* ones(1, length(1: numIterRK4)),...
    'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', sz1,...
    'Color', clr4, 'MarkerFaceColor', clr4)
plot(1:numIterRK4A, abs(hListRK4A) * (180 / pi),...
    'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', sz1,...
    'Color', clr5, 'MarkerFaceColor', clr5)
% plot(1:numIterRK4AN, abs(hListRK4AN),...
%     'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', sz1,...
%     'Color', clr6, 'MarkerFaceColor', clr6)


%% Plot Settings
    xlabel('Number of Iterations [-]', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    ylabel(['h (Step Size) [', char(176), ']'], 'FontSize', sizeLabel, 'FontWeight', 'bold')
    
    set(gca,'fontsize', sizeGCA)
    % title('Channel 0')

% ----------------------------------------------------------------------------------------------------------- %
    % Legend Settings
    [~, objh] = legend({['<R-K>: 4th'], ['<R-K>: 4th, Adaptive']}, 'Location', 'northeast', 'FontSize', sizeLegend, 'Orientation', 'vertical');
    % title(lgd, 'Kd')
    
    % Change Legend Marker Size
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12);          %// set marker size as desired
    set(objhl, 'linewidth', 12); 
    
% ----------------------------------------------------------------------------------------------------------- %
    % Minor Ticks
    set(gca,'XMinorTick','on')
    set(gca,'yMinorTick','on')
    % xtickformat('%1.0f')
    % ytickformat('%.4f')
    
    set(gca,'linewidth', 1);
    
    % Axes Thickness
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.XAxis.Color = 'k';
    ax.YAxis.LineWidth = 2;
    ax.YAxis.Color = 'k';
    
    % Set Axes Limits
%     axis([-inf +inf -inf -10^(0.5)])
    
    % Log Scale
%     set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    
    grid on 
    box  on
    hold off
        
%% Helpers
function V_ = MtoV_(M, pu)
    
    V_ = (1 / (pu .* M.^2) + 1).^(- 1 / 2);

end

function M = V_toM(V_, pu)
    
    M  = sqrt(1 ./ (pu .* (V_.^(- 2) - 1)));

end

function [phiCone, xList, vList, numIter] = rk4(h, v, f, a, v_a)

    % Iteration
    xList = [];
    vList = [];
    
    I = 0;
    while true
            
        if I == 0
            x = a;
            v = v_a;
            
        else
            k1 = f(x, v);
            k2 = f(x + (h / 2), v + (h / 2) .* k1);
            k3 = f(x + (h / 2), v + (h / 2) .* k2);
            k4 = f(x + h, v + h .* k3);

            k = (1 / 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4);
            v = v + h .* k;
            
            % march
            x = x + h;
            
        end
        
        % V_phi
        V_phi = v(2);
        
        if V_phi > 0   
            break;      
        end
        
        % log
        xList = [xList x];
        vList = [vList v];
        
        % Count Iterations
        I = I + 1;

    end
    
    % Transpose
    xList = xList' * (180 / pi);
    vList = vList';
    
    % Cone Angle & Unit Conversion
    phiCone = x * (180 / pi);
    
    % Number of Iterations
    numIter = I;

end

function [phiCone, xList, vList, hList, numIter] = rk4Adapt(h, v, f, a, v_a, epsTol)

    % Iteration
    xList = [];
    vList = [];
    hList = [];
    
    I = 0;
    while true

        if I == 0
            x = a;
            v = v_a;

        else
            % 1. Step Size: h1 = h
            h1 = h;
            v1 = v;
            
            k1 = f(x, v1);
            k2 = f(x + (h1 / 2), v1 + (h1 / 2) .* k1);
            k3 = f(x + (h1 / 2), v1 + (h1 / 2) .* k2);
            k4 = f(x + h1, v1 + h1 .* k3);

            k = (1 / 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4);
            v1 = v1 + h1 .* k;
            
            % 2. Step Size: h2 = h1 / 2
            h2 = h1 / 2;
            v2 = v;
            
            for J = 1 : 2
                k1 = f(x, v2);
                k2 = f(x + (h2 / 2), v2 + (h2 / 2) .* k1);
                k3 = f(x + (h2 / 2), v2 + (h2 / 2) .* k2);
                k4 = f(x + h2, v2 + h2 .* k3);

                k = (1 / 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4);
                v2 = v2 + h2 .* k;
                
            end
            
            % Evaluation
            v = v1;
            
            % Update Step Size
            h = h1 * nthroot((15 / 16) * epsTol / max(v2 - v1), 5);
            h = h1 * (abs((15 / 16) * epsTol / max(v2 - v1))) ^ (1 / 5);    % abs(): is it right to use hear?
            
            % march
            x = x + h;
            
        end
        
        % V_phi
        V_phi = v(2);
        
        if V_phi > 0   
            break;      
        end

        % log
        xList = [xList x];
        vList = [vList v];
        hList = [hList h];
        
        % Count Iterations
        I = I + 1;

    end
    
    % Transpose
    xList = xList'  * (180 / pi);
    vList = vList';
    hList = hList';
    
    % Cone Angle & Unit Conversion
    phiCone = x * (180 / pi);
    
    % Number of Iterations
    numIter = I;

end

function [phiCone, xList, vList, hList, numIter] = rk4AdaptNew(h_start, v, f, a, v_a)
% Why this Not converge?

    % Iteration
    xList = [];
    vList = [];
    hList = [];
    
    I = 0;
    while true

        if I == 0
            x = a;
            v = v_a;
            h = h_start;

        else
            k1 = f(x, v);
            k2 = f(x + (h / 2), v + (h / 2) .* k1);
            k3 = f(x + (h / 2), v + (h / 2) .* k2);
            k4 = f(x + 1, v + h .* k3);

            k = (1 / 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4);
            v = v + h .* k;
            
            % Update Step Size
            k_h = 10 ^ - 5;
            if h + k_h < 0
                h = h + k_h;
            end
            
            % march
            x = x + h;
            
        end
        
        % V_phi
        V_phi = v(2);
        
        if V_phi > 0   
            break;      
        end

        % log
        xList = [xList x];
        vList = [vList v];
        hList = [hList h];
        
        % Count Iterations
        I = I + 1;

    end
    
    % Transpose
    xList = xList'  * (180 / pi);
    vList = vList';
    hList = hList';
    
    % Cone Angle & Unit Conversion
    phiCone = x * (180 / pi);
    
    % Number of Iterations
    numIter = I;

end

function [xList, vList] = abm3(h, numIter, v, f, a, v_a, epsTol)
% needs verification & improvements

xList = [];
vList = [];
% <R-4>: 4th for the first some steps
    for I = 0 : 2

        if I == 0
            x = a;
            v = v_a;

        else
            k1 = f(x, v);
            k2 = f(x + (h / 2), v + (h / 2) .* k1);
            k3 = f(x + (h / 2), v + (h / 2) .* k2);
            k4 = f(x + h, v + h .* k3);

            k = (1 / 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4);
            v = v + h .* k;

            % march
            x = x + h;

        end
        
        % log
        xList = [xList x];
        vList = [vList v];
        
    end
    
% <A-B> for the rest steps
    vListTemp = [];
    for I = 3 : numIter
        
        k = (23 .* f(x, vList(:, I)) - 16 .* f(x, vList(:, I - 1)) + 5 .* f(x, vList(:, I - 2))) ./ 12;
        v = v + h .* k;
        vList = [vList v];
        vListTemp = [vListTemp v];
        
        diffRatioList = [];
%         while true
%             k = (5 .* f(x, vList(:, I + 1)) + 8 .* f(x, vList(:, I)) - f(x, vList(:, I - 1))) ./ 12;
%             v = v + h .* k;
%             vList(:, I + 1) = v;
%             vListTemp = [vListTemp v];
% 
%             diffRatio = max(abs((vListTemp(:, end) ./ vListTemp(:, end - 1)) - 1));
%             diffRatioList = [diffRatioList diffRatio];
%             if size(vListTemp, 2) >= 2
%                 if (diffRatio < epsTol)
%                     break
%                     
%                 elseif (length(diffRatioList) >= 2) && (diffRatioList(end) > diffRatioList(end - 1))
%                     break
%                     
%                 end
%             end
%         end 

        % march
        x = x + h;

        % log
        xList = [xList x];
        
    end
    vListTemp = [];
  
    xList = xList';
    vList = vList';
    
end


%% Prototypes
% delta = atan((2 * cot(phiShock) * (M1^2 * (sin(phiShock))^2 - 1) / (M1^2 * (gma + cos(2 * phiShock))) + 2));

%         % V_phi
%         if I == 0
%             V_phi = V_phi2;
%             
%         elseif I == 1
%             V_phi = (vList(1, I + 1) - vList(1, I)) / h; % 2 Point FWD Difference
%             
%         else
%             V_phi = (vList(1, I + 1) - vList(1, I - 1)) / (2 * h); % 2 Point CNT Difference
%             
%         end

% f = @(x, v) [v(2)
%              pu.*(2*((v(1)).^3-v(1))+(v(2)-v(2).*v(1).^2-v(2).^3).*cot(x)-2*v(1).*v(2).^2)-(v(1).^2.*v(2).^2)./...
%              pu.*((v(2)).^2+(v(1)).^2 - 1) + (v(2)).^2];
%% Cubic Splines
% apply cubic splines to discrete data set.
% in this case, altitude data from KSP.

close all
clear all
clc


%% Plot Settings
% =============================================================
% Plot Font Sizes
sizeLabel = 38;
sizeLegend = 32;
sizeGCA = 32;
  
sz1 = 7;
sz2 = 45;
sz3 = 25;
clr0 = [0.0000 0.0000 0.0000];
clr1 = [0.0000 0.4470 0.7410];
clr2 = [0.8500 0.3250 0.0980];
clr3 = [0.9290 0.6940 0.1250];
clr4 = [0.4940 0.1840 0.5560];

%% Q.8
% Predictors & Responses
y = [71.39725771255326
     71.397855296731
     71.39839147555176
     71.39741565159056
     71.3974290312035
     71.36871784948744
     71.23688948049676
     70.97120439854916
     71.20402172242757
    103.34535335097462
    161.07866703846958
    239.0367215400329
    334.667786396225
    448.0987208780134
    727.2068122352939
   1061.507182889036
   1283.178696324234
   1293.8158878617687
   1246.8578533106484
   1117.750316832331
    741.628542227787
    440.0477964589372
    291.3427893407643
    169.56280187994707
    117.05500181450043
     93.23358271701727
     75.59569252317306
     69.88938237878028
     74.3204794091871
     71.33144801703747];

x = [0 : 3 : 3 * (length(y) - 1)]';

% Number of Data Points
N = length(x);

% symbolic
syms x_
f = @(x_) [x_ ^ 3, x_ ^ 2, x_ ^ 1, x_ ^ 0];
Df1 = diff(f, x_, 1);
Df2 = diff(f, x_, 2);


%% Construct the X & Y Arrays
% <cond1> Penetration
for k = 1 : N - 1
    X1(2 * (k - 1) + 1 : 2 * k, 4 * (k - 1) + 1 : 4 * k) = [f(x(k)    )
                                                            f(x(k + 1))];
end

Y = [y(1); repelem(y(2 : end - 1), 2); y(end)];

% <cond2> Continuity: 1st Derivative
for k = 1 : N - 2
    X2(k, 4 * (k - 1) + 1 : 4 * (k    )) = + double(subs(Df1, x_, x(k + 1)));
    X2(k, 4 * (k    ) + 1 : 4 * (k + 1)) = - double(subs(Df1, x_, x(k + 1)));
end

Y = [Y; zeros(N - 2, 1)];

% <cond3> Continuity: 2nd Derivative
for k = 1 : N - 2
    X3(k, 4 * (k - 1) + 1 : 4 * (k    )) = + double(subs(Df2, x_, x(k + 1)));
    X3(k, 4 * (k    ) + 1 : 4 * (k + 1)) = - double(subs(Df2, x_, x(k + 1)));
end

Y = [Y; zeros(N - 2, 1)];

% <cond4> Natural Cubic Splines
    X4(1,               1 : 4          ) = double(subs(Df2, x_, x(1  )));
    X4(2, 4 * (N - 2) + 1 : 4 * (N - 1)) = double(subs(Df2, x_, x(end)));

Y = [Y; zeros(2, 1)];

% Final Matrix X
X = [X1; X2; X3; X4];

%% Solve for Coefficients a
a = X \ Y;

%% Retrieve Polynomials
for k = 1 : N - 1
    F(k) = [x_ ^ 3, x_ ^ 2, x_ ^ 1, x_ ^ 0] * a(4 * (k - 1) + 1 : 4 * k);
end


%% Estimate Value at X0
x0 = 12.7;
for k = 1 : N - 1
    if x0 >= x(k) && x0 <= x(k + 1)
        k0 = k;
    end
end
double(subs(F(k0), x_, x0))


%% Plot
hold on
plot(x, y, 'Color', clr0, 'LineStyle', '--', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', sz1, 'MarkerFaceColor', clr0)
for k = 1 : N - 1
    fplot(F(k), [x(k), x(k + 1)], 'LineStyle', '-', 'Marker', 'none', 'LineWidth', 2, 'MarkerSize', sz1)
    pause(0.025)
end

%% Plot Settings
    xlabel('x', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    ylabel('y', 'FontSize', sizeLabel, 'FontWeight', 'bold')
    
    set(gca,'fontsize', sizeGCA)
    % title('Channel 0')

% ----------------------------------------------------------------------------------------------------------- %
%     % Legend Settings
%     [~, objh] = legend({['x1'], ['x2'], ['x3']}, 'Location', 'southeast', 'FontSize', sizeLegend, 'Orientation','horizontal');
%     % title(lgd, 'Kd')
%     
%     % Change Legend Marker Size
%     objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
%     set(objhl, 'Markersize', 12);          %// set marker size as desired
%     set(objhl, 'linewidth', 12); 
    
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
%     axis([x(1) x(end) -inf +inf])
    
    % Log Scale
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    
    grid on 
    box  on
    hold off
    
%% Helpers

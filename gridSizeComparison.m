% Grid size comparison %

close all
set(gcf,'units','centimeters','position',[2 2 30 20])
set(0,'defaultTextInterpreter','tex');
set(0,'defaultfigurecolor',[1 1 1]);
% Adapt path to where you keep simulation data
path = "C:\Users\tonge\OneDrive\Documents\Durham\Engineering\Year Four\Capstone Project\Simulations\Precursor Simulations\Matlab Files\";
% Case variables, enter the name of the case folder 
case1 = "U08"; 
case2 = "U10"; 
case3 = "U20"; 


[k, Uhub, zhub, heights, startAvg, endAvg] = readConst();

% Specify the surface roughness you entered into SOWFA
z01 = 0.1;               
z02 = 0.1;               
z03 = 0.1;

% Some variables for the legends 
Z0 = "$$z_{0} = $$" + num2str(z01) + "$$m$$";



% Read in the data 

[U1, UU1, T1] = readData(path + case1 + "/");
[U2, UU2, T2] = readData(path + case2 + "/");
[U3, UU3, T3] = readData(path + case3 + "/");

[U1, UU1, T1] = removeData(U1, UU1, T1, startAvg, endAvg);
[U2, UU2, T2] = removeData(U2, UU2, T2, startAvg, endAvg);
[U3, UU3, T3] = removeData(U3, UU3, T3, startAvg, endAvg);

% Directly specify heights

heights1 = [4 12 20 28 36 44 52 60 68 76 84 92 100 108 116 124 132 140 148 156 164 172 180 188 196 204 212 220 228 236 244 252 260 268 276 284 292 300 308 316 324 332 340 348 356 364 372 380 388 396 404 412 420 428 436 444 452 460 468 476 484 492 500 508 516 524 532 540 548 556 564 572 580 588 596 604 612 620 628 636 644 652 660 668 676 684 692 700 708 716 724 732 740 748 756 764 772 780 788 796 804 812 820 828 836 844 852 860 868 876 884 892 900 908 916 924 932 940 948 956 964 972 980 988 996 1004 1012 1020 1028 1036 1044 1052 1060 1068 1076 1084 1092 1100 1108 1116 1124 1132 1140 1148 1156 1164 1172 1180 1188 1196 1204 1212 1220 1228 1236 1244 1252 1260 1268 1276];
heights2 = heights;
heights3 = [10 30 50 70 90 110 130 150 170 190 210 230 250 270 290 310 330 350 370 390 410 430 450 470 490 510 530 550 570 590 610 630 650 670 690 710 730 750 770 790 810 830 850 870 890 910 930 950 970 990 1010 1030 1050 1070 1090 1110 1130 1150 1170 1190 1210 1230 1250 1270];


% Write data to individual arrays ready for plotting 
[u1, v1, w1, M1, uu1, uv1, uw1, vv1, vw1, ww1, t1] = createArrays(U1, UU1, T1, heights1);
[u2, v2, w2, M2, uu2, uv2, uw2, vv2, vw2, ww2, t2] = createArrays(U2, UU2, T2, heights2);
[u3, v3, w3, M3, uu3, uv3, uw3, vv3, vw3, ww3, t3] = createArrays(U3, UU3, T3, heights3);

%%%%%%%%%%%%
close all
set(gcf,'units','centimeters','position',[2 2 11 9])

% 
% 
figure 
X1 = (heights1/z01);
fit1 = polyfit(log(X1(5:7)),M1(5:7),1);
plot1 = polyval(fit1,log(X1));
semilogx(X1,M1,'o')
hold on
semilogx(X1,plot1)
grid on 
ustar1 = fit1(1) * k;
figure
plot(M1,heights1,M2,heights2,M3,heights3);
legend('8m','10m','20m');

figure 
X1 = (heights1/z01);
fit1 = polyfit(log(X1(5:7)),M1(5:7),1);
plot1 = polyval(fit1,log(X1));
semilogx(X1,M1,'o')
hold on
semilogx(X1,plot1)

ustar1 = fit1(1) * k;

X2 = (heights2/z02);
X3 = (heights3/z03);


fit2 = polyfit(log(X2(5:7)),M2(5:7),1);
fit3 = polyfit(log(X3(5:7)),M3(5:7),1);


plot2 = polyval(fit2,log(X2));
plot3 = polyval(fit3,log(X3));



ustar2 = fit2(1) * k;
ustar3 = fit3(1) * k;



% Dimensionless mean shear


dudz1 = gradient(M1,8);
dudz2 = gradient(M2,10);
dudz3 = gradient(M3,20);


phi1 = (k*heights1/ustar1) .* dudz1;
phi2 = (k*heights2/ustar2) .* dudz2;
phi3 = (k*heights3/ustar3) .* dudz3;

figure
plot(phi1,heights1, '-xk')
hold on
plot(phi2,heights2, '-xb')
hold on
plot(phi3,heights3, '-xr')

axis([0 2.4 0 100])

legend('8m','10m','20m','AutoUpdate','off')
yline(75,'--')
xline(1,'--')


% f1 = fit(M1',heights1','smoothingspline');
% f2 = fit(M2',heights2','smoothingspline');
% f3 = fit(M3',heights3','smoothingspline');

f1 = fit(heights1',M1','smoothingspline');
f2 = fit(heights2',M2','smoothingspline');
f3 = fit(heights3',M3','smoothingspline');

fitted1 = f1(heights1(1:160))';
fitted2 = f2(heights1(1:160))';
fitted3 = f3(heights1(1:160))';

figure
plot(fitted1,heights1,...
    fitted2,heights1,...
    fitted3,heights1)

pc_fine   = abs((fitted1-fitted2)./fitted1)*100;
pc_course = abs((fitted1-fitted3)./fitted1)*100;

pc_data1 = (pc_course(1:158));
pc_data2 = (pc_fine(1:158));

close all
set(gcf,'units','centimeters','position',[2 2 11 9])

h = heights1/60;
h_data = (h(1:158));
plot(pc_data1,h_data,'--k','LineWidth',1);
hold on
plot(pc_data2,h_data,':r','LineWidth',1);
yline(11.7,'--');
yline(13.2,'--');
axis([-1 12 0 16])
% legend('\it \fontname{Times New Roman} \fontsize{12} S20 to S08', '\it \fontname{Times New Roman} \fontsize{12} S10 to S08');
% xlabel('\it \fontname{Times New Roman} \fontsize{14} Difference (%)')
% ylabel('\it \fontname{Times New Roman} \fontsize{14} z/z_h')
l = legend('U20 to U08','U10 to U08','AutoUpdate','off');
set(l,'Interpreter', 'latex', 'FontName', 'Times New Roman', 'fontsize', 11, 'FontAngle', 'italic');
d = ylabel("$$z/z_{h}$$");
set(d,'Interpreter', 'latex', 'fontsize',11);
d = xlabel("$$\Delta\overline{U}_\%$$"+" "+"$$(\%)$$");
set(d,'Interpreter', 'latex', 'fontsize',11);

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 11;

ax = gca;
exportgraphics(ax,"pcdiff.png","Resolution",300)

ave_fine = mean(pc_fine);
ave_course = mean(pc_course);



clc;

t = out.tout;
TS_POS_I = out.TS_POS_I;
ts_q_b = out.TS_q;  % target body -> ECI, scalar-first [qw qx qy qz]

idx0 = 1;
idxf = length(t);

%% STL file
stlFile = 'ISS_stationary.stl';
TR = stlread(stlFile);

scale = 5000;
TR_big = triangulation(TR.ConnectivityList, TR.Points * scale);

figure(111); clf
set(gcf, 'Color', 'w');
ax = axes;
hold on;
grid off;
axis off;

%% Draw Earth
n = 50;
Re = 6378137;
[x, y, z] = sphere(n);
X = Re * x;
Y = Re * y;
Z = Re * z;

texfile = 'earth_texture.jpg';
img = imread(texfile);
img = flipud(img);

h_earth = surf(X, Y, Z, img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.08, ...
    'DiffuseStrength', 0.9, ...
    'AmbientStrength', 0.35, ...
    'FaceAlpha', 0.4);

%% ECI coordinate axes
eci_axis_len = Re * 1.15;
quiver3(0,0,0,eci_axis_len,0,0,'k','LineWidth',2,'MaxHeadSize',0.35);
quiver3(0,0,0,0,eci_axis_len,0,'k','LineWidth',2,'MaxHeadSize',0.35);
quiver3(0,0,0,0,0,eci_axis_len,'k','LineWidth',2,'MaxHeadSize',0.35);
text(eci_axis_len,0,0, 'ECI X', 'Color', 'k', 'FontSize', 12);
text(0,eci_axis_len,0, 'ECI Y', 'Color', 'k', 'FontSize', 12);
text(0,0,eci_axis_len, 'ECI Z', 'Color', 'k', 'FontSize', 12);

%% Target orbit path
h_orbit = plot3(TS_POS_I(:,1), TS_POS_I(:,2), TS_POS_I(:,3), ...
    'k:', 'LineWidth', 3);

%% STL-to-target-body rotation
% R_BS maps the STL/model coordinate frame into the target body frame.
R_BS = [ 0  0  1;
         0 -1  0;
         1  0  0 ];

%% Initial and final ISS
p0 = TS_POS_I(idx0, :).';
pf = TS_POS_I(idxf, :).';
R0_TI = quat2rotm(ts_q_b(idx0, :)); % target body -> ECI
Rf_TI = quat2rotm(ts_q_b(idxf, :)); % target body -> ECI

h0 = hgtransform('Parent', ax);
hf = hgtransform('Parent', ax);

h_initial = trisurf(TR_big, ...
    'FaceColor',  [0.95 0.25 0.10], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.55, ...
    'Parent', h0);

h_final = trisurf(TR_big, ...
    'FaceColor', [0.95 0.25 0.10], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.55, ...
    'Parent', hf);

H0 = eye(4);
H0(1:3, 1:3) = R0_TI * R_BS;
H0(1:3, 4) = p0;
h0.Matrix = H0;

Hf = eye(4);
Hf(1:3, 1:3) = Rf_TI * R_BS;
Hf(1:3, 4) = pf;
hf.Matrix = Hf;

plot3(p0(1), p0(2), p0(3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
plot3(pf(1), pf(2), pf(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);

 
%% Target body axes at each ISS attitude
body_axis_len = 2.0e6;
draw_body_axes(p0, R0_TI, body_axis_len, '0');
draw_body_axes(pf, Rf_TI, body_axis_len, 'f');

%% View options
camlight('right');
lighting gouraud;
material dull;
axis equal;
axis off;
axis(7.5e6 * [-1 1 -1 1 -1 1]);
view(-15, 40);
grid off;
axis off;
legend([h_earth, h_orbit, h_initial, h_final], ...
    {'Earth', 'Target orbit', 'Initial ISS', 'Final ISS'}, ...
    'Location', 'northeast');

%% 

figure(112);clf
set(gcf, 'Color', 'w');
ax = axes;
grid on; hold on;
xlabel('ECI X [m]');
ylabel('ECI Y [m]');
zlabel('ECI Z [m]');

% Transform 객체 생성
ht = hgtransform('Parent', ax);

% STL patch 생성 (Parent를 ht로 지정)
h = trisurf(TR_big, ...
    'FaceColor', [1 0 0], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1.0, ...
    'Parent', ht);

% 보기 옵션 
title('STL Viewer');

% camlight('headlight');
camlight('right');
lighting gouraud;
material dull;

% draw earth
%%%%%%%%%%%%%%%%%%%%%%%%
% rotate3d on;
n = 50;
Re = 6378137;

% sphere
[x, y, z] = sphere(n);
X = Re * x;
Y = Re * y;
Z = Re * z;

% texture image
texfile = 'earth_texture.jpg';
img = imread(texfile);

% orientation fix
img = flipud(img);

% 변환 객체 생성
ht_earth = hgtransform;

% Earth
hs = surf(X, Y, Z, img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.08, ...
    'DiffuseStrength', 0.9, ...
    'AmbientStrength', 0.35, ...
    'FaceAlpha', 0.5, ...
    'Parent', ht_earth);   % <-- 여기 중요
view(-17,30)
camlight('right');
axis equal;
axis(8e+06*[-1 1 -1 1 -1 1])

% initial plot
R = quat2rotm(ts_q_b(1,:));
ex = R * [1;0;0] * 2e+6;
ey = R * [0;0;1] * 2e+6;
ez = R * [0;1;0] * 2e+6;

% ECI 좌표축
quiver3(0,0,0,Re*1.4,0,0,'k','LineWidth',2);
quiver3(0,0,0,0,Re*1.4,0,'k','LineWidth',2);
quiver3(0,0,0,0,0,Re*1.4,'k','LineWidth',2);
text(Re*1.3,0,0, 'ECI X [m]', ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');
text(0,Re*1.3,0, 'ECI Y [m]', ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');
text(0,0,Re*1.3, 'ECI Z [m]', ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

% Target 경로
plot3(TS_POS_I(:,1),TS_POS_I(:,2),TS_POS_I(:,3),'k:',"LineWidth",3)
% Target 좌표축
qx = quiver3(TS_POS_I(1,1),TS_POS_I(1,2),TS_POS_I(1,3),ex(1),ex(2),ex(3),1,'r',LineWidth=1.5);
qy = quiver3(TS_POS_I(1,1),TS_POS_I(1,2),TS_POS_I(1,3),ey(1),ey(2),ey(3),1,'g',LineWidth=1.5);
qz = quiver3(TS_POS_I(1,1),TS_POS_I(1,2),TS_POS_I(1,3),ez(1),ez(2),ez(3),1,'b',LineWidth=1.5);

% legend([qx, qy, qz],{'x','y','z'})
R_BS = [ 0  0  1;
         0 -1  0;
         1  0  0 ];

H_BS = eye(4);
H_BS(1:3,1:3) = R_BS;

for i=1:300:length(t)
    eul = quat2eul(ts_q_b(i,:), 'ZYX');

    yaw   = eul(1);
    pitch = eul(2);
    roll  = eul(3);

    % Rz = makehgtform('zrotate', yaw); % rad
    % Ry = makehgtform('yrotate', pitch);
    % Rx = makehgtform('xrotate', roll);
    % T = makehgtform('translate',TS_POS_I(i,:));
    % 적용 순서: 오른쪽부터 먼저 적용
    % ht.Matrix = T * Rz * Ry * Rx;

    T = makehgtform('translate',TS_POS_I(i,:),'zrotate', yaw,'yrotate', pitch,'xrotate', roll);
    ht.Matrix = T*H_BS;

    R = quat2rotm(ts_q_b(i,:));
    ex = R * [1;0;0] * 2e+6;
    ey = R * [0;1;0] * 2e+6;
    ez = R * [0;0;1] * 2e+6; 
    set(qx, 'XData', TS_POS_I(i,1), 'YData', TS_POS_I(i,2), 'ZData', TS_POS_I(i,3), ...
    'Udata',ex(1),'Vdata',ex(2),'Wdata',ex(3));
    set(qy, 'XData', TS_POS_I(i,1), 'YData', TS_POS_I(i,2), 'ZData', TS_POS_I(i,3), ...
    'Udata',ey(1),'Vdata',ey(2),'Wdata',ey(3));
    set(qz, 'XData', TS_POS_I(i,1), 'YData', TS_POS_I(i,2), 'ZData', TS_POS_I(i,3), ...
    'Udata',ez(1),'Vdata',ez(2),'Wdata',ez(3));

    drawnow;
end

 






%% Local functions
function draw_body_axes(p, R_TI, axis_len, suffix)
    ex = R_TI * [1;0;0] * axis_len*1.1;
    ey = R_TI * [0;1;0] * axis_len*1.1;
    ez = R_TI * [0;0;1] * axis_len*1.1;

    quiver3(p(1), p(2), p(3), ex(1)/1.1, ex(2)/1.1, ex(3)/1.1, ...
        'r', 'LineWidth', 1.6, 'MaxHeadSize', 0.35);
    quiver3(p(1), p(2), p(3), ey(1)/1.1, ey(2)/1.1, ey(3)/1.1, ...
        'g', 'LineWidth', 1.6, 'MaxHeadSize', 0.35);
    quiver3(p(1), p(2), p(3), ez(1), ez(2), ez(3), ...
        'b', 'LineWidth', 1.6, 'MaxHeadSize', 0.35);

    text(p(1)+ex(1), p(2)+ex(2), p(3)+ex(3), ['X_B_', suffix], ...
        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
    text(p(1)+ey(1), p(2)+ey(2), p(3)+ey(3), ['Y_B_', suffix], ...
        'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
    text(p(1)+ez(1), p(2)+ez(2), p(3)+ez(3), ['Z_B_', suffix], ...
        'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');
end
% close all;
clc;

sim_t = out.tout;
ss_pos_i = out.SS_POS_I;
ss_vel_i = out.SS_VEL_I;
ss_acc_i = out.SS_ACC_I;
ss_angle_i = out.SS_angle_I; % deg
ss_q_b = out.SS_q_b; 
ss_rate_b = out.SS_w_b; % deg/s
ss_rated_b = out.SS_wd_b; %

ts_pos_i = out.TS_POS_I;
ts_vel_i = out.TS_VEL_I;
ts_acc_i = out.TS_ACC_I;
ts_rated_b = out.TS_wd; 
ts_rate_b = out.TS_wb; % deg
ts_angle_i = out.TS_angle_I; % deg
ts_q_b = out.TS_q;

rel_ts_pos = out.REL_TS_POS_I; % eci 기준 port
rel_ts_vel = out.REL_TS_VEL_I; % eci 기준 port
rel_ts_acc = out.REL_TS_ACC_I;
rel_pos_cg = ts_pos_i - ss_pos_i;

rel_pos_b = out.rel_pos_b;
rel_vel_b = out.rel_vel_b;
rel_angle_b = out.rel_angle_b;
rel_w_b = out.rel_w_b;

% los_angle = out.LoS_angle;
wheel_rate = out.wheel_rate; % rpm
wheel_rated = out.wheel_rated;

qp_F_cmd = out.RCS_QP_F_cmd;
qp_tau_cmd = out.RCS_QP_tau_cmd;
qp_rw_cmd = out.RW_QP_tau_cmd;
s_tau = out.S_tau;
s_dump = out.S_dump;
s_F = out.S_F;
q_err = quatmultiply(quatinv(ts_q_b), ss_q_b);

RCS_onoff = out.RCS_onoff;

%% plot
% 상대 위치 오차
figure(1);clf
sgtitle('Target-Chaser relative state in chaser body frame')
subplot(2,1,1)
hold on;grid on;ylabel('Pos [m]');ylim padded
plot(sim_t, rel_pos_b(:,1))
plot(sim_t, rel_pos_b(:,2))
plot(sim_t, rel_pos_b(:,3))
legend('x','y','z')
subplot(2,1,2)
hold on;grid on;ylabel('Vel [m/s]');xlabel('time [s]');ylim padded
plot(sim_t, rel_vel_b(:,1))
plot(sim_t, rel_vel_b(:,2))
plot(sim_t, rel_vel_b(:,3))
legend('x','y','z')
 
% 상대 자세 오차
figure(2);clf
sgtitle('Relative attitude error')
subplot(2,1,1)
hold on;grid on;ylabel('ang_{err} [deg]');ylim padded
plot(sim_t, rel_angle_b)
legend('roll','pitch','yaw')
subplot(2,1,2)
hold on;grid on;ylabel('\omega_{err} [deg/s]');xlabel('time [s]');ylim padded
plot(sim_t, rel_w_b)
legend('x','y','z')

% RW 상태
figure(3);clf
title('Wheel speed')
hold on;grid on;ylabel('[rpm]');xlabel('time [s]');ylim([-1000 4500]);xlim tight
plot(sim_t,wheel_rate)
yline(Omega_ref(1)*rad2rpm,'--','LineWidth',1)
legend('w_1','w_2','w_3','w_4','reference','Location','northeast')



% qp cmd
figure(4);clf
sgtitle('QP command')
subplot(3,1,1)
hold on;grid on;ylabel('RCS F cmd [N]');ylim padded;xlim tight
plot(sim_t,qp_F_cmd)
legend('x','y','z')
subplot(3,1,2)
hold on;grid on;ylabel('RCS torque cmd [N m]');ylim padded;xlim tight
plot(sim_t,qp_tau_cmd)
legend('x','y','z')
subplot(3,1,3)
hold on;grid on;ylabel('RW torque cmd [N m]');xlabel('time [s]');ylim padded;xlim tight
plot(sim_t,qp_rw_cmd)
legend('w_1','w_2','w_3','w_4')
  
%slack 변수
figure(5);clf
sgtitle('QP slack variable')
subplot(3,1,1)
hold on;grid on;ylabel('s_{\tau} [N m]');ylim padded;xlim tight
plot(sim_t,s_tau)
legend('x','y','z')
subplot(3,1,2)
hold on;grid on;ylabel('s_{dump} [N m]');ylim padded;xlim tight
plot(sim_t,s_dump)
legend('x','y','z')
subplot(3,1,3)
hold on;grid on;ylabel('s_F [N]');xlabel('time [s]');ylim padded;xlim tight
plot(sim_t,s_F)
legend('x','y','z')
 

% RCS duty
figure(7);clf
sgtitle('RCS on/off cmd')
for i = 1:16
    subplot(4,4,i)
    plot(sim_t,RCS_onoff(:,i))
    grid on; ylim([-0.2 1.2])
end

% RCS on cmd count
figure(8);clf
sgtitle('RCS on command count'); hold on; grid on;
bar(1:16,sum(RCS_onoff)); xlabel('Thruster Num');xticks = 1:16;


% figure('Name','SS Inertial POS Data Log [m]');
% plot(sim_t, ss_pos_i(:,1));
% hold on; grid on;
% plot(sim_t, ss_pos_i(:,2));
% plot(sim_t, ss_pos_i(:,3));
% 
% figure('Name','SS Inertial POS Data Log 3D [m]');
% plot3(ss_pos_i(:,1), ss_pos_i(:,2), ss_pos_i(:,3),'k.');
% grid on;
% axis equal;
% xlabel('POS I X [m]');
% ylabel('POS I Y [m]');
% zlabel('POS I Z [m]');
%
% figure('Name','SS Inertial VEL Data Log [m/s]');
% plot(sim_t, ss_vel_i(:,1));
% hold on; grid on;
% plot(sim_t, ss_vel_i(:,2));
% plot(sim_t, ss_vel_i(:,3));
% 
% figure('Name','SS Inertial ACC Data Log [m/s^2]');
% plot(sim_t, ss_acc_i(:,1));
% hold on; grid on;
% plot(sim_t, ss_acc_i(:,2));
% plot(sim_t, ss_acc_i(:,3));
% 
% figure('Name','SS Inertial Angle Data Log [deg]');
% plot(sim_t, ss_angle_i(:,1));
% hold on; grid on;
% plot(sim_t, ss_angle_i(:,2));
% plot(sim_t, ss_angle_i(:,3));
% 
% figure('Name','SS Inertial Rate Data Log [deg/s]');
% plot(sim_t, ss_rate_b(:,1));
% hold on; grid on;
% plot(sim_t, ss_rate_b(:,2));
% plot(sim_t, ss_rate_b(:,3));
% 
% figure('Name','SS Inertial Rate d Data Log [deg/s^2]');
% title('SS Inertial Rate d Data Log [deg/s^2]')
% plot(sim_t, ss_rated_b(:,1));
% hold on; grid on;
% plot(sim_t, ss_rated_b(:,2));
% plot(sim_t, ss_rated_b(:,3));
% 
%
% figure('Name','TS Relative Position wrt SS [m]');
% plot3(rel_ts_pos(:,1),rel_ts_pos(:,2),rel_ts_pos(:,3),'k-');
% grid on;
% xlabel('Rel POS I X [m]');
% ylabel('Rel POS I Y [m]');
% zlabel('Rel POS I Z [m]');

 
% simulation result 
x_ref = 1.0; % [m]
idx_end = length(sim_t);

pos_end = rel_pos_b(idx_end,:);
vel_end = rel_vel_b(idx_end,:);
ang_end = rel_angle_b(idx_end,:);
w_end   = rel_w_b(idx_end,:);

axial_error = abs(pos_end(1) - x_ref);
lateral_error = norm(pos_end(2:3));
closing_rate = -vel_end(1);
lateral_rate = norm(vel_end(2:3));
att_error_norm = norm(ang_end);
ang_rate_norm = norm(w_end);
pitch_yaw_error = norm(rel_angle_b(idx_end,2:3));
pitch_yaw_rate_error = norm(rel_w_b(idx_end,2:3));

fprintf('\n===== Initial Contact Condition Metrics =====\n');
fprintf('Time : %.3f s | Axial separation : %.6f m | Axial error : %.6f m\n', ...
    sim_t(idx_end), pos_end(1), axial_error);
fprintf('Position : x = %.6f m | y = %.6f m | z = %.6f m | lateral = %.6f m\n', ...
    pos_end(1), pos_end(2), pos_end(3), lateral_error);
fprintf('Velocity : vx = %.6f m/s | vy = %.6f m/s | vz = %.6f m/s | closing = %.6f m/s | lateral = %.6f m/s\n', ...
    vel_end(1), vel_end(2), vel_end(3), closing_rate, lateral_rate);
fprintf('Angle : roll = %.6f deg | pitch = %.6f deg | yaw = %.6f deg | norm = %.6f deg\n', ...
    ang_end(1), ang_end(2), ang_end(3), pitch_yaw_error);
fprintf('Angular rate : wx = %.6f deg/s | wy = %.6f deg/s | wz = %.6f deg/s | norm = %.6f deg/s\n', ...
    w_end(1), w_end(2), w_end(3), pitch_yaw_rate_error);
fprintf('=============================================\n\n');


return

%% 3d motion
figure(100);clf
set(gcf, 'Color', 'w');
ax = axes;
title('Target-Chaser rel motion')
hold on; grid on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
axis equal
axis([-33 9 -4 4 -4 4])
view(-40,30)

% Target 모델
TS_stl = stlread('ISS_stationary.stl');
TS_stl = triangulation(TS_stl.ConnectivityList, TS_stl.Points * 0.36);
ht = hgtransform('Parent', ax);
h1 = trisurf(TS_stl, ...
    'FaceColor', [0.5 0.5 0.8], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.7, ...
    'Parent', ht);
camlight('right') 
% lighting gouraud

T = eye(4);
T(1:3,1:3) = [ 0  0  1;
               0 -1  0;
               1  0  0];
T(1:3,4) = [31.3 0.8 0]';
ht.Matrix = T;



% chaser 모델
SS_stl = stlread('geo_mod.stl');

V = SS_stl.Points*1;  %   확대


[F,V2] = reducepatch(SS_stl.ConnectivityList, V, 0.5); % 20%만 유지
SS_stl = triangulation(F, V2);


% Transform 객체 생성
hs = hgtransform('Parent', ax);

% STL patch 생성 (Parent를 ht로 지정)
h2 = trisurf(SS_stl, ...
    'FaceColor', [0.9 0.9 0.9], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.6, ...
    'Parent', hs);

 
 
% ss
s_p1 = plot3(0,0,0,'k.',"MarkerSize",20);  % chaser CG
s_qvx1 = quiver3(0,0,0,0,0,0,3,'r','LineWidth',1);
s_qvy1 = quiver3(0,0,0,0,0,0,3,'g','LineWidth',1);
s_qvz1 = quiver3(0,0,0,0,0,0,3,'b','LineWidth',1);

s_p2 = plot3(0,0,0,'r.',"MarkerSize",20);  % chaser docking port
s_qvx2 = quiver3(0,0,0,0,0,0,3,'r','LineWidth',1);
s_qvy2 = quiver3(0,0,0,0,0,0,3,'g','LineWidth',1);
s_qvz2 = quiver3(0,0,0,0,0,0,3,'b','LineWidth',1);

[cyl_x, cyl_y, cyl_z] = cylinder(0.07, 32); % rad, n
h_port = hgtransform('Parent', ax);
s_qvp = surf(cyl_x, cyl_y, cyl_z, ...
    'FaceColor', [0.5 0.5 0.5], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1, ...
    'Parent', h_port);


% ts
Rs = quat2rotm(ts_q_b(1,:));
ex = Rs * [1;0;0] * 0.5; 
ey = Rs * [0;1;0] * 0.5;
ez = Rs * [0;0;1] * 0.5;
t_p = plot3(0,0,0,'b.',"MarkerSize",20);
t_qvx = quiver3(0,0,0,ex(1),ex(2),ex(3),3,'r','LineWidth',1);
t_qvy = quiver3(0,0,0,ey(1),ey(2),ey(3),3,'g','LineWidth',1);
t_qvz = quiver3(0,0,0,ez(1),ez(2),ez(3),3,'b','LineWidth',1);

% legend([t_p s_p s_qvx s_qvy s_qvz], {'target','chaser','x','y','z'},'Location','southeast');
 
for i = 1:100:length(sim_t) 
 
    % ss 
    Rs = quat2rotm(q_err(i,:));
    r_s = Rs*0.5;
    rel_cg = quat2dcm(ts_q_b(i,:)) * rel_pos_cg(i,:)';
    set(s_p1,"XData",-rel_cg(1),"YData",-rel_cg(2),"ZData",-rel_cg(3))

    set(s_qvx1,"XData",-rel_cg(1),"YData",-rel_cg(2),"ZData",-rel_cg(3), ...
        "UData", r_s(1,1),"VData",r_s(2,1),"WData",r_s(3,1))
    set(s_qvy1,"XData",-rel_cg(1),"YData",-rel_cg(2),"ZData",-rel_cg(3), ...
        "UData", r_s(1,2),"VData",r_s(2,2),"WData",r_s(3,2))
    set(s_qvz1,"XData",-rel_cg(1),"YData",-rel_cg(2),"ZData",-rel_cg(3), ...
        "UData", r_s(1,3),"VData",r_s(2,3),"WData",r_s(3,3))

    % ss port
    rel_port = quat2dcm(ts_q_b(i,:)) * rel_ts_pos(i,:)';
    
    set(s_p2,"XData",-rel_port(1),"YData",-rel_port(2),"ZData",-rel_port(3))

    set(s_qvx2,"XData",-rel_port(1),"YData",-rel_port(2),"ZData",-rel_port(3), ...
        "UData", r_s(1,1),"VData",r_s(2,1),"WData",r_s(3,1))
    set(s_qvy2,"XData",-rel_port(1),"YData",-rel_port(2),"ZData",-rel_port(3), ...
        "UData", r_s(1,2),"VData",r_s(2,2),"WData",r_s(3,2))
    set(s_qvz2,"XData",-rel_port(1),"YData",-rel_port(2),"ZData",-rel_port(3), ...
        "UData", r_s(1,3),"VData",r_s(2,3),"WData",r_s(3,3))
    
    % port cylinder
 
    T_port = eye(4);
    T_port(1:3,1:3) = Rs*[ 0  0 -1;
                           0  1  0;
                           1  0  0];
    T_port(1:3,4)   = -rel_port(:);        % 원통 시작점, chaser docking port 위치

    h_port.Matrix = T_port;

    % spacecraft
    T = eye(4);
    T(1:3,1:3) = Rs;
    T(1:3,4) = -rel_cg(:)' ;
    hs.Matrix = T;

    
    drawnow;
    pause(0.003);
end

 
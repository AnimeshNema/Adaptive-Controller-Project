clc;
clear all;

vrep=remApi('remoteApi'); % using the prototype file (remoteApiProto.m)
vrep.simxFinish(-1); % just in case, close all opened connections

clientID=vrep.simxStart('127.0.0.1',19997,true,true,5000,5);

if (clientID>-1) % If no connection is established, client id will be -1.
    disp('Connected to remote API server');
    vrep.simxSynchronous(clientID,true);

    % joint_names = ['Joint1'; 'Joint2'];
    % ------- joint target velocities discussed below
    % joint_target_velocities = ones(size(joint_names(:,1))) * 10000.0;

    %-------- get the handles for each joint and set up streaming
    [returnCode, joint_handle_1] = vrep.simxGetObjectHandle(clientID, 'shoulder', vrep.simx_opmode_blocking);
    [returnCode, joint_handle_2] = vrep.simxGetObjectHandle(clientID, 'elbow', vrep.simx_opmode_blocking);

    %--------- get handle for target and set up streaming
    % [returnCode, target_handle] = vrep.simxGetObjectHandle(clientID,'target', vrep.simx_opmode_blocking);

    dt = .001; % specify a simulation time step
    % [returnCode] = vrep.simxSetFloatingParameter(clientID, vrep.sim_floatparam_simulation_time_step, dt, vrep.simx_opmode_oneshot);

    % --------------------- Start the simulation

    %-------- start our simulation in lockstep with our code
    [returnCode] = vrep.simxStartSimulation(clientID, vrep.simx_opmode_blocking);
%     [returnCode] = vrep.simxSetJointTargetPosition(clientID, joint_handle_1,-90, vrep.simx_opmode_blocking);
%     [returnCode] = vrep.simxSetJointTargetPosition(clientID, joint_handle_2,-10, vrep.simx_opmode_blocking);
    
%     [returnCode] = vrep.simxSetJointForce(clientID, joint_handle_1,2.5, vrep.simx_opmode_blocking);
%     [returnCode] = vrep.simxSetJointForce(clientID, joint_handle_2,0.5, vrep.simx_opmode_blocking);
    count = 0;
%     track_hand = []
%     track_target = []
    MBar = 0.0;
    while count < 1 % run for 1 simulated second

        [returnCode, q1] = vrep.simxGetJointPosition(clientID, joint_handle_1, vrep.simx_opmode_blocking);
        [returnCode, q2] = vrep.simxGetJointPosition(clientID, joint_handle_2, vrep.simx_opmode_blocking);

        [returnCode, dq1] = vrep.simxGetObjectFloatParameter(clientID, joint_handle_1, 2012, vrep.simx_opmode_blocking);
        [returnCode, dq2] = vrep.simxGetObjectFloatParameter(clientID, joint_handle_2, 2012, vrep.simx_opmode_blocking);

        % Execute Controller
        dx = Controller(q1, q2, dq1, dq2, MBar,count);
        q1 = q1 + dx(1);
        q2 = q2 + dx(2);
        dq1 = dq1 + dx(3);
        dq2 = dq2 + dx(4);
        MBar = dx(5);
        torque1 = dx(6);
        torque2 = dx(7);

        [returnCode] = vrep.simxSetJointTargetVelocity(clientID,joint_handle_1,10000,vrep.simx_opmode_blocking);
        [returnCode] = vrep.simxSetJointTargetVelocity(clientID,joint_handle_2,10000,vrep.simx_opmode_blocking);

        [returnCode] = vrep.simxSetJointForce(clientID, joint_handle_1, torque1, vrep.simx_opmode_blocking);
        [returnCode] = vrep.simxSetJointForce(clientID, joint_handle_2, torque2, vrep.simx_opmode_blocking);


        vrep.simxSynchronousTrigger(clientID);
        count = count + dt;
    end
     
    %--stop the simulation
    [returnCode] = vrep.simxStopSimulation(clientID, vrep.simx_opmode_blocking);

    %----------- Before closing the connection to V-REP,
    %---------make sure that the last command sent out had time to arrive.
    [returnCode] = vrep.simxGetPingTime(clientID);

    %-------- Now close the connection to V-REP:
    vrep.simxFinish(-1); % just in case, close all opened connections
    disp("Connection Closed....... ");
end

vrep.delete(); % This is the destructor




function [dx] = Controller(q1,q2,dq1,dq2,MBar,count)

    payload = 0.5236;
    % desired trajectories
    theta_d = [90; 90];
    dtheta_d = [10000; 10000];
    ddtheta_d = [0; 0];

    % given trajectories
    theta = [q1; q2];
    dtheta= [dq1; dq2];

    % errors
    global lambda e de a v r m1 l1 l2 g
    lambda = 0.99;
    e = theta - theta_d;
    de = dtheta - dtheta_d;
    a = ddtheta_d - (lambda*de);
    v = dtheta_d - (lambda*e);
    r = de + (lambda*e);

    % a positive definite matrix (to be used later for W_update)
    P = 0.2*eye(11); 

    % True model
    % global M1 G1 M2 C2 G2 M C G PM PG PC
    % actual dynamic model of the system is characterized by M and C
    % for link 1
    M1 = [(1/3)*2.5*(0.1)^2 , 0 ; 0, 0];
    G1 = [(1/2)*2.5*9.8*0.1*cos(q1); 0];

    % for link 2
    PM = [0.1^2 + (1/3)*0.1^2 + 0.1*0.1*cos(q2) , (1/3)*0.1^2 + (1/2)*0.1*0.1*cos(q2); (1/3)*0.1^2 + (1/2)*0.1*0.1*cos(q2), (1/3)*0.1^2];
    PC = [0, -((1/2)*0.1*0.1*sin(q2)*dq2 + 0.1*0.1*sin(q2)*dq1) ; (1/2)*0.1*0.1*sin(q2)*dq1, 0];
    PG = [(1/2)*9.8*0.1*cos(q1) + (1/2)*9.8*0.1*cos(q1 + q2) ; (1/2)*9.8*0.1*cos(q1+ q2)];
    M2 = (2.5 + payload)*PM;
    G2 = (2.5 + payload)*PG;
    C2 = (2.5 + payload)*PC;

    %actual model
    M = M1 + M2;
    C = C2;
    G = G1 + G2;
    invM = inv(M);
    invMC = inv(M)*C;
    % disp('Reached')
    % Fourier Series
    Z = [(1/2); cos((pi*(count))/5);sin((pi*(count))/5);cos((2*pi*(count))/5);sin((2*pi*(count))/5);cos((3*pi*(count))/5);sin((3*pi*(count))/5);cos((4*pi*(count))/5);sin((4*pi*(count))/5);cos((5*pi*(count))/5);sin((5*pi*(count))/5)];
    % disp('Not Reached')
    m2t_bar = MBar;

    % Estimated model
    global M_bar C_bar G_bar M2_bar C2_bar G2_bar
    M2_bar = m2t_bar*[0.1^2 + (1/3)*0.1^2 + 0.1*0.1*cos(q2) , (1/3)*0.1^2 + (1/2)*0.1*0.1*cos(q2); (1/3)*0.1^2 + (1/2)*0.1*0.1*cos(q2), (1/3)*0.1^2];
    C2_bar = m2t_bar*[0, -((1/2)*0.1*0.1*sin(q2)*dq2 + 0.1*0.1*sin(q2)*dq1) ; (1/2)*0.1*0.1*sin(q2)*dq1, 0];
    G2_bar = m2t_bar*[(1/2)*9.8*0.1*cos(q1) + (1/2)*9.8*0.1*cos(q1+ q2) ; (1/2)*9.8*0.1*cos(q1+ q2)];

    M_bar = M1 + M2_bar;
    C_bar = C2_bar;
    G_bar = G1 + G2_bar;
    
    % Torque
    tau = adaptive_ctrl(theta_d, dtheta_d, ddtheta_d, theta, dtheta);
    %global torque
    %torque = [torque, tau];

    %update the system state, compute dx
    dx=zeros(7,1);
    dx(1) = dq1;
    dx(2) = dq2;
    dx(3:4) = -invMC* [dq1; dq2] - invM*G + invM*tau; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
    %update law - W_update
    W_update = -inv(P)*[Z*transpose(e)*PM*a + Z*transpose(e)*PC*v + Z*transpose(e)*PG];
    dx(5) = transpose(W_update)*Z;
    dx(6) = tau(1,1);
    dx(7) = tau(2,1);

end

% function to calculate torque
function tau = adaptive_ctrl(theta_d, dtheta_d, ddtheta_d, theta, dtheta)
    global M C M_bar C_bar G_bar lambda e de a v r
    %Kp = 100*eye(1);
    Kv = 500*eye(2);
    tau = (M_bar*a)+ (C_bar*v) + (G_bar) - Kv*r;
end
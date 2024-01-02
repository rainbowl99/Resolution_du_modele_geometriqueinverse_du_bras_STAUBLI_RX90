function [conf] = IK_RX90(p,R,q_k,IndSolLoc)

    % p : End-effector position
    % R : End-effector orientation
    % q_k : Current configuration
    % IndSolLoc : Solution number
    % Get RX90 data (load length and DH params)
    [L2,L3,L6,dh] = RX90data;
    % There are 8 solutions to the IK of an RX90
    % We will store them in 8x6 array (one solution per line) called q
    q = zeros(8,6);
    % Compute 04 position in R0
    p_4 = p - R * [0;0;L6];
    %-------------------------------------------------------------

    %Q5.6:Traiter les poses non-atteignablese.
    %Quand c'est la pose non-atteignable, toutes les thetas sont 0.
    r = L2+L3+L6;
    if p(1)^2 + p(2)^2 + p(3)^2 >= r^2
        fprintf('La pose non-atteignable');
        conf = q(IndSolLoc,:);
        return
    end


    %-------------------------------------------------------------
    % Pour la situation que cos(pi/2) (sin(pi)) n'est pas egale a 0
    c = 10^(-15);
    %resolution of q1
    for i = 1:4
        q(i,1) = atan2(p_4(2),p_4(1));
    end
    for i = 5:8
        q(i,1) = atan2(p_4(2),p_4(1))+pi;
    end    
    for i = 1:8
        if rem(i,2) == 1
            epsilon = -1;
        else
            epsilon = 1;
        end
        %resolution of q2 and q3
        lambda = p_4(1) * cos(q(i,1)) + p_4(2) * sin(q(i,1));
        A = 2*L2*lambda;
        B = 2*L2*p_4(3);
        C = lambda^2 + p_4(3)^2 + L2^2 - L3^2;
        q(i,2) = atan2(B*C - epsilon*A*sqrt(A^2+B^2-C^2),A*C + epsilon*B*sqrt(A^2+B^2-C^2));
        q(i,3) = atan2(lambda * cos(q(i,2)) + p_4(3) * sin(q(i,2)) - L2, lambda * sin(q(i,2)) - p_4(3)* cos(q(i,2)));
        T_03 = TH(q(i,1),dh(1,:)) * TH(q(i,2),dh(2,:)) * TH(q(i,3),dh(3,:));
        T_06 = [R(1,1) R(1,2) R(1,3) p(1);R(2,1) R(2,2) R(2,3) p(2);R(3,1) R(3,2) R(3,3) p(3);0 0 0 1];
        T_36 = inv(T_03)*T_06;
        %resolution of q4,q5,q6
        if i == 1 || i == 2
            q(i,5) = atan2(sqrt(T_36(1,3)^2 + T_36(2,3)^2),T_36(3,3));
            if sin(q(i,5)) < c
                q(i,4) = -pi/2;
            else
                q(i,4) = atan2(T_36(2,4),T_36(1,4));
            end
            T_05 = TH(q(i,1),dh(1,:)) * TH(q(i,2),dh(2,:)) * TH(q(i,3),dh(3,:)) * TH(q(i,4),dh(4,:)) * TH(q(i,5),dh(5,:));
            T_56 = inv(T_05)*T_06;
            q(i,6) = atan2(T_56(2,1),T_56(1,1));
        elseif  i == 5 || i == 6
            q(i,5) = atan2(sqrt(T_36(1,3)^2 + T_36(2,3)^2),T_36(3,3));
            if sin(q(i,5)) < c
                q(i,4) = pi/2;
            else
                q(i,4) = atan2(T_36(2,4),T_36(1,4));
            end
            T_05 = TH(q(i,1),dh(1,:)) * TH(q(i,2),dh(2,:)) * TH(q(i,3),dh(3,:)) * TH(q(i,4),dh(4,:)) * TH(q(i,5),dh(5,:));
            T_56 = inv(T_05)*T_06;
            q(i,6) = atan2(T_56(2,1),T_56(1,1));    
        else
            q(i,5) = -1 * q(i-2,5);
            q(i,4) = q(i-2,4) + pi;
            q(i,6) = q(i-2,6) + pi;
        end
    end

    conf = q(IndSolLoc,:);
    disp('Q = ');
    disp(q);  
    

    %-------------------------------------------------------------
    %Q5.5:L'existence de singularités doit être prise en compte.
    q_cor = [];
    for i = 1:8
        if abs(q(i,2) - q(i,3)) == pi || abs(q(i,2) - q(i,3)) < c || abs(q(i,3) - q(i,5)) == pi || abs(q(i,3) - q(i,5)) < c || abs(q(i,2) - q(i,5)) == pi || abs(q(i,2) - q(i,5)) < c
            q_cor = [q_cor;[]];
        else
            q_cor = [q_cor;q(i,:)];
        end
    end
    conf = q_cor(IndSolLoc,:);
    disp('Q(sans singularité) = ');
    disp(q_cor);  
end

% For example :
% For p_E = [-0.5;0;0.6] and R_E = [-1 0 0;0 -1 0;0 0 1];
% 
% q =
% 
%     3.1416    1.3495         0   -0.0000    1.7921   -3.1416
%     3.1416   -0.0000   -3.1416   -1.5708    0.0000   -1.5708
%     3.1416    1.3495         0    3.1416   -1.7921         0
%     3.1416   -0.0000   -3.1416    1.5708   -0.0000    1.5708
%     6.2832   -3.1416   -0.0000    1.5708    0.0000   -1.5708
%     6.2832    1.7921    3.1416    3.1416    1.7921   -3.1416
%     6.2832   -3.1416   -0.0000    4.7124   -0.0000    1.5708
%     6.2832    1.7921    3.1416    6.2832   -1.7921         0
 




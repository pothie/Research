% x: grid of distance (vector)
% T: Total time
% ux0: initial condition [U1(0,t),U2(0,t)]
% v: speed (kTotal)
% dv: derivative of speed
% q: flow
% pce: pce values of different classes
function [U,U1,U2,tgrid] = CU3(x,T,ux0,v,dv,q,xT)
    %pt = [0.25 10;0.5 50;1 50;1.25 10];
    %CFL = 0.7;
    dx = x(2)-x(1);
    % preallocate U,U1,U2
    U = zeros(length(x),ceil(T/2));
    U1 = U;
    U2 = U;
    U1(:,1) = ux0(:,1);
    U2(:,1) = ux0(:,2);
    U(:,1) = xT(U1(:,1),U2(:,1));
    
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0%length(tgrid)<100
         
        Up1 = U1(:,tstep);
        Up2 = U2(:,tstep);
        
        [uav1,uav2,dt] = CUscheme2(Up1,Up2,dx,q,dv,v,xT);
        
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        % RK2 k1
        k1_1 = dt*(uav1);
        k1_2 = dt*(uav2);
        
        % calculating k2 = un+k1/2
        % periodic BC
%         Up1(2:end) = Up1(2:end)+k1_1;
%         Up2(2:end) = Up2(2:end)+k1_2;
%         Up1(1) = Up1(end);
%         Up2(1) = Up2(end);
        Up1(2:end-1) = Up1(2:end-1)+k1_1;
        Up2(2:end-1) = Up2(2:end-1)+k1_2;
        Up1(end) = Up1(end-1);
        Up2(end) = Up2(end-1);
        Up1(1) = Up1(2);
        Up2(1) = Up2(2);
        
        [uav1,uav2,~] = CUscheme2(Up1,Up2,dx,q,dv,v,xT);
        k2_1 = dt*(uav1);
        k2_2 = dt*(uav2);
        Up1= U1(:,tstep)*0.75+Up1*0.25;
%         Up2= U2(:,tstep)*0.75+Up2*0.25;
%         Up1(2:end)= Up1(2:end)+k2_1/4;
%         Up2(2:end)= Up2(2:end)+k2_2/4;
%         Up1(1) = Up1(end);
%         Up2(1) = Up2(end);
        Up1(2:end-1)= Up1(2:end-1)+k2_1/4;
        Up2(2:end-1)= Up2(2:end-1)+k2_2/4;
        Up1(end) = Up1(end-1);
        Up2(end) = Up2(end-1);
        Up1(1) = Up1(2);
        Up2(1) = Up2(2);

        [uav1,uav2,~] = CUscheme2(Up1,Up2,dx,q,dv,v,xT);
        k3_1 = dt*(uav1);
        k3_2 = dt*(uav2);
        Up1= U1(:,tstep)/3+Up1*2/3;
        Up2= U2(:,tstep)/3+Up2*2/3;
        Up1(2:end-1)= Up1(2:end-1)+k3_1*2/3;
        Up2(2:end-1)= Up2(2:end-1)+k3_2*2/3;
%         Up1(2:end)= Up1(2:end)+k3_1*2/3;
%         Up2(2:end)= Up2(2:end)+k3_2*2/3;
%         Up1(1) = Up1(end);
%         Up2(1) = Up2(end);

        %Boundary points
        Up1(end) = Up1(end-1);
        Up2(end) = Up2(end-1);
        Up1(1) = Up1(2);
        Up2(1) = Up2(2);

        U1(:,tstep+1) = Up1;
        U2(:,tstep+1) = Up2;

        U(:,tstep+1) = xT(U1(:,tstep+1),U2(:,tstep+1)); 
        
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
end
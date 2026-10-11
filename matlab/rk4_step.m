function next=rk4_step(x,u0,u1,h)
% RK4 with linearly interpolated endpoint inputs.
um=(u0+u1)/2;k1=fw_motion(x,u0);k2=fw_motion(x+h/2*k1,um);
k3=fw_motion(x+h/2*k2,um);k4=fw_motion(x+h*k3,u1);next=x+h/6*(k1+2*k2+2*k3+k4);
end

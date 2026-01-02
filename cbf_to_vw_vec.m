function vw = cbf_to_vw_vec(x,y,th,v,goal,obs_c,obs_D,T,u_min,u_max,gamma,epsM,W)
[v_cmd, w_cmd] = cbf_to_vw_wrapper(x,y,th,v,goal,obs_c,obs_D,T,u_min,u_max,gamma,epsM,W);
vw = [v_cmd; w_cmd];
end

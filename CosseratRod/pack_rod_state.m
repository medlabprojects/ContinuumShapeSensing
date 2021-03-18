function state = pack_rod_state(p, R, n, m)

R = orthogonalize_rotation_matrix(R);

state = [p; R(:,1); R(:,2); R(:,3); n; m];

end
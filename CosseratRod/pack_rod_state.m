function x = pack_rod_state(p, q, n, m)

if size(q) == [3,3] % If q is a rotation matrix
    q = rotm2quat(q)'; % Convert it to a quaternion
end

% Unitize the quaternion
q = q./(norm(q));

x = [p; q; n; m];
end
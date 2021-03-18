function R = orthogonalize_rotation_matrix(R0)
    % A simple, computationally cheap and stable procedure to 
    % re-orthogonalize the rotation matrix. This orthogonalization 
    % procedure is taken from Direction Cosine Matrix IMU: Theory by 
    % William Premerlani and Paul Bizard; equations 19-21.
    
    % Row vectors of slightly messed up R
    x = R0(1,:);
    y = R0(2,:);
    z = R0(3,:);
    
    % Let error=dot(x,y) where dot() is the dot product. If the matrix was 
    % orthogonal, the error would be zero. The error is spread across x and 
    % y equally: x_ort=x-(error/2)*y and y_ort=y-(error/2)*x. The third row 
    % z_ort=cross(x_ort, y_ort), which is, by definition orthogonal to 
    % x_ort and y_ort. Now, you still need to normalize x_ort, y_ort and 
    % z_ort as these vectors are supposed to be unit vectors.
    
    error = dot(x, y);
    x_ort = x - (error/2)*y;
    y_ort = y - (error/2)*x;
    z_ort = cross(x_ort, y_ort);
    
    x_new = 0.5*(3 - dot(x_ort, x_ort))*x_ort;
    y_new = 0.5*(3 - dot(y_ort, y_ort))*y_ort;
    z_new = 0.5*(3 - dot(z_ort, z_ort))*z_ort;

    R = [x_new; y_new; z_new];
end




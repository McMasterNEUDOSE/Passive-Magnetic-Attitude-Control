function Beta = calcBeta(t, q, B)
z=[0, 0, 1];                    % CHANGE THIS to reflect alignment axis 
Beta = zeros(1, length(t));

for  i = 1:length(t)
    mag         = quatrotate(q(i,:),B(i,:)./norm(B(i,:)));
    Beta(:,i)   = atan2d(norm(cross(mag,z)),dot(mag,z));
end
end

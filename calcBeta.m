function Beta = calcBeta(t, q, B)
y=[0, 1, 0];                    % CHANGE THIS to reflect alignment axis 
Beta = zeros(1, length(t));

for  i = 1:length(t)
    mag         = quatrotate(q(i,:),B(i,:)./norm(B(i,:)));
    Beta(:,i)   = atan2d(norm(cross(mag,y)),dot(mag,y));
end
end

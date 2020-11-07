function Nad = calcNad(t, q, Earth_vector)
% Calculate angle between Z axis and nadir 
z=[0, 0, 1];                    % Long axis 
Nad = zeros(1, length(t));

for  i = 1:length(t)
    unit_evec         = quatrotate(q(i,:),Earth_vector(i,:)./norm(Earth_vector(i,:)));
    Nad(:,i) = abs(atand(norm(cross(unit_evec,z))/dot(unit_evec,z)));
end
end

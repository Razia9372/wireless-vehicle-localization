function a = wrap180(a)
%WRAP180  Wrap angle(s) in degrees to (-180, 180].
%   Used for azimuth/elevation innovations before an EKF/NLS update.
    a = mod(a + 180, 360) - 180;
    a(a == -180) = 180;   % map -180 to +180 for the open-closed convention
end

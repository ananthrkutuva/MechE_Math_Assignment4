%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
%   V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
%   dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
    %unpack values from struct
    m_sun = orbit_params.m_sun;
    m_planet = orbit_params.m_planet;
    G = orbit_params.G;

    r = V(1:2);     %position of the planet [x_p; y_p]
    v = V(3:4);     %velocity of the planet [dxdt_p; dydt_p]

    %distance from the sun (sun is stuck at the origin)
    r_mag = norm(r);

    %newton's law of gravitation (eqn 17), points from the planet back
    %towards the sun so it has the minus sign
    F = -m_planet*m_sun*G*r/r_mag^3;

    %F = m*a, so divide by the planet mass to get the acceleration
    %(m_planet cancels out, the orbit doesn't depend on it)
    a = F/m_planet;

    %rate of change of position is the velocity,
    %rate of change of velocity is the acceleration
    dVdt = [v; a];
end

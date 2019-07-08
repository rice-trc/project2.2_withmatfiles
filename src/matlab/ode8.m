function Y = ode8(odefun,tspan,y0,varargin)
%ODE5  Solve differential equations with a non-adaptive method of order 5.
%   Y = ODE5(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE8(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the Dormand-Prince method of order 8 in a general 
%   framework of explicit Runge-Kutta methods.
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode5(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);

% Method coefficients -- Butcher's tableau
%  
%   C | A
%   --+---
%     | B

C = [0.0; 1 / 18; 1 / 12; 1 / 8; 5 / 16; 3 /8; 59 / 400; 93 /200; ...
    5490023248 / 9719169821; 13 /20; 1201146811 /1299019798; 1; 1];

A = [ 1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
      1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
      1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
      5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
      3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0; 
      29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
      16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
      39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
      246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
     -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
      185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
      403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0];

% RK8
B = [14005451 /335480064, 0, 0, 0, 0, - 59238493 / 1068277825, ...
    181606767 / 758867731, 561292985 / 797845732, ...
    - 1041891430 / 1371343529, 760417239 / 1151165299, ...
    118820643 / 751138087, - 528747749 / 2220607170, 1 / 4];

% RK7
%B = [13451932 / 455176623, 0, 0, 0, 0, - 808719846 / 976000145, ...
%    1757004468 / 5645159321, 656045339 / 265891186, ...
%    - 3867574721 / 1518517206, 465885868 / 322736535, ...
%    53011238 / 667516719, 2 / 45, 0];

% More convenient storage
A = A.'; 
B = B(:);      

nstages = length(B);
F = zeros(neq,nstages);

Y(:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y(:,i-1);  
  
  % General explicit Runge-Kutta framework
  F(:,1) = feval(odefun,ti,yi,varargin{:});  
  for stage = 2:nstages
    tstage = ti + C(stage-1)*hi;
    ystage = yi + F(:,1:stage-1)*(hi*A(1:stage-1,stage-1));
    F(:,stage) = feval(odefun,tstage,ystage,varargin{:});
  end  
  Y(:,i) = yi + F*(hi*B);
  
end
Y = Y.';
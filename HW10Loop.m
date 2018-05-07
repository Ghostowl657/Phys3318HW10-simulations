% Phys 3318 HW #10

% Driven Damped Oscillator (Pendulum):
% d(phi)^2/dt^2+2B*d(phi)/dt+w0^2sin(phi)=(gamma)w0^2cos(wt)
% d(phi)/dt(0) = 0, phi(0)=0

clear all % prevent weird behavior


% defining parameters
fineness = 100000;
finalT = 100;
stepsize = finalT/fineness;
peakDeviation = 0.06;
G = 1.06;
periodicity = [];
precision = 10000;
j=1;
currentPeriod = 1;
gammaNs = [];
k=1;
while G<=1.09
    j=j+1
    % setting up the diff-eq
    syms y(t) g
    eqn = diff(y,t,2)+(3*pi/2)*diff(y,t)+(3*pi)^2*sin(y) == (3*pi)^2*g*cos(2*pi*t);
    V = odeToVectorField(eqn);
    V = subs(V,g,G);
    M = matlabFunction(V,'Vars',{'t','Y'});
    
    % solving it
    tspan = linspace(0,finalT,fineness);
    y0 = [0 0];
    sol = ode45(M,tspan,y0);
    
    % plotting the solution
    yVal = deval(sol,tspan,1);
%     plot(tspan,yVal)
%     grid on
%     
%     xlabel('t')
%     ylabel('y(t)')
%     title('ODE Solution')
    
    
    tSide = 40;
    actualT = tSide/stepsize;
    [yMin,tMin] = min(yVal(actualT:length(tspan)));
    tMin = tMin+actualT-1;
    
    %plot(tspan(actualT:length(tspan)),yVal(actualT:length(tspan)))
    %findpeaks(yVal(actualT:length(tspan)))
    [pks, tpks] = findpeaks(yVal(actualT:length(tspan)));
    for i=2:length(pks)
        if pks(i)<pks(1)+peakDeviation && pks(i)>pks(1)-peakDeviation
            firstLoc = tpks(1);
            secondLoc = tpks(i);
            period = (secondLoc-firstLoc)*stepsize;
            periodicity(j) = period;
            if period==currentPeriod*2
                currentPeriod = currentPeriod*2
                %precision = precision*2;
                gammaNs(k)=G
                k=k+1;
            end
            break
        end
    end
    G=G+1/precision
end
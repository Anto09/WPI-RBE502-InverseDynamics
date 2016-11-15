%{
    HW4: function for calculating all outputs needed for homework.
    Input args: 
                l1 l2 l3 = lengths, 
                t1 t2 t3 = joint variables q
                td1 td2 td3 = joint velocities q'
                tdd1 tdd2 tdd3 = joint accelerations q''
                m1 m2 m3 = masses
                fx fy fz = external forces
                tx ty tz = external torques
                ml = load mass
                sym = bool var if any joint vars are symbolic
%}
function output = hw4(l1, l2, l3, t1, t2, t3, td1, td2, td3, tdd1, tdd2, tdd3, m1, m2, m3, fx, fy, fz, tx, ty, tz, ml, sym)

    %express everything in meters first
    if (sym == 0)
        l1 = l1/100;
        l2 = l2/100;
        l3 = l3/100;
    end
    
    %convert degrees to radians for angles if operating in numerical mode
    if (sym == 0)
        t1 = t1/180 * pi;
        t2 = t2/180 * pi;
        t3 = t3/180 * pi;
    end
    
    t = [t1; t2; t3];
    
    %same goes for angular velocities
    if (sym == 0)
        td1 = td1/180 * pi;
        td2 = td2/180 * pi;
        td3 = td3/180 * pi;
    end
    td = [td1; td2; td3];
    
    %same goes for angular accelerations
    if (sym == 0)
        tdd1 = tdd1/180 * pi;
        tdd2 = tdd2/180 * pi;
        tdd3 = tdd3/180 * pi;
    end
    tdd = [tdd1; tdd2; tdd3];
        
    %{
        Each lci is just half the length of the given links, since the
        center of mass is assumed to be halfway through the links.
    %}  
    lc1 = l1/2;
    lc2 = l2/2;
    lc3 = l3/2;

    %the local positions of the center of masses
    lc = [lc1; lc2; lc3; l3];
    l = [l1; l2; l3];
    
    %{
        Compute the position kinematics for the tip lengths
    %}
    
    x1 = l1*cos(t1);
    x2 = x1 + l2*cos(t1+t2);
    x3 = x2 + l3*cos(t1+t2+t3);
    
    y1 = l1*sin(t1);
    y2 = y1 + l2*sin(t1+t2);
    y3 = y2 + l3*sin(t1+t2+t3);
    
    'link tip positions'
    x_vec = [x1; x2; x3]
    y_vec = [y1; y2; y3]
    
    %{
        Compute the position kinematics for the center of masses
    %}
    
    xc1 = lc1*cos(t1);
    xc2 = xc1 + lc2*cos(t1+t2);
    xc3 = xc2 + lc3*cos(t1+t2+t3);
    xc4 = x3;
    
    yc1 = lc1*sin(t1);
    yc2 = yc1 + lc2*sin(t1+t2);
    yc3 = yc2 + lc3*sin(t1+t2+t3);
    yc4 = y3;

    'link point mass positions'
    xc_vec = [xc1; xc2; xc3; xc4]
    yc_vec = [yc1; yc2; yc3; yc4]
    
    %{
        Compute closed form solution. Each equation corresponds to the
        position derivative of the tip of link i with respect to its
        configuration and the configuration of the links before it
    %}
    xd1 = -l1*sin(t1)*td1;
    xd2 = xd1 - l2*sin(t1+t2)*(td1+td2);
    xd3 = xd2 - l3*sin(t1+t2+t3)*(td1+td2+td3);
    
    yd1 = l1*cos(t1)*td1;
    yd2 = yd1 + l2*cos(t1+t2)*(td1+td2);
    yd3 = yd2 + l3*cos(t1+t2+t3)*(td1+td2+td3);

    %compute 6-DOF Jacobian
    J = [-l1*sin(t1) - l2*sin(t1+t2) - l3*sin(t1+t2+t3), -l2*sin(t1+t2) - l3*sin(t1+t2+t3), -l3*sin(t1+t2+t3);
         l1*cos(t1) + l2*cos(t1+t2) + l3*cos(t1+t2+t3),   l2*cos(t1+t2) + l3*cos(t1+t2+t3),  l3*cos(t1+t2+t3);
         0,                                               0,                                 0;
         0,                                               0,                                 0;
         0,                                               0,                                 0;
         1,                                               1,                                 1];    

    %print 6-DOF Jacobian
    '3-DOF Jacobian'
    J
    
    %closed form solution
    'Closed form velocities'
    [xd3; yd3; 0; 0; 0; (td1+td2+td3)]
    
    %jacobian derived velocities
    'Jacobian velocities'
    velocities = J*[td1; td2; td3]
    
    %jacobian derived forces
    torques = J'*[fx;fy;fz;tx;ty;tz]
    
    m = [m1; m2; m3; ml];
    %Lagrange calculation
    if (m1 ~= 0 && m2 ~= 0 && m3 ~= 0)
        %{
            Each Jwi corresponds to the matrix which, if multiplied to the
            vector of joint derivatives, results in the angular velocity
            for center of mass i
        %}
        
        Jw1 = [0 0 0;
               0 0 0;
               1 0 0];
           
        Jw2 = [0 0 0;
               0 0 0;
               1 1 0];
        
        Jw3 = [0 0 0;
               0 0 0;
               1 1 1];
           
        Jw4 = Jw3;
        
        %{
            Each Jvci is the matrix which, if multiplied to the vector of
            joint derivatives, results in the position derivative of the
            center of mass i (its linear velocity)
        %}
        
        Jvc1 = [-lc1*sin(t1) 0 0;
                 lc1*cos(t1) 0 0;
                 0           0 0];
             
        Jvc2 = [-l1*sin(t1)-lc2*sin(t1+t2)  -lc2*sin(t1+t2) 0;
                 l1*cos(t1)+lc2*cos(t1+t2)   lc2*cos(t1+t2) 0;
                 0                           0              0];
             
        Jvc3 = [-l1*sin(t1)-l2*sin(t1+t2)-lc3*sin(t1+t2+t3)  -l2*sin(t1+t2)-lc3*sin(t1+t2+t3) -lc3*sin(t1+t2+t3);
                 l1*cos(t1)+l2*cos(t1+t2)+lc3*sin(t1+t2+t3)   l2*cos(t1+t2)+lc3*cos(t1+t2+t3)  lc3*cos(t1+t2+t3);
                 0                                            0                                0];
             
             
        Jvc4 = [-l1*sin(t1)-l2*sin(t1+t2)-l3*sin(t1+t2+t3)  -l2*sin(t1+t2)-l3*sin(t1+t2+t3) -l3*sin(t1+t2+t3);
                 l1*cos(t1)+l2*cos(t1+t2)+l3*sin(t1+t2+t3)   l2*cos(t1+t2)+l3*cos(t1+t2+t3)  l3*cos(t1+t2+t3);
                 0                                           0                               0];     
             
        
        I = 1/12*[m*l1*l1; m*l2*l2; m*l3*l3; 0]; %use inertia tensor of infinitely thin rod with center of mass in the middle
        
        %{
            Each Ri is just the identity matrix since all joints lie on the
            same plane, therefore we do not include them in the calculation
            for the arm dynamics.
        
            The calculation here is representative of the derivation in the
            scanned part of the pdf.
        
            D is the total kinetic energy of the system. Kl is the linear
            term of D and Ka is the angular term.
            
            P is the total potential energy. Here, each rci is actually the
            position of the center of mass.
        %}
        D = zeros(3,3);
        Kl = zeros(3,3);
        Ka = zeros(3,3);
        P = 0;
        
        %gravity 
        g = -9.821;
        
        for i=1:4
            if (i == 1)
                Kl = Kl + m(i)*(Jvc1.'*Jvc1);
                Ka = Ka + I(i)*(Jw1.'*Jw1);
                D = D + m(i)*(Jvc1.'*Jvc1) + I(i)*(Jw1.'*Jw1);
            elseif (i == 2)
                Kl = Kl + m(i)*(Jvc2.'*Jvc2);
                Ka = Ka + I(i)*(Jw2.'*Jw2);
                D = D + m(i)*(Jvc2.'*Jvc2) + I(i)*(Jw2.'*Jw2);
            elseif (i == 3)
                Kl = Kl + m(i)*(Jvc3.'*Jvc3);
                Ka = Ka + I(i)*(Jw3.'*Jw3);
                D = D + m(i)*(Jvc3.'*Jvc3) + I(i)*(Jw3.'*Jw3);
            elseif (i == 4)
                Kl = Kl + m(i)*(Jvc4.'*Jvc4);
                Ka = Ka + I(i)*(Jw4.'*Jw4);
                D = D + m(i)*(Jvc4.'*Jvc4) + I(i)*(Jw4.'*Jw4);
            end
            P = P + m(i)*g*yc_vec(i);
        end
        
        'D matrix'
        D
        'Linear terms'
        Kl
        'Angular terms'
        Ka
        'Kinetic Energy'
        0.5*td'*D*td
        'Potential Energy'
        P
        
        %{
            Use euler-lagrange to solve for torque
            lagrange_inertias represents the first summation
            lagrange_coriolis represents the second summation
            %lagrange_gravitational for potential energy terms
        %}
        lagrange_torques = [];
        lagrange_inertias = [];
        lagrange_coriolis = [];
        lagrange_gravitational = [];
        if (sym == 1)
            for k=1:3
                lagrange_gravitational = [lagrange_gravitational; diff(P,t(k))];
                ck = 0;
                dk = 0;
                for i=1:3
                    %cijk = cjik, so just go in increasing order
                    for j = 1:i
                        %Christoffel symbols method
                        cijk = 0.5 * (diff(D(k,j),t(i)) + diff(D(k,i),t(j)) - diff(D(i,j)/t(k)));
                        ck = ck + cijk*td(i)*td(j);
                    end
                    dk = dk + D(k,i)*tdd(i);
                end
                lagrange_coriolis = [lagrange_coriolis; ck];
                lagrange_inertias = [lagrange_inertias; dk];
            end
        else
            I1 = I(1);
            I2 = I(2);
            I3 = I(3);
            I4 = I(4);
            %{
                potential energy partial derivatives obtained with symbolic
                differentiation
            %}
            lagrange_gravitational = [- (9821*m3*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/200 + (l3*cos(t1 + t2 + t3))/200))/1000 - (9821*ml*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100))/1000 - (9821*m2*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/200))/1000 - (9821*l1*m1*cos(t1))/200000;
                                                                                                       - (9821*m3*((l2*cos(t1 + t2))/200 + (l3*cos(t1 + t2 + t3))/200))/1000 - (9821*ml*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/1000 - (9821*l2*m2*cos(t1 + t2))/200000;
                                                                                                                                                                                                                - (9821*l3*m3*cos(t1 + t2 + t3))/200000 - (9821*l3*ml*cos(t1 + t2 + t3))/100000];
                                                                                                                                                                                                            
            
            %{
                coriolis and centrifugal partial derivatives 
                obtained with symbolic differentiation
            %}                                                                                                                                                                                                
            lagrange_coriolis = [ (m3*(2*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) - 2*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/200)) + ml*(2*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) - 2*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/100)) - (m1*((l1*cos(t1)^2)/20000 + (l1*sin(t1)^2)/20000) + m2*((cos(t1)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/50 + (sin(t1)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/50) + m3*((cos(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50) + ml*((cos(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50))/(2*t1))*td1^2 + (- (m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/200) - ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/200) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)))/2 - (ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/100) - ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)))/2 - (m3*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/100) + ml*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100) + m2*((l2*cos(t1 + t2)*cos(t1))/20000 + (l2*sin(t1 + t2)*sin(t1))/20000))/(2*t1) - (m2*((l2*sin(t1 + t2)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/100 - (l2*cos(t1 + t2)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/100))/2 - (m3*(2*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100) - 2*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100)))/2 - (ml*(2*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100) - 2*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100)))/2)*td1*td2 + ((m3*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100))/2 + (ml*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50))/2 + (m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/200))/200 - (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/200))/200 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200))/2 + (ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100))/100 - (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/100))/100 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100))/2 - (m3*((l3*cos(t1 + t2 + t3)*cos(t1))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1))/10000))/(2*t1))*td1*td3 + (- m2*((l2*sin(t1 + t2)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/200 - (l2*cos(t1 + t2)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/200) - (m3*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/50) + ml*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/50) + m2*((l2*cos(t1 + t2)^2)/20000 + (l2*sin(t1 + t2)^2)/20000))/(2*t1) - m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/200) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)) - ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)))*td2^2 + ((m3*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200))/2 - (m3*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/10000))/(2*t1) - (m3*((l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/200))/200 - (l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200))/2 + (ml*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100))/2 - (ml*((l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/100))/100 - (l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100))/2)*td2*td3 + (ml*((l3^2*cos(t1 + t2 + t3)^2)/10000 + (l3^2*cos(t1 + t2 + t3)*sin(t1 + t2 + t3))/10000 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100) - (ml*((l3*cos(t1 + t2 + t3)^2)/5000 + (l3*sin(t1 + t2 + t3)^2)/5000) + m3*((l3*cos(t1 + t2 + t3)^2)/20000 + (l3*sin(t1 + t2 + t3)^2)/20000))/(2*t1) + m3*((l3^2*cos(t1 + t2 + t3)^2)/40000 + (l3^2*cos(t1 + t2 + t3)*sin(t1 + t2 + t3))/40000 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200))*td3^2;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (- m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/200) - ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/200) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)) - ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/100) - ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)) - (m1*((l1*cos(t1)^2)/20000 + (l1*sin(t1)^2)/20000) + m2*((cos(t1)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/50 + (sin(t1)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/50) + m3*((cos(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50) + ml*((cos(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50))/(2*t2))*td1^2 + ((m3*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200))/2 - (m3*((l3*cos(t1 + t2 + t3)*cos(t1))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1))/10000))/(2*t2) + (ml*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100))/2)*td1*td3 - td2*((m3*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/100) + ml*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100) + m2*((l2*cos(t1 + t2)*cos(t1))/20000 + (l2*sin(t1 + t2)*sin(t1))/20000))/(2*t2) + (m2*((l2*sin(t1 + t2)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/200 - (l2*cos(t1 + t2)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/200))/2 + (m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/200) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)))/2 + (ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100) - ((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100) + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)))/2)*td1 + (- (ml*((l3*cos(t1 + t2 + t3)^2)/5000 + (l3*sin(t1 + t2 + t3)^2)/5000) + m3*((l3*cos(t1 + t2 + t3)^2)/20000 + (l3*sin(t1 + t2 + t3)^2)/20000))/(2*t2) - m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200) - ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100))*td3^2 - td2*((m3*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/10000))/(2*t2) + (m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/100))/2 + (ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/50 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/50))/2)*td3 - (td2^2*(m3*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/50) + ml*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/50) + m2*((l2*cos(t1 + t2)^2)/20000 + (l2*sin(t1 + t2)^2)/20000)))/(2*t2);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           td1^2*(m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/200))/200 - (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/200))/200 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200) + ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100))/100 - (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/100))/100 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100) - (m1*((l1*cos(t1)^2)/20000 + (l1*sin(t1)^2)/20000) + m2*((cos(t1)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/50 + (sin(t1)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/50) + m3*((cos(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50) + ml*((cos(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50))/(2*t3)) - (td3^2*(ml*((l3*cos(t1 + t2 + t3)^2)/5000 + (l3*sin(t1 + t2 + t3)^2)/5000) + m3*((l3*cos(t1 + t2 + t3)^2)/20000 + (l3*sin(t1 + t2 + t3)^2)/20000)))/(2*t3) - (td2^2*(m3*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/50) + ml*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/50) + m2*((l2*cos(t1 + t2)^2)/20000 + (l2*sin(t1 + t2)^2)/20000)))/(2*t3) - td1*td2*((m3*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/100) + ml*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100) + m2*((l2*cos(t1 + t2)*cos(t1))/20000 + (l2*sin(t1 + t2)*sin(t1))/20000))/(2*t3) + (m3*((l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/200))/200 - (l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200))/2 + (ml*((l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/100))/100 - (l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100))/2) + td1*td3*((ml*((l3^2*cos(t1 + t2 + t3)^2)/10000 + (l3^2*cos(t1 + t2 + t3)*sin(t1 + t2 + t3))/10000 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100))/2 - (m3*((l3*cos(t1 + t2 + t3)*cos(t1))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1))/10000))/(2*t3) + (m3*((l3^2*cos(t1 + t2 + t3)^2)/40000 + (l3^2*cos(t1 + t2 + t3)*sin(t1 + t2 + t3))/40000 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200))/2) - td2*td3*((m3*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/10000))/(2*t3) + (m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200))/2 + (ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100))/2)];
 
                                                                                                                                                  
                                                                                                                                                                                                            
            %{
                inertial partial derivatives obtained with symbolic 
                differentiation
            %}                                                                                                                                                                                                  
            lagrange_inertias = [tdd2*(I2 + I3 + I4 + m2*((l2*cos(t1 + t2)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/200 + (l2*sin(t1 + t2)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/200) + m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100) + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100)) + ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100) + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))) + tdd3*(I3 + I4 + m3*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200) + ml*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100)) + tdd1*(I1 + I2 + I3 + I4 + m1*((l1^2*cos(t1)^2)/40000 + (l1^2*sin(t1)^2)/40000) + m2*(((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100)^2 + ((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100)^2) + m3*(((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100)^2 + ((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)^2) + ml*(((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100)^2 + ((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100)^2));
                                                                                                                                                                                                                                                                                tdd2*(I2 + I3 + I4 + m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)^2 + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)^2) + ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)^2 + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)^2) + m2*((l2^2*cos(t1 + t2)^2)/40000 + (l2^2*sin(t1 + t2)^2)/40000)) + tdd1*(I2 + I3 + I4 + m2*((l2*cos(t1 + t2)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/200 + (l2*sin(t1 + t2)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/200) + m3*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100) + ((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100)) + ml*(((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100) + ((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))) + tdd3*(I3 + I4 + m3*((l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200) + ml*((l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      tdd3*(I3 + I4 + ml*((l3^2*cos(t1 + t2 + t3)^2)/10000 + (l3^2*sin(t1 + t2 + t3)^2)/10000) + m3*((l3^2*cos(t1 + t2 + t3)^2)/40000 + (l3^2*sin(t1 + t2 + t3)^2)/40000)) + tdd1*(I3 + I4 + m3*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200) + ml*((l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100)) + tdd2*(I3 + I4 + m3*((l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200) + ml*((l3*cos(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100))];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/200))/200 - (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/200))/200 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200) + ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l1*cos(t1))/100 + (l3*cos(t1 + t2 + t3))/100))/100 - (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 + (l1*sin(t1))/100 - (l3*cos(t1 + t2 + t3))/100))/100 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100) - (m3*((l3*cos(t1 + t2 + t3)*cos(t1))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1))/10000))/t3 - (ml*((l3*cos(t1 + t2 + t3)^2)/5000 + (l3*sin(t1 + t2 + t3)^2)/5000) + m3*((l3*cos(t1 + t2 + t3)^2)/20000 + (l3*sin(t1 + t2 + t3)^2)/20000))/(2*t3) - (m3*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/100) + ml*((cos(t1)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100) + m2*((l2*cos(t1 + t2)*cos(t1))/20000 + (l2*sin(t1 + t2)*sin(t1))/20000))/t3 - (m1*((l1*cos(t1)^2)/20000 + (l1*sin(t1)^2)/20000) + m2*((cos(t1)*((l2*cos(t1 + t2))/200 + (l1*cos(t1))/100))/50 + (sin(t1)*((l2*sin(t1 + t2))/200 + (l1*sin(t1))/100))/50) + m3*((cos(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50) + ml*((cos(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/50 + (sin(t1)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/50))/(2*t3) + ml*((l3^2*cos(t1 + t2 + t3)^2)/10000 + (l3^2*cos(t1 + t2 + t3)*sin(t1 + t2 + t3))/10000 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100) + m3*((l3^2*cos(t1 + t2 + t3)^2)/40000 + (l3^2*cos(t1 + t2 + t3)*sin(t1 + t2 + t3))/40000 - (l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200) - (m3*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/20000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/20000) + ml*((l3*cos(t1 + t2 + t3)*cos(t1 + t2))/10000 + (l3*sin(t1 + t2 + t3)*sin(t1 + t2))/10000))/t3 - (m3*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/50) + ml*((sin(t1 + t2)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/50 + (cos(t1 + t2)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/50) + m2*((l2*cos(t1 + t2)^2)/20000 + (l2*sin(t1 + t2)^2)/20000))/(2*t3) - m3*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100))/200) - ml*((l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100))/100) - m3*((l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/200 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/200 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/200 + (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/200))/200 - (l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/200))/200) - ml*((l3*sin(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*cos(t1 + t2))/100 + (l1*cos(t1))/100))/100 - (l3*cos(t1 + t2 + t3)*((l3*sin(t1 + t2 + t3))/100 + (l2*sin(t1 + t2))/100 + (l1*sin(t1))/100))/100 + (l3*cos(t1 + t2 + t3)*((l2*sin(t1 + t2))/100 - (l3*cos(t1 + t2 + t3))/100))/100 - (l3*sin(t1 + t2 + t3)*((l2*cos(t1 + t2))/100 + (l3*cos(t1 + t2 + t3))/100))/100);
        end
        lagrange_torques = lagrange_inertias + lagrange_coriolis + lagrange_gravitational;
        'M(q)qdd'
        lagrange_inertias
        'C(q,qd)'
        lagrange_coriolis
        'g(q)'
        lagrange_gravitational 
        'M(q)qdd + C(q,qd) + g(q)'
        lagrange_torques                  
        
    end
    
    %NE calculation
    if (m1 ~= 0 && m2 ~= 0 && m3 ~= 0)
        %set the rotation matrices
        r01 = rotation(t1);
        r12 = rotation(t2);
        r23 = rotation(t3);  
        r02 = rotation(t1+t2);
        r03 = rotation(t1+t2+t3);
        %alternatively, but producing longer symbolic matrix terms
        %r02 = r01*r12;
        %r03 = r02*r23;
        
        %set the vectors from and to the center of masses
        i = [1; 0; 0];
        r = [l1*i l2*i l3*i];
        rc = [lc1*i lc2*i lc3*i];
        rce = [(l1-lc1)*i (l2-lc2)*i (l3-lc3)*i];
        
        %set the angular velocities and acceleraitons; no need to recurse
        k = [0; 0; 1];
        %start at 0, end at 3
        w = [zeros(3,1)];
        wa = [zeros(3,1)];
        wm = [zeros(3,1)];
        
        %initialize linear components
        ac = [zeros(3,1)];
        ae = [zeros(3,1)];
        
        %know components for final link
        ne_f = [];
        ne_t = [];
        
        %gravitational vectors
        j = [0; 1; 0];
        gv = [-r01*g*j -r02*g*j -r03*g*j];
        
        %{
            Forward recursion. Each ac,i in class is actually ac,i+1 in the
            array due to there being an ac,0. Same goes for ae,i. Each rc,i
            and re,i is actually rc,i-1 and re,i-1 due to the lack a fourth
            term in the matrices
        %}
        for i=2:4
            R = zeros(3,3);
            R0 = zeros(3,3);
            if (i == 2)
                R = r01;
                R0 = r01;
            elseif (i == 3)
                R = r12;
                R0 = r02;
            elseif (i==4)
                R = r23;
                R0 = r03;
            end
            wi = R.'*w(:,i-1) + td(i-1)*R0.'*k;
            wai = R.'*(wa(:,i-1) + R0.'*k*tdd(i-1) + cross(wi, R0.'*k*td(i-1)));
            
            w = [w wi];
            wa = [wa wai];
            
            
            aei = R.'*ae(:,i-1) + cross(wa(:,i), r(:,i-1)) + cross(w(:,i), cross(w(:,i), r(:,i-1)));
            aci = R.'*ac(:,i-1) + cross(wa(:,i), rc(:,i-1)) + cross(w(:,i), cross(w(:,i), rc(:,i-1)));
            
            ac = [ac aci];
            ae = [ae aei];
        end
        gv
        
        %to access elements of force and torque matrices
        index = 1;
        for i = 4:-1:1
            if i == 4 %for the point mass at the end, add an external force
                %ne_fi = m(i)*ae(:,i) - m(i)*gv(:,i-1);
                ne_fi = - m(i)*[0;g;0];
                ne_f = [ne_f ne_fi];
                ne_t = [ne_t zeros(3,1)];
            else
                R = zeros(3,3);
                if (i == 1)
                    R = r12;
                elseif (i == 2)
                    R = r23;
                elseif (i==3)
                    R = eye(3);
                end
                
                ne_fi = m(i)*ac(:,i+1) - m(i)*gv(:,i) + R*ne_f(:,index-1);
                ne_ti = R*ne_t(:,index-1) - cross(ne_fi, rc(:,i)) + cross(R*ne_f(:,index-1), rce(:,i)) + wa(:,i) + cross(w(:,i), I(i)*w(:,i));
                
                ne_f = [ne_f ne_fi];
                ne_t = [ne_t ne_ti];
            end
            index = index + 1;
        end
        %flip them because they were ordered in reverse to allow for better
        %addition of elements
        ne_f = fliplr(ne_f);
        ne_t = fliplr(ne_t);
        
        %print
        'Final forces'
         ne_f(:,1:3)
        'Final torques'
         ne_t(:,1:3)
    end
end
% Suboptimal Quaternion of Accelerometer-Magnetometer combination (OQAM)
% Proposed by Jin Wu
% e-mail: jin_wu_uestc@hotmail.com


function q = SQAM(Ab, Mb, Mr)

    ax = Ab(1);       ay = Ab(2);       az = Ab(3);
    mx = Mb(1);       my = Mb(2);       mz = Mb(3);

    mN = Mr(1);
    mD = Mr(3);
    
    alpha = ax * mx + ay * my + az * mz;
    beta = sqrt(1 - alpha * alpha);
    t2 = sqrt(2) * sqrt(1 + alpha * mD + beta * mN);

    q = [
             ay * (beta + mN + t2 * mx) + (beta * mD - alpha * mN - ax * t2) * my;
    
           - (ax * mD) + beta * t2 + mx + beta * mN * mx - ...
             az * (beta + mN + t2 * mx) - beta * mD * mz + ...
             ax * t2 * mz - alpha * (ax - mD * mx - mN * mz);
   
    
           - (alpha * ay) - ay * mD + my + alpha * mD * my + beta * mN * my - ...
             az * t2 * my + ay * t2 * mz;
             
           - (az * mD) + ax * (beta + mN) + beta * mD * mx + mz + beta * mN * mz - ...
             alpha * (az + mN * mx - mD * mz)
    ];

    q = q ./ norm(q);
end
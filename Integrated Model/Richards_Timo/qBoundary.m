function qTop = qBoundary(t)    
    
    if t<25,
        qTop = 0;
    elseif t>=25 && t<200
        qTop = -0.002;
    elseif t>=200
        qTop = 0;
    end
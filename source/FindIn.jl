# this function finds the indium content fraction given fixed aluminum to reach a target bandgap
function FindIn(TargetBandgap, AlFraction)
    x = AlFraction
    a = 0.0
    b = 1.0
    NMax = 20
    n = 1
    cprev = AlFraction
    c = (a+b)*0.5
    while abs(c-cprev) > 0.000001
        y = a
        z = 1 - AlFraction - a
        f_a =  TargetBandgap - (x*EgAlN + y*Eg + z*EgInN - bAlInN*x*(1-x) - bInGaN*y*(1-y) - bAlGaN*x*y + (bInGaN+bAlInN)*x*y - Beg*x*y*z)
        y = c
        z = 1 - AlFraction - c
        f_c =  TargetBandgap - (x*EgAlN + y*Eg + z*EgInN - bAlInN*x*(1-x) - bInGaN*y*(1-y) - bAlGaN*x*y + (bInGaN+bAlInN)*x*y - Beg*x*y*z)
        if f_a*f_c<0
            b = c
        else
            a = c
        end
        n = n+1
        cprev = c
        c = (a+b)*0.5
    end
    return 1-c-x
end

function FindAl(TargetBandgap)
    a = 0.0
    b = 1
    NMax = 20
    n = 1
    cprev = 1
    c = (a+b)*0.5
    while abs(c-cprev) > 0.0001
        x = a
        y = 1-a
        f_a =  TargetBandgap - (x*EgAlN + y*Eg - bAlGaN*x*y)
        x = c
        y = 1-c
        f_c =  TargetBandgap - (x*EgAlN + y*Eg - bAlGaN*x*y)
        if f_a*f_c<0
            b = c
        else
            a = c
        end
        n = n+1
        cprev = c
        c = (a+b)*0.5
    end
    return c
end
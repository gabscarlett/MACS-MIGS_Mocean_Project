function [Fex] = waveChanDamping(Comp, CompP, cdPlate)

% damping matrix for plate
Dplate = zeros(3, 3);

% apply damping to the heave element
Dplate(1,1) = cdPlate;
 
CompP.SetDpar(Dplate); % set the damping

xi = CompP.Motions; % plate motions

omega = 2*pi./Comp.T; % wave frequency

Fex0 = Comp.Fex; % wave channel forced without damping plate

% pre allocation
Xi = zeros(size(Fex0));
Frad = zeros(size(Fex0));

inds = 7:9; % indices of the plate heave, roll and pitch
Xi(:,:,inds) = xi;

for n = 1:length(Comp.T)
    
    a_ = zeros(Comp.DoF);
    b_ = zeros(Comp.DoF);
    
    for p = 1:Comp.DoF
        for q = 1:Comp.DoF
            a_(p,q) = Comp.A(n,p,q);
            b_(p,q) = Comp.B(n,p,q);
        end
    end
    
    for m = 1:1
        Frad(n,m,:) = (-omega(n)^2.*a_ + 1i*omega(n)*b_)*squeeze(Xi(n,m,:));
    end
end

Fex = Fex0 - Frad;
Fex = squeeze(Fex);

end
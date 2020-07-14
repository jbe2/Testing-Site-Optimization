clear all
sz = 5; #size of mesh
beta = 0.6; #infection coeff.
gama = 0.2;  #recovery coeff.
tau = 0.2   #movement b/t nodes coeff.
I = zeros(sz);
I(1,1) = 1;   #initial condition for infection
S = floor(100.*(rand(sz)+1));  #initial suceptible people
N = I+S;  # population across grid

#For plotting
N_tot = sum(sum(N));  
I_tot = sum(sum(I));
S_tot = sum(sum(S));

dt = 1;
final_time = 50;
numsteps = final_time/dt


I_neibs = zeros(sz);
for t = 1:numsteps
    #Creating the Matrix for infected neighbors
    I_neibs(1,1) = I(1,2) + I(2,1);
    I_neibs(1,sz) = I(1,sz-1) + I(2,sz);
    I_neibs(sz,1) = I(sz,2) + I(sz-1,1);
    I_neibs(sz,sz) = I(sz-1,sz) + I(sz,sz-1);
    I_neibs(1,2:sz-1) = I(1,1:sz-2)+I(1,3:sz) + I(2,2:sz-1); 
    I_neibs(sz,2:sz-1) = I(sz,1:sz-2)+I(sz,3:sz) + I(sz-1,2:sz-1);
    I_neibs(2:sz-1,1) = I(1:sz-2,1)+I(3:sz,1) + I(2:sz-1,2);
    I_neibs(2:sz-1,sz) = I(1:sz-2,sz)+I(3:sz,sz) + I(2:sz-1,sz-1);
    for i = 2:sz-1
        for j = 2:sz-1
            I_neibs(i,j) = I(i+1,j) + I(i-1,j) + I(i,j+1) + I(i,j-1);
        end
    end
    
    #Solving the ODE
    dsdt = -beta.*(I.*S./N) + gama.*I - tau.*I_neibs;
    didt = beta.*(I.*S./N) - gama.*I +  tau.*I_neibs;
    I = I + didt.*dt
    S = N-I;
    
    #Make sure the system is clsoed

    
    #For plotting
    I_tot = [I_tot, sum(sum(I))];
    S_tot = [S_tot, sum(sum(S))];
end

plot(I_tot)
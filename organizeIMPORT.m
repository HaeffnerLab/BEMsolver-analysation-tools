function data = organizeIMPORT(elstruct)

data.X = elstruct.X;
data.Y = elstruct.Y;
data.Z = elstruct.Z;
data.grid = elstruct.grid;
data.date_imported = elstruct.date_imported;

data.Vrf = elstruct.EL_phi0 - elstruct.EL_phi1 ; #takes care of out of phase drive
data.W1 = elstruct.EL_phi2;

%for compatibility:
NUM_AXIS = size(data.X,1);
data.N1 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N2 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N3 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N4 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N5 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N6 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N7 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N8 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N9 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.N10 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W2 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W3 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W4 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W5 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W6 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W7 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W8 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W9 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.W10 = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
data.Vc = zeros(NUM_AXIS ,NUM_AXIS ,NUM_AXIS );
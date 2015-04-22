function G4x = G4coeff(n)
%G4COEFF  Evaluate coefficients for C_4 for great ellipse
%
%   G4x = G4COEFF(n) evaluates the coefficients of epsilon^l in expansion
%   of the greate ellipse area (expressed in terms of n and epsi).  n is a
%   scalar.  G4x is a 1 x 21 array.

  nG4 = 6;
  nG4x = (nG4 * (nG4 + 1)) / 2;
  G4x = zeros(1, nG4x);
  G4x(0+1) = (n*(n*(n*(n*(200*n+416)+1144)+6864)+21021)+15015)/90090;
  G4x(1+1) = (n*(n*((-117944*n-110552)*n-84227)-41184)-9009)/120120;
  G4x(2+1) = (n*(n*(6417449*n+3013374)+1012583)+172458)/720720;
  G4x(3+1) = ((-135037988*n-32774196)*n-4232371)/5765760;
  G4x(4+1) = (138833443*n+13938873)/5765760;
  G4x(5+1) = -13200233/1537536;
  G4x(6+1) = (n*(n*(n*(117944*n+110552)+84227)+41184)+9009)/1081080;
  G4x(7+1) = (n*((-5975241*n-2676466)*n-847847)-136422)/4324320;
  G4x(8+1) = (n*(71379996*n+16424252)+1987557)/17297280;
  G4x(9+1) = (-39452953*n-3753828)/8648640;
  G4x(10+1) = 2625577/1537536;
  G4x(11+1) = (n*(n*(203633*n+80106)+20735)+2574)/1441440;
  G4x(12+1) = ((-3634676*n-741988)*n-76219)/5765760;
  G4x(13+1) = (2443153*n+208182)/2882880;
  G4x(14+1) = -5512967/15375360;
  G4x(15+1) = (n*(48020*n+8372)+715)/1153152;
  G4x(16+1) = (-71477*n-5317)/768768;
  G4x(17+1) = 22397/439296;
  G4x(18+1) = (1407*n+91)/329472;
  G4x(19+1) = -5453/1317888;
  G4x(20+1) = 21/146432;
end

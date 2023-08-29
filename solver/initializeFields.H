
#include "clock.H"
Random randNumber(clock::getTime());


Info<< "Initialising fields... " << endl;
double xMax = max(mesh.C().component(vector::X)).value();
double yMax = max(mesh.C().component(vector::Y)).value();
double zMax = max(mesh.C().component(vector::Z)).value();
double xCenter, yCenter, zCenter, randTheta[3], randRadius, maxvol[int(phases)-1], curvol[int(phases)-1];
//int NS;
if (phases == 2)
{
//! For 2D simulations
if (dimensions == 2)
{
//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01*vector(0,0,1);
    cEigenStrain1 = 0*cEigenStrain1;
    //NS = 1;
    
    //dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    //Info << "Probability of nucleation: " << prob.value() << endl;
    
//for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/2);
    Info<< "random theta: " << randTheta[2] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

    //if (prob.value() > randNumber.globalScalar01()){
    theta = theta + randTheta[2]*gaussianSeed*vector(0,0,1);
    phi_1 = phi_1 + gaussianSeed;
    //}
//}
}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01*vector(0,0,1);
    cEigenStrain1 = 0*cEigenStrain1;
    
for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/2);
    Info<< "random theta: " << randTheta[2] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

    theta = theta + randTheta[2]*gaussianSeed*vector(0,0,1);
    phi_1 = phi_1 + gaussianSeed;
}
}

//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;
    
    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
maxvol[0] = volumeFraction[0][0]*xMax*yMax;
curvol[0] = 0;

if (numSeeds != 1){
while(curvol[0] <= maxvol[0]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    curvol[0] = curvol[0] + pi.value()*randRadius*randRadius;
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(randRadius*randRadius));
    
    //volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi_1 = phi_1 + gaussianSeed;
    //}
}
}

if (numSeeds == 1){
xCenter = centerX[0][0];
yCenter = centerY[0][0];
phi_1 = phi_1 + (1-phi_1)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) - seedRadius[0][0]*seedRadius[0][0]));//(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = centerX[0][0];
yCenter = centerY[0][0];
phi_1 = phi_1 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

theta = theta_z.value()*vector(0,0,1);
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;
}

//! For 3D simulations
if (dimensions == 3)
{

//! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    theta = 0.01*vector(1,1,1);
    cEigenStrain1 = 0*cEigenStrain1;
    //NS = 1;

    /*dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;

for(int i=0; i<numSeeds; i++){*/
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetaz: " << randTheta[2] << endl;
    randTheta[0] = randNumber.globalScalar01()*(pi.value()/4);
    randTheta[1] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetax: " << randTheta[0] << endl;
    Info<< "random thetay: " << randTheta[1] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

            //if (prob.value() > randNumber.globalScalar01()){
     theta = theta + gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
     phi_1 = phi_1 + gaussianSeed;

}

//! Initial conditions for multigrain simulation
if ((swcool == 0)&&(swch == 1))
{
    theta = 0.01*vector(1,1,1);
    cEigenStrain1 = 0*cEigenStrain1;

for(int i=0; i<numSeeds; i++){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    randTheta[2] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetaz: " << randTheta[2] << endl;
    randTheta[0] = randNumber.globalScalar01()*(pi.value()/4);
    randTheta[1] = randNumber.globalScalar01()*(pi.value()/4);
    Info<< "random thetax: " << randTheta[0] << endl;
    Info<< "random thetay: " << randTheta[1] << endl;

    Info<< "Filling phi and theta fields in seeds" << endl;
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

    theta = theta + gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
    phi_1 = phi_1 + gaussianSeed;
}
}

//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
maxvol[0] = volumeFraction[0][0]*xMax*yMax;
curvol[0] = 0;

if (numSeeds != 1){
while(curvol[0] <= maxvol[0]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    randRadius = (1 + spread[0][0]*(2*randNumber.globalScalar01() - 1))*seedRadius[0][0];
    curvol[0] = curvol[0] + pi.value()*(4/3)*randRadius*randRadius;
    volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(randRadius*randRadius));

    //volScalarField gaussianSeed = (1-phi)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter) - randRadius*randRadius));

    //if (prob.value() > randNumber.globalScalar01()){
    //theta = theta + randTheta*gaussianSeed;
    phi_1 = phi_1 + gaussianSeed;
    //}
}
}

if (numSeeds == 1){
xCenter = centerX[0][0];
yCenter = centerY[0][0];
zCenter = centerZ[0][0];
phi_1 = phi_1 + (1-phi_1)*0.5*(1 - sign((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter) - seedRadius[0][0]*seedRadius[0][0]));//(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));
}
}

//! Initial conditions for single dendrite solidification simulation
if (swch == 0)
{
xCenter = centerX[0][0];
yCenter = centerY[0][0];
zCenter = centerZ[0][0];
phi_1 = phi_1 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

theta = theta_x.value()*vector(1,0,0) + theta_y.value()*vector(0,1,0) + theta_z.value()*vector(0,0,1);
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
}

if (phases == 3)
{
if (dimensions == 2)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax;
curvol[i_phase] = 0;
}

if (numSeeds != 1){
while((curvol[0]+curvol[1]) <= (maxvol[0]+maxvol[1])){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    if (curvol[0] <= maxvol[0]){
    curvol[0] = curvol[0] + pi.value()*seedRadius[0][0]*seedRadius[0][0];
    
phi_1 = phi_1 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));
}

if (curvol[1] <= maxvol[1]){
    curvol[1] = curvol[1] + pi.value()*seedRadius[1][0]*seedRadius[1][0];

phi_2 = phi_2 + (1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));
}

}
/*
while(curvol[1] <= maxvol[1]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[1] = curvol[1] + pi.value()*seedRadius[1][0]*seedRadius[1][0];

phi_2 = phi_2 + (1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));

}
*/
phi_3 = 1 - (phi_1 + phi_2);

theta = theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;
}

if (dimensions == 3)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax;
curvol[i_phase] = 0;
}

if (numSeeds != 1){
while(curvol[0] <= maxvol[0]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[0] = curvol[0] + pi.value()*(4/3)*seedRadius[0][0]*seedRadius[0][0];
    
phi_1 = phi_1 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

}

while(curvol[1] <= maxvol[1]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[1] = curvol[1] + pi.value()*(4/3)*seedRadius[1][0]*seedRadius[1][0];

phi_2 = phi_2 + (1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[1][0]*seedRadius[1][0]));

}

phi_3 = 1 - (phi_1 + phi_2);

theta = theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
}

if (phases == 4)
{
if (dimensions == 2)
{
//! Initial conditions for precipitate growth simulation
if (swch == 1)
{
    cEigenStrain1 = 0*cEigenStrain1;
    cEigenStrain2 = 0*cEigenStrain2;
    cEigenStrain3 = 0*cEigenStrain3;
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax;
curvol[i_phase] = 0;
}

if (numSeeds != 1){
while((curvol[0]+curvol[1]+curvol[2]) <= (maxvol[0]+maxvol[1]+maxvol[2])){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    if (curvol[0] <= maxvol[0]){
    curvol[0] = curvol[0] + pi.value()*seedRadius[0][0]*seedRadius[0][0];
    
phi_1 = phi_1 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));
}

else if (curvol[1] <= maxvol[1]){
    curvol[1] = curvol[1] + pi.value()*seedRadius[1][0]*seedRadius[1][0];

phi_2 = phi_2 + (1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));
}

else if (curvol[2] <= maxvol[2]){
    curvol[2] = curvol[2] + pi.value()*seedRadius[2][0]*seedRadius[2][0];

phi_3 = phi_3 + (1-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[2][0]*seedRadius[2][0]));
}

//phi_2 = phi_2*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

//phi_3 = phi_3*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

//phi_4 = phi_4*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));
}

//phi_4 = 1 - phi_1;
/*
while(curvol[1] <= maxvol[1]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[1] = curvol[1] + pi.value()*seedRadius[1][0]*seedRadius[1][0];

//phi_1 = phi_1*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));

phi_2 = phi_2 + (1-phi_1-phi_2)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));

//phi_3 = phi_3*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));

//phi_4 = phi_4*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0]));
}

//phi_4 = 1 - (phi_1 + phi_2);

while(curvol[2] <= maxvol[2]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[2] = curvol[2] + pi.value()*seedRadius[2][0]*seedRadius[2][0];

//phi_1 = phi_1*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[2][0]*seedRadius[2][0]));

//phi_2 = phi_2*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[2][0]*seedRadius[2][0]));

phi_3 = phi_3 + (1-phi_1-phi_2-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[2][0]*seedRadius[2][0]));

//phi_4 = phi_4*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[2][0]*seedRadius[2][0]));
}

phi_4 = 1 - phi_1*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0])) - phi_2*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[1][0]*seedRadius[1][0])) - phi_3*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[2][0]*seedRadius[2][0]));
*/
phi_4 = 1 - (phi_1 + phi_2 + phi_3);

phase_index = phi_1 + 2*phi_2 + 3*phi_3 + 4*phi_4;

theta = theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max phi_4:" << min(phi_4).value() << ' ' << max(phi_4).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;
}

if (dimensions == 3)
{
//! Initial conditions for precipitate growth simulation
if (swch == 2)
{
    /*theta = 0.01;
    //NS = 1;

    dimensionedScalar prob = Foam::exp(-0.5*(log(T0-initial)/0.5)*(log(T0-initial)/0.5))/(0.5*Foam::sqrt(2*pi)*(T0-initial));
    Info << "Probability of nucleation: " << prob.value() << endl;*/
for (int i_phase = 0; i_phase < (phases-1); i_phase++){   
maxvol[i_phase] = volumeFraction[i_phase][0]*xMax*yMax;
curvol[i_phase] = 0;
}

if (numSeeds != 1){
while(curvol[0] <= maxvol[0]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[0] = curvol[0] + pi.value()*seedRadius[0][0]*seedRadius[0][0];
    
phi_1 = phi_1 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));

}

while(curvol[1] <= maxvol[1]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[1] = curvol[1] + pi.value()*seedRadius[1][0]*seedRadius[1][0];

phi_2 = phi_2 + (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[1][0]*seedRadius[1][0]));

}

while(curvol[2] <= maxvol[2]){
    Info<< "creating a new seed" << endl;
    xCenter = randNumber.globalScalar01()*xMax;
    yCenter = randNumber.globalScalar01()*yMax;
    zCenter = randNumber.globalScalar01()*zMax;
    //xCenter = 0.5*xMax;
    //yCenter = 0.5*yMax;
    Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
    //randTheta = randNumber.globalScalar01()*(pi.value()/2);
    //Info<< "random theta: " << randTheta << endl;

    Info<< "Filling phi field in seeds" << endl;
    //randRadius = (1 + spread*(2*randNumber.globalScalar01() - 1))*seedRadius[0];
    curvol[2] = curvol[2] + pi.value()*seedRadius[2][0]*seedRadius[2][0];

phi_3 = phi_3 + (1-phi_3)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[2][0]*seedRadius[2][0]));

}

phi_4 = 1 - (phi_1 + phi_2 + phi_3);

theta = theta_x.value()*vector(1,0,0) + theta_y.value()*vector(0,1,0) + theta_z.value()*vector(0,0,1);
}
}

Info<< "Initialised fields." << endl;

Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
Info<< "Min/max phi_4:" << min(phi_4).value() << ' ' << max(phi_4).value() << endl;
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
}


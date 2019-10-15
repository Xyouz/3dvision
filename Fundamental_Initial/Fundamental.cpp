// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

// Marius Dufraisse MVA

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

// Select 8 points and return the corresponding fundamental matrix
FMatrix<float,3,3> randF(vector<Match>& matches){
    FMatrix<float,9,9> A(-1.);
    int n_match = matches.size();
    float x1, y1, x2, y2;
    for (int i=0; i<8; i++){
        Match match = matches[rand()%n_match];
        // Normalization
        x1 = 0.001 * match.x1; y1 = 0.001 * match.y1;
        x2 = 0.001 * match.x2; y2 = 0.001 * match.y2;
        // Building the linear system
        A(i,0) = x1 * x2;
        A(i,1) = x1 * y2;
        A(i,2) = x1;
        
        A(i,3) = y1 * x2;
        A(i,4) = y1 * y2;
        A(i,5) = y1;
        
        A(i,6) = x2;
        A(i,7) = y2;
        A(i,8) = 1;
        
        
        A(8,i) = 0;
    }
    A(8,8) = 0;

    // Solve linear system using svd
    FVector<float,9> S;
    FMatrix<float,9,9> U, V;
    svd(A,U,S,V);
    // Extracting F out of V
    FMatrix<float,3,3> F;
    for (int l = 0; l < 3; l++){
        for (int k = 0; k < 3; k++){
            F(l,k) = V(8,3*l+k);
        }
    }
    // Enforcing det(F) = 0
    FVector<float,3> Sf;
    FMatrix<float,3,3> Uf, Vf;
    svd(F,Uf,Sf,Vf);
    Sf[2] = 0;
    F = Uf * Diagonal(Sf) * Vf;

    // Normalization
    FMatrix<float,3,3> N(0.0);
    N(0,0) = 0.001; N(1,1) = 0.001; N(2,2) = 1;
    F = N * F * N;
    return F;
}

// Given a fundamental matrix F and a vector of matches, return
// the vector containing the indices of matches that are inliers for F
vector<int> findInliers(FMatrix<float,3,3>& F, vector<Match>& matches, float distMax){
    vector<int> inliers;
    distMax = distMax;
    float x1, x2, y1, y2;
    for (int i = 0; i < matches.size(); i++){
        x1 = matches[i].x1; y1 = matches[i].y1;
        x2 = matches[i].x2; y2 = matches[i].y2;
        FloatPoint3 point, point1;
        point[0]=x1;
        point[1]=y1;
        point[2]=1;

        point1[0]=x2;
        point1[1]=y2;
        point1[2]=1;

        // Compute distance between point1  and the epipolar line of point
        point = transpose(F) * point;
        float dist = abs(point1*point);
        dist = dist/sqrt(point[0]*point[0]+point[1]*point[1]);
        if (dist <= distMax){
            inliers.push_back(i);
        }
    }
    return inliers;
}


// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;

    //FMatrix<float,3,3> F;
    vector<int> Inliers;
    
    int n = matches.size();
    int m = 0; // Estimated number of inliers
    int c = 0;
    cout << "Computing F...\n";
    while (c < Niter){
        FMatrix<float,3,3> F = randF(matches);

        Inliers = findInliers(F, matches, distMax);

        if (Inliers.size() > m){ // More inliers than the previous best
            bestF = F;
            bestInliers = Inliers;
            m = Inliers.size();
            Niter = (int)(log(BETA)/log(1-0.0001-pow(float(m)/float(n),8)));
        }
        c++;
    }
    cout << "Number of iterations " << c << "\n";
    cout << "Number of inliers " << m << "\n";
    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    // F not refined using the inliers as I wasn't able to compute the svd
    // of a FMatrix without fixing it's size during the compilation
    cout << "Done computing F\n"; 
    return bestF;
}

// Draw a cross at position (x,y)
void drawCross(int x, int y, Color color, int size = 5, int penWidth=1){
    drawLine(x-size, y-size, x+size, y+size, color, penWidth);
    drawLine(x-size, y+size, x+size, y-size, color, penWidth);
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    int c = 0;
    while(true) {
        int x,y;
        DoublePoint3 point1, point2;
        int middle = I1.width();
        if(getMouse(x,y) == 3)
            break;
        c ++;
        Color color = Color(156+rand()%100,86+rand()%170,86+rand()%170);
        drawCross(x,y,color);

        if (x < middle){ // Click in left image
            point1[0] = (float)x;
            point1[1] = (float)y;
            point1[2] = 1.0;
            point2 = transpose(F) * point1; // Line equation
            float yMiddle = -point2[2]/point2[1]; // Applied in x=0
            float yMax = -(point2[0]*(float)(I2.width()) + point2[2])/point2[1]; // Applied in x=I1.width
            drawLine(middle,(int)yMiddle,I1.width()+I2.width(),(int)yMax,color);
        }
        else { // Click in right image
            point1[0] = (float)(x - middle);
            point1[1] = (float)y;
            point1[2] = 1.0;
            point2 = F * point1;
            float y0 = -point2[2]/point2[1];
            float yMiddle = -(point2[0]*(float)middle + point2[2])/point2[1];
            drawLine(0,(int)y0,middle,(int)yMiddle,color);
        }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}

#include <Imagine/Images.h>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "maxflow/graph.h"

using namespace std;
using namespace Imagine;

typedef Image<byte> byteImage;
typedef Image<double> doubleImage;

/// 4-neighbors
static const int dx[]={+1,  0, -1,  0};
static const int dy[]={ 0, -1,  0, +1};

static const int nNeighbors = (sizeof(dx)/sizeof(*dx));

/*
 * Return image of mean intensity value over (2n+1)x(2n+1) patch
 */
doubleImage meanImage(const doubleImage&I, int n)
{
    // Create image for mean values
    int w = I.width(), h = I.height();
    doubleImage IM(w,h);
    // Compute patch area
    double area = (2*n+1)*(2*n+1);
    // For each pixel
    for (int i=0;i<w;i++)
	for (int j=0;j<h;j++) {
	    // If pixel is close to border (less than n pixels)
	    if (j<n || j>=h-n || i<n ||i>=w-n) {
		// Mean is meaningless
		IM(i,j)=0;
		continue;
	    }
	    // Otherwise, initialize sum
	    double sum = 0;
	    // For each pixel displacement in patch
	    for(int x = i-n; x <= i+n; x++)
		for (int y = j-n; y <= j+n; y++)
		    // Sum pixel intensity
		    sum += I(x,y);
	    // Remember average intensity
	    IM(i,j) = sum / area;
	}
    return IM;
}


/*
 * Compute correlation between two pixels in images 1 and 2
 */
double correl(
    const doubleImage& I1,  // Image 1
    const doubleImage& I1M, // Image of mean intensity value over patch
    const doubleImage& I2,  // Image2
    const doubleImage& I2M, // Image of mean intensity value over patch
    int u1, int v1,         // Pixel of interest in image 1
    int u2, int v2,         // Pixel of interest in image 2
    int n)                  // Half patch size
{
    // Initialize correlation
    double c = 0;
    // For each pixel displacement in patch
    for (int x = -n; x <= n; x++) 
	for(int y = -n; y <= n; y++)
	    // Sum correllation
	    c += (I1(u1+x,v1+y) - I1M(u1,v1)) * (I2(u2+x,v2+y) - I2M(u2,v2));
    // Return correlation normalized by patch size
    return c / ((2*n+1)*(2*n+1));
}


/*
 * Compute ZNCC between two patches in images 1 and 2
 */
double zncc(
    const doubleImage& I1,  // Image 1
    const doubleImage& I1M, // Image of mean intensity value over patch
    const doubleImage& I2,  // Image2
    const doubleImage& I2M, // Image of mean intensity value over patch
    int u1, int v1,         // Pixel of interest in image 1
    int u2, int v2,         // Pixel of interest in image 2
    int n)                  // Half patch size
{
    // Compute variance of patch in image 1
    double var1 = correl(I1, I1M, I1, I1M, u1, v1, u1, v1, n);
    if (var1 == 0)
	return 0;
    // Compute variance of patch in image 2 
    double var2 = correl(I2, I2M, I2, I2M, u2, v2, u2, v2, n);
    if (var2 == 0)
	return 0;
    // Return normalized cross correlation
    return correl(I1, I1M, I2, I2M, u1, v1, u2, v2, n) / sqrt(var1 * var2);
}

double rho(double c){
    return (c < 0.0) ? 1 : sqrt(1-c);
}

/*
 * Load two rectified images.
 * Compute the disparity of image 2 w.r.t. image 1.
 * Display disparity map.
 * Display 3D mesh of corresponding depth map.
 */
int main()
{
    ///// Load and crop images
    cout << "Loading images... " << flush;
    byteImage I;
    doubleImage I1,I2; // First and second image
    // Load Thierry's face (seen from bottom)
    load(I, srcPath("face00R.png"));
    // Crop it
    I1 = I.getSubImage(IntPoint2(20,30), IntPoint2(430,420));
    // Load Thierry's face (seen from top)
    load(I, srcPath("face01R.png"));
    // Crop it
    I2 = I.getSubImage(IntPoint2(20,30), IntPoint2(480,420));
    // Done
    cout << "done" << endl;


    ///// Set parameters
    cout << "Setting parameters... " << flush;
    // Generic parameters
    int zoom = 2;      // Zoom factor (to speedup computations)
    int n = 3;         // Consider correlation patches of size (2n+1)*(2n+1)
    float lambdaf = 0.1; // Weight of regularization (smoothing) term
    int wcc = max(1+int(1/lambdaf),20); // Energy discretization precision [as we build a graph with 'int' weights]
    int lambda = lambdaf * wcc; // Weight of regularization (smoothing) term [must be >= 1]
    float sigma = 3;   // Gaussian blur parameter for disparity
    // Image-specific, hard-coded parameters for approximate 3D reconstruction,
    // as real geometry before rectification is not known
    int dmin = 10;     // Minmum disparity
    int dmax = 55;     // Maximum disparity
    float fB = 40000;  // Depth factor 
    float db = 100;    // Disparity base
    // Done
    cout << "done" << endl;


    ////// Display images
    cout << "Displaying images... " << flush;
    // Open window large enough for both images
    int w1 = I1.width(), w2 = I2.width(), h = I1.height();
    openWindow(w1+w2, h);
    // Display first image
    display(grey(I1));
    // Display second image
    display(grey(I2), w1, 0);
    // Done
    cout << "done" << endl;


    ///// Construct graph
    cout << "Constructing graph (be patient)... " << flush;
    // Precompute images of mean intensity value over patch
    doubleImage I1M = meanImage(I1,n), I2M = meanImage(I2,n);
    // Zoomed image dimension, disregarding borders (strips of width equal to patch half-size)
    int nx = (w1 - 2*n) / zoom, ny = (h - 2*n) / zoom;
    // Disparity range
    int nd = dmax-dmin;
    // "Infinite" value
    int INF=1000000;
    // Value of the K penalizer
    double K = 1 + (nd-1) * nNeighbors * lambda;
    // Create graph
    /////------------------------------------------------------------
    auto node = [nx, ny, nd](int x, int y, int d){
        return d + nd * (x + nx * y);
    };
    
    auto noded = [nd](int n){
        return n % nd;
    };

    auto nodex = [nx, nd](int n){
        return (n / nd) % nx;
    };

    auto nodey = [nx, ny, nd](int n){
        return (n / (nd* nx));
    };

    int nPixel = nx * ny;
    int nNode = nPixel * nd;
    int nEdge = nPixel * (1+ nd) - nx -ny;
    Graph<double,double,double> G(nNode,nEdge);
    G.add_node(nNode);
    
    cout << "\n!!! FIND GOOD CONDITIONS FOR X IN WINDOW !!\n";
    cout << "Add % meter, check cap computation\n";
    for (int x = 5; x < nx - 10; x++){
        for (int y = 5; y < ny - 5; y++){
            for (int d = 0; d < nd; d++){ // Peut-être nd+1 mais il faudrait changer nnode
                // Edges to neighbors
                for(int i=0; i<nNeighbors; i++) { // Bourrin (arrete bidir pourrait être créée en une fois)
                    int xn=x+dx[i], yn=y+dy[i];
                    int noden = node(xn, yn, d);
                    if ((noden < 0) || (noden >= nNode)){continue;}
                    G.add_edge(node(x,y,d),noden,lambda,0);
                }
                // Edges between levels of disparity
                double cap = K;
                cap += wcc * rho(zncc(I1, I1M, I2, I2M, x, y, x + d + dmin, y, n));
                if (d == 0){
                    G.add_tweights(node(x,y,0), cap, 0);
                }
                else if (d == (nd - 1)){
                    G.add_tweights(node(x,y,d), 0, cap);
                }
                else {
                    G.add_edge(node(x,y,d-1), node(x,y,d), cap, 0);
                }                
            }
        }
    }
    /////------------------------------------------------------------
    // Done
    cout << "done" << endl;


    ///// Compute cut
    cout << "Computing minimum cut... " << flush;
    int f = G.maxflow();
    // Done
    cout << "done" << endl
	 << "  max flow = " << f << endl;


    ///// Extract disparity map from minimum cut
    cout << "Extracting disparity map from minimum cut... " << flush;
    doubleImage D(nx,ny);
    // For each pixel
    for (int i=0;i<nx;i++) {
	for (int j=0;j<ny;j++) {
	    ///// Extract disparity from minimum cut
	    /////------------------------------------------------------------
	    /////  BEGIN CODE TO BE COMPLETED: define disparity map D from graph G and minimum cut
	    /////
	    /* WARNING: dummy code to replace */ D(i,j) = 40+2*I1(n+i*zoom,n+j*zoom)/(1+(i>j));
	    /* WARNING: dummy code to replace */ if (i+j==0) cout<<endl<<"***\n*** Dummy computation: code has to be completed!\n***   (3D generated by default does not make sense)\n***"<<endl;
	    /////
	    /////  END CODE TO BE COMPLETED
	    /////------------------------------------------------------------
	}
    }
    // Done
    cout << "done" << endl;


    ///// Display disparity map
    cout << "Displaying disparity map... " << flush;
    // Display disparity map
    display(enlarge(grey(D),zoom), n, n);
    // Done
    cout << "done" << endl;


    ///// Compute and display blured disparity map
    cout << "Click to compute and display blured disparity map... " << flush;
    click();
    D=blur(D,sigma);
    display(enlarge(grey(D),zoom),n,n);
    // Done
    cout << "done" << endl;


    ///// Compute depth map and 3D mesh renderings
    cout << "Click to compute depth map and 3D mesh renderings... " << flush;
    click();
    // Open new window and make it active
    Window W = openWindow3D(512, 512, "3D");
    setActiveWindow(W);
    ///// Compute 3D points
    Array<FloatPoint3> p(nx*ny);
    Array<Color> pcol(nx*ny);
    // For each pixel in image 1
    for (int i = 0;i < nx; i++)
	for (int j = 0; j < ny; j++) {
	    // Compute depth: magic constants depending on camera pose
	    float depth = fB / (db + D(i,j)-dmin);
	    // Create 3D point
	    p[i+nx*j] = FloatPoint3(float(i), float(j), -depth);
	    // Get point color in original image
	    byte g = byte( I1(n+i*zoom,n+j*zoom) );
	    pcol[i+nx*j] = Color(g,g,g);
	}
    ///// Create mesh from 3D points
    Array<Triangle> t(2*(nx-1)*(ny-1));
    Array<Color> tcol(2*(nx-1)*(ny-1));
    // For each pixel in image 1 (but last line/column)
    for (int i = 0; i < nx-1; i++)
	for (int j = 0; j < ny-1; j++) {
	    // Create triangles with next pixels in line/column
	    t[2*(i+j*(nx-1))] = Triangle(i+nx*j, i+1+nx*j, i+nx*(j+1));
	    t[2*(i+j*(nx-1))+1] = Triangle(i+1+nx*j, i+1+nx*(j+1), i+nx*(j+1));
	    tcol[2*(i+j*(nx-1))] = pcol[i+nx*j];
	    tcol[2*(i+j*(nx-1))+1] = pcol[i+nx*j];
	}
    // Create first mesh as textured with colors taken from original image
    Mesh Mt(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0, FACE_COLOR);
    Mt.setColors(TRIANGLE, tcol.data());
    // Create second mesh with artificial light
    Mesh Mg(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0,
	    CONSTANT_COLOR, SMOOTH_SHADING);
    // Done
    cout << "done" << endl;


    ///// Display 3D mesh renderings
    cout << "***** 3D mesh renderings *****" << endl;
    cout << "- Button 1: toggle textured or gray rendering" << endl;
    cout << "- SHIFT+Button 1: rotate" << endl;
    cout << "- SHIFT+Button 3: translate" << endl;
    cout << "- Mouse wheel: zoom" << endl;
    cout << "- SHIFT+a: zoom out" << endl;
    cout << "- SHIFT+z: zoom in" << endl;
    cout << "- SHIFT+r: recenter camera" << endl;
    cout << "- SHIFT+m: toggle solid/wire/points mode" << endl;
    cout << "- Button 3: exit" << endl;
    // Display textured mesh
    showMesh(Mt);
    bool textured = true;
    // Wait for mouse events
    Event evt;
    while (true) {
	getEvent(5, evt);
	// On mouse button 1
	if (evt.type == EVT_BUT_ON && evt.button == 1) {
	// Toggle textured rendering and gray rendering
	    if (textured) {
		hideMesh(Mt,false);
		showMesh(Mg,false);
	    } else {
		hideMesh(Mg,false);
		showMesh(Mt,false);
	    }
	    textured = !textured;
	}
	// On mouse button 3
	if (evt.type == EVT_BUT_ON && evt.button == 3)
	    // Exit
	    break;
    }
    // That's all folks!
    endGraphics();
    return 0;
}

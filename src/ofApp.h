#pragma once
#include "ofMain.h"
#include "ofxBox2d.h"
#include "ofxGui.h"
#include "ofxCsv.h"
#include "particleClass.hpp"

// -------------------------------------------------


class ofApp : public ofBaseApp {
	
public:
	
	void setup();
	void update();
	void draw();
	
	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void resized(int w, int h);
	
	ofxBox2d box2d;   // the box2d world
    ParticleSystem myParticleSystem;
    ParticleSystem particleArrayBorderFixed;  // the cavity
    ParticleSystem particleArrayBorderMobile; // the mobile part of the "cavity" (piston for instance)
   
    void addParticle(const ofVec2f& _pos, const ofVec2f& _vel, const float _mass,  const float _charge, const float _size, const float _damp, ParticleSystem &  _particleArray);
    uint16_t addSpiralParticles(const ofVec2f& _center, const float _radius, const float sizeBot, ParticleSystem& _particleArray);
    void deleteLastParticle(ParticleSystem &  _particleArray);
    
    
    ofxFloatSlider gravity;
    void gravityChange(float &_gravity);
    
    float scaleFactor; // modelview scale
   
    vector <shared_ptr<ofxBox2dEdge> >   edges;
    
    ofRectangle confinementRectangle;
    enum modeContainer {disk, rect}; // torus?
    modeContainer modeCavity;
    void createCavity();
    void makeContainer();
    void makeContainerRect(const ofRectangle & _rect);
    void makeContainerDisk(float cx, float cy, float radius, int numPoints);
    
    
    
    // Container with fixed masses:
    void makeContainerDiskFixedMasses(const ofVec2f& _center, float radius, uint16_t _numMasses, ParticleSystem& _particleArrayBorderFixed);
    void makeContainerFixedRect(const ofRectangle & _rect, uint16_t _numMasses, ParticleSystem& _particleArrayBorderFixed); // return num masses vertical (for oscillation)
    void makeMobileWall(const ofVec2f _A, const ofVec2f _B,  uint16_t _numMasses, ParticleSystem& _particleArrayBorderMobile);
     
    void translateParticleSystem(const ofVec2f& _translation, ParticleSystem& _particleArrayBorderMObile);
    ofxFloatSlider impulseStrength;
    
    
    ofxToggle oscillateSide;
    bool oscillateSidePrevious = true; // need to start as true to reset the offset timer for the phase
    float timeStop, timeOffset =0;
    ofxFloatSlider freqOscillateSide;
    ofxFloatSlider ampOscillateSide;
    void resetMobileWall();
    float phaseStart = 0;
    
    enum colorSpeedMode {colorTangent, colorRadial, colorNorm};
    colorSpeedMode  modeColor = colorTangent;
   
    ofxToggle runStop;
    
    //Interaction Force Type:
    // ... there is no good widget for that!
    
    // Drawing modes:
    ofxToggle drawParticles;
    ofxToggle drawDisplacements;
    ofxToggle drawSpeeds;
   // ofxToggle colorTangent;
    
    ofxFloatSlider sizeParticle;
    ofxToggle sizeMode; // if true, resize from speed norm
    
    // bool addForcesRepulsion = false;
    ofxToggle addForcesRepulsion;
    ofxFloatSlider amplitudeRepulsion;
    
 
    ofxFloatSlider charge;
    void chargeUpdateChange(float & _charge);
    
    ofxFloatSlider massParticle;
    void massUpdateChange(float & _mass);
    bool flagChangeMass = false;
    
    ofxToggle addForceCentrale;
    ofxFloatSlider amplitudeForceCentrale;
    
    ofxToggle blackHole;
    ofxFloatSlider blackHoleRadius;
    ofxToggle rebirth;
    
    //bool addForcesRandom = true;
    ofxToggle addForcesRandom;
    ofxFloatSlider amplitude;
     // float amplitude = 150;
    
    
    //float scatterAngle = 10;
    ofxFloatSlider scatterAngle;
    
    //float dampingFactor = 0.0000000001;
    ofxFloatSlider dampingFactor;
    void changeDampingFactor(float & _dampingFactor);
    
    uint64_t timerExcitation;
    //uint64_t periodMillisExcitation = 10;
    ofxFloatSlider frequencyExcitation;
    
    //bool freeJumpMode = false;
    ofxToggle parabolicJump;
 
   // uint64_t periodMillisStopMotion = 10;
   // ofxIntSlider periodMillisStopMotion;
    
     ofxFloatSlider frequencyUpdate;
    void frequenyUpdateChange(float &_freqUpdate);
    
    ofxButton resetSpeed;
    void resetSpeedAll();
    
    // Time reversal ("push button"):
    ofxButton timeReversal;
    void timeReversalAll();
    
    ofxButton restart;
    void restartSystem();
    
    
    ofxButton resetToEquilibriumPositions;
    void resetToEquilibriumPos();
    ofxButton saveEquilibriumPositions;
    void saveEquilibriumPos();
    
    
    // =========== GUI PANEL & screen info ============
    ofxPanel gui;
    ofxToggle bHide;
    
    
    // LOGGING DATA:
    ofxToggle logPhaseSpaceData;
    ofxButton clearPhaseSpaceData;
    ofxButton savePhaseSpacDataFile;
    
    // Saving a "snapshot" of the system, including speeds (note: the type of force, etc, is not saved yet):
    
    // (a) in memory
    ofxButton  saveSnapshotStateMemory;
    ofxButton  loadSnapshotStateMemory;
    void loadConfigurationMemory();
    void saveConfigurationMemory();
    
    // (b) in file
    ofxCsv csv;
    ofxButton  saveSnapshotState;
    ofxButton  loadSnapshotState;
    void loadConfiguration();
    void saveConfiguration();
    
    void savePhaseSpaceTrajectory();
    void clearLogData();
    void toggleLogging();
    
    // ==================================
    
    ofPolyline waveLine;
    
    //===========================
    // CONFINEMENT SHAPE:
    ofVec2f posCenter;// = ofVec2f(ofGetWidth()/2+100, ofGetHeight()/2);
        float radiusContainer;
    bool continueAddingCenter = false;
       bool continueAddingLines = false;
       vector <ofPolyline> lines;
    void deleteLastBorder();
    float displacement = -125;
       int addCount = 0;
    
    
    ofEasyCam cam; // add mouse controls for camera movement
    bool cam3D;
};


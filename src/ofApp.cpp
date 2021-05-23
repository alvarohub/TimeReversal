#include "ofApp.h"


#define FACTOR_MASS 10000.0

//--------------------------------------------------------------
void ofApp::setup() {

    modeCavity = rect; // disk
    
    modeColor = colorNorm;
    
    box2d.init();
    //box2d.createGround();
    box2d.enableGrabbing();
    box2d.registerGrabbing();
    
    myParticleSystem.init(&box2d);
    particleArrayBorderFixed.init(&box2d);
    particleArrayBorderMobile.init(&box2d);
    
    radiusContainer =ofGetHeight()/2.3;
    posCenter = ofVec2f(ofGetWidth()/2+200, ofGetHeight()/2);
    //confinementRectangle.setFromCenter(posCenter, ofGetWidth()-500, 1.6*radiusContainer/15);
    confinementRectangle.setFromCenter(posCenter, 2*radiusContainer, 1.6*radiusContainer);///15);
    
    // Setting GUI:
    massParticle.addListener(this, &ofApp::massUpdateChange);
    charge.addListener(this, &ofApp::chargeUpdateChange);
    frequencyUpdate.addListener(this, &ofApp::frequenyUpdateChange);
    dampingFactor.addListener(this, &ofApp::changeDampingFactor);
    gravity.addListener(this, &ofApp::gravityChange);
    resetSpeed.addListener(this, &ofApp::resetSpeedAll);
    timeReversal.addListener(this, &ofApp::timeReversalAll);
    restart.addListener(this, &ofApp::restartSystem);
    //oscillateSide.addListener(this, &ofApp::resetMobileWall); // Arg, listeners cannot be added to toggles?
    
    saveSnapshotStateMemory.addListener(this, &ofApp::saveConfigurationMemory);
    loadSnapshotStateMemory.addListener(this, &ofApp::loadConfigurationMemory);
    
    saveSnapshotState.addListener(this, &ofApp::saveConfiguration);
    loadSnapshotState.addListener(this, &ofApp::loadConfiguration);
    
    savePhaseSpacDataFile.addListener(this, &ofApp::savePhaseSpaceTrajectory);
    clearPhaseSpaceData.addListener(this, &ofApp::clearLogData);
    
    resetToEquilibriumPositions.addListener(this, &ofApp::resetToEquilibriumPos);
    saveEquilibriumPositions.addListener(this, &ofApp::saveEquilibriumPos);
    
    
    gui.setup();
    //gui.enableHiDpi();
    gui.setUseTTF(true);// loadFont(const std::string& filename, int fontsize, bool _bAntiAliased = true, bool _bFullCharacterSet = false, int dpi = 0);
    
    gui.add(runStop.setup("Run/Stop", true));
    gui.add(restart.setup("RESTART system"));
    
    // TODO: SAVE ALSO CURRENT SPEEDS!!!
    gui.add(saveEquilibriumPositions.setup("Save as EQUILIBRIUM"));
    gui.add(resetToEquilibriumPositions.setup("Reset to EQUILIBRIUM (v=0)"));
    gui.add(saveSnapshotStateMemory.setup("Snapshot state in MEM"));
    gui.add(loadSnapshotStateMemory.setup("Retrieve state from MEM"));
    gui.add(saveSnapshotState.setup("Snapshot state in FILE"));
    gui.add(loadSnapshotState.setup("Load state from FILE"));
    gui.add(logPhaseSpaceData.setup("Start/Stop LOG", false));
    gui.add(clearPhaseSpaceData.setup("Clear LOG"));
    gui.add(savePhaseSpacDataFile.setup("Save LOG in FILE"));
    
    gui.add(resetSpeed.setup("RESET SPEED"));
    gui.add(timeReversal.setup("REVERSE SPEED"));
    
    gui.add(frequencyUpdate.setup("Freq. Update",60,0,500));
    
    gui.add(massParticle.setup("Mass Particle", 1,0,2));
    flagChangeMass = false;
    gui.add(charge.setup("Charge", 1,-1.0,1.0));
    
    gui.add(dampingFactor.setup("Damping", 0.3, 0, 20));
    gui.add(gravity.setup("Gravity", 0,0,20));
    
    gui.add(addForceCentrale.setup("Force Centrale", false));
    gui.add(amplitudeForceCentrale.setup("Amplitude Force Centrale", 0, -20, 20));
    
    gui.add(blackHole.setup("Black Hole", false));
    gui.add(blackHoleRadius.setup("Radius Black Hole", 5, 0, 50));
    gui.add(rebirth.setup("Particle Rebirth", false));
    
    gui.add(addForcesRepulsion.setup("Repulsion", true));
    gui.add(amplitudeRepulsion.setup("Amplitude Repulsion", .1, 0, 30.0));
    
    gui.add(addForcesRandom.setup("Vibration", false));
    gui.add(frequencyExcitation.setup("Frequency Vibration",50,0,200));
    gui.add(amplitude.setup("Amplitude Vibration", 0, 20, 1000));
    gui.add(scatterAngle.setup("Scatter angle", 40, 0, 360));
    gui.add(parabolicJump.setup("Parabolic Jump", false));
    
    
    //gui.add(colorTangent.setup("Color Tangent/Speed", false));
    gui.add(drawParticles.setup("Draw particles", true));
    gui.add(drawDisplacements.setup("Draw displacements", true));
    //gui.add(drawSpeeds.setup("Draw particles", true));
    
    gui.add(sizeParticle.setup("Size particle", 10,0,50));
    gui.add(sizeMode.setup("Resize from Speed", false));
    
    gui.add(oscillateSide.setup("Side Oscilaltion", false));
    gui.add(freqOscillateSide.setup("Freq. Oscillation (Hz)", 0.5, 0, 10));
    gui.add(ampOscillateSide.setup("Amp. Oscillation", 10, 0, 500));
    gui.add(impulseStrength.setup("Impulse (horizontal) on particle", 0, -1, 1));
    
    gui.add(bHide.setup("Show/Hide keyboard commands & data", false));
    
    
    addCount = 1;
    
    scaleFactor = 1;
    
    ofSetVerticalSync(true);
    ofBackground(0);
    ofSetLogLevel(OF_LOG_NOTICE);
    ofSetVerticalSync(true);
    ofDisableAntiAliasing();
    
    createCavity();
    
    // NOTE: set as default the size of the particules in the system equal to the size of the confinement particlues:
    sizeParticle = particleArrayBorderMobile.back()->getRadius();
    
}

void ofApp:: createCavity() {
    int numParticlesMobileWall;
    numParticlesMobileWall = 200.0/confinementRectangle.getPerimeter()*confinementRectangle.getHeight();
    switch(modeCavity) {
        case disk:
            makeContainerDiskFixedMasses(posCenter,radiusContainer, 200, particleArrayBorderFixed);
            makeMobileWall(confinementRectangle.getBottomLeft(), confinementRectangle.getTopLeft(), numParticlesMobileWall, particleArrayBorderMobile);
            //makeContainerDisk(posCenter.x, posCenter.y,radiusContainer, 60);
            //makeContainerDisk(posCenter.x, posCenter.y,radiusContainer-sizeParticle*3.5, 60);
            break;
            //        case torus:
            //            makeContainerDisk(posCenter.x, posCenter.y, ofGetHeight()/4, 30);
            //            makeContainerDisk(posCenter.x, posCenter.y, ofGetHeight()/5, 30);
            //            break;
        case rect:
            makeContainerFixedRect(confinementRectangle, 200, particleArrayBorderFixed);
            makeMobileWall(confinementRectangle.getBottomLeft(), confinementRectangle.getTopLeft(), numParticlesMobileWall, particleArrayBorderMobile);
            
            //  makeContainerRect(confinementRectangle);
            break;
        default:
            break;
    }
    
    // timerStopMotion = ofGetElapsedTimeMillis();
    timerExcitation = ofGetElapsedTimeMicros();
}

void ofApp::addParticle(const ofVec2f& _pos, const ofVec2f& _vel, const float _mass,  const float _charge, const float _size, const float _damp, ParticleSystem &  _particleArray) {
    // NOTE: if _mass = 0, the body is FIXED.
    ofVec2f T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
    ofVec2f X = (_pos - T)/scaleFactor ;
    _particleArray.addParticle(X, X, _vel, _mass, _charge, _size, _damp);
}



// Add bots on an archimedian spiral (equal *steph* length, not constant angle step!)
uint16_t ofApp::addSpiralParticles(const ofVec2f& _center, const float _radius, const float _radiusBot, ParticleSystem & _particleArray) {
    //   const float _radiusArm, // r = _radiusArm * theta
    //   const float _numTours,
    //   const uint16_t _numPoints,
    
    float radiusArm = 2*_radiusBot;
    float numTours = 1.0*floor(_radius / radiusArm);
    
    float phi = 2.0 * PI * numTours;
    float length = radiusArm / (4.0 * PI) * (phi * sqrt(1 + phi * phi) + log(phi + sqrt(1 + phi * phi)));
    
    uint16_t numAddedBots = floor (length /_radiusBot/2);
    
    float stepLength = 7 * length / (numAddedBots - 1);
    float theta = 0, stepTheta = 0;
    
    uint16_t addCounterBots = 0;
    
    // 1) Go outwards:
    while (theta <= phi)
    {
        float r = radiusArm * theta / 2 / PI;
        
        cout << "adding bot: "<<addCounterBots<<endl;
        
        ofVec2f pos(_center.x + r * cos(theta), _center.y + r * sin(theta));
        addParticle(pos, ofVec2f(0,0), massParticle/FACTOR_MASS, charge, sizeParticle, dampingFactor, _particleArray);
        
        // Use dicotomy to find the stephTheta such that the length increase is equal to stepLength:
        // float stepTheta = 1.0*phi/_numPoints; // the step should be smaller than that
        // ... OR, for large number of points, we have the approximation:
        stepTheta = stepLength / (radiusArm * sqrt(1 + theta * theta));
        
        theta += stepTheta;
        addCounterBots++;
    }
    
    return(addCounterBots);
}

void ofApp::deleteLastParticle(ParticleSystem &  _particleArray) {
    _particleArray.deleteLastParticle();
}

//--------------------------------------------------------------
void ofApp::update() {
    if (runStop) {
        
        bool updateExcitation = false, updateStopMotion = false;
        
        // oscillation left side container:
        if (!oscillateSidePrevious && oscillateSide) {
            oscillateSidePrevious = true;
            //resetMobileWall();
            timeOffset+= ofGetElapsedTimef() - timeStop;
        }
        if (oscillateSidePrevious && !oscillateSide) {
            float phase = 2.0*PI*freqOscillateSide*(ofGetElapsedTimef() - timeOffset);
            float cycles = phase/2.0/PI;
            if (cycles -  floor(cycles) < 0.1) {
                oscillateSidePrevious = false;
                timeStop = ofGetElapsedTimef();
            }
            // resetMobileWall();
        }
        
        if (oscillateSide || oscillateSidePrevious) {
            float phase = 2.0*PI*freqOscillateSide*(ofGetElapsedTimef() - timeOffset);
            ofVec2f translation = ofVec2f( ampOscillateSide*(sin(phase)), 0);
            translateParticleSystem(translation, particleArrayBorderMobile);
            
        }
        
        
        //    if (  ofGetElapsedTimeMicros() - timerUpdate > 1000000.0/frequencyUpdate ) {
        //              updateDynamics = true;
        //              timerUpdate = ofGetElapsedTimeMicros();
        //          }
        
        
        //ofVec2f mouse(ofGetMouseX(), ofGetMouseY());
        //float minDis = ofGetMousePressed() ? 300 : 100;
        
        ofSetHexColor(0x444342);
        ofNoFill();
        for (int i=0; i<lines.size(); i++) {
            lines[i].draw();
        }
        
        /*
         for (int i=0; i<edges.size(); i++) {
         edges[i].get()->draw();
         }
         */
        
        
        // Add repulsion?
        if (addForcesRepulsion) {
            
            
            for(int i=0; i<myParticleSystem.size(); i++) {
                
                Particle* ptr_particleA = myParticleSystem.getParticle(i);
                ofVec2f posA = ptr_particleA->getPosition();
                
                //(a) between confined particles:
                for(int j=i+1; j<myParticleSystem.size(); j++) {
                    Particle* ptr_particleB = myParticleSystem.getParticle(j);
                    ofVec2f posB = ptr_particleB->getPosition();
                    
                    //
                    //                    ptr_particleA->addPointForce_InvR(posB, amplitudeRepulsion/100000000.0);
                    //                    ptr_particleB->addPointForce_InvR(posA, amplitudeRepulsion/100000000.0);
                    //
                    ptr_particleA->addPointForce_InvR3(posB, amplitudeRepulsion/100000000.0);
                    ptr_particleB->addPointForce_InvR3(posA, amplitudeRepulsion/100000000.0);
                    
                    //                    ptr_particleA->addPointForce_InvR4(posB, amplitudeRepulsion/100000000.0);
                    //                    ptr_particleB->addPointForce_InvR4(posA, amplitudeRepulsion/100000000.0);
                    //
                    //                    ptr_particleA->addPointForceSpring(posB, .0001, 3*sizeParticle);
                    //                    ptr_particleB->addPointForceSpring(posA, .0001, 3*sizeParticle);
                    
                    
                }
                
                //(b) between confined particles and the confinement box and mobile wall:
                // (the mobile wall, and the confinement box particles dont interact between themselves)
                for(int j=0; j<particleArrayBorderMobile.size(); j++) {
                    Particle* ptr_particleB = particleArrayBorderMobile.getParticle(j);
                    ofVec2f posB = ptr_particleB->getPosition();
                    ptr_particleA->addPointForce_InvR4(posB, amplitudeRepulsion/100000000.0);
                    
                    //   ptr_particleA->addPointForceSpring(posB, .000001, 2*sizeParticle);
                    
                }
                for(int j=0; j<particleArrayBorderFixed.size(); j++) {
                    Particle* ptr_particleB = particleArrayBorderFixed.getParticle(j);
                    ofVec2f posB = ptr_particleB->getPosition();
                    ptr_particleA->addPointForce_InvR4(posB, amplitudeRepulsion/100000000.0);
                    //    ptr_particleA->addPointForceSpring(posB, .000001, 2*sizeParticle);
                }
                
                
            }
        }
        
        if (addForceCentrale) {
            for(int i= 0; i<myParticleSystem.size(); i++) {
                Particle* ptr_particleA = myParticleSystem.getParticle(i);
                ofVec2f posA = ptr_particleA->getPosition();
                // ofVec2f posCenter = ofVec2f(ofGetWidth()/2+100, ofGetHeight()/2);
                
                if ((posA-posCenter).length() > blackHoleRadius*2/3)
                {
                    //particleArray[i].get()->addRepulsionForceInv(posCenter, -amplitudeForceCentrale/FACTOR_MASS);
                    ptr_particleA->addPointForceConst(posCenter, -amplitudeForceCentrale/FACTOR_MASS);
                    
                } else { // make it repel from the exact center (stronger)
                    ptr_particleA->addPointForceConst(posCenter, amplitudeForceCentrale/1000.0);
                    
                }
                
                
            }
        }
        
        
        if (blackHole) {
            for(int i= 0; i<myParticleSystem.size(); i++) {
                ofVec3f pos = myParticleSystem.getParticle(i)->getPosition();
                ofVec3f OP = pos - posCenter;
                float distCenter = OP.length();
                
                if (distCenter<=blackHoleRadius)
                {
                    myParticleSystem.erase(i);
                    
                    // Make it appear again on the border?
                    if (rebirth) {
                        ofVec2f posBirth = posCenter;
                        float angle = ofRandom(0, 2*PI);
                        posBirth = posBirth + (radiusContainer-2*sizeParticle)*ofVec2f(cos(angle), sin(angle));
                        
                        // add radial velocity?
                        //ofVec3f ur = (posBirth-posCenter).normalize();
                        //ofVec2f utheta(-ur.y, ur.x);
                        //addParticle(posBirth.x, posBirth.y, massParticle/FACTOR_MASS, sizeParticle, 5*utheta.rotate(20), 1);
                        ofVec2f vit(0,0);
                        addParticle(ofVec2f(posBirth.x, posBirth.y), vit, massParticle/FACTOR_MASS, charge, sizeParticle, dampingFactor, myParticleSystem);
                        
                    }
                    else
                        i=i-1;
                    
                    
                }
            }
        }
        
        
        if (  ofGetElapsedTimeMicros() - timerExcitation > 1000000.0/frequencyExcitation ) {
            updateExcitation = true;
            timerExcitation = ofGetElapsedTimeMicros();
        }
        for(int i= 0; i<myParticleSystem.size(); i++) {
            //float dis = mouse.distance(particleArray[i].get()->getPosition());
            //if(dis < minDis) particleArray[i].get()->addRepulsionForce(mouse, 9);
            //else particleArray[i].get()->addAttractionPoint(mouse, 4.0);
            
            Particle* ptr_particle = myParticleSystem.getParticle(i);
            
            ptr_particle->myEventData.onFloor = true;
            
            if(parabolicJump) {
                if (ofGetElapsedTimeMicros() -  ptr_particle->myEventData.lastTimeExcitation > 0.5*ptr_particle->myEventData.timeFreeFall*100000) {
                    ptr_particle->setVelocity(0, 0);
                } else  ptr_particle->myEventData.onFloor = false;
            }
            
            if (updateExcitation &&  ptr_particle->myEventData.onFloor) {
                
                ptr_particle->myEventData.lastTimeExcitation = ofGetElapsedTimeMicros();
                
                
                if (addForcesRandom) {
                    
                    float headAngle = ptr_particle->getRotation();
                    ofVec2f direction(1.0,0.0);
                    direction.rotate(headAngle);
                    
                    // 1) Impulse:
                    //virtual void addImpulseForce(ofVec2f pt, ofVec2f amt);
                    ptr_particle->addImpulseForce(ptr_particle->getPosition(), 0.000001*amplitude*direction);
                    
                    // Force:
                    //   particleArray[i].get()->addForce(direction, amplitude/1000);
                    
                    // Using speed setting is not realistic: (problem with collisions)
                    // particleArray[i].get()->setVelocity(0.01*amplitude*direction);//0, 0);
                    
                    // 2)  Thumble:
                    ptr_particle->setRotation(headAngle+ofRandom(-scatterAngle/2,scatterAngle/2));
                    
                }
            }
        }
        
        
        if (flagChangeMass) {
            for(int i= 0; i<myParticleSystem.size(); i++) {
                Particle* ptr_particle = myParticleSystem.getParticle(i);
                ptr_particle->setMass(massParticle/FACTOR_MASS);
                ptr_particle->setFixedRotation(true);
            }
            flagChangeMass = false;
        }
        
        box2d.update();
        
        if (logPhaseSpaceData) myParticleSystem.updateLog();
        
        
        for(int i= 0; i<myParticleSystem.size(); i++) {
            Particle* ptr_particle = myParticleSystem.getParticle(i);
            if (updateExcitation && ptr_particle->myEventData.onFloor) ptr_particle->myEventData.timeFreeFall = ptr_particle->getVelocity().length();
        }
        
        updateExcitation = false;
        updateStopMotion = false;
        
    }
    
}


//--------------------------------------------------------------
void ofApp::draw() {
    
    cam.begin();
    ofScale(1, -1, 1); // flip the y axis and zoom in a bit
    ofTranslate(-ofGetWidth()/2,-ofGetHeight()/2,0);
    
    ofPushMatrix();
    // scaleFactor=1.1;
    //ofScale(scaleFactor);
    // ofTranslate(-ofGetWidth()/2/scaleFactor, -ofGetHeight()/2/scaleFactor);
    //ofTranslate(100, 0);
    
    
    // Draw the confinement box:
    for(int i=0; i< particleArrayBorderFixed.size(); i++) {
        particleArrayBorderFixed.getParticle(i)->draw(false);
    }
    
    // Draw the mobile wall (if any):
    for(int i=0; i< particleArrayBorderMobile.size(); i++) {
        particleArrayBorderMobile.getParticle(i)->draw(false);
    }
    
    
    for(int i= 0; i<myParticleSystem.size(); i++) {
        Particle* ptr_particle = myParticleSystem.getParticle(i);
        ofVec2f pos = ptr_particle->getPosition();
        ofVec2f vel = ptr_particle->getVelocity();
        
        // ============= Draw the particle system =============
        if (drawParticles) {
            // (a) Choose color mode:
            ofFill();
            switch(modeColor) {
                case colorTangent:
                    ptr_particle->setColorTangentVelocity(posCenter);
                    break;
                case colorRadial:
                    ptr_particle->setColorRadialVelocity(posCenter);
                    break;
                case colorNorm:
                    ptr_particle->setColorNormVelocity();
                    break;
            }
            
            
            // Draw the particles, with or without forward arrow:
            if (addForcesRandom) ptr_particle->draw(true); else ptr_particle->draw(false);
            
            //  APPARENT size in order to represent something (like speed)
            if (sizeMode) {
                ofVec2f speed = ptr_particle->getVelocity();
                float newRadius = sizeParticle*ofMap( speed.length(), 0.0, 1.0, 1.0, 2.5, true);
                ofDrawCircle(pos, newRadius);
            }
        }
        
        
        // ============= Draw the DISPLACEMENT =============
        // **** TODO: make it using a polyline (and saving in file?)
        if (drawDisplacements) {
            
            ofVec2f equilPosition = ptr_particle -> getEquilibriumPosition();
            ofVec2f displacementVector = ptr_particle->getDisplacement();
            float angle = displacementVector.angle(ofVec2f(0,1));
            
            float dispNorm = displacementVector.length()*4;
            ofColor c = ofColor::fromHsb(dispNorm, 95+160.0*(170-dispNorm)/170, 255);//, speedNorm);
            ofSetColor(c);
            
            ofPushMatrix();
            ofVec3f equil3D = ofVec3f(equilPosition);
            ofTranslate(equil3D.x, equil3D.y,0);
            ofRotateYDeg(90);
            ofRotateZDeg(angle);
            ofDrawLine(0,0,dispNorm,0);//dispNorm*ur);
            ofTranslate(dispNorm,0,0);
            ofDrawSphere(0,0, sizeParticle);
            ofPopMatrix();
            
            //ofSetLineWidth(dispNorm/2);
            // ofDrawLine(equilPosition, pos);//dispNorm*ur);
            
            
            
            //            ofPushMatrix();
            //            ofTranslate(equil3D.x, equil3D.y, dispNorm);
            //            ofRotateXDeg(90);
            //            ofDrawCylinder(sizeParticle,2*dispNorm);
            //            ofPopMatrix();
            
            // as Z component:
            //             ofPushMatrix();
            //             ofTranslate(equil3D.x, equil3D.y, dispNorm);
            //             ofDrawCircle(0,0,sizeParticle);
            //             ofPopMatrix();
            
            
            //ofVec3f displacementZ3D= equil3D+dispNorm*ofVec3f(0,0,1);
            //ofDrawLine(equil3D, displacementZ3D);//dispNorm*ur);
            //ofDrawSphere(displacementZ3D, sizeParticle);
            
            
            //            if (modeCavity == disk) {
            //                // Represent the deviation from the initial position as an arrow towards the center (to see the wave):
            //                ofVec3f ur = 1.0*(equilPosition-posCenter).normalize();
            //                //ofVec3f velocity3D = ofVec3f(vel.x, vel.y, 0);
            //                //float tangentSpeed = abs(1000*ur.getCrossed( velocity3D ).z); // can be negative or positive
            //
            //                ofDrawLine(pos, pos-dispNorm*ur);
            //            }
            //            else if (modeCavity == rect) {
            //                // Represent the deviation from the initial position as an upward arrow t(to see the wave):
            //                ofVec2f ur = ofVec2f(0,-1.0);
            //                //ofDrawLine(equilPosition, equilPosition + (equilPosition-pos).length()*ur);
            //                ofDrawLine(equilPosition, equilPosition+dispNorm*ur);
            //                // ofDrawLine(pos, equilPosition);//dispNorm*ur);
            //            }
            
            
        }
        
        
        //        if(addForcesRandom) {
        //            particleArray[i].get()->draw(true);
        //        }
        //        else {
        //           particleArray[i].get()->draw(false);
        //        }
        
        
        
        
    }
    
    
    //ofSetHexColor(0x444342);
    ofSetColor(255,255,255);
    ofNoFill();
    ofSetLineWidth(2);
    for (int i=0; i<lines.size(); i++) {
        lines[i].draw();
    }
    
    ofSetLineWidth(1);
    ofNoFill();
    ofSetColor(255,255,255,255);
    ofDrawCircle(posCenter, blackHoleRadius);
    if (blackHole) {
        ofFill();
        ofSetColor(255,0,0,100);
        ofDrawCircle(posCenter, blackHoleRadius);
    }
    
    // show circle around current mouse position:
    ofSetLineWidth(.5);
    ofNoFill();
    ofSetColor(255,255,255,255);
    ofDrawCircle(mouseX, mouseY, blackHoleRadius);
    
    /*
     for (int i=0; i<edges.size(); i++) {
     edges[i].get()->draw();
     }
     */
    
    // draw the ground
    //box2d.drawGround();
    
    ofPopMatrix();
    
    
    cam.end();
    
    // Draw the GUI:
    gui.draw();
    // Draw additional data
    if(!bHide){
        string info = "";
        info += "[.] - add one particle under cursor\n";
        info += "[c] - add one particle at center\n";
        info += "[n] -  particle under cursor\n";
        info += "[/] - add 10 particles under cursor\n";
        info += "[,] - delete last particle\n\n";
        
        info += "[i] - Impulse on particle under cursor\n\n";
        
        info += "[o] - save equilibrium position\n";
        info += "[p] - reset to equilibrium position\n";
        info += "[w] - show displacement from equilibrium \n\n";
        
        info += "[1] - color mode norm velocity\n";
        info += "[2] - color mode norm tengential vel\n";
        info += "[3] - color mode norm radial vel\n\n";
        //        info += "BOUNDARIES:\n";
        //        info += "* Use mouse drag to create boundaries\n";
        //        info += "[0] - delete last boundary\n\n";
        //        info += "VISUALIZATION:\n";
        //        info += "[</>] - scale image\n";
        info += "Total particles: "+ofToString(box2d.getBodyCount())+"\n";
        // info += "Total Joints: "+ofToString(box2d.getJointCount())+"\n\n";
        info += "FPS: "+ofToString(ofGetFrameRate(), 1)+"\n";
        ofSetColor(200, 200, 0);
        ofDrawBitmapString(info, ofGetWidth()-350, ofGetHeight()-250);
        
    }
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
    
    if(key == 'h'){bHide = !bHide;}
    else if(key == 's'){ gui.saveToFile("settings.xml");}
    else if(key == 'l'){ gui.loadFromFile("settings.xml");}
    
    else if(key == '1') modeColor = colorNorm;
    else if(key == '2') modeColor = colorTangent;
    else if(key == '3') modeColor = colorRadial;
    
    else if (key == 'f') {addForcesRepulsion=!addForcesRepulsion;}
    else if (key == 'd') {addForcesRandom=!addForcesRandom;}
    
    else if(key == 'c') addParticle(posCenter, ofVec2f(0,0), massParticle/FACTOR_MASS, charge, sizeParticle, dampingFactor, myParticleSystem);
    
    else if(key == '.') addParticle(ofVec2f(mouseX, mouseY), ofVec2f(0,0), massParticle/FACTOR_MASS, charge, sizeParticle, dampingFactor, myParticleSystem);
    //addParticle(ofVec2f(mouseX+ofRandom(0.5*sizeParticle), mouseY+ofRandom(0.5*sizeParticle)), ofVec2f(0,0), massParticle/FACTOR_MASS, charge, sizeParticle, dampingFactor, myParticleSystem);
    
    else if(key == '/') { // add several particles at once:
        for (int i=0; i<5; i++) addParticle(ofVec2f(mouseX, mouseY), ofVec2f(0,0), massParticle/FACTOR_MASS, charge, sizeParticle, dampingFactor, myParticleSystem);
    }
    else if(key == ',') deleteLastParticle(myParticleSystem);
    
    // Delete particle near mouse:
    else if (key=='n') {
        for(int i=0; i<myParticleSystem.size(); i++) {
            ofVec3f pos = myParticleSystem.getParticle(i)->getPosition();
            ofPoint T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
            ofPoint MousePos = (ofPoint(mouseX, mouseY)- T)/scaleFactor ;
            float distMouse = (MousePos-pos).length();
            if (distMouse<=sizeParticle)
            {
                myParticleSystem.erase(i);
                i=i-1;
            }
        }
    }
    
    else if (key == 'm') {
        // make the particle fixed (mass = 0):
        for(int i=0; i<myParticleSystem.size(); i++) {
            ofVec2f impulse = ofVec2f(10,0);
            Particle* ptr_particle = myParticleSystem.getParticle(i);
            ofVec3f pos = ptr_particle->getPosition();
            ofPoint T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
            ofPoint MousePos = (ofPoint(mouseX, mouseY)- T)/scaleFactor ;
            float distMouse = (MousePos-pos).length();
            if (distMouse<=sizeParticle)
            {
                ptr_particle->setMass(0);
            }
        }
    }
    else if (key == 'M') {
        // make the particle fixed (mass = 0):
        for(int i=0; i<myParticleSystem.size(); i++) {
            ofVec2f impulse = ofVec2f(10,0);
            Particle* ptr_particle = myParticleSystem.getParticle(i);
            ofVec3f pos = ptr_particle->getPosition();
            ofPoint T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
            ofPoint MousePos = (ofPoint(mouseX, mouseY)- T)/scaleFactor ;
            float distMouse = (MousePos-pos).length();
            if (distMouse<=sizeParticle)
            {
                ptr_particle->setMass(massParticle/FACTOR_MASS);
            }
        }
    }
    
    
    else if (key=='i') {
        // Give a kick to the right to a particle near the mouse:
        for(int i=0; i<myParticleSystem.size(); i++) {
            ofVec2f impulse = ofVec2f(10,0);
            Particle* ptr_particle = myParticleSystem.getParticle(i);
            ofVec3f pos = ptr_particle->getPosition();
            ofPoint T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
            ofPoint MousePos = (ofPoint(mouseX, mouseY)- T)/scaleFactor ;
            float distMouse = (MousePos-pos).length();
            if (distMouse<=sizeParticle)
            {
                ptr_particle->addImpulseForce(ofVec2f(0,0), impulse/1000);
            }
        }
    }
    
    else if (key == 'o') resetToEquilibriumPos();
    else if (key == 'p') saveEquilibriumPos();
    
    
    else if(key == 't') ofToggleFullscreen();
    else if(key == ' ') continueAddingLines=true; //!continueAdding;
    // if(key == '.') continueAddingCenter=true; //!continueAdding;
    
    //     else if (key == 'a') {gravity += 10; box2d.setGravity(0, gravity);}
    //     else if (key == 'q') {gravity -= 10; box2d.setGravity(0, gravity);}
    //     else if (key == 'w') amplitude *= 1.5;
    //     else if (key == 's') {amplitude /= 1.5; if (amplitude<0) amplitude = 0;}
    //
    
    
    
    else if (key == '>') scaleFactor *= 1.1;
    else if (key == '<') scaleFactor /= 1.1;
    
    
    //     else if (key==OF_KEY_RIGHT){ scatterAngle+=1;if (scatterAngle>180) scatterAngle=180;}
    //     else if (key==OF_KEY_LEFT) {scatterAngle-=1; if (scatterAngle<0) scatterAngle=0;}
    //
    //     else if (key==OF_KEY_UP) periodMillisExcitation+=1;
    //     else if (key==OF_KEY_DOWN) {periodMillisExcitation-=1; if (periodMillisExcitation<0) periodMillisExcitation=0;}
    //
    //     else if (key=='1') {
    //        dampingFactor*=1.2; if (dampingFactor>100) dampingFactor =100;
    //          for(int i=0; i<particleArray.size(); i++)
    //              particleArray[i].get()->setDamping(dampingFactor);
    //    }
    //     else if (key=='2') {
    //        dampingFactor/=1.2; if (dampingFactor<0.1) dampingFactor =0.1;
    //        for(int i=0; i<particleArray.size(); i++)
    //            particleArray[i].get()->setDamping(dampingFactor);
    //    }
    
    //     else if (key=='5') {freeJumpMode=!freeJumpMode;}
    //     else if (key=='3') periodMillisStopMotion+=1;
    //     else if (key=='4') {periodMillisStopMotion-=1; if (periodMillisStopMotion<0) periodMillisStopMotion=0;}
    //
    
    else if (key=='\\'){ //OF_KEY_BACKSPACE) {
        deleteLastBorder();
    }
    
    else if (key==';') {
        radiusContainer-=sizeParticle/2;
        deleteLastBorder();
        makeContainerDisk(ofGetWidth()/2, ofGetHeight()/2,radiusContainer, 30);
    }
    else if (key==':') {
        radiusContainer+=sizeParticle/2;
        deleteLastBorder();
        makeContainerDisk(ofGetWidth()/2, ofGetHeight()/2,radiusContainer, 30);
    }
    
    else if (key=='w') {
        drawDisplacements=!drawDisplacements;
    }
    
    //    else if (key == 'p') displacement += 10;
    //    else if (key == 'o') displacement -= 10;
    
    
}


void ofApp::clearLogData() {
    myParticleSystem.clearLog();
}

void ofApp::savePhaseSpaceTrajectory() {
    
    //ofFileDialogResult result = ofSystemLoadDialog("Select file to save", true);
    ofFileDialogResult saveFileResult = ofSystemSaveDialog("Simulation" + ofGetTimestampString() + ".txt", "Saving position & time stamp for all particles");
    if (saveFileResult.bSuccess) {
        myParticleSystem.savePhaseSpaceTrajectory(saveFileResult.filePath);
        ofLogVerbose("File saved");
    }
    else ofLogVerbose("Saving failed!");
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {
    
}

void ofApp::makeContainerDiskFixedMasses(const ofVec2f& _center, float _radius, uint16_t _numMasses, ParticleSystem& _particleArrayBorderFixed) {
    float sectionAngle = 2.0*PI/_numMasses;
    for (int k=0; k<_numMasses; k++) {
        //float x = cx + 0.9*radius * cos(k*sectionAngle);
        //float y = cy + 1.1*radius * sin(k*sectionAngle);
        
        ofVec2f pos = _center + _radius*ofVec2f(cos(k*sectionAngle), sin(k*sectionAngle));
        
        //  void addParticle(const float x, const float y, const float _mass, const ofVec2f _vel, const float _damp, const float _size);
        // mass = 0 => fixed body! Also, charge should be positive (for repulsion, see update)
        addParticle(pos, ofVec2f(0,0), 0, 1.0, 5, dampingFactor, _particleArrayBorderFixed);
        _particleArrayBorderFixed.back()->setColor(ofColor(255,255,255));
    }
}



void ofApp::makeContainerDisk(float cx, float cy, float radius, int numPoints) {
    
    lines.push_back(ofPolyline());
    
    float sectionAngle = 2*PI/numPoints;
    for (int k=0; k<numPoints+1; k++) {
        //float x = cx + 0.9*radius * cos(k*sectionAngle);
        //float y = cy + 1.1*radius * sin(k*sectionAngle);
        float x = cx + radius * cos(k*sectionAngle);
        float y = cy + radius * sin(k*sectionAngle);
        lines.back().addVertex(x, y);
    }
    
    auto edge = std::make_shared<ofxBox2dEdge>();
    
    lines.back().simplify();
    
    for (auto i=0; i<lines.back().size(); i++) {
        edge.get()->addVertex(lines.back()[i]);
    }
    
    edge.get()->create(box2d.getWorld());
    //setPhysics(float density, float bounce, float friction)
    edge.get()->setPhysics(-1, 1.0, 0.0);
    
    edges.push_back(edge);
    
}

void ofApp::deleteLastBorder() {
    if (lines.size()>0) {
        lines.pop_back();
        edges.pop_back();
    }
}

void ofApp::makeContainerRect(const ofRectangle & _rect) {
    
    
    lines.push_back(ofPolyline::fromRectangle(_rect));
    
    auto edge = std::make_shared<ofxBox2dEdge>();
    
    // lines.back().simplify();
    
    for (auto i=0; i<lines.back().size(); i++) {
        edge.get()->addVertex(lines.back()[i]);
    }
    // close the rectangle (strange):
    edge.get()->addVertex(lines.back()[0]);
    
    
    edge.get()->create(box2d.getWorld());
    
    //edge.get()->setPhysics(0.0011, 1.0, 0.0);  // uncomment this to see it fall!
    edges.push_back(edge);
    
}

void ofApp::makeContainerFixedRect(const ofRectangle & _rect, uint16_t _numMasses, ParticleSystem& _particleArrayBorderFixed) {
    
    float perimeter = _rect.getPerimeter();
    float stepMass = perimeter/_numMasses;
    
    // Create horizontal sides:
    int numMassWidth = glm::ceil(_rect.getWidth()/stepMass);
    float numMassHeight = glm::ceil(_rect.getHeight()/stepMass);
    float correctedRight = _rect.getLeft() + stepMass * numMassWidth;
    
    //Create vertical sides:
    //    // IMPORTANT: the FIRST index (0) is for the vertical side on the left, for all the masses!!!
    //    for (int i=0; i<numMassHeight; i++) { // 1 and -1 (not to duplicate the corner masses)
    //        //  void addParticle(const float x, const float y, const float _mass, const ofVec2f _vel, const float _damp, const float _size);
    //        // mass = 0 => fixed body!
    //        addParticle(_rect.getLeft(), _rect.getBottom() - i*stepMass, 0, ofVec2f(0,0), dampingFactor, stepMass/2.2, _particleArrayBorderFixed);
    //        _particleArrayBorderFixed.back()->setColor(ofColor(255,255,255));
    //    }
    
    //Mass should be 0 and charge should be positive (for repulsion, see update)
    for (int i=0; i<numMassHeight; i++) {
        addParticle(ofVec2f(correctedRight, _rect.getBottom()  - i*stepMass), ofVec2f(0,0), 0, 1.0, stepMass/2.2, dampingFactor, _particleArrayBorderFixed); // mass = 0 => fixed body
        _particleArrayBorderFixed.back()->setColor(ofColor(255,255,255));
    }
    
    for (int i=1; i<numMassWidth; i++) {
        addParticle(ofVec2f(correctedRight -i*stepMass , _rect.getTop()) , ofVec2f(0,0), 0, 1.0, stepMass/2.2, dampingFactor, _particleArrayBorderFixed);// mass = 0 => fixed body
        _particleArrayBorderFixed.back()->setColor(ofColor(255,255,255));
        
        addParticle(ofVec2f(correctedRight-i*stepMass , _rect.getBottom()), ofVec2f(0,0), 0, 1.0, stepMass/2.2, dampingFactor, _particleArrayBorderFixed); // mass = 0 => fixed body
        _particleArrayBorderFixed.back()->setColor(ofColor(255,255,255));
        
    }
    
    //return(numMassHeight);
}


void ofApp::makeMobileWall(const ofVec2f _A, const ofVec2f _B,  uint16_t _numMasses, ParticleSystem& _mobileWallParticles) {
    ofVec2f u =(_B-_A);
    float length = u.length();
    u/=length;
    float step = length/_numMasses;
    ofVec2f newParticlePos = _A;
    for (float dist = -step; dist <= length ; dist += step) {
        //Mass should be 0 and charge should be positive (for repulsion, see update)
        addParticle(newParticlePos, ofVec2f(0,0), 0, 1.0, step/2.2, dampingFactor, _mobileWallParticles);
        _mobileWallParticles.back()->setColor(ofColor(255,0,255));
        newParticlePos+=step*u;
    }
}


void ofApp::translateParticleSystem(const ofVec2f& _translation, ParticleSystem& _mobileWallParticles) {
    // IMPORTANT: the FIRST index (0) is for the vertical side on the left, for all the masses!!!
    for (int i=0; i<_mobileWallParticles.size(); i++) {
        Particle* particle = _mobileWallParticles.getParticle(i);
        ofVec2f newPos = particle->myData.posInit + _translation;
        particle->setPosition(newPos.x, newPos.y);
    }
    
}


void ofApp::makeContainer() {
    
    lines.push_back(ofPolyline());
    lines.back().addVertex(ofGetWidth()-10, ofGetHeight()-10);
    lines.back().addVertex(10, ofGetHeight()-10);
    lines.back().addVertex(10,10);
    lines.back().addVertex(ofGetWidth()-10, 10);
    lines.back().addVertex(ofGetWidth()-10, ofGetHeight()-10);
    
    auto edge = std::make_shared<ofxBox2dEdge>();
    
    lines.back().simplify();
    
    for (auto i=0; i<lines.back().size(); i++) {
        edge.get()->addVertex(lines.back()[i]);
    }
    
    edge.get()->create(box2d.getWorld());
    edge.get()->setPhysics(1, 1.0,0.0);  // uncomment this to see it fall!
    edges.push_back(edge);
    
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ) {}
void ofApp::mousePressed(int x, int y, int button) {}
void ofApp::mouseDragged(int x, int y, int button) {}
void ofApp::mouseReleased(int x, int y, int button) {}

//void ofApp::mousePressed(int x, int y, int button) {
//addParticle(x, y, massParticle/FACTOR_MASS, sizeParticle, 1);

//    ofPoint T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
//    ofPoint X = (ofPoint(x,y)- T)/scaleFactor ;
//
//    lines.push_back(ofPolyline());
//    lines.back().addVertex(X);//x*scaleFactor, y*scaleFactor);
//}

//--------------------------------------------------------------
//void ofApp::mouseDragged(int x, int y, int button) {
//    ofPoint T(ofGetWidth()/2*(1.0-scaleFactor), ofGetHeight()/2*(1.0-scaleFactor));
//    ofPoint X = (ofPoint(x,y)- T)/scaleFactor ;
//    lines.back().addVertex(X);
//}

//--------------------------------------------------------------
//void ofApp::mouseReleased(int x, int y, int button) {
//
//    auto edge = std::make_shared<ofxBox2dEdge>();
//    lines.back().simplify();
//
//    for (auto i=0; i<lines.back().size(); i++) {
//        edge.get()->addVertex(lines.back()[i]);
//    }
//
//    //poly.setPhysics(1, .2, 1);  // uncomment this to see it fall!
//    edge.get()->create(box2d.getWorld());
//    edges.push_back(edge);
//
//    //lines.clear();
//}



//--------------------------------------------------------------
void ofApp::resized(int w, int h){
    
}


void ofApp::frequenyUpdateChange(float &_freqUpdate){
    box2d.setFPS(_freqUpdate);
}

void ofApp::gravityChange(float &_gravity) {
    box2d.setGravity(-_gravity, 0);
}

void ofApp::massUpdateChange(float & _mass) {
    // note: to change the mass, the world should be unlocked (not in the middle of a time step)
    //massParticle = _mass;
    flagChangeMass = true;
}

void ofApp::chargeUpdateChange(float & _charge) {
    for(int i=0; i<myParticleSystem.size(); i++)
        myParticleSystem.getParticle(i)->myData.charge=_charge;
}

void ofApp::changeDampingFactor(float & _dampingFactor) {
    for(int i=0; i<myParticleSystem.size(); i++)
        myParticleSystem.getParticle(i)->setDamping(_dampingFactor);
    //dampingFactor = _dampingFactor;
}

void ofApp::restartSystem() {
    myParticleSystem.clearSystem();
    particleArrayBorderFixed.clearSystem();
    particleArrayBorderMobile.clearSystem();
    createCavity();
    //    float phase = 2.0*PI*freqOscillateSide*(ofGetElapsedTimef() - timeOffset);
    //    ofVec2f translation = ofVec2f( ampOscillateSide*(sin(phase)), 0);
    //    translateParticleSystem(translation, particleArrayBorderMobile);
}

void ofApp::resetSpeedAll() {
    for(int i=0; i<myParticleSystem.size(); i++)
        myParticleSystem.getParticle(i)->setVelocity(0,0);
}


void ofApp::timeReversalAll() {
    for(int i=0; i<myParticleSystem.size(); i++) {
        Particle* particle = myParticleSystem.getParticle(i);
        ofVec2f vel =  -1.0*particle->getVelocity();
        particle->setVelocity(vel);
    }
}


void ofApp::resetMobileWall() {
    particleArrayBorderMobile.resetToEquilibriumPos();
}


void ofApp::resetToEquilibriumPos() {
    // Store current position of particles as "equilibrium" positions:
    myParticleSystem.resetToEquilibriumPos();
}

void ofApp::saveEquilibriumPos() {
    // Store current position of particles as "equilibrium" positions:
    myParticleSystem.saveEquilibriumPos();
}


void ofApp::saveConfigurationMemory() {
    myParticleSystem.saveSnapshotStateMemory();
}

void ofApp::loadConfigurationMemory(){
    myParticleSystem.loadSnapshotStateMemory();
}


void ofApp::loadConfiguration() {
    ofFileDialogResult loadFileState = ofSystemLoadDialog("Load state snapshot of the system");
    if (loadFileState.bSuccess) {
        myParticleSystem.loadSnapshotState(loadFileState.filePath);
        ofLogVerbose("State load");
    }
    else ofLogVerbose("Loading failed!");
}

void ofApp::saveConfiguration() {
    ofFileDialogResult saveFileResult = ofSystemSaveDialog("Snapshot" + ofGetTimestampString() + ".txt", "Saving a snapshot of the system");
    if (saveFileResult.bSuccess) {
        myParticleSystem.saveSnapshotState(saveFileResult.filePath);
        ofLogVerbose("File saved");
    }
    else ofLogVerbose("Saving failed!");
}

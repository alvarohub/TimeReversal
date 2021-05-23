//
//  particleClass.hpp
//  Physics_janusDaemons
//
//  Created by Alvaro Cassinelli on 7/3/2021.
//

#pragma once
#include "ofMain.h"
#include "ofxBox2d.h"
#include "ofxCsv.h"

#define DEFAULT_PERIOD_LOG 0.001

class Particle : public ofxBox2dCircle {
    
public:
    
    struct ParticleData {
        //NOTE: data may be duplicated here (because Box2d particle also has these fields, but here
        // we can integrate all of them for simpler access (saving, loading).
        // ... TODO
        ofVec2f posInit;// = ofVec2f(0,0);
        ofVec2f velInit = ofVec2f(0,0);
        ofVec2f equilibriumPosition = posInit;
        
        float mass = 1.0;
        float charge = 1.0;
        float size = 10;
        float damp = 10;
    };
    
    
    struct EventTimeData {
        float lastTimeExcitation; // to compute the time in free fall in parabolic jump
        float timeFreeFall;
        bool onFloor;
    };
    
    struct PhaseSpaceData {
        // float time; // evolution time is common to all particles
        ofVec2f pos;
        ofVec2f velocity;
        ofVec2f displacement; // displacement from equilibrium
        // ofVec2f acc;
    };
    
    Particle() { };
    Particle(const ofVec2f& _posInit, const ofVec2f& _posEquilibrium, const ofVec2f& _velInit, const float _mass, const float _charge, const float _size, const float _damp ) {
        myData.posInit = _posInit;
        myData.equilibriumPosition = _posEquilibrium;
        myData.velInit = _velInit;
        myData.mass = _mass;
        myData.size = _charge;
        myData.damp = _damp;
        myData.size = _size;
    };
    // Creation of a particle from data with type ParticleData:
    // NOTE: this is called when doing instantiation like this: Particel newParticle = oldParticle(data)
    Particle(const ParticleData &_data ) {
        myData = _data;
    };
    
    //Copy constructor:
    // NOTE:  it's provided automatically by the compiler and it's fine when there isn't pointers pointing to
    // content needed to be (deep) copied. And we are not interested in copying the history of the phase space.
    Particle(const Particle &_old) { // note: this is called when doing instantiation like this: Partice newParticle = oldParticle
        myData = _old.myData;
        myEventData = _old.myEventData;
    } // this copy constructor will be called with assigning ex: Particle newP = oldP;
    
    // Assignement operator:
    // note: called as: particleA = particleB. Actually same than particleA.operator=(particleB);
    // If not declared, the compiler gives you one implicitly (read copy constructor remarks above)
    Particle& operator = (const Particle &_otherParticle)
    {
        //cout<<"Assignment operator called "<<endl;
        myData = _otherParticle.myData;
        myEventData = _otherParticle.myEventData;
        return *this;
    }
    
    // FORCES ------------------------------
    void addPointForceConst(ofVec2f pt, float _force);
    // Note: the charge will change signs on the following forces:
    void addPointForceSpring(ofVec2f pt, float _k, float _relaxDist);
    void addPointForce_InvR(ofVec2f pt,  float _amplitude);
    void addPointForce_InvR2(ofVec2f pt, float _amplitude);
    void addPointForce_InvR3(ofVec2f pt, float _amplitude);
    void addPointForce_InvR4(ofVec2f pt, float _amplitude);
    void addMagneticForce(ofVec3f _magneticField);
    
    void setColorNormVelocity();
    void setColorTangentVelocity(ofVec2f _center);
    void setColorRadialVelocity(ofVec2f _center);
    
    // ====== Phase Space Log  ==================
    
    PhaseSpaceData getPhaseSpacePoint(uint32_t _index) {return phaseSpaceTrajectory[_index];}
    const PhaseSpaceData& getLastSpacePoint() {return phaseSpaceTrajectory.back();}
    
    void resetLog() {
        phaseSpaceTrajectory.clear();
    }
    
    void phaseSpaceLog() { //float _time) {
        PhaseSpaceData _newPoint;
        _newPoint.pos = getPosition();
        _newPoint.velocity = getVelocity();
        _newPoint.displacement = _newPoint.pos - myData.equilibriumPosition;
        phaseSpaceTrajectory.push_back(_newPoint);
    }
    
    void setEventData(float _lastTimeExcitation, float _timeFreeFall, bool _onFloor) {
        myEventData.lastTimeExcitation = _lastTimeExcitation;
        myEventData.timeFreeFall = _timeFreeFall;
        myEventData.onFloor = _onFloor;
    }
    void setEventData(const EventTimeData& _newEventData) {myEventData = _newEventData;};
    const EventTimeData& getEventData() {return myEventData;}
    
    ParticleData& getParticleData() {return myData;}
    
    const ofVec2f& getEquilibriumPosition() {return myData.equilibriumPosition;}
    void setEquilibriumPosition(const ofVec2f& _pos) {
        myData.equilibriumPosition = _pos;
    }
    
    const ofVec2f getDisplacement() {return (getPosition() - myData.equilibriumPosition);}
    
    
    ParticleData myData;
    EventTimeData myEventData;
    
private:
    
    vector<PhaseSpaceData> phaseSpaceTrajectory;
    
};

// ======================================================================

class ParticleSystem {
    
public:
    
    ParticleSystem() {}
    
    void init(ofxBox2d* _ptr_box2d) {
        ptr_box2d = _ptr_box2d;
        resetTimeStamp();
        timeStamp.clear();
        lastLogTime = ofGetElapsedTimef();
        periodLog = DEFAULT_PERIOD_LOG; // in seconds, note: this is the minimum, calculation can take longer but the time stamp will be correct
    }
    
    void setPeriodLog(float _periodLog) {periodLog = _periodLog;}
    
    void setSystemWorld(ofxBox2d* _box2d) { ptr_box2d = _box2d; };
    
    void addParticle(Particle& _particle);
    void addParticle(const Particle::ParticleData& _particleData);
    void addParticle(const ofVec2f& _pos, const ofVec2f& _equil, const ofVec2f& _vel, float _mass, float _charge, float _size, float _damp);
    void addParticle(float x, float y, float vx, float vy, float _ex, float _ey, float _mass, float _charge, float _size, float _damp);
    
    Particle* getParticle(uint32_t _index) {return particleArray[_index].get();}
    Particle* back() {return particleArray.back().get();}
    
    void deleteLastParticle() {if (!particleArray.empty()) particleArray.pop_back();}
    void erase(uint16_t _index) {if (!particleArray.empty()) particleArray.erase(particleArray.begin()+_index);}
    
    uint16_t size() {return particleArray.size();}
    
    ofVec2f getPosition(uint32_t _index) { return particleArray[_index].get()->getPosition();}
    
    void clearSystem() {
        particleArray.clear();
        timeStamp.clear();
    }
    
    void clearLog() {
        for (auto particle : particleArray) particle->resetLog();
        timeStamp.clear();
        resetTimeStamp();
    }
    
    void resetTimeStamp() {ofResetElapsedTimeCounter();}
    
    void updateLog() { // only if active, and with the correct sampling rate:
        if (ofGetElapsedTimef() - lastLogTime > periodLog) {
            for (auto particle : particleArray) particle->phaseSpaceLog();
            timeStamp.push_back(lastLogTime);
            lastLogTime = ofGetElapsedTimef();
        }
    }
    
    void savePhaseSpaceTrajectory(const string& _nameFile) {
        csvRecorder.clear();
        ofxCsvRow row;
        
        for (uint32_t i = 0 ; i< timeStamp.size() ; i++) {
            row.clear();
            row.setFloat(0, timeStamp[i]);
            for (uint16_t j=0; j < particleArray.size(); j++) {
                //  cout << j << endl;
                Particle::PhaseSpaceData data = particleArray[j].get()->getPhaseSpacePoint(i);
                row.setFloat(7*j+1, j); // particle index
                row.setFloat(7*j+2, data.pos.x);
                row.setFloat(7*j+3, data.pos.y);
                row.setFloat(7*j+4, data.velocity.x);
                row.setFloat(7*j+5, data.velocity.y);
                row.setFloat(7*j+6, data.displacement.x);
                row.setFloat(7*j+7, data.displacement.y);
            }
            csvRecorder.addRow(row);
        }
        
        // finally, save to file:
        csvRecorder.save(_nameFile);
    }
    
    void saveSnapshotStateMemory() {
        particleDataArray.clear(); // note: the phase space trajectory is not recorded
        for (auto particle : particleArray) {
            Particle::ParticleData dataParticle = particle.get()->myData;
            // NOTE: the position should be the current one, not the equilibrium one nor the initial one (which we can either keep, or
            // reset here to the present position. I will do the last, if not one has to set the position explicitly afterwards (using setPosition(...))
            dataParticle.posInit = particle->getPosition();
            particleDataArray.push_back(dataParticle);
        }
    }
    
    void loadSnapshotStateMemory() {
        clearSystem();
        for (auto dataParticle : particleDataArray) {
            Particle newParticle(dataParticle);
            addParticle(newParticle);
        }
    }
    
    void loadSnapshotState(const string& _nameFile) {
        // Load a CSV File.
        if(csv.load(_nameFile)) {
            
            // ofLog() << "Print out the first value";
            // ofLog() << csv.getRow(0).getString(0);
            // Print the table to the console.
            // ofLog() << "Print the table";
            //csv.print(); // Uses default separator ",".
            // ... or do it manually
            //    for(auto row : csv) {
            //        ofLog() << ofJoinString(row, "|");
            //    }
            
            // (a) First clear everything:
            clearSystem();
            
            // (b) Instantiate and initialize state of all the new particles.
            // NOTE: each row in the state csv file correspond to data for one particle, with an index for clarity.
            // The data correspond to the parameters of the method addParticle:
            //      addParticle(float x, float y, float _mass, ofVec2f _vel, float _damp, float _size)
            // Example:
            //          numrow (0...), x, y, mass, vx, vy, damping, size
            
            cout << "Loading " + ofToString(csv.getNumRows()) + " particles." << endl;
            
            for(auto row : csv) {
                // ofLog() << ofJoinString(row, "|");
                // for (int i=0; i<row.getNumCols();i++) ...
                cout << " Loading particle " + ofToString(row.getInt(0)) << endl;
                ofVec2f pos = ofVec2f(row.getFloat(1), row.getFloat(2));
                ofVec2f equil = ofVec2f(row.getFloat(3), row.getFloat(4));
                ofVec2f vel = ofVec2f(row.getFloat(5), row.getFloat(6));
                float mass = row.getFloat(7);
                float charge = row.getFloat(8);
                float size = row.getFloat(9);
                float damp = row.getFloat(10);
                Particle newParticle = Particle(pos, equil, vel, mass, charge, size, damp);
                addParticle(newParticle);
            }
            
        }
    }
    
    
    void saveSnapshotState(const string& _nameFile) {
        
        csvRecorder.clear();
        ofxCsvRow row;
        for (uint32_t i = 0 ; i< particleArray.size() ; i++) {
            row.clear();
            
            Particle* ptr_particle = particleArray[i].get();
            
            //note: pos and vel are not from the particle data (that is the init conditions)
            ofVec2f pos = ptr_particle->getPosition();
            ofVec2f vel = ptr_particle->getVelocity();
            
            ofVec2f equil = ptr_particle->myData.equilibriumPosition;
            
            row.setInt(0, i); // particle index (not necessary, but for clarity)
            row.setFloat(1, pos.x);  row.setFloat(2, pos.y);
            row.setFloat(3, equil.x);  row.setFloat(4, equil.y);
            row.setFloat(5, vel.x);  row.setFloat(6, vel.y);
            row.setFloat(7, ptr_particle->myData.mass);
            row.setFloat(8, ptr_particle->myData.charge);
            row.setFloat(9, ptr_particle->myData.size);
            row.setFloat(10, ptr_particle->myData.damp);
            
            csvRecorder.addRow(row);
        }
        
        // finally, save to file:
        csvRecorder.save(_nameFile);
        
    }
    
    void saveEquilibriumPos() {
        for (uint16_t j=0; j < particleArray.size(); j++) {
            Particle* ptr_particle = particleArray[j].get();
            ptr_particle->setEquilibriumPosition(ptr_particle->getPosition());
        }
    }
    
    void resetToEquilibriumPos() { // note: we reset AND set the speeds to ZERO!
        for (uint16_t j=0; j < particleArray.size(); j++) {
            Particle* ptr_particle = particleArray[j].get();
            ptr_particle->setPosition(ptr_particle->getEquilibriumPosition());
        }
    }
    
    
    ofxBox2d* ptr_box2d;
    vector< ofPtr<Particle> > particleArray;
    vector<Particle::ParticleData> particleDataArray;
    
private:
    float lastLogTime, periodLog;
    vector <float> timeStamp;
    uint32_t numLogs;
    ofxCsv csv;
    ofxCsv csvRecorder;
};

//
//  particleClass.cpp
//  Physics_janusDaemons
//
//  Created by Alvaro Cassinelli on 7/3/2021.
//

#include "particleClass.hpp"

void Particle::addPointForceConst(ofVec2f pt, float _force) {
    if(body != NULL) {
        b2Vec2 P(pt.x/OFX_BOX2D_SCALE, pt.y/OFX_BOX2D_SCALE);
        b2Vec2 D = P - body->GetPosition();
        D.Normalize();
        b2Vec2 F = _force*D;
   //     F.y =0;
        body->ApplyForceToCenter(-F,true);
        
    }
}

void Particle:: addPointForceSpring(ofVec2f pt, float _k, float _relaxDist) {
    if(body != NULL) {
        float scaledRelaxDist = _relaxDist / OFX_BOX2D_SCALE;
        b2Vec2 P(pt.x/OFX_BOX2D_SCALE, pt.y/OFX_BOX2D_SCALE);
        b2Vec2 D = P - body->GetPosition();
        float dist = D.Length();
        if (dist<2*scaledRelaxDist) {
            D*=1.0/dist; // D.Normalize();
            b2Vec2 F = _k*(_relaxDist - dist)*D;
    //        F.y =0;
            body->ApplyForceToCenter(-F,true);
        }
    }
}

void Particle::addPointForce_InvR(ofVec2f pt,  float _amplitude) {
    if(body != NULL) {
        b2Vec2 P(pt.x/OFX_BOX2D_SCALE, pt.y/OFX_BOX2D_SCALE);
        b2Vec2 D = P - body->GetPosition();
        float distSq = D.LengthSquared();
        //D*=1.0/dist;
        //D.Normalize();
        b2Vec2 F = myData.charge*_amplitude/(distSq+0.00001)*D;
        body->ApplyForceToCenter(-F,true);
        
    }
}

void Particle::addPointForce_InvR2(ofVec2f pt, float _amplitude) {
    if(body != NULL) {
        b2Vec2 P(pt.x/OFX_BOX2D_SCALE, pt.y/OFX_BOX2D_SCALE);
        b2Vec2 D = P - body->GetPosition();
        float dist = D.Length();
        //D.Normalize();
        b2Vec2 F = myData.charge*_amplitude/(dist*dist*dist+0.0000001)*D;
        body->ApplyForceToCenter(-F,true);
        
    }
}

void Particle::addPointForce_InvR3(ofVec2f pt, float _amplitude) {
    
    if(body != NULL) {
        b2Vec2 P(pt.x/OFX_BOX2D_SCALE, pt.y/OFX_BOX2D_SCALE);
        b2Vec2 D = P - body->GetPosition();
        float distSq = D.LengthSquared();
        //D*=1.0/dist;
        //D.Normalize();
        b2Vec2 F = myData.charge*_amplitude/(distSq*distSq+0.0000001)*D;
        
     //   F.y =0;
        body->ApplyForceToCenter(-F,true);
        
    }
    
}


void Particle::addPointForce_InvR4(ofVec2f pt, float _amplitude) {
    if(body != NULL) {
        b2Vec2 P(pt.x/OFX_BOX2D_SCALE, pt.y/OFX_BOX2D_SCALE);
        b2Vec2 D = P - body->GetPosition();
        float dist = D.Length(); // avoid computing the squareroot if possible...
        // if (dist<20*getRadius()) {
        //D.Normalize();
        b2Vec2 F = myData.charge*_amplitude/(dist*dist*dist*dist*dist+0.00000001)*D;
    //    F.y =0;
        body->ApplyForceToCenter(-F,true);
        //}
    }
}

void Particle::addMagneticForce(ofVec3f _magneticField) {
    if(body != NULL) {
        // Compute lorentz force (assuming charge constant and positive):
        b2Vec2 vel = body->GetLinearVelocity();
        ofVec2f lorenzForceOf = myData.charge*ofVec3f(vel.x, vel.y, 0).cross(_magneticField);
        b2Vec2 lorentzForce(lorenzForceOf.x,lorenzForceOf.y);
        body->ApplyForceToCenter(lorentzForce, true);
    }
}


void Particle::setColorNormVelocity() {
    ofVec2f force = getForce();
    ofVec2f speed = getVelocity();
    float speedNorm = ofMap( speed.length(), 0, 0.3, 170,0, true);
    //cout <<speed.length()<< endl;
    ofColor c = ofColor::fromHsb(speedNorm, 95+160.0*(170-speedNorm)/170, 255);//, speedNorm);
    setColor(c);
}

void Particle::setColorTangentVelocity(ofVec2f _center) {
    
    ofVec3f pos = getPosition();
    ofVec3f ur = (pos-_center).normalize();
    ofVec2f velocity2D = getVelocity();
    ofVec3f velocity3D = ofVec3f(velocity2D.x, velocity2D.y, 0);
    float tangentSpeed = ur.getCrossed( velocity3D ).z; // can be negative or positive
    // cout << tangentSpeed << endl; // between 0.001 and 1
    float maxSpeed = 0.7, minSpeed = -0.7;
    
    // 1) Linear map:
    float colormaping  = ofMap(tangentSpeed, minSpeed, maxSpeed, 0, 255, true);
    ofColor c = ofColor::fromHsb(colormaping, 255.0, 255.0);
    
    //2)  log map:
    //    float eps = 1.0; //color unit
    //    float sigma = (maxSpeed-minSpeed)/log(255.0/eps);
    //    float logColorMap = 255.0*(1-exp(-(tangentSpeed-minSpeed)/sigma));
    //    cout << int(logColorMap) << endl;
    //    ofColor c = ofColor::fromHsb(logColorMap, 255.0, 255.0);
    
    setColor(c);
}

void Particle::setColorRadialVelocity(ofVec2f _center) {
    ofVec3f pos = getPosition();
    ofVec3f ur = (pos-_center).normalize();
    ofVec2f velocity2D = getVelocity();
    float radialSpeed = ur.dot(velocity2D);
    // cout <<radialSpeed << endl;
    ofColor c = ofColor::fromHsb(20.0*radialSpeed+127,255.0, 255.0);
    setColor(c);
}



// =======================================================================

void ParticleSystem::addParticle(Particle& _particle) {
    addParticle(_particle.getParticleData());
}

void ParticleSystem::addParticle(const Particle::ParticleData& _particleData) {
    addParticle(_particleData.posInit, _particleData.equilibriumPosition, _particleData.velInit, _particleData.mass,_particleData.charge, _particleData.size, _particleData.damp);
}

void ParticleSystem::addParticle(const ofVec2f& _pos, const ofVec2f& _equil, const ofVec2f& _vel, float _mass, float _charge, float _size, float _damp) {
    
    //auto particle = ofPtr<Particle>(); //std::make_shared<Particle>();
    auto particle = std::make_shared<Particle>();
    
    particle.get()->myData.posInit = _pos;
    particle.get()->myData.equilibriumPosition = _equil;
    particle.get()->myData.velInit = _vel;
    particle.get()->myData.mass = _mass;
    particle.get()->myData.charge = _charge;
    particle.get()->myData.size = _size;
    particle.get()->myData.damp = _damp;
    
    //virtual void setPhysics(float density, float bounce, float friction);
    // note: I modified this: now size and mass are independent: "density" is just mass
    particle.get()->setPhysics(_mass, 1.0, 0.0); //0.005, 0.1, 0.0); // note: friction is between objects, not over the "world surface"
    // note: first set physics, then setup (it seems)
    
    particle.get()->setup( ptr_box2d->getWorld(), _pos.x, _pos.y, _size); // note: I modified this: now size and mass are independent
    particle.get()->setVelocity(_vel); //)(0,0);
    
    particle.get()->setDamping(_damp);
    particle.get()->setMass(_mass);
    
    particle.get()->setFixedRotation(true); //true => no inertia moment
    particle.get()->setRotation(ofRandom(0,360.0));
    
    particle.get()->setBullet(true); //more precise calculation, but slower...
    
    particle.get()->setEventData(ofGetElapsedTimeMicros(), 0, true);
    
    particleArray.push_back(particle);
}

void ParticleSystem::addParticle(float x, float y, float ex, float ey, float vx, float vy,  float _mass, float _charge, float _size, float _damp) {
    addParticle(ofVec2f(x,y), ofVec2f(vx, vy), ofVec2f(ex, ey), _mass, _charge, _size, _damp);
}

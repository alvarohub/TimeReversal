#include "ofMain.h"
#include "ofApp.h"

int main() {
	ofSetupOpenGL(3000, 1000, OF_WINDOW);
    ofRunApp( new ofApp());
	//return ofRunApp(std::make_shared<ofApp>());
}

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "World.h"

using namespace std;

int main(int argc, char * argv[]) {
    // Check command line arguments
    if (argc < 7) {
        cout << "ERROR: Too few command line arguments." << endl;
        cout << "USAGE: > " << argv[0] << " recordfile world-x-width world-y-width world-z-width objectsInWorld numberOfFrames" << endl;
        return 1;
    }
    // Get parameters from command line
    char * filename = argv[1];
    double worldSize[3];
    worldSize[0] = atof(argv[2]);
    worldSize[1] = atof(argv[3]);
    worldSize[2] = atof(argv[4]);
    unsigned totalObjects = atoi(argv[5]);
    unsigned nFrames = atoi(argv[6]);

    // seed randomizer
    srand (time(NULL));
 
    // Create file
    ofstream fout;
    fout.open(filename, std::ios::out|std::ios::binary);
 
    // Create and initialize the objects in the world
    World world(&fout, worldSize);

    // Add objects to the world
    world.createObjects(totalObjects);
 
  
    // Write the file header
    short headerLengthInBytes = 6;
    fout.write((char*)&headerLengthInBytes, UNSIGNEDSHORT);
    fout.write((char*)&nFrames, UNSIGNED);

    // Initial report
    clock_t startTime, endTime;
    cout << "----------------------------------" << endl;
    cout << "Producing animation data using serial implementation" << endl;
    cout << "  Generating " << nFrames << " frames" << endl;
    cout << "  World contains " << totalObjects << " objects" << endl;
    cout << "  World size: " << worldSize[0] << "x" << worldSize[1] << "x" << worldSize[2] << endl; 
    startTime = clock()   
    // Initial state
    world.record_frame();
    for (int i=1; i<nFrames; i++) {
        // Advance position and velocity one full step
        // -- Simple explicit Euler steps
        world.advance_full_step(); 

        // Collision detection
        world.detect_collisions();

        // Record the new frame
        world.record_frame();
    }
    
    // End report
    endTime = clock();
    double elapsedTime = double(endTime - startTime)/CLOCKS_PER_SEC;
    cout << "Finished producing timeseries data in file: " << filename << "." << endl;
    cout << "  Elapsed time: " << elapsedTime << " seconds" << endl;
    
    // Close file
    fout.close();

    return 0;
}

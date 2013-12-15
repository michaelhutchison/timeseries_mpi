#include <iostream>
#include <cstdlib>
#include "Slice_serial.h"

using namespace std;

int main(int argc, char * argv[]) {
    // Check command line arguments
    if (argc < 10) {
        cout << "ERROR: Too few command line arguments." << endl;
        cout << "USAGE: > " << argv[0] << " recordfile \\"
                                       << " world-x-width world-y-width world-z-width \\"
                                       << " x-rows y-rows z-rows \\"
                                       << " objectsInWorld numberOfFrames" << endl;
        return 1;
    }
    // Get parameters from command line
    char * filename = argv[1];
    double worldSize[3];
    worldSize[0] = atof(argv[2]);
    worldSize[1] = atof(argv[3]);
    worldSize[2] = atof(argv[4]);
    unsigned sliceRows[3];
    sliceRows[0] = atof(argv[5]);
    sliceRows[1] = atof(argv[6]);
    sliceRows[2] = atof(argv[7]);
    unsigned totalObjects = atoi(argv[8]);
    unsigned nFrames = atoi(argv[9]);

    // seed randomizer
    srand (time(NULL));
 
    // Create file
    ofstream fout;
    bool writeOn = true;
    if (writeOn) fout.open(filename, std::ios::out|std::ios::binary);
 
    // Create and initialize Slices
    unsigned size = sliceRows[0] * sliceRows[1] * sliceRows[2];
    unsigned objectsPerSlice = totalObjects / size;
    vector<Slice_serial *> slices;
    for (unsigned rank=0; rank<size; rank++) {
        // Create the slice
        Slice_serial * slice = new Slice_serial(rank, size, &fout, worldSize, sliceRows); 
        slices.push_back(slice);
        // Add objects to the slice
        slice->createObjects(objectsPerSlice);
        slice->setTotalObjects(size * objectsPerSlice); 
    }
    for (unsigned rank=0; rank<size; rank++) {
        slices[rank]->set_neighbors(slices);
    }
    // Write the file header
    if (writeOn) {
        short headerLengthInBytes = 6;
        fout.write((char*)&headerLengthInBytes, TYPEUNSIGNEDSHORT);
        fout.write((char*)&nFrames, TYPEUNSIGNED);
    }
    // Initial report
    clock_t startTime, endTime;
    cout << "----------------------------------" << endl;
    cout << "Producing animation data using serial implementation" << endl;
    cout << "  Generating " << nFrames << " frames" << endl;
    cout << "  World contains " << totalObjects << " objects" << endl;
    cout << "  World size: " << worldSize[0] << "x" << worldSize[1] << "x" << worldSize[2] << endl; 
    cout << "  Slice array: " << sliceRows[0] << "x" << sliceRows[1] << "x" << sliceRows[2] << endl;
    startTime = clock();   

    // Initial state
    if (writeOn) {
        for (unsigned i=0; i<slices.size(); i++)
            slices[i]->record_frame(); 
    }
    for (unsigned i=1; i<nFrames; i++) {
       
        // Exchange ghost objects
        for (unsigned i=0; i<slices.size(); i++)
            slices[i]->clear_ghost_objects();
        for (unsigned i=0; i<slices.size(); i++)
            slices[i]->exchange_ghost_objects();

        // Advance position and velocity one full step
        // -- Simple explicit Euler steps
        for (unsigned i=0; i<slices.size(); i++)
            slices[i]->advance_full_step(); 

        // Collision detection
        for (unsigned i=0; i<slices.size(); i++)
            slices[i]->detect_collisions();

        // MPI Communication step 2: update of remote and notification of new bodies
        for (unsigned i=0; i<slices.size(); i++)
            slices[i]->exchange_objects();

        // Record each slice's part of the new frame
        if (writeOn) {
            for (unsigned i=0; i<slices.size(); i++)
                slices[i]->record_frame();
        }
    }

    // Final report
    endTime = clock();
    double elapsedTime = double(endTime - startTime)/CLOCKS_PER_SEC;
    cout << "Finished producing timeseries data in file: " << filename << "." << endl;
    cout << "  Elapsed time: " << elapsedTime << " seconds" << endl;
    
    // Close file
    if (writeOn) fout.close();

    // Clean up dynamic memory
    for (unsigned i=0; i<slices.size(); i++)
        delete slices[i];

    return 0;
}

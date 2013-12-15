#include <iostream>
#include <cstdlib>
#include "Slice.h"
#include "mpi.h"

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
    sliceRows[0] = atoi(argv[5]);
    sliceRows[1] = atoi(argv[6]);
    sliceRows[2] = atoi(argv[7]);
    unsigned totalObjects = atoi(argv[8]);
    unsigned nFrames = atoi(argv[9]);

    // Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // seed randomizer
    srand (time(NULL) + rank);
 
    // Create file
    MPI_File fh;
    bool writeOn = true;
    if (writeOn) {
        MPI_File_open(MPI_COMM_WORLD,
                      filename,
                      MPI_MODE_WRONLY | MPI_MODE_CREATE,
                      MPI_INFO_NULL,
                      &fh);
    }
    // Create and initialize a Slice
    Slice slice(rank, size, &fh, worldSize, sliceRows); 
    // Add objects to the slice
    unsigned objectsPerSlice = totalObjects / size;
    slice.createObjects(objectsPerSlice);
    slice.setTotalObjects(size * objectsPerSlice); 
    
    // Write the file header
    if (writeOn && rank == 0) {
        MPI_Status status;
        // Write header length in bytes
        short headerLengthInBytes = 6;
        MPI_File_write( fh,
                        &headerLengthInBytes,
                        1,
                        MPI_SHORT,
                        &status);
        // Write number of frames
        MPI_File_write( fh,
                                &nFrames,
                                1,
                                MPI_UNSIGNED,
                                &status);
    }
    
    /* Parallel FFD Algorithm */
    double startTime, endTime;
    if (rank == 0) {
        cout << "----------------------------------" << endl;
        cout << "Producing animation data using mpi" << endl;
        cout << "  Using " << size << " processors" << endl;
        cout << "  Generating " << nFrames << " frames" << endl;
        cout << "  World contains " << size * objectsPerSlice << " objects" << endl;
        cout << "  World size: " << worldSize[0] << "x" << worldSize[1] << "x" << worldSize[2] << endl; 
        cout << "  Slice array: " << sliceRows[0] << "x" << sliceRows[1] << "x" << sliceRows[2] << endl;
            
        startTime = MPI_Wtime();   
    }
    // Initial state
    if (writeOn) slice.record_frame();
    for (int i=1; i<nFrames; i++) {
        // MPI Communication step 1: force synchronization
        // -- Not implemented. This is where one would update
        //    a force field, such as a fluid flow, and coordinate
        //    forces at boundaries between processes.
       
        // Exchange ghost objects
        slice.exchange_ghost_objects();

        // Advance position and velocity one full step
        // -- Simple explicit Euler steps
        slice.advance_full_step(); 

        // Collision detection
        slice.detect_collisions();

        // MPI Communication step 2: update of remote and notification of new bodies
        slice.exchange_objects();

        // Record this slice's part of the new frame
        if (writeOn) slice.record_frame();
    }
    if (rank == 0) {
        endTime = MPI_Wtime();
        double elapsedTime = endTime - startTime;
        cout << "Finished producing timeseries data in file: " << filename << "." << endl;
        cout << "  Elapsed time: " << elapsedTime << " seconds" << endl;
    }

    // Close file
    if (writeOn) MPI_File_close(&fh);

    MPI_Finalize();

    return 0;
}

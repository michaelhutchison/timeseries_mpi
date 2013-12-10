#include <iostream>
#include "Slice.h"
#include "mpi.h"

#define OBJECTSPERSLICE 200
#define FRAMES 400

using namespace std;

int main(int argc, char * argv[]) {
    // Check command line arguments
    if (argc < 2) {
        cout << "ERROR: Too few command line arguments." << endl;
        cout << "USAGE: > " << argv[0] << " filename " << endl;
        return 1;
    }
    // Get filename from command line
    char * filename = argv[1];

    // Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // seed randomizer
    srand (time(NULL) + rank);
 
    // Create file
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD,
                  filename,
                  MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL,
                  &fh);

    // Create a slice
    Slice slice(rank, size, &fh); 
    
    // Add objects to the slice
    slice.createObjects(OBJECTSPERSLICE);

    // Write the file header
    if (rank == 0) {
        MPI_Status status;
        // Write header length in bytes
        short headerLengthInBytes = 6;
        MPI_File_write( fh,
                        &headerLengthInBytes,
                        1,
                        MPI_SHORT,
                        &status);
        // Write number of frames
        unsigned nFrames = FRAMES;
        MPI_File_write( fh,
                                &nFrames,
                                1,
                                MPI_UNSIGNED,
                                &status);
    }


    /* Parallel FFD Algorithm */
    double startTime, endTime;
    if (rank == 0) {
        cout << "Producing animation data" << endl;
        cout << "Generating " << FRAMES << " frames" << endl;
        cout << "World contains " << size * OBJECTSPERSLICE << " objects" << endl;
        startTime = MPI_Wtime();   
    }
    // Initial state
    slice.record_frame();
    for (int i=1; i<FRAMES; i++) {
        // MPI Communication step 1: force synchronization
        // -- Not implemented. This is where one would update
        //    a force field, such as a fluid flow, and coordinate
        //    forces at boundaries between processes.
        
        // Advance position and velocity one full step
        // -- Simple explicit Euler steps
        slice.advance_full_step(); 

        // MPI Communication step 2: update of remote and notification of new bodies
        slice.exchange_objects();

        // Collision detection
        slice.handle_collisions();

        // Record this slice's part of the new frame
        slice.record_frame();
    }
    if (rank == 0) {
        endTime = MPI_Wtime();
        double elapsedTime = endTime - startTime;
        cout << "Finished producing timeseries data in " << filename << "." << endl;
        cout << "Elapsed time: " << elapsedTime << " seconds" << endl;
    }

    // Close file
    MPI_File_close(&fh);

    MPI_Finalize();

    return 0;
}

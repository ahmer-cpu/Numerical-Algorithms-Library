#include <iostream>

using namespace std;

int main()
{
    // This C++ code is only to demonstrate the algorithms employed and expose the mathematical structure.
    // We also use this for initial testing and validation (verifying code and algorithmic coherence).
    // Comprehensive testing drivers and reporting will be built in Python using NumPy and Matplotlib functionalities.
    // No exception handling is done in this assigment, e.g. validating user input!

    // User input for array initialization
    int CDIM = 4;
    int RDIM = 4;
 
    // We can add RDIM and CDIM variables here, but this is bad practice in C++.
    // Change this manually everytime you change the value for the variables above!
    int L[4][4];

    // Printing the array to show the effect of undefined initialization
    cout << "\nUninitialized array:\n\n";
    for (int i=0; i < RDIM; i++)
    {
        for (int j=0; j < CDIM; j++) {
            cout << L[i][j] << " ";
        }
        cout << endl;
    }

    // Garbage values are clearly stored automatically for local static arrays
    // We will consider garbage values "undefined"


    // User initialization for a lower unit triangular matrix (Toy example)
    // Indexing start from 1, i.e. the second row has the first non-trivial storage
    for (int i=1; i < RDIM; i++)
    {
        for (int j=0; j <= i-1; j++) {
            L[i][j] = i+j+1;
        }
    }
    // We only consider these S = [n(n-1)]/2 storage locations as used

    // Printing the array to show the effect of undefined initialization and our particular initializations
    cout << "\nAfter our values are stored:\n\n";
    for (int i=0; i < RDIM; i++)
    {
        for (int j=0; j < CDIM; j++) {
            cout << L[i][j] << " ";
        }
        cout << endl;
    }

    // Printing array to validate input
    cout << endl << "Relevant storage values and position:\n\n" <<  "L =\n[";
    for (int i=1; i < RDIM; i++)
    {
        for (int j = 0; j <= i-1; j++)
        {
            cout << L[i][j] << " ";
        }
        if (i == RDIM-1)
        {
            cout << "\b]";
        }
        cout << endl;
        cout << " ";
    }


    // Sub=routine 1 (Unit Lower Triangular Matrix - Vector (2D)) //

    // Toy Example setup (L is as above)
    int v[4] = {21,33,41,5};

    int w[4] = {0};                                     // Initializes a zeros array of size RDIM
    
    w[0] = v[0];
    for (int i=1; i < RDIM; i++)                        // Loop runs as many times as size(v)
    {
        for (int j=0; j < i; j++)
        {
            w[i] = w[i] + ((L[i][j])*v[j]);             // Inner-product until the diagonal of L (exclusive of the diagonal value)
        }
        w[i] = w[i] + v[i];                             // Add back the vector because of the 1s diagonal
    }

    // Printing the output vector
    cout << endl << "\nOutput vector after Lv multiplication (2d):\n" << "w = [";
    for (int i=0; i<RDIM; i++)
    {
        cout << w[i] << " ";
    }
    cout << "\b]";
    cout << endl;
    // Analytic Answer: [21 75 236 500]


    // Sub-routine 2 (Unit Lower Triangular Matrix - Vector (1D)) //

    // We convert the same L matrix to a 1d array in the compressed row form
    int L_size = ((RDIM)*(RDIM-1))/2;
    // Manully input this value for the dimension of the 1d array
    int L_1d[6] = {0};

    int i = 0;
    for (int j = 1; j < RDIM; j++)
    {
        for (int k = 0; k <= j-1; k++)
        {
            L_1d[i] = L[j][k];                          // Assigning the values to the 1d array
            i++;
        }
    }

    // Printing L in the new compressed row form
    cout << endl << "1D array storage:" << endl << "L' = [";
    for (int i=0; i < L_size; i++)
    {
        cout << L_1d[i] << " ";
    }
    cout << "\b]";

    int w_2[4] = {0};                                   // Initializing the output vector

    w_2[0] = v[0];
    int counter = 0;                                    // We need this because the indexing of the 1d array starts at how many values have been used
    for (int i=1; i < RDIM; i++)
    {
        int index = counter;                            // Where we start indexing the 1d array
        for (int j=0; j < i; j++)
        {
            w_2[i] = w_2[i] + (L_1d[index])*(v[j]);     // The indexing of the vector remains 0:i-1 
            index++;
            counter++;
        }
       w_2[i] = w_2[i] + v[i];                          // Add back the vector because of the 1s diagonal
    }

    // Printing output vector
    cout << endl << "\nOutput vector after Lv multiplication (1d):\n" << "w_2 = [";
    for (int i=0; i < RDIM; i++)
    {
        cout << w_2[i] << " ";
    }
    cout << "\b]";
    cout << endl;


    // Sub-routine 3 (Unit Lower Triangular Banded Matrix - Vector) //

    int NEWRDIM = 5;
    int NEWCDIM = 5;
                                
    // We assume that we have a row compressed form of the banded matrix 
    // i.e. two values for each row except row 1 that has one non-trivial storage
    // There are always 2n-3 storage locations for a unit banded matrix.

    int L_banded[7] = {3,4,5,6,7,8,9};                              // Values of \lambda_{21,31,32,42,43,53,54},  for example and n=5
    int v_banded[5] = {2,2,3,4,5};

    int SIZE_banded = 2*(NEWRDIM) - 3;

    // Printing L_banded in the compressed row form
    cout << endl << "1D array storage of banded matrix:" << endl << "L_banded  = [";
    for (int i=0; i < SIZE_banded; i++)
    {
        cout << L_banded[i] << " ";
    }
    cout << "\b]\n";

    int w_banded[5] = {0};

    w_banded[0] = v_banded[0];                                      // Simple assignment due to row 1 [1 0 0 0 ...]
    w_banded[1] = (L_banded[0])*(v_banded[0]);                      // Only inner product with just one multiplication (ignoring diagonal 1s for now)
    int index = 1;
    for (int i=2; i < NEWRDIM; i++)
    {
        w_banded[i] = w_banded[i] + ((L_banded[index])*(v_banded[i-2]));
        w_banded[i] = w_banded[i] + ((L_banded[index+1])*(v_banded[i-1]));
        index = index + 2;
    }

    for( int i=1; i < NEWRDIM; i++)
    {
        w_banded[i] = w_banded[i] + v_banded[i];                    // Adding back the vector just as before
    }

    // Printing output vector
    cout << "\nOutput vector after Lv multiplication (banded):\n" << "w_banded = [";
    for (int i=0; i < NEWRDIM; i++)
    {
        cout << w_banded[i] << " ";
    }
    cout << "\b]";
    // Analytic solution: w = [2 8 21 37 65]


    // Sub-routine 4 (Lower Matrix - Upper Matrix (Single 2D)) //
    int n = 4;                                                              // Dimension of LU combined
    int LU[4][4] = {                                                        // LU is stored in a single 2D array (L is unit)
        {1,2,3,4},
        {4,5,6,7},
        {8,9,10,11},
        {12,13,14,15}
        };


    int M[4][4] = {};                                                       // Zeros array (n^2 non-trivial storage)            
    for (int i=0; i<n; i++)                                                 // The first row and column of M is pre-determined
    {
        M[0][i] = LU[0][i];
    }
    for(int i=1; i<n; i++)
    {
        M[i][0] = (LU[i][0])*(LU[0][0]);                                    // The first column of M is pre-determined
    }
    for (int i=1; i < n; i++)                                           
    {
        for (int j=1; j<n; j++)
        {
            int min = j+1;
            if((min-i)>0)
            {
                min = i;
            }
            int index = 0;
            for (int counter = 0; counter < min; counter++)             // min tells us how many times we need to run this loop exactly
            {   
                M[i][j] = M[i][j] + (LU[i][counter])*(LU[counter][j]);  // Inner product until (exclusive) we encounter 1 (this is of course implicit)
                index++;
            }
            if(i <= j)
            {
                M[i][j] = M[i][j] + LU[index][j];                             // Implicit multiplication with 1 of the remaining term in the inner product
            }
        }
    }

    // Printing output Matrix
    cout << endl << "\nOutput Matrix of L-U multiplication:\n";
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            cout << M[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}